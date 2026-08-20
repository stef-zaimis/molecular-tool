"""
End-to-end tests for the `project.*` RPC surface.

These drive the service exactly as Electron does — a child process, one JSON
request per line on stdin, one JSON response per line on stdout — rather than
calling `ProjectService` in-process. That is the point: it is the only way to
catch a method that is missing from the dispatch table, a payload that is not
JSON-serialisable, or a result whose shape does not match the TypeScript
contract the renderer is compiled against.

Everything on stdout must be protocol. Diagnostics belong on stderr, and a
stray `print` in the scientific code would break these tests, which is
intentional.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

RECORDS = [
    ("focal_1|AU|Target", "AACCGGTT"),
    ("focal_2|GB|Target", "AACCGGTA"),
    ("other_1|FR|Contrast", "AGCTGGTT"),
    ("other_2|DE|Contrast", "CGTTGGTT"),
]

RUN_OPTIONS = {
    "ignoreGaps": False,
    "giveBenefitOfDoubtToAmbiguousBases": False,
    "minCandidateSize": 1,
    "maxCandidateSize": 2,
}


def write_fasta(path: Path, records) -> Path:
    path.write_text(
        "".join(f">{header}\n{sequence}\n" for header, sequence in records), encoding="utf-8"
    )
    return path


class ServiceClient:
    """A live service subprocess speaking the newline-delimited JSON protocol."""

    def __init__(self) -> None:
        self.proc = subprocess.Popen(
            [sys.executable, "-m", "molecular_diagnosis.service"],
            cwd=str(REPO_ROOT),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            bufsize=1,
        )
        self._counter = 0
        ready = json.loads(self.proc.stdout.readline())
        assert ready["id"] == "ready" and ready["ok"], ready

    def call(self, method: str, params: dict | None = None) -> dict:
        self._counter += 1
        request = {"id": str(self._counter), "method": method, "params": params or {}}
        self.proc.stdin.write(json.dumps(request) + "\n")
        self.proc.stdin.flush()
        line = self.proc.stdout.readline()
        if not line:
            raise AssertionError("the service closed stdout:\n" + self.proc.stderr.read())
        return json.loads(line)

    def ok(self, method: str, params: dict | None = None) -> dict:
        response = self.call(method, params)
        assert response["ok"], response.get("error")
        return response["result"]

    def close(self) -> None:
        try:
            self.proc.stdin.close()
            self.proc.wait(timeout=10)
        except Exception:
            self.proc.kill()


@pytest.fixture()
def service():
    client = ServiceClient()
    yield client
    client.close()


@pytest.fixture()
def fasta(tmp_path):
    return write_fasta(tmp_path / "beetles.fasta", RECORDS)


@pytest.fixture()
def project_dir(tmp_path):
    return tmp_path / "project"


def open_with_focal_set(service, project_dir, fasta):
    """Open a project, link the FASTA, and build a two-entry focal set."""
    service.ok("project.create", {"projectDir": str(project_dir), "title": "Test"})
    file_id = service.ok("project.linkFasta", {"path": str(fasta)})["fastaFileId"]
    set_id = service.ok("project.createFocalSet", {"title": "Targets"})["focalSetId"]
    service.ok("project.addFocalEntries", {"focalSetId": set_id, "query": "Target"})
    return file_id, set_id


# ---------------------------------------------------------------------------
# Lifecycle
# ---------------------------------------------------------------------------


def test_creating_a_project_reports_capabilities_and_no_sources(service, project_dir):
    result = service.ok("project.create", {"projectDir": str(project_dir), "title": "Fresh"})

    assert result["metadata"]["title"] == "Fresh"
    assert result["sources"] == []
    capabilities = result["capabilities"]
    assert capabilities["sqliteVersion"]
    # Whether FTS is available depends on the build; what must hold is that the
    # reported capability and the accelerator agree.
    assert capabilities["acceleratedSearch"] == (
        capabilities["fts5"] and capabilities["trigram"]
    )


def test_the_database_file_lives_in_the_project_directory(service, project_dir):
    service.ok("project.create", {"projectDir": str(project_dir)})
    assert (project_dir / "project.sqlite").is_file()
    assert (project_dir / "outputs").is_dir()


def test_calling_a_project_method_before_opening_fails_cleanly(service):
    response = service.call("project.listFocalSets")
    assert response["ok"] is False
    assert response["error"]["code"] == "NO_PROJECT_OPEN"


def test_opening_a_folder_without_a_project_does_not_create_one(service, project_dir):
    """
    "Open Existing Project" must not manufacture one.

    Creating an empty project.sqlite in whatever folder the user picked makes a
    mis-click look exactly like "all my work disappeared".
    """
    project_dir.mkdir(parents=True, exist_ok=True)

    response = service.call("project.open", {"projectDir": str(project_dir)})

    assert response["ok"] is False
    assert response["error"]["code"] == "PROJECT_NOT_FOUND"
    assert not (project_dir / "project.sqlite").exists()


def test_creating_a_project_where_one_already_exists_is_refused(service, project_dir):
    service.ok("project.create", {"projectDir": str(project_dir), "title": "First"})
    service.ok("project.close")

    response = service.call("project.create", {"projectDir": str(project_dir), "title": "Second"})

    assert response["ok"] is False
    assert response["error"]["code"] == "PROJECT_ALREADY_EXISTS"
    # The original project is untouched and still openable under its own title.
    assert service.ok("project.open", {"projectDir": str(project_dir)})["metadata"]["title"] == (
        "First"
    )


def test_opening_an_existing_project_never_renames_it(service, project_dir):
    service.ok("project.create", {"projectDir": str(project_dir), "title": "Original"})
    service.ok("project.close")

    reopened = service.ok("project.open", {"projectDir": str(project_dir)})
    assert reopened["metadata"]["title"] == "Original"

    renamed = service.ok("project.setTitle", {"title": "Renamed"})
    assert renamed["metadata"]["title"] == "Renamed"

    service.ok("project.close")
    assert service.ok("project.open", {"projectDir": str(project_dir)})["metadata"]["title"] == (
        "Renamed"
    )


# ---------------------------------------------------------------------------
# Persistence across sessions
# ---------------------------------------------------------------------------


def test_a_project_reopens_with_its_links_and_focal_sets(service, project_dir, fasta):
    _file_id, set_id = open_with_focal_set(service, project_dir, fasta)
    service.ok("project.close")

    result = service.ok("project.open", {"projectDir": str(project_dir)})

    # `unverified`, not `current`: opening is deliberately stat-only, so a
    # project linking a multi-gigabyte alignment opens instantly. Size and
    # mtime match, which is evidence but not proof, and the UI shows it as
    # provisional rather than confidently current.
    assert [source["state"] for source in result["sources"]] == ["unverified"]
    assert result["sources"][0]["available"] is True

    # Asking for proof promotes it, and costs one hash for the whole session.
    strong = service.ok("project.refreshSources", {"strong": True})["sources"]
    assert strong[0]["state"] == "current"

    focal_sets = service.ok("project.listFocalSets")["focalSets"]
    assert len(focal_sets) == 1
    assert focal_sets[0]["id"] == set_id
    assert [entry["header"] for entry in focal_sets[0]["entries"]] == [
        "focal_1|AU|Target",
        "focal_2|GB|Target",
    ]


def test_a_missing_file_is_surfaced_on_open_without_losing_the_focal_set(
    service, project_dir, fasta, tmp_path
):
    """The case the whole source-status UI exists for."""
    open_with_focal_set(service, project_dir, fasta)
    service.ok("project.close")
    os.replace(fasta, tmp_path / "moved.fasta")

    result = service.ok("project.open", {"projectDir": str(project_dir)})

    source = result["sources"][0]
    assert source["state"] == "missing"
    assert source["available"] is False
    assert source["indexUsable"] is False
    # A message the UI can show verbatim, not a code it has to translate.
    assert source["message"]

    # The link and the focal set survive: the file is gone, the project is not.
    focal_sets = service.ok("project.listFocalSets")["focalSets"]
    assert len(focal_sets[0]["entries"]) == 2


def test_relinking_a_moved_file_rebuilds_nothing(service, project_dir, fasta, tmp_path):
    file_id, _set_id = open_with_focal_set(service, project_dir, fasta)
    moved = tmp_path / "moved.fasta"
    os.replace(fasta, moved)

    result = service.ok("project.relinkFasta", {"fastaFileId": file_id, "path": str(moved)})

    assert result["identical"] is True
    assert result["reindexed"] is False
    assert result["source"]["state"] == "current"


# ---------------------------------------------------------------------------
# Change detection and recovery
# ---------------------------------------------------------------------------


def test_editing_the_file_is_detected_and_recovered_by_reindexing(
    service, project_dir, fasta
):
    file_id, set_id = open_with_focal_set(service, project_dir, fasta)

    write_fasta(fasta, [RECORDS[0], RECORDS[2], RECORDS[3]])

    sources = service.ok("project.refreshSources", {"strong": True})["sources"]
    assert sources[0]["state"] == "stale"

    source = service.ok("project.reindexFasta", {"fastaFileId": file_id})["source"]
    assert source["state"] == "current"
    assert source["sequenceCount"] == 3

    # The focal entry whose record was deleted is now reported missing; the
    # other one is still located.
    presence = service.ok(
        "project.focalPresence", {"focalSetId": set_id, "selectedFastaFileId": file_id}
    )["entries"]
    assert sorted(entry["state"] for entry in presence) == ["missing", "present_current"]


def test_a_cheap_refresh_reports_state_without_hashing(service, project_dir, fasta):
    """The call the UI makes on window focus and on every screen change."""
    open_with_focal_set(service, project_dir, fasta)

    sources = service.ok("project.refreshSources", {"strong": False})["sources"]
    # Already proven during linking, so it stays confidently current.
    assert sources[0]["state"] == "current"


# ---------------------------------------------------------------------------
# Search and focal sets
# ---------------------------------------------------------------------------


def test_search_returns_exact_headers_in_file_order(service, project_dir, fasta):
    service.ok("project.create", {"projectDir": str(project_dir)})
    service.ok("project.linkFasta", {"path": str(fasta)})

    hits = service.ok("project.searchHeaders", {"query": "Target"})["hits"]
    assert [hit["header"] for hit in hits] == ["focal_1|AU|Target", "focal_2|GB|Target"]
    assert [hit["ordinal"] for hit in hits] == [0, 1]


def test_adding_and_removing_focal_entries_by_query(service, project_dir, fasta):
    _file_id, set_id = open_with_focal_set(service, project_dir, fasta)

    removed = service.ok(
        "project.removeFocalEntries", {"focalSetId": set_id, "query": "focal_1"}
    )["removed"]
    assert removed == ["focal_1|AU|Target"]

    focal_sets = service.ok("project.listFocalSets")["focalSets"]
    assert [entry["header"] for entry in focal_sets[0]["entries"]] == ["focal_2|GB|Target"]


def test_the_whole_focal_set_api_is_reachable_over_the_wire(service, project_dir, fasta):
    """
    Every method the focal UI needs, exercised through the real dispatch table.

    A missing entry here is exactly the failure this file exists to catch: the
    renderer compiles fine against a method the service does not have.
    """
    _file_id, set_id = open_with_focal_set(service, project_dir, fasta)

    assert service.ok("project.getFocalSet", {"focalSetId": set_id})["focalSet"]["title"] == (
        "Targets"
    )

    renamed = service.ok(
        "project.renameFocalSet", {"focalSetId": set_id, "title": "Renamed"}
    )["focalSet"]
    assert renamed["title"] == "Renamed"

    replaced = service.ok(
        "project.replaceFocalEntries",
        {"focalSetId": set_id, "text": "focal_1|AU|Target ; other_1|FR|Contrast"},
    )
    assert replaced["headers"] == ["focal_1|AU|Target", "other_1|FR|Contrast"]

    locked = service.ok(
        "project.setFocalSetLocked", {"focalSetId": set_id, "locked": True}
    )["focalSet"]
    assert locked["locked"] is True

    unlocked = service.ok(
        "project.setFocalSetLocked", {"focalSetId": set_id, "locked": False}
    )["focalSet"]
    assert unlocked["locked"] is False

    assert service.ok("project.deleteFocalSet", {"focalSetId": set_id})["deleted"] == set_id
    assert service.ok("project.listFocalSets")["focalSets"] == []


def test_a_locked_focal_set_is_enforced_by_the_backend(service, project_dir, fasta):
    """Not by a disabled button: the refusal has to survive a direct call."""
    _file_id, set_id = open_with_focal_set(service, project_dir, fasta)
    service.ok("project.setFocalSetLocked", {"focalSetId": set_id, "locked": True})

    for method, params in (
        ("project.renameFocalSet", {"focalSetId": set_id, "title": "x"}),
        ("project.deleteFocalSet", {"focalSetId": set_id}),
        ("project.addFocalEntries", {"focalSetId": set_id, "query": "Contrast"}),
        ("project.removeFocalEntries", {"focalSetId": set_id, "query": "focal_1"}),
        ("project.replaceFocalEntries", {"focalSetId": set_id, "text": "anything"}),
    ):
        response = service.call(method, params)
        assert response["ok"] is False, method
        assert response["error"]["code"] == "FOCAL_SET_LOCKED", method

    entries = service.ok("project.listFocalSets")["focalSets"][0]["entries"]
    assert [entry["header"] for entry in entries] == [
        "focal_1|AU|Target",
        "focal_2|GB|Target",
    ]


def test_an_unknown_focal_set_id_is_a_coded_refusal(service, project_dir, fasta):
    open_with_focal_set(service, project_dir, fasta)

    response = service.call("project.getFocalSet", {"focalSetId": "not-a-real-id"})

    assert response["ok"] is False
    assert response["error"]["code"] == "UNKNOWN_FOCAL_SET"
    # A SQLite integrity error would have arrived with a traceback attached.
    assert "traceback" not in response["error"]


def test_an_empty_mutation_query_cannot_erase_the_focal_set(service, project_dir, fasta):
    _file_id, set_id = open_with_focal_set(service, project_dir, fasta)

    for method in ("project.addFocalEntries", "project.removeFocalEntries"):
        response = service.call(method, {"focalSetId": set_id, "query": "   "})
        assert response["ok"] is False, method
        assert response["error"]["code"] == "FOCAL_QUERY_EMPTY", method

    entries = service.ok("project.listFocalSets")["focalSets"][0]["entries"]
    assert len(entries) == 2


def test_search_reports_files_it_could_not_search(service, project_dir, fasta, tmp_path):
    service.ok("project.create", {"projectDir": str(project_dir)})
    service.ok("project.linkFasta", {"path": str(fasta)})
    gone = write_fasta(tmp_path / "gone.fasta", [("g1|Target", "AACCGGTT")])
    service.ok("project.linkFasta", {"path": str(gone)})
    os.remove(gone)

    result = service.ok("project.searchHeaders", {"query": "Target"})

    assert [hit["header"] for hit in result["hits"]] == [
        "focal_1|AU|Target",
        "focal_2|GB|Target",
    ]
    assert [entry["displayName"] for entry in result["unavailable"]] == ["gone.fasta"]


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------


def test_a_single_file_run_writes_into_the_project_outputs_directory(
    service, project_dir, fasta
):
    file_id, set_id = open_with_focal_set(service, project_dir, fasta)

    result = service.ok(
        "project.runMolecularDiagnosis",
        {
            "focalSetId": set_id,
            "fastaFileIds": [file_id],
            "singleFile": True,
            "options": RUN_OPTIONS,
        },
    )

    # Exact headers, not the substring that found them.
    assert result["focalHeaders"] == ["focal_1|AU|Target", "focal_2|GB|Target"]
    assert result["sequenceCount"] == 4
    for path in (result["outputs"]["reportTxt"], result["outputs"]["workbookXlsx"]):
        assert Path(path).is_file()
        assert Path(path).parent == project_dir / "outputs"
    # The same dmc block the non-project run returns, so the UI renders both
    # kinds of run through one code path.
    assert result["dmc"]["singleSites"]
    assert "diagnostics" in result["dmc"]


def test_a_single_file_run_is_refused_when_a_focal_entry_is_absent(
    service, project_dir, fasta
):
    file_id, set_id = open_with_focal_set(service, project_dir, fasta)
    write_fasta(fasta, [RECORDS[0], RECORDS[2], RECORDS[3]])
    service.ok("project.reindexFasta", {"fastaFileId": file_id})

    response = service.call(
        "project.runMolecularDiagnosis",
        {
            "focalSetId": set_id,
            "fastaFileIds": [file_id],
            "singleFile": True,
            "options": RUN_OPTIONS,
        },
    )

    assert response["ok"] is False
    assert response["error"]["code"] == "FOCAL_ENTRIES_NOT_IN_FILE"
    assert "focal_2|GB|Target" in response["error"]["detail"]


def test_an_expected_refusal_carries_no_traceback(service, project_dir, fasta):
    """
    A coded refusal is a decision, not a crash.

    Shipping a traceback for one trains the UI (and the reader) to treat
    ordinary outcomes as bugs.
    """
    _file_id, set_id = open_with_focal_set(service, project_dir, fasta)

    response = service.call(
        "project.runMolecularDiagnosis",
        {
            "focalSetId": set_id,
            "fastaFileIds": [],
            "singleFile": False,
            "options": RUN_OPTIONS,
        },
    )

    assert response["ok"] is False
    assert response["error"]["code"] == "NO_FILES_SELECTED"
    assert "traceback" not in response["error"]


def test_a_multi_file_run_combines_two_alignments(service, project_dir, fasta, tmp_path):
    file_id, set_id = open_with_focal_set(service, project_dir, fasta)
    second = write_fasta(
        tmp_path / "more.fasta",
        [("extra_1|ES|Contrast", "AACCGGTC"), ("extra_2|IT|Contrast", "AGCCGGTT")],
    )
    second_id = service.ok("project.linkFasta", {"path": str(second)})["fastaFileId"]

    result = service.ok(
        "project.runMolecularDiagnosis",
        {
            "focalSetId": set_id,
            "fastaFileIds": [file_id, second_id],
            "singleFile": False,
            "options": RUN_OPTIONS,
        },
    )

    # Merged in memory: no combined FASTA is written anywhere.
    assert result["sequenceCount"] == 6
    assert result["alignmentLength"] == 8
    assert sorted(result["fastaFileIds"]) == sorted([file_id, second_id])


def test_incompatible_alignment_lengths_are_refused_with_a_useful_detail(
    service, project_dir, fasta, tmp_path
):
    file_id, set_id = open_with_focal_set(service, project_dir, fasta)
    short = write_fasta(tmp_path / "short.fasta", [("short_1|X", "AACC")])
    short_id = service.ok("project.linkFasta", {"path": str(short)})["fastaFileId"]

    response = service.call(
        "project.runMolecularDiagnosis",
        {
            "focalSetId": set_id,
            "fastaFileIds": [file_id, short_id],
            "singleFile": False,
            "options": RUN_OPTIONS,
        },
    )

    assert response["ok"] is False
    assert response["error"]["code"] == "INCOMPATIBLE_ALIGNMENT_LENGTHS"
    # Names and lengths, so the user can tell which file is the odd one out.
    assert "short.fasta" in response["error"]["detail"]


def test_a_run_is_refused_when_a_selected_file_has_duplicate_headers(
    service, project_dir, fasta, tmp_path
):
    file_id, set_id = open_with_focal_set(service, project_dir, fasta)
    dupes = write_fasta(
        tmp_path / "dupes.fasta",
        [("d1|X", "AACCGGTT"), ("d1|X", "TTTTGGTT"), ("d2|X", "GGGGGGTT")],
    )
    dupes_id = service.ok("project.linkFasta", {"path": str(dupes)})["fastaFileId"]

    response = service.call(
        "project.runMolecularDiagnosis",
        {
            "focalSetId": set_id,
            "fastaFileIds": [file_id, dupes_id],
            "singleFile": False,
            "options": RUN_OPTIONS,
        },
    )

    assert response["ok"] is False
    assert response["error"]["code"] == "DUPLICATE_HEADER_IN_FILE"
    assert "dupes.fasta" in response["error"]["detail"]
    # Refusing to analyse it is not refusing to keep it: the link survives.
    assert dupes_id in [
        source["fastaFileId"] for source in service.ok("project.refreshSources")["sources"]
    ]


def test_an_all_files_run_is_refused_when_an_entry_is_outside_the_scope(
    service, project_dir, fasta, tmp_path
):
    """
    The multi-file half of the completeness rule.

    `elsewhere|QQ|Target` is a member of the focal set and lives in a project
    file that is NOT part of this run, so the run would have silently used a
    two-sequence focal group where the user specified three.
    """
    file_id, set_id = open_with_focal_set(service, project_dir, fasta)
    inside = write_fasta(tmp_path / "inside.fasta", [("in_1|X", "AAAAAAAA")])
    inside_id = service.ok("project.linkFasta", {"path": str(inside)})["fastaFileId"]
    outside = write_fasta(tmp_path / "outside.fasta", [("elsewhere|QQ|Target", "AACCGGTT")])
    service.ok("project.linkFasta", {"path": str(outside)})
    service.ok("project.addFocalEntries", {"focalSetId": set_id, "query": "Target"})

    response = service.call(
        "project.runMolecularDiagnosis",
        {
            "focalSetId": set_id,
            "fastaFileIds": [file_id, inside_id],
            "singleFile": False,
            "options": RUN_OPTIONS,
        },
    )

    assert response["ok"] is False
    assert response["error"]["code"] == "FOCAL_ENTRIES_NOT_IN_SCOPE"
    assert "elsewhere|QQ|Target" in response["error"]["detail"]


def test_focal_presence_can_be_scoped_to_the_files_a_run_will_read(
    service, project_dir, fasta, tmp_path
):
    file_id, set_id = open_with_focal_set(service, project_dir, fasta)
    outside = write_fasta(tmp_path / "outside.fasta", [("elsewhere|QQ|Target", "AACCGGTT")])
    service.ok("project.linkFasta", {"path": str(outside)})
    service.ok("project.addFocalEntries", {"focalSetId": set_id, "query": "Target"})

    scoped = service.ok(
        "project.focalPresence", {"focalSetId": set_id, "fastaFileIds": [file_id]}
    )["entries"]
    states = {entry["header"]: entry["state"] for entry in scoped}
    assert states["elsewhere|QQ|Target"] == "missing"

    project_wide = service.ok("project.focalPresence", {"focalSetId": set_id})["entries"]
    wide_states = {entry["header"]: entry["state"] for entry in project_wide}
    assert wide_states["elsewhere|QQ|Target"] == "present_current"


def test_a_failed_open_leaves_the_current_project_open(service, project_dir, fasta, tmp_path):
    """Picking the wrong folder must not close the project already in use."""
    _file_id, set_id = open_with_focal_set(service, project_dir, fasta)
    empty = tmp_path / "not-a-project"
    empty.mkdir()

    response = service.call("project.open", {"projectDir": str(empty)})
    assert response["ok"] is False
    assert response["error"]["code"] == "PROJECT_NOT_FOUND"

    focal_sets = service.ok("project.listFocalSets")["focalSets"]
    assert [entry["id"] for entry in focal_sets] == [set_id]
