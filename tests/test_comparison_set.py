"""
The optional explicit Comparison Set.

Semantics under test:

* EMPTY comparison: the run is exactly what it was before comparison sets
  existed — every non-focal sequence in the selected scope is the contrast.
* NON-EMPTY comparison: the analysis sees exactly the focal sequences plus the
  listed comparison sequences, and the scientific core is called unchanged.
* Focal and comparison never share an exact header; the service refuses it
  on its own, whatever the UI did.
* Comparison entries are saved with, loaded with, and deleted with their focal
  set. There is no separate comparison library.
"""

from __future__ import annotations

import sqlite3

import pytest

from molecular_diagnosis.project import service as service_module
from molecular_diagnosis.project.db import (
    MIGRATIONS_DIR,
    SCHEMA_VERSION,
    connect,
    open_project_db,
)
from molecular_diagnosis.project.service import ProjectError, ProjectService

# Focal F1/F2. The contrast specimens are chosen so the answer DEPENDS on which
# of them are compared against: C1/C2 differ from the focal at site 2, C3 and
# C4 only elsewhere.
RECORDS = [
    ("C1_other", "AGCCGGTT"),
    ("F1_target", "AACCGGTT"),
    ("C2_other", "AGCCGGTA"),
    ("F2_target", "AACCGGTT"),
    ("C3_other", "TACCGGTA"),
    ("C4_other", "AACCGGTA"),
]

FOCAL = ["F1_target", "F2_target"]


def write_fasta(path, records):
    path.write_text(
        "".join(f">{header}\n{sequence}\n" for header, sequence in records), encoding="utf-8"
    )
    return path


@pytest.fixture()
def project(tmp_path):
    service = ProjectService(tmp_path / "proj")
    service.open_project("Comparison")
    yield service
    service.close()


def link(project, tmp_path, name, records=RECORDS):
    row, _status = project.link_fasta(str(write_fasta(tmp_path / name, records)))
    return row.id


def save(project, *, comparison=None, focal=FOCAL, set_id=None, title="targets"):
    return project.save_focal_set(
        focal_set_id=set_id, title=title, headers=focal, comparison_headers=comparison
    )


def run(project, set_id, file_ids, *, single_file=None):
    return project.run_molecular_diagnosis(
        focal_set_id=set_id,
        fasta_file_ids=list(file_ids),
        single_file=len(file_ids) == 1 if single_file is None else single_file,
        options={"minCandidateSize": 1, "maxCandidateSize": 2},
    )


@pytest.fixture()
def captured(monkeypatch):
    """Record the sequences the scientific pipeline is actually handed."""
    from molecular_diagnosis import pipeline

    calls: list[dict] = []
    real = pipeline.run_pipeline_on_sequences

    def spy(**kwargs):
        calls.append(kwargs)
        return real(**kwargs)

    monkeypatch.setattr(pipeline, "run_pipeline_on_sequences", spy)
    return calls


# ---------------------------------------------------------------------------
# Blank comparison: unchanged behaviour
# ---------------------------------------------------------------------------


def test_blank_comparison_passes_the_full_scope_unchanged(project, tmp_path, captured):
    file_id = link(project, tmp_path, "a.fasta")
    set_id = save(project)["id"]

    result = run(project, set_id, [file_id])

    (call,) = captured
    assert list(call["sequences"]) == [header for header, _ in RECORDS]
    assert call["non_focal_headers"] == ["C1_other", "C2_other", "C3_other", "C4_other"]
    assert call["comparison_summary"] is None
    assert result["comparisonHeaders"] == []
    assert result["sequenceCount"] == result["scopeSequenceCount"] == len(RECORDS)


def test_blank_comparison_returns_the_same_dictionary_object():
    sequences = {"a": "AC", "b": "AG"}
    from molecular_diagnosis.focal import ExactHeaders

    out = ProjectService.apply_comparison_set(
        sequences, ExactHeaders(["a"]), [], single_file=True
    )
    assert out is sequences


def test_blank_comparison_report_has_no_comparison_line(project, tmp_path):
    file_id = link(project, tmp_path, "a.fasta")
    set_id = save(project)["id"]
    result = run(project, set_id, [file_id])
    report = open(result["outputs"]["reportTxt"], encoding="utf-8").read()
    assert "Comparison set:" not in report
    assert "Sequences analysed" not in report
    assert "Total sequences read: 6" + "\n" in report
    assert "Non-focal sequences: 4" in report


def test_a_save_without_comparison_matches_an_explicit_empty_one(project, tmp_path):
    file_id = link(project, tmp_path, "a.fasta")
    legacy = save(project, title="legacy")  # comparison_headers=None
    explicit_empty = save(project, title="empty", comparison=[])
    assert legacy["comparisonHeaders"] == explicit_empty["comparisonHeaders"] == []
    assert (
        run(project, legacy["id"], [file_id])["dmc"]
        == run(project, explicit_empty["id"], [file_id])["dmc"]
    )


# ---------------------------------------------------------------------------
# Explicit comparison: the analysis is limited to focal + comparison
# ---------------------------------------------------------------------------


def test_one_comparison_entry_limits_the_analysis(project, tmp_path, captured):
    file_id = link(project, tmp_path, "a.fasta")
    set_id = save(project, comparison=["C4_other"])["id"]

    result = run(project, set_id, [file_id])

    (call,) = captured
    # Scope order is kept, so the reference (first focal) is unchanged.
    assert list(call["sequences"]) == ["F1_target", "F2_target", "C4_other"]
    assert call["focal_headers"] == FOCAL
    assert call["non_focal_headers"] == ["C4_other"]
    assert result["comparisonHeaders"] == ["C4_other"]
    assert result["sequenceCount"] == 3
    assert result["scopeSequenceCount"] == len(RECORDS)


def test_multiple_comparison_entries_keep_scope_order(project, tmp_path, captured):
    file_id = link(project, tmp_path, "a.fasta")
    # Saved in a different order from the FASTA; analysis follows the FASTA.
    set_id = save(project, comparison=["C2_other", "C1_other"])["id"]

    run(project, set_id, [file_id])

    (call,) = captured
    assert list(call["sequences"]) == ["C1_other", "F1_target", "C2_other", "F2_target"]
    assert call["non_focal_headers"] == ["C1_other", "C2_other"]
    assert call["focal_headers"] == FOCAL


def test_explicit_comparison_equals_analysing_only_those_records(tmp_path):
    """
    The scientific core is not changed, only fed: an explicit comparison must
    give the same science as a FASTA that holds only focal + comparison.
    """
    chosen = {"F1_target", "F2_target", "C1_other", "C2_other"}

    narrowed = ProjectService(tmp_path / "narrowed")
    narrowed.open_project("narrowed")
    a = link(narrowed, tmp_path, "full.fasta")
    narrowed_set = save(narrowed, comparison=["C1_other", "C2_other"])["id"]
    with_comparison = run(narrowed, narrowed_set, [a])
    narrowed.close()

    plain = ProjectService(tmp_path / "plain")
    plain.open_project("plain")
    b = link(plain, tmp_path, "subset.fasta", [r for r in RECORDS if r[0] in chosen])
    plain_set = save(plain)["id"]
    without = run(plain, plain_set, [b])
    plain.close()

    assert with_comparison["dmc"] == without["dmc"]


def test_explicit_comparison_changes_the_answer_when_it_should(project, tmp_path):
    file_id = link(project, tmp_path, "a.fasta")
    everything = run(project, save(project, title="all")["id"], [file_id])
    near = run(
        project,
        save(project, title="near", comparison=["C1_other", "C2_other"])["id"],
        [file_id],
    )
    # Site 2 (index 1) separates focal from C1/C2 on its own; it cannot
    # separate them from C4, which shares the focal state there.
    assert [1] in near["dmc"]["combinationsByLength"].get("1", [])
    assert [1] not in everything["dmc"]["combinationsByLength"].get("1", [])


def test_explicit_comparison_is_described_in_the_report(project, tmp_path):
    file_id = link(project, tmp_path, "a.fasta")
    set_id = save(project, comparison=["C1_other", "C2_other"])["id"]
    report = open(run(project, set_id, [file_id])["outputs"]["reportTxt"], encoding="utf-8").read()
    assert "Total sequences read: 6" + "\n" in report
    assert "Sequences analysed (focal + comparison set): 4" + "\n" in report
    assert "Non-focal sequences: 2" in report
    assert (
        "Comparison set: explicit, 2 selected specimen(s); other non-focal sequences in the "
        "selected FASTA scope were not analysed" in report
    )


def test_comparison_entries_from_any_selected_file_are_accepted(project, tmp_path, captured):
    first = link(project, tmp_path, "a.fasta", RECORDS[:4])
    second = link(project, tmp_path, "b.fasta", RECORDS[4:])
    set_id = save(project, comparison=["C1_other", "C4_other"])["id"]

    run(project, set_id, [first, second], single_file=False)

    (call,) = captured
    assert call["non_focal_headers"] == ["C1_other", "C4_other"]


# ---------------------------------------------------------------------------
# Refusals
# ---------------------------------------------------------------------------


def test_overlap_between_focal_and_comparison_is_refused(project, tmp_path, captured):
    file_id = link(project, tmp_path, "a.fasta")
    set_id = save(project, comparison=["C1_other", "F2_target"])["id"]

    with pytest.raises(ProjectError) as caught:
        run(project, set_id, [file_id])

    assert caught.value.code == "COMPARISON_OVERLAPS_FOCAL"
    assert caught.value.detail == "F2_target"
    assert captured == []


def test_overlap_is_exact_header_equality_not_substring(project, tmp_path):
    records = RECORDS + [("F1_target_extra", "AGCCGGTT")]
    file_id = link(project, tmp_path, "a.fasta", records)
    set_id = save(project, comparison=["F1_target_extra"])["id"]
    result = run(project, set_id, [file_id])
    assert result["comparisonHeaders"] == ["F1_target_extra"]


def test_comparison_header_missing_from_the_selected_file_is_refused(project, tmp_path):
    first = link(project, tmp_path, "a.fasta", RECORDS[:4])
    link(project, tmp_path, "b.fasta", RECORDS[4:])
    # C3 exists in the project, but not in the file being analysed.
    set_id = save(project, comparison=["C1_other", "C3_other"])["id"]

    with pytest.raises(ProjectError) as caught:
        run(project, set_id, [first])

    assert caught.value.code == "COMPARISON_ENTRIES_NOT_IN_FILE"
    assert caught.value.detail == "C3_other"


def test_comparison_header_missing_from_every_selected_file_is_refused(project, tmp_path):
    first = link(project, tmp_path, "a.fasta", RECORDS[:4])
    second = link(project, tmp_path, "b.fasta", RECORDS[4:])
    set_id = save(project, comparison=["C1_other", "Nowhere"])["id"]

    with pytest.raises(ProjectError) as caught:
        run(project, set_id, [first, second], single_file=False)

    assert caught.value.code == "COMPARISON_ENTRIES_NOT_IN_SCOPE"


def test_existing_focal_validation_still_applies(project, tmp_path):
    file_id = link(project, tmp_path, "a.fasta")
    set_id = save(project, focal=["F1_target", "Missing_focal"], comparison=["C1_other"])["id"]
    with pytest.raises(ProjectError) as caught:
        run(project, set_id, [file_id])
    assert caught.value.code == "FOCAL_ENTRIES_NOT_IN_FILE"


def test_refusals_reach_the_rpc_boundary_with_their_own_codes(project, tmp_path):
    from molecular_diagnosis.service.projects import _as_service_error

    error = ProjectError("COMPARISON_OVERLAPS_FOCAL", "msg", detail="x")
    assert _as_service_error(error).code == "COMPARISON_OVERLAPS_FOCAL"


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------


def test_save_and_load_round_trip(project, tmp_path):
    link(project, tmp_path, "a.fasta")
    saved = save(project, comparison=[" C2_other ", "C1_other", "C2_other", ""])
    assert saved["comparisonHeaders"] == ["C2_other", "C1_other"]

    loaded = project.get_focal_set(saved["id"])
    assert [entry["header"] for entry in loaded["entries"]] == FOCAL
    assert loaded["comparisonHeaders"] == ["C2_other", "C1_other"]
    (listed,) = project.list_focal_sets()
    assert listed["comparisonHeaders"] == ["C2_other", "C1_other"]


def test_round_trip_survives_reopening_the_project(tmp_path):
    first = ProjectService(tmp_path / "proj")
    first.open_project("p")
    set_id = save(first, comparison=["C1_other"])["id"]
    first.close()

    again = ProjectService(tmp_path / "proj", create_if_missing=False)
    assert again.get_focal_set(set_id)["comparisonHeaders"] == ["C1_other"]
    again.close()


def test_editing_and_deleting_comparison_entries(project):
    set_id = save(project, comparison=["C1_other", "C2_other"])["id"]

    edited = save(project, set_id=set_id, comparison=["C2_other", "C3_other"])
    assert edited["comparisonHeaders"] == ["C2_other", "C3_other"]

    cleared = save(project, set_id=set_id, comparison=[])
    assert cleared["comparisonHeaders"] == []
    # Focal membership was not disturbed by either edit.
    assert [entry["header"] for entry in cleared["entries"]] == FOCAL


def test_a_save_without_comparison_leaves_the_stored_list_alone(project):
    set_id = save(project, comparison=["C1_other"])["id"]
    renamed = save(project, set_id=set_id, title="renamed")  # comparison=None
    assert renamed["comparisonHeaders"] == ["C1_other"]


def test_unchanged_comparison_entries_keep_their_rows(project):
    set_id = save(project, comparison=["C1_other", "C2_other"])["id"]
    before = project.connection.execute(
        "SELECT header, id FROM focal_set_comparison_entry ORDER BY header"
    ).fetchall()
    save(project, set_id=set_id, comparison=["C1_other", "C2_other"])
    after = project.connection.execute(
        "SELECT header, id FROM focal_set_comparison_entry ORDER BY header"
    ).fetchall()
    assert [tuple(row) for row in before] == [tuple(row) for row in after]


def test_overlap_may_be_saved_as_a_draft_state(project):
    """Saving never drops one side of an overlap; the run is what refuses it."""
    saved = save(project, comparison=["F1_target"])
    assert saved["comparisonHeaders"] == ["F1_target"]
    assert [entry["header"] for entry in saved["entries"]] == FOCAL


def test_focal_and_comparison_save_atomically(project, monkeypatch):
    set_id = save(project, comparison=["C1_other"])["id"]

    def boom(*_args, **_kwargs):
        raise sqlite3.OperationalError("disk on fire")

    monkeypatch.setattr(project.repository, "replace_comparison_entries", boom)
    with pytest.raises(sqlite3.OperationalError):
        save(project, set_id=set_id, focal=["F1_target"], comparison=["C2_other"], title="x")

    unchanged = project.get_focal_set(set_id)
    assert unchanged["title"] == "targets"
    assert [entry["header"] for entry in unchanged["entries"]] == FOCAL
    assert unchanged["comparisonHeaders"] == ["C1_other"]


def test_a_locked_set_refuses_comparison_edits(project):
    set_id = save(project, comparison=["C1_other"])["id"]
    project.set_focal_set_locked(set_id, True)
    with pytest.raises(ProjectError) as caught:
        save(project, set_id=set_id, comparison=[])
    assert caught.value.code == "FOCAL_SET_LOCKED"
    assert project.get_focal_set(set_id)["comparisonHeaders"] == ["C1_other"]


def test_deleting_a_focal_set_deletes_its_comparison_entries(project):
    set_id = save(project, comparison=["C1_other", "C2_other"])["id"]
    project.delete_focal_set(set_id)
    remaining = project.connection.execute(
        "SELECT COUNT(*) FROM focal_set_comparison_entry"
    ).fetchone()[0]
    assert remaining == 0


def test_save_rpc_accepts_and_omits_comparison_headers(project, monkeypatch):
    from molecular_diagnosis.service import projects as rpc

    monkeypatch.setitem(rpc._open, "current", project)
    created = rpc.save_focal_set(
        {"focalSetId": None, "title": "t", "headers": FOCAL, "comparisonHeaders": ["C1_other"]}
    )["focalSet"]
    assert created["comparisonHeaders"] == ["C1_other"]

    kept = rpc.save_focal_set({"focalSetId": created["id"], "title": "t2", "headers": FOCAL})
    assert kept["focalSet"]["comparisonHeaders"] == ["C1_other"]
    monkeypatch.delitem(rpc._open, "current")


# ---------------------------------------------------------------------------
# Migration: existing projects
# ---------------------------------------------------------------------------


def test_a_version_2_project_migrates_with_empty_comparison_sets(tmp_path):
    path = tmp_path / "old.sqlite"
    con = connect(path)
    for script in ("001_initial.sql", "002_fasta_file_locked.sql"):
        con.executescript((MIGRATIONS_DIR / script).read_text(encoding="utf-8"))
    assert con.execute("PRAGMA user_version").fetchone()[0] == 2
    con.execute(
        "INSERT INTO focal_set(id, title, locked, sort_order, created_at_ms, updated_at_ms)"
        " VALUES ('old-set', 'Old', 0, 0, 1, 1)"
    )
    con.execute(
        "INSERT INTO focal_set_entry(id, focal_set_id, header, sort_order,"
        " created_at_ms, updated_at_ms) VALUES ('e1', 'old-set', 'F1_target', 0, 1, 1)"
    )
    con.close()

    con, _caps, _fts = open_project_db(path)
    assert con.execute("PRAGMA user_version").fetchone()[0] == SCHEMA_VERSION == 3
    assert con.execute("SELECT COUNT(*) FROM focal_set_comparison_entry").fetchone()[0] == 0
    assert con.execute("SELECT header FROM focal_set_entry").fetchone()[0] == "F1_target"
    con.close()


def test_an_old_saved_set_loads_empty_and_runs_as_before(tmp_path, captured):
    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    con = connect(project_dir / service_module.PROJECT_DB_NAME)
    for script in ("001_initial.sql", "002_fasta_file_locked.sql"):
        con.executescript((MIGRATIONS_DIR / script).read_text(encoding="utf-8"))
    con.close()

    project = ProjectService(project_dir, create_if_missing=False)
    project.open_project("old")
    file_id = link(project, tmp_path, "a.fasta")
    # A set written the old way: focal entries only.
    set_id = project.save_focal_set(focal_set_id=None, title="old", headers=FOCAL)["id"]
    assert project.get_focal_set(set_id)["comparisonHeaders"] == []

    run(project, set_id, [file_id])
    (call,) = captured
    assert len(call["sequences"]) == len(RECORDS)
    project.close()
