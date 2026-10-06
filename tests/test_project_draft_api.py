"""
The three operations the project-backed workspace edits against.

`saveFocalSet` is the ONLY thing here that writes. `headerPresence` and
`matchFocalHeaders` exist so a renderer can edit an unsaved working copy —
colouring it, and running `-` over it — without the act of typing committing
anything to the database or touching the filesystem.
"""

from __future__ import annotations

import pytest

from molecular_diagnosis.project.service import ProjectError, ProjectService

RECORDS = [
    ("focal_1|AU|Target", "AACC"),
    ("focal_2|GB|Target", "AACC"),
    ("other_1|FR|Contrast", "AGCT"),
    ("other_2|DE|Contrast", "CGTT"),
]


def write_fasta(path, records):
    path.write_text(
        "".join(f">{header}\n{sequence}\n" for header, sequence in records), encoding="utf-8"
    )
    return path


@pytest.fixture()
def project(tmp_path):
    service = ProjectService(tmp_path / "proj")
    service.open_project("Test")
    yield service
    service.close()


@pytest.fixture()
def fasta(tmp_path):
    return write_fasta(tmp_path / "a.fasta", RECORDS)


@pytest.fixture()
def linked(project, fasta):
    row, _status = project.link_fasta(str(fasta))
    return row


# ---------------------------------------------------------------------------
# saveFocalSet — the explicit Save boundary
# ---------------------------------------------------------------------------


def test_saving_a_new_focal_set_creates_it_and_its_entries_together(project, linked):
    saved = project.save_focal_set(
        focal_set_id=None,
        title="Targets",
        headers=["focal_1|AU|Target", "focal_2|GB|Target"],
    )

    assert saved["title"] == "Targets"
    assert saved["locked"] is False
    assert [entry["header"] for entry in saved["entries"]] == [
        "focal_1|AU|Target",
        "focal_2|GB|Target",
    ]
    # One set — not a titled empty one left behind by a half-done create.
    assert len(project.list_focal_sets()) == 1


def test_saving_an_existing_set_updates_title_and_membership_in_place(project, linked):
    first = project.save_focal_set(
        focal_set_id=None, title="Draft", headers=["focal_1|AU|Target", "focal_2|GB|Target"]
    )
    kept_id = next(e["id"] for e in first["entries"] if e["header"] == "focal_1|AU|Target")
    locations_before = project.repository.locations_for_entry(kept_id)
    assert locations_before

    second = project.save_focal_set(
        focal_set_id=first["id"],
        title="Renamed",
        headers=["focal_1|AU|Target", "other_1|FR|Contrast"],
    )

    assert second["id"] == first["id"]
    assert second["title"] == "Renamed"
    assert [entry["header"] for entry in second["entries"]] == [
        "focal_1|AU|Target",
        "other_1|FR|Contrast",
    ]
    # The surviving header keeps its identity, and with it its location cache.
    assert next(e["id"] for e in second["entries"] if e["header"] == "focal_1|AU|Target") == kept_id
    assert project.repository.locations_for_entry(kept_id) == locations_before
    assert len(project.list_focal_sets()) == 1


def test_saving_trims_and_deduplicates_preserving_first_occurrence(project, linked):
    saved = project.save_focal_set(
        focal_set_id=None,
        title="  Spaced  ",
        headers=[" focal_2|GB|Target ", "focal_1|AU|Target", "focal_2|GB|Target", "  ", ""],
    )

    assert saved["title"] == "Spaced"
    assert [entry["header"] for entry in saved["entries"]] == [
        "focal_2|GB|Target",
        "focal_1|AU|Target",
    ]


def test_saving_keeps_headers_that_match_no_fasta(project, linked):
    saved = project.save_focal_set(
        focal_set_id=None, title="Mixed", headers=["focal_1|AU|Target", "typed_by_hand"]
    )

    assert [entry["header"] for entry in saved["entries"]] == [
        "focal_1|AU|Target",
        "typed_by_hand",
    ]


def test_saving_a_locked_set_is_refused_and_changes_nothing(project, linked):
    saved = project.save_focal_set(
        focal_set_id=None, title="Locked", headers=["focal_1|AU|Target"]
    )
    project.set_focal_set_locked(saved["id"], True)

    with pytest.raises(ProjectError) as caught:
        project.save_focal_set(
            focal_set_id=saved["id"], title="Changed", headers=["other_1|FR|Contrast"]
        )

    assert caught.value.code == "FOCAL_SET_LOCKED"
    current = project.get_focal_set(saved["id"])
    assert current["title"] == "Locked"
    assert [entry["header"] for entry in current["entries"]] == ["focal_1|AU|Target"]


def test_saving_an_unknown_id_is_a_coded_refusal(project, linked):
    with pytest.raises(ProjectError) as caught:
        project.save_focal_set(focal_set_id="nope", title="T", headers=["h"])
    assert caught.value.code == "UNKNOWN_FOCAL_SET"


def test_saving_requires_a_title(project, linked):
    with pytest.raises(ProjectError) as caught:
        project.save_focal_set(focal_set_id=None, title="   ", headers=["focal_1|AU|Target"])
    assert caught.value.code == "INVALID_PARAMETER"


def test_an_empty_membership_is_savable(project, linked):
    """A set the user emptied is a state, not an error; Run gates on it instead."""
    saved = project.save_focal_set(focal_set_id=None, title="Empty", headers=[])
    assert saved["entries"] == []


# ---------------------------------------------------------------------------
# headerPresence — colours for an unsaved draft
# ---------------------------------------------------------------------------


def test_header_presence_answers_for_unsaved_headers(project, fasta, linked, tmp_path):
    other = write_fasta(tmp_path / "b.fasta", [("b1|X|Target", "AACC"), ("b2|X", "TTTT")])
    other_row, _status = project.link_fasta(str(other))

    entries = project.header_presence(
        ["focal_1|AU|Target", "b1|X|Target", "typed_by_hand"],
        selected_file_id=linked.id,
    )
    states = {entry["header"]: entry["state"] for entry in entries}

    assert states["focal_1|AU|Target"] == "present_current"
    # In the project, but not in the selected file.
    assert states["b1|X|Target"] == "present_other"
    assert states["typed_by_hand"] == "missing"
    assert entries[0]["occurrences"] == {linked.id: 1}
    assert entries[1]["occurrences"] == {other_row.id: 1}


def test_header_presence_writes_nothing(project, linked):
    project.header_presence(["focal_1|AU|Target"], selected_file_id=linked.id)

    # A draft header belongs to no focal set, so nothing may be persisted for it.
    rows = project.connection.execute("SELECT COUNT(*) FROM focal_entry_location").fetchone()[0]
    assert rows == 0
    assert project.repository.list_focal_sets() == []


def test_header_presence_never_hashes_or_reindexes(project, fasta, linked, monkeypatch):
    """Typing must not cost filesystem work, so the expensive paths are barred."""
    from molecular_diagnosis.project import service as service_module

    def explode(*_args, **_kwargs):
        raise AssertionError("presence must not read or hash the FASTA")

    monkeypatch.setattr(service_module, "scan_fasta", explode)
    monkeypatch.setattr(service_module, "hash_file", explode)
    monkeypatch.setattr(service_module.ProjectService, "reindex", explode)

    entries = project.header_presence(["focal_1|AU|Target"], selected_file_id=linked.id)
    assert entries[0]["state"] == "present_current"


def test_header_presence_for_all_files_has_no_orange(project, fasta, linked, tmp_path):
    other = write_fasta(tmp_path / "b.fasta", [("b1|X|Target", "AACC")])
    project.link_fasta(str(other))

    entries = project.header_presence(
        ["focal_1|AU|Target", "b1|X|Target", "nowhere"], selected_file_id=None
    )
    states = {entry["header"]: entry["state"] for entry in entries}

    assert states["focal_1|AU|Target"] == "present_current"
    assert states["b1|X|Target"] == "present_current"
    assert states["nowhere"] == "missing"
    assert "present_other" not in set(states.values())


def test_header_presence_reports_unknown_when_a_source_is_unavailable(
    project, fasta, linked, tmp_path
):
    gone = write_fasta(tmp_path / "b.fasta", [("b1|X|Target", "AACC")])
    gone_row, _status = project.link_fasta(str(gone))
    gone.unlink()

    # Not found anywhere, and one file could not be consulted: unproven, not absent.
    assert project.header_presence(["nowhere_at_all"], selected_file_id=None)[0]["state"] == (
        "unknown"
    )

    # Selecting the vanished file makes every answer about it unknown.
    assert project.header_presence(["focal_1|AU|Target"], selected_file_id=gone_row.id)[0][
        "state"
    ] == "unknown"

    # A header proven present somewhere readable stays green: an unavailable
    # file cannot make it less present.
    assert project.header_presence(["focal_1|AU|Target"], selected_file_id=None)[0]["state"] == (
        "present_current"
    )


def test_header_presence_is_unknown_while_a_source_is_stale(project, fasta, linked):
    """A file edited since indexing cannot answer from its stale index."""
    write_fasta(fasta, [("completely|different", "AACC")])

    entries = project.header_presence(["focal_1|AU|Target"], selected_file_id=linked.id)
    assert entries[0]["state"] == "unknown"


def test_header_presence_can_be_scoped(project, fasta, linked, tmp_path):
    other = write_fasta(tmp_path / "b.fasta", [("b1|X|Target", "AACC")])
    other_row, _status = project.link_fasta(str(other))

    scoped = project.header_presence(
        ["b1|X|Target"], selected_file_id=None, fasta_file_ids=[linked.id]
    )
    assert scoped[0]["state"] == "missing"

    wide = project.header_presence(
        ["b1|X|Target"], selected_file_id=None, fasta_file_ids=[linked.id, other_row.id]
    )
    assert wide[0]["state"] == "present_current"


def test_header_presence_deduplicates_and_ignores_blanks(project, linked):
    entries = project.header_presence(
        ["focal_1|AU|Target", " ", "focal_1|AU|Target", ""], selected_file_id=linked.id
    )
    assert [entry["header"] for entry in entries] == ["focal_1|AU|Target"]
    assert project.header_presence([], selected_file_id=linked.id) == []


def test_header_presence_handles_more_headers_than_sqlite_binds(project, linked):
    """The batched lookup chunks, so a large draft stays one query per chunk."""
    headers = [f"synthetic_{index}" for index in range(1500)] + ["focal_1|AU|Target"]

    entries = project.header_presence(headers, selected_file_id=linked.id)

    assert len(entries) == 1501
    states = {entry["header"]: entry["state"] for entry in entries}
    assert states["focal_1|AU|Target"] == "present_current"
    assert states["synthetic_0"] == "missing"


def test_header_presence_counts_duplicate_records(project, tmp_path):
    path = write_fasta(
        tmp_path / "dupes.fasta", [("dupe|X", "AACC"), ("dupe|X", "TTTT"), ("solo|X", "GGGG")]
    )
    row, _status = project.link_fasta(str(path))

    entries = project.header_presence(["dupe|X"], selected_file_id=row.id)
    assert entries[0]["occurrences"] == {row.id: 2}


# ---------------------------------------------------------------------------
# matchFocalHeaders — `-` over a working copy
# ---------------------------------------------------------------------------


def test_match_focal_headers_uses_casefolded_substring_semantics(project):
    headers = ["focal_1|AU|Target", "focal_2|GB|Target", "other_1|FR|Contrast"]

    assert project.match_focal_headers("target", headers)["matched"] == [
        "focal_1|AU|Target",
        "focal_2|GB|Target",
    ]
    assert project.match_focal_headers("FOCAL_1", headers)["matched"] == ["focal_1|AU|Target"]
    assert project.match_focal_headers("zzz", headers)["matched"] == []


def test_match_focal_headers_folds_case_the_way_python_does(project):
    # "ß" folds to "ss" under casefold(); JavaScript's toLowerCase() leaves it
    # alone. That divergence is the whole reason this lives in Python.
    assert project.match_focal_headers("STRASSE", ["Straße_01"])["matched"] == ["Straße_01"]


def test_match_focal_headers_refuses_an_empty_query(project):
    for query in ("", "   "):
        with pytest.raises(ProjectError) as caught:
            project.match_focal_headers(query, ["anything"])
        assert caught.value.code == "FOCAL_QUERY_EMPTY"


def test_match_focal_headers_does_not_search_the_fasta(project, linked):
    """`-` operates on the working copy only; the file is irrelevant to it."""
    # focal_2 exists in the FASTA but is not in the supplied working copy.
    assert project.match_focal_headers("focal_2", ["focal_1|AU|Target"])["matched"] == []


def test_match_focal_headers_agrees_with_the_persistent_remove(project, linked):
    """The working-copy helper and the saved-set mutation must not diverge."""
    saved = project.save_focal_set(
        focal_set_id=None, title="S", headers=["focal_1|AU|Target", "focal_2|GB|Target"]
    )
    stored = [entry["header"] for entry in saved["entries"]]

    matched = project.match_focal_headers("focal_2", stored)["matched"]
    removed = project.remove_focal_entries_by_query(saved["id"], "focal_2")["removed"]

    assert matched == removed


# ---------------------------------------------------------------------------
# resolveFocalAddQuery — what `+` would add, uncapped
# ---------------------------------------------------------------------------


def test_resolve_add_query_returns_every_match_with_no_cap(project, tmp_path):
    """
    The regression this exists for.

    The renderer used to expand `+` through the capped `search_headers`
    preview, so a query matching more headers than the cap silently produced a
    smaller focal set than the user asked for.
    """
    records = [(f"BIG{index:04d}|AU|Target", "AACC") for index in range(240)]
    records.append(("outgroup|XX|Contrast", "AGCT"))
    path = write_fasta(tmp_path / "many.fasta", records)
    project.link_fasta(str(path))

    result = project.resolve_focal_add_query("Target")

    assert len(result["headers"]) == 240
    assert result["headers"][0] == "BIG0000|AU|Target"
    assert result["headers"][-1] == "BIG0239|AU|Target"
    # The capped preview is a different operation and still caps.
    assert len(project.search_headers("Target", limit=200)["hits"]) == 200


def test_resolve_add_query_writes_nothing(project, linked):
    project.resolve_focal_add_query("Target")

    assert project.repository.list_focal_sets() == []
    rows = project.connection.execute("SELECT COUNT(*) FROM focal_entry_location").fetchone()[0]
    assert rows == 0


def test_resolve_add_query_matches_the_mutating_plus_exactly(project, fasta, linked, tmp_path):
    """Both go through `_expand_add_query`, so they cannot drift apart."""
    other = write_fasta(tmp_path / "b.fasta", [("b1|ES|Target", "AACC"), ("b2|X", "TTTT")])
    other_row, _status = project.link_fasta(str(other))

    resolved = project.resolve_focal_add_query("target")["headers"]

    focal = project.create_focal_set("S")
    added = project.add_focal_entries_from_query(focal["id"], "target")["matched"]

    assert resolved == added
    assert resolved == ["focal_1|AU|Target", "focal_2|GB|Target", "b1|ES|Target"]
    assert other_row.id


def test_resolve_add_query_honours_the_requested_scope(project, fasta, linked, tmp_path):
    other = write_fasta(tmp_path / "b.fasta", [("b1|ES|Target", "AACC")])
    other_row, _status = project.link_fasta(str(other))

    only_a = project.resolve_focal_add_query("Target", fasta_file_ids=[linked.id])["headers"]
    assert only_a == ["focal_1|AU|Target", "focal_2|GB|Target"]

    only_b = project.resolve_focal_add_query("Target", fasta_file_ids=[other_row.id])["headers"]
    assert only_b == ["b1|ES|Target"]


def test_resolve_add_query_refuses_an_unsearchable_scope(project, fasta, linked, tmp_path):
    gone = write_fasta(tmp_path / "b.fasta", [("b1|ES|Target", "AACC")])
    project.link_fasta(str(gone))
    gone.unlink()

    with pytest.raises(ProjectError) as caught:
        project.resolve_focal_add_query("Target")

    assert caught.value.code == "SEARCH_SCOPE_UNAVAILABLE"
    assert "b.fasta" in (caught.value.detail or "")


def test_resolve_add_query_refuses_an_empty_query(project, linked):
    for query in ("", "   "):
        with pytest.raises(ProjectError) as caught:
            project.resolve_focal_add_query(query)
        assert caught.value.code == "FOCAL_QUERY_EMPTY"
