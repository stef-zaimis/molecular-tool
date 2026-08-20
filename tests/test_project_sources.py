"""
Source verification, reindexing, relinking, focal locations, search, the
alignment cache and multi-file scopes.

These are the behavioural tests for the invariant that the linked FASTA is
always authoritative and everything in SQLite is derived cache.
"""

from __future__ import annotations

import os
import shutil

import pytest

from molecular_diagnosis.project.indexing import SourceChangedDuringRead, read_header_at, scan_fasta
from molecular_diagnosis.project.locations import PresenceState
from molecular_diagnosis.project.search import MIN_TRIGRAM_LENGTH, quote_fts_query, search_headers
from molecular_diagnosis.project.service import ProjectError, ProjectService
from molecular_diagnosis.project.sources import SourceState

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


def set_mtime(path, mtime_ns: int) -> None:
    os.utime(path, ns=(mtime_ns, mtime_ns))


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
# Project open / missing files
# ---------------------------------------------------------------------------


def test_project_open_reports_a_missing_linked_file(project, fasta, linked, tmp_path):
    fasta.unlink()
    payload = project.open_project("Test")
    source = payload["sources"][0]
    assert source["state"] == SourceState.MISSING.value
    assert source["available"] is False
    assert source["indexUsable"] is False
    # The registration itself survives so the user can relocate it.
    assert project.repository.get_fasta_file(linked.id) is not None


def test_a_missing_file_keeps_its_rows_but_blocks_stale_use(project, fasta, linked):
    assert project.repository.record_count(linked.id) == 4
    fasta.unlink()

    status = project.verifier.cheap_status(project.repository.get_fasta_file(linked.id))
    assert status.state is SourceState.MISSING
    assert status.index_usable is False
    # Rows are still there — they are cache, not garbage — but unusable.
    assert project.repository.record_count(linked.id) == 4
    with pytest.raises(ProjectError, match="not available"):
        project.load_alignment(linked.id)


def test_missing_file_locations_are_preserved(project, fasta, linked):
    focal_set = project.repository.create_focal_set("S")
    project.add_focal_entries_from_query(focal_set.id, "Target")
    assert len(project.repository.locations_for_set(focal_set.id)) == 2

    fasta.unlink()
    presence = project.focal_presence(focal_set.id, selected_file_id=linked.id)

    assert {entry.state for entry in presence} == {PresenceState.UNKNOWN}
    # Persisted locations are not discarded merely because the file vanished.
    assert len(project.repository.locations_for_set(focal_set.id)) == 2


# ---------------------------------------------------------------------------
# Cheap vs strong verification
# ---------------------------------------------------------------------------


def test_unchanged_file_verifies_without_reindexing(project, fasta, linked):
    revision = linked.index_revision
    project.verifier.forget_all()

    status = project.verifier.ensure_current(project.repository.get_fasta_file(linked.id))
    assert status.state is SourceState.CURRENT
    assert project.repository.get_fasta_file(linked.id).index_revision == revision


def test_first_use_in_a_session_strongly_verifies_then_caches(project, fasta, linked, monkeypatch):
    project.verifier.forget_all()
    calls = {"n": 0}
    import molecular_diagnosis.project.sources as sources

    real = sources.hash_file

    def counted(path):
        calls["n"] += 1
        return real(path)

    monkeypatch.setattr(sources, "hash_file", counted)

    row = project.repository.get_fasta_file(linked.id)
    assert project.verifier.ensure_current(row).state is SourceState.CURRENT
    assert calls["n"] == 1

    # Repeated checks — the focal-box keystroke case — must not re-hash.
    for _ in range(5):
        assert project.verifier.ensure_current(row).state is SourceState.CURRENT
    assert calls["n"] == 1


def test_metadata_changed_but_bytes_identical_keeps_the_index(project, fasta, linked):
    original = project.repository.get_fasta_file(linked.id)
    project.verifier.forget_all()

    # Same content, new mtime: the cheap detector fires, the strong one clears it.
    set_mtime(fasta, 1_500_000_000_000_000_000)

    status = project.verifier.ensure_current(project.repository.get_fasta_file(linked.id))
    updated = project.repository.get_fasta_file(linked.id)

    assert status.state is SourceState.CURRENT
    assert updated.index_revision == original.index_revision  # no reindex
    assert updated.indexed_sha256 == original.indexed_sha256
    assert updated.indexed_mtime_ns == os.stat(fasta).st_mtime_ns  # cheap metadata refreshed


def test_real_content_change_is_detected_and_reindexes(project, fasta, linked):
    original = project.repository.get_fasta_file(linked.id)
    write_fasta(fasta, RECORDS + [("extra|NL|Contrast", "TTTT")])
    project.verifier.forget_all()

    status = project.verifier.ensure_current(project.repository.get_fasta_file(linked.id))
    assert status.state is SourceState.STALE

    project.reindex(linked.id)
    updated = project.repository.get_fasta_file(linked.id)
    assert updated.index_revision == original.index_revision + 1
    assert updated.sequence_count == 5
    assert project.repository.record_count(linked.id) == 5


def test_content_changed_while_size_and_mtime_are_restored(project, fasta, linked):
    """
    The case the cheap detector cannot see: an edit of identical length with
    the original mtime put back. Only the hash catches it.
    """
    original = project.repository.get_fasta_file(linked.id)
    stat_before = os.stat(fasta)

    edited = list(RECORDS)
    edited[3] = ("other_2|DE|Contrast", "CGTA")  # same length, different bytes
    write_fasta(fasta, edited)
    os.truncate(fasta, stat_before.st_size)
    set_mtime(fasta, stat_before.st_mtime_ns)

    stat_after = os.stat(fasta)
    assert (stat_after.st_size, stat_after.st_mtime_ns) == (
        stat_before.st_size, stat_before.st_mtime_ns
    ), "the fixture failed to restore the cheap metadata"

    project.verifier.forget_all()
    status = project.verifier.ensure_current(project.repository.get_fasta_file(linked.id))

    assert status.state is SourceState.STALE
    assert project.repository.get_fasta_file(linked.id).index_revision == original.index_revision


def test_stat_change_invalidates_a_session_verification(project, fasta, linked):
    row = project.repository.get_fasta_file(linked.id)
    assert project.verifier.ensure_current(row).state is SourceState.CURRENT
    assert project.verifier.session_verified(linked.id)

    write_fasta(fasta, RECORDS[:3])
    status = project.verifier.cheap_status(project.repository.get_fasta_file(linked.id))

    assert status.state is SourceState.STALE
    assert not project.verifier.session_verified(linked.id)


# ---------------------------------------------------------------------------
# Reading race
# ---------------------------------------------------------------------------


def test_a_file_changing_during_the_read_is_rejected(tmp_path, monkeypatch):
    path = write_fasta(tmp_path / "racy.fasta", RECORDS)
    import molecular_diagnosis.project.indexing as indexing

    real_fstat = os.fstat
    state = {"calls": 0}

    def shifting_fstat(fd):
        stat = real_fstat(fd)
        state["calls"] += 1
        if state["calls"] > 1:
            # Pretend a writer touched the file between the two stats.
            return os.stat_result(
                (stat.st_mode, stat.st_ino, stat.st_dev, stat.st_nlink, stat.st_uid,
                 stat.st_gid, stat.st_size + 7, stat.st_atime, stat.st_mtime, stat.st_ctime)
            )
        return stat

    monkeypatch.setattr(indexing.os, "fstat", shifting_fstat)
    with pytest.raises(SourceChangedDuringRead):
        scan_fasta(path)


def test_a_racing_read_does_not_commit_an_index(project, fasta, linked, monkeypatch):
    original = project.repository.get_fasta_file(linked.id)
    import molecular_diagnosis.project.service as service_module

    def always_racing(_path):
        raise SourceChangedDuringRead("simulated")

    monkeypatch.setattr(service_module, "scan_fasta", always_racing)
    with pytest.raises(ProjectError) as caught:
        project.reindex(linked.id)

    assert caught.value.code == "SOURCE_CHANGED_DURING_READ"
    unchanged = project.repository.get_fasta_file(linked.id)
    assert unchanged.index_revision == original.index_revision
    assert unchanged.indexed_sha256 == original.indexed_sha256


def test_an_invalid_new_version_does_not_replace_the_index(project, fasta, linked):
    original = project.repository.get_fasta_file(linked.id)
    fasta.write_text(">ragged_1\nACGT\n>ragged_2\nACG\n", encoding="utf-8")

    with pytest.raises(ProjectError) as caught:
        project.reindex(linked.id)
    assert caught.value.code == "FASTA_NOT_ALIGNED"

    kept = project.repository.get_fasta_file(linked.id)
    assert kept.index_revision == original.index_revision
    assert kept.indexed_sha256 == original.indexed_sha256
    assert project.repository.record_count(linked.id) == 4
    # ...but the source is not usable, because it no longer matches the index.
    project.verifier.forget_all()
    assert project.verifier.ensure_current(kept).index_usable is False


# ---------------------------------------------------------------------------
# Byte offsets
# ---------------------------------------------------------------------------


def test_stored_offsets_resolve_to_their_headers(project, fasta, linked):
    for ordinal, header in project.repository.headers_for_file(linked.id):
        rows = project.repository.locate_header(linked.id, header)
        row = next(r for r in rows if r["ordinal"] == ordinal)
        assert read_header_at(fasta, row["record_start_byte"]) == header


def test_a_stale_offset_is_detected_rather_than_shown(project, fasta, linked):
    rows = project.repository.locate_header(linked.id, "other_1|FR|Contrast")
    offset = rows[0]["record_start_byte"]

    # Prepend a record: every later offset now points at the wrong place.
    write_fasta(fasta, [("inserted|XX|Contrast", "GGGG"), *RECORDS])
    assert read_header_at(fasta, offset) != "other_1|FR|Contrast"


# ---------------------------------------------------------------------------
# Focal locations
# ---------------------------------------------------------------------------


def test_focal_location_survives_an_edit_elsewhere_in_the_file(project, fasta, linked):
    focal_set = project.repository.create_focal_set("S")
    project.add_focal_entries_from_query(focal_set.id, "focal_1")
    before = project.repository.locations_for_set(focal_set.id)
    assert len(before) == 1

    # Change a LATER record only: focal_1 keeps its ordinal and byte offset.
    edited = list(RECORDS)
    edited[3] = ("other_2|DE|Contrast", "CGTA")
    write_fasta(fasta, edited)

    presence = project.focal_presence(focal_set.id, selected_file_id=linked.id)
    assert presence[0].state is PresenceState.PRESENT_CURRENT
    after = project.repository.locations_for_set(focal_set.id)
    assert after[0].ordinal == before[0].ordinal
    assert after[0].record_start_byte == before[0].record_start_byte
    assert after[0].last_verified_at_ms is not None


def test_header_moving_to_a_different_ordinal_is_repaired(project, fasta, linked):
    focal_set = project.repository.create_focal_set("S")
    project.add_focal_entries_from_query(focal_set.id, "focal_1")
    assert project.repository.locations_for_set(focal_set.id)[0].ordinal == 0

    write_fasta(fasta, [("inserted|XX|Contrast", "GGGG"), *RECORDS])
    presence = project.focal_presence(focal_set.id, selected_file_id=linked.id)

    assert presence[0].state is PresenceState.PRESENT_CURRENT
    assert project.repository.locations_for_set(focal_set.id)[0].ordinal == 1


def test_focal_entry_disappearing_after_an_edit(project, fasta, linked):
    focal_set = project.repository.create_focal_set("S")
    project.add_focal_entries_from_query(focal_set.id, "focal_1")

    write_fasta(fasta, RECORDS[1:])
    presence = project.focal_presence(focal_set.id, selected_file_id=linked.id)

    assert presence[0].state is PresenceState.MISSING
    assert project.repository.locations_for_set(focal_set.id) == []


def test_focal_entry_appearing_after_an_edit(project, fasta, linked):
    focal_set = project.repository.create_focal_set("S")
    project.repository.begin()
    project.repository.add_focal_entries(focal_set.id, ["late|ZZ|Target"])
    project.repository.commit()

    assert project.focal_presence(focal_set.id, selected_file_id=linked.id)[0].state is (
        PresenceState.MISSING
    )

    write_fasta(fasta, [*RECORDS, ("late|ZZ|Target", "AACC")])
    presence = project.focal_presence(focal_set.id, selected_file_id=linked.id)

    assert presence[0].state is PresenceState.PRESENT_CURRENT
    assert project.repository.locations_for_set(focal_set.id)[0].ordinal == 4


def test_sequence_changing_while_the_header_stays_keeps_the_location(project, fasta, linked):
    focal_set = project.repository.create_focal_set("S")
    project.add_focal_entries_from_query(focal_set.id, "focal_1")

    edited = list(RECORDS)
    edited[0] = ("focal_1|AU|Target", "TTTT")
    write_fasta(fasta, edited)

    presence = project.focal_presence(focal_set.id, selected_file_id=linked.id)
    assert presence[0].state is PresenceState.PRESENT_CURRENT
    # But the alignment cache must not hand back the old sequence.
    alignment = project.load_alignment(linked.id)
    assert alignment.sequences["focal_1|AU|Target"] == "TTTT"


def test_duplicate_headers_produce_multiple_locations(project, tmp_path):
    path = write_fasta(
        tmp_path / "dup.fasta",
        [("dup|H", "AAAA"), ("dup|H", "CCCC"), ("solo|H", "GGGG")],
    )
    row, _status = project.link_fasta(str(path))
    focal_set = project.repository.create_focal_set("S")
    project.add_focal_entries_from_query(focal_set.id, "dup")

    locations = project.repository.locations_for_set(focal_set.id)
    assert sorted(location.ordinal for location in locations) == [0, 1]
    presence = project.focal_presence(focal_set.id, selected_file_id=row.id)
    assert presence[0].occurrences[row.id] == 2


def test_present_in_another_file_is_orange_only_for_a_selected_file(project, tmp_path, fasta, linked):
    other = write_fasta(tmp_path / "b.fasta", [("only_here|QQ|Target", "AACC"), ("filler|X", "GGGG")])
    other_row, _ = project.link_fasta(str(other))

    focal_set = project.repository.create_focal_set("S")
    project.repository.begin()
    project.repository.add_focal_entries(focal_set.id, ["only_here|QQ|Target"])
    project.repository.commit()

    selected = project.focal_presence(focal_set.id, selected_file_id=linked.id)
    assert selected[0].state is PresenceState.PRESENT_OTHER

    all_files = project.focal_presence(focal_set.id, selected_file_id=None)
    assert all_files[0].state is PresenceState.PRESENT_CURRENT  # green, never orange


# ---------------------------------------------------------------------------
# Relinking
# ---------------------------------------------------------------------------


def test_relink_to_a_byte_identical_move_avoids_reindexing(project, fasta, linked, tmp_path):
    focal_set = project.repository.create_focal_set("S")
    project.add_focal_entries_from_query(focal_set.id, "Target")
    before_locations = project.repository.locations_for_set(focal_set.id)
    before = project.repository.get_fasta_file(linked.id)

    moved = tmp_path / "moved" / "a.fasta"
    moved.parent.mkdir()
    shutil.move(str(fasta), str(moved))

    result = project.relink_fasta(linked.id, str(moved))
    after = project.repository.get_fasta_file(linked.id)

    assert result["identical"] is True
    assert result["reindexed"] is False
    assert after.index_revision == before.index_revision
    assert after.source_path == str(moved)
    assert project.repository.record_count(linked.id) == 4
    assert project.repository.locations_for_set(focal_set.id) == before_locations


def test_relink_to_changed_valid_fasta_reindexes(project, fasta, linked, tmp_path):
    before = project.repository.get_fasta_file(linked.id)
    replacement = write_fasta(
        tmp_path / "replacement.fasta", [*RECORDS, ("extra|NL|Contrast", "TTTT")]
    )

    result = project.relink_fasta(linked.id, str(replacement))
    after = project.repository.get_fasta_file(linked.id)

    assert result["identical"] is False
    assert result["reindexed"] is True
    assert after.index_revision == before.index_revision + 1
    assert after.sequence_count == 5


def test_relink_to_an_invalid_file_is_rejected_cleanly(project, fasta, linked, tmp_path):
    before = project.repository.get_fasta_file(linked.id)
    bad = tmp_path / "bad.fasta"
    bad.write_text(">a\nACGT\n>b\nACG\n", encoding="utf-8")

    with pytest.raises(ProjectError) as caught:
        project.relink_fasta(linked.id, str(bad))
    assert caught.value.code == "FASTA_NOT_ALIGNED"

    after = project.repository.get_fasta_file(linked.id)
    assert after.source_path == before.source_path
    assert after.index_revision == before.index_revision
    assert project.repository.record_count(linked.id) == 4


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------


def test_duplicate_headers_are_preserved_by_the_index(project, tmp_path):
    path = write_fasta(
        tmp_path / "dup.fasta", [("dup|H", "AAAA"), ("dup|H", "CCCC"), ("solo|H", "GGGG")]
    )
    row, _status = project.link_fasta(str(path))
    assert project.repository.record_count(row.id) == 3
    assert project.repository.get_fasta_file(row.id).sequence_count == 2
    assert project.repository.get_fasta_file(row.id).duplicate_header_count == 1


def test_fts_is_used_when_supported(project, linked):
    if not project.repository.fts_enabled:
        pytest.skip("this SQLite has no FTS5 trigram support")
    hits = search_headers(project.repository, "target")
    assert [hit.header for hit in hits] == ["focal_1|AU|Target", "focal_2|GB|Target"]


def test_fallback_matches_the_accelerated_path(project, fasta, linked, tmp_path):
    extra = write_fasta(
        tmp_path / "c.fasta",
        [("MiXeD|Case|Target", "AACC"), ("weird;semi|X", "GGGG"), ("a-b|OR|NEAR", "TTTT")],
    )
    project.link_fasta(str(extra))

    queries = [
        "target", "TARGET", "Target", "contrast", "|AU|", "semi", "a-b", "OR",
        "focal_1|AU|Target", "zzz", ";", "MiXeD",
    ]
    for query in queries:
        accelerated = search_headers(project.repository, query)
        fallback = search_headers(project.repository, query, force_fallback=True)
        assert [hit.header for hit in accelerated] == [hit.header for hit in fallback], query


def test_short_queries_take_the_fallback_and_still_match(project, linked):
    for query in ["A", "AU", "|"]:
        assert len(query) < MIN_TRIGRAM_LENGTH
        hits = search_headers(project.repository, query)
        expected = search_headers(project.repository, query, force_fallback=True)
        assert [hit.header for hit in hits] == [hit.header for hit in expected]
    assert search_headers(project.repository, "AU")  # would be empty via trigram alone


def test_fts_syntax_in_user_input_is_treated_as_text(project, tmp_path):
    path = write_fasta(
        tmp_path / "syntax.fasta",
        [("alpha OR beta|X", "AACC"), ("gamma|X", "GGGG"), ('quo"ted|X', "TTTT")],
    )
    project.link_fasta(str(path))

    assert [hit.header for hit in search_headers(project.repository, "alpha OR beta")] == [
        "alpha OR beta|X"
    ]
    # A bare OR must not behave as a boolean operator.
    assert search_headers(project.repository, "alpha OR gamma") == []
    assert [hit.header for hit in search_headers(project.repository, 'quo"ted')] == ['quo"ted|X']


def test_quote_fts_query_escapes_embedded_quotes():
    assert quote_fts_query('a"b') == '"a""b"'


# ---------------------------------------------------------------------------
# Alignment cache
# ---------------------------------------------------------------------------


def test_alignment_cache_is_reused_for_an_unchanged_source(project, linked):
    first = project.load_alignment(linked.id)
    before = project.alignments.hits
    second = project.load_alignment(linked.id)

    assert second is first
    assert project.alignments.hits == before + 1


def test_alignment_cache_is_invalidated_after_a_change(project, fasta, linked):
    first = project.load_alignment(linked.id)
    assert "extra|NL|Contrast" not in first.sequences

    write_fasta(fasta, [*RECORDS, ("extra|NL|Contrast", "TTTT")])
    project.verifier.forget_all()
    second = project.load_alignment(linked.id)

    assert second is not first
    assert "extra|NL|Contrast" in second.sequences


def test_alignment_cache_is_not_used_when_the_source_is_missing(project, fasta, linked):
    project.load_alignment(linked.id)
    fasta.unlink()
    with pytest.raises(ProjectError):
        project.load_alignment(linked.id)


# ---------------------------------------------------------------------------
# Focal-set semantics
# ---------------------------------------------------------------------------


def test_add_expands_a_query_to_exact_headers(project, linked):
    focal_set = project.repository.create_focal_set("S")
    result = project.add_focal_entries_from_query(focal_set.id, "target")

    assert result["added"] == ["focal_1|AU|Target", "focal_2|GB|Target"]
    stored = [entry.header for entry in project.repository.list_focal_entries(focal_set.id)]
    # The query itself never becomes a member.
    assert "target" not in stored
    assert stored == ["focal_1|AU|Target", "focal_2|GB|Target"]


def test_add_deduplicates_across_queries(project, linked):
    focal_set = project.repository.create_focal_set("S")
    project.add_focal_entries_from_query(focal_set.id, "target")
    again = project.add_focal_entries_from_query(focal_set.id, "focal_1")

    assert again["added"] == []
    assert len(project.repository.list_focal_entries(focal_set.id)) == 2


def test_remove_drops_members_matching_the_query(project, linked):
    focal_set = project.repository.create_focal_set("S")
    project.add_focal_entries_from_query(focal_set.id, "target")
    result = project.remove_focal_entries_by_query(focal_set.id, "focal_2")

    assert result["removed"] == ["focal_2|GB|Target"]
    assert [e.header for e in project.repository.list_focal_entries(focal_set.id)] == [
        "focal_1|AU|Target"
    ]


def test_membership_is_exact_not_substring(project, tmp_path):
    path = write_fasta(
        tmp_path / "prefix.fasta",
        [("ABC123", "AACC"), ("ABC123_extra", "AACC"), ("other|X", "AGCT"), ("other2|X", "CGTT")],
    )
    row, _status = project.link_fasta(str(path))
    focal_set = project.repository.create_focal_set("S")
    project.repository.begin()
    project.repository.add_focal_entries(focal_set.id, ["ABC123"])
    project.repository.commit()

    selector = project.focal_selector(focal_set.id)
    scope = project.build_scope([row.id])
    from molecular_diagnosis.focal import partition_headers

    focal, non_focal = partition_headers(list(scope.sequences), selector)
    assert focal == ["ABC123"]
    assert "ABC123_extra" in non_focal


# ---------------------------------------------------------------------------
# Multi-file scopes
# ---------------------------------------------------------------------------


def test_compatible_multi_file_scope_combines_in_memory(project, fasta, linked, tmp_path):
    other = write_fasta(tmp_path / "b.fasta", [("b1|X", "AAAA"), ("b2|X", "TTTT")])
    other_row, _status = project.link_fasta(str(other))

    scope = project.build_scope([linked.id, other_row.id])
    assert scope.ok
    assert scope.alignment_length == 4
    assert len(scope.sequences) == 6


def test_incompatible_alignment_lengths_are_rejected(project, fasta, linked, tmp_path):
    other = write_fasta(tmp_path / "b.fasta", [("b1|X", "AAAAAA"), ("b2|X", "TTTTTT")])
    other_row, _status = project.link_fasta(str(other))

    scope = project.build_scope([linked.id, other_row.id])
    assert not scope.ok
    assert any(p.code == "INCOMPATIBLE_ALIGNMENT_LENGTHS" for p in scope.problems)


def test_cross_file_duplicate_headers_are_rejected(project, fasta, linked, tmp_path):
    clash = write_fasta(tmp_path / "b.fasta", [("focal_1|AU|Target", "TTTT"), ("b2|X", "AAAA")])
    clash_row, _status = project.link_fasta(str(clash))

    scope = project.build_scope([linked.id, clash_row.id])
    assert not scope.ok
    problem = next(p for p in scope.problems if p.code == "DUPLICATE_HEADER_ACROSS_FILES")
    assert "focal_1|AU|Target" in problem.message
    assert "a.fasta" in (problem.detail or "") and "b.fasta" in (problem.detail or "")


def test_single_file_run_blocked_when_an_entry_is_absent(project, fasta, linked, tmp_path):
    other = write_fasta(tmp_path / "b.fasta", [("elsewhere|QQ|Target", "AACC"), ("b2|X", "GGGG")])
    project.link_fasta(str(other))

    focal_set = project.repository.create_focal_set("S")
    project.add_focal_entries_from_query(focal_set.id, "Target")

    with pytest.raises(ProjectError) as caught:
        project.run_molecular_diagnosis(
            focal_set_id=focal_set.id, fasta_file_ids=[linked.id], single_file=True, options={}
        )
    assert caught.value.code == "FOCAL_ENTRIES_NOT_IN_FILE"
    assert "elsewhere|QQ|Target" in (caught.value.detail or "")


def test_single_file_run_writes_into_the_project_outputs_dir(project, fasta, linked):
    focal_set = project.repository.create_focal_set("S")
    project.add_focal_entries_from_query(focal_set.id, "Target")

    result = project.run_molecular_diagnosis(
        focal_set_id=focal_set.id, fasta_file_ids=[linked.id], single_file=True, options={}
    )
    for value in result["outputs"].values():
        if value:
            assert str(project.outputs_dir) in value


def test_multi_file_run_uses_exact_membership(project, fasta, linked, tmp_path):
    other = write_fasta(tmp_path / "b.fasta", [("b1|X", "AAAA"), ("b2|X", "TTTT")])
    other_row, _status = project.link_fasta(str(other))

    focal_set = project.repository.create_focal_set("S")
    project.add_focal_entries_from_query(focal_set.id, "Target")

    result = project.run_molecular_diagnosis(
        focal_set_id=focal_set.id,
        fasta_file_ids=[linked.id, other_row.id],
        single_file=False,
        options={},
    )
    assert result["sequenceCount"] == 6
    assert result["focalHeaders"] == ["focal_1|AU|Target", "focal_2|GB|Target"]


def test_a_build_without_fts5_still_searches_correctly(tmp_path, monkeypatch):
    """
    The accelerator is optional. On a SQLite without FTS5 (or without the
    trigram tokenizer) the project must still open, index and search — just
    without the index. This simulates that build rather than assuming it.
    """
    from molecular_diagnosis.project import db as db_module

    real = db_module.detect_capabilities
    monkeypatch.setattr(
        db_module,
        "detect_capabilities",
        lambda con=None: db_module.Capabilities(real(con).sqlite_version, False, False),
    )

    service = ProjectService(tmp_path / "nofts")
    try:
        service.open_project("No FTS")
        assert service.repository.fts_enabled is False
        assert not db_module.fts_available(service.repository.con)

        write_fasta(tmp_path / "a.fasta", RECORDS)
        service.link_fasta(str(tmp_path / "a.fasta"))

        hits = search_headers(service.repository, "target")
        assert [hit.header for hit in hits] == ["focal_1|AU|Target", "focal_2|GB|Target"]

        # And the focal workflow that depends on search still works end to end.
        focal_set = service.repository.create_focal_set("S")
        result = service.add_focal_entries_from_query(focal_set.id, "Target")
        assert result["added"] == ["focal_1|AU|Target", "focal_2|GB|Target"]
        assert len(service.repository.locations_for_set(focal_set.id)) == 2
    finally:
        service.close()


def test_capability_reporting_matches_the_installed_fts_table(project):
    """`fts_enabled` must describe reality, not intent."""
    from molecular_diagnosis.project import db as db_module

    assert project.repository.fts_enabled == db_module.fts_available(project.repository.con)
