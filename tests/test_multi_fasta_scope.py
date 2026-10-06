"""
What "All files" ACTUALLY does today, written down.

This file documents current behaviour; it does not endorse it. The overlapping
multi-FASTA rules (what a repeated header across files should mean, whether an
identical sequence is a duplicate or a confirmation, how a full alignment plus
a focal-only subset should combine) are an open design question. These tests
exist so that design can start from what the code does rather than from what we
assume it does, and so a future change to the policy is a deliberate edit to
named expectations rather than a silent behaviour change.

The reported case: File A is a full alignment, File B is the focal subset of
it. Selecting All files refuses the run.
"""

from __future__ import annotations

import pytest

from molecular_diagnosis.project.service import ProjectError, ProjectService

# A full alignment: two focal sequences plus a contrast set.
FULL = [
    ("Leptacis_tipulae_1", "AACCGGTT"),
    ("Leptacis_tipulae_2", "AACCGGTT"),
    ("Leptacis_other_1", "AGCCGGTT"),
    ("Leptacis_other_2", "AGCCGGTA"),
]

# The focal-only subset of it: the same records, byte for byte.
SUBSET = [
    ("Leptacis_tipulae_1", "AACCGGTT"),
    ("Leptacis_tipulae_2", "AACCGGTT"),
]

# The same headers, DIFFERENT residues — an edited or realigned copy.
SUBSET_EDITED = [
    ("Leptacis_tipulae_1", "TTCCGGTT"),
    ("Leptacis_tipulae_2", "TTCCGGTT"),
]

# No headers in common with FULL.
DISJOINT = [
    ("Outgroup_1", "AGGGGGTT"),
    ("Outgroup_2", "AGGGGGTA"),
]

FOCAL_HEADERS = ["Leptacis_tipulae_1", "Leptacis_tipulae_2"]


def write_fasta(path, records):
    path.write_text(
        "".join(f">{header}\n{sequence}\n" for header, sequence in records), encoding="utf-8"
    )
    return path


@pytest.fixture()
def project(tmp_path):
    service = ProjectService(tmp_path / "proj")
    service.open_project("Multi FASTA")
    yield service
    service.close()


def link(project, tmp_path, name, records):
    row, _status = project.link_fasta(str(write_fasta(tmp_path / name, records)))
    return row.id


def focal_set(project, headers=FOCAL_HEADERS):
    saved = project.save_focal_set(focal_set_id=None, title="tipulae", headers=headers)
    return saved["id"]


# ---------------------------------------------------------------------------
# The reported failure
# ---------------------------------------------------------------------------


def test_full_alignment_plus_its_focal_subset_is_refused_under_all_files(project, tmp_path):
    """
    THE reported case, reproduced.

    Every header in the subset is also in the full file, so `build_scope`
    records one DUPLICATE_HEADER_ACROSS_FILES problem per shared header and the
    run never starts. Nothing about the sequences is examined: the refusal is
    about the header appearing twice, not about disagreement.
    """
    full = link(project, tmp_path, "full.fasta", FULL)
    subset = link(project, tmp_path, "subset.fasta", SUBSET)
    set_id = focal_set(project)

    with pytest.raises(ProjectError) as raised:
        project.run_molecular_diagnosis(
            focal_set_id=set_id,
            fasta_file_ids=[full, subset],
            single_file=False,
            options={},
        )

    assert raised.value.code == "DUPLICATE_HEADER_ACROSS_FILES"
    assert "full.fasta and subset.fasta" in (raised.value.detail or "")


def test_the_refusal_is_the_same_when_the_repeated_records_disagree(project, tmp_path):
    """
    Identical or contradictory, the answer today is the same.

    No sequence comparison happens, so "the same record twice" and "two
    different records under one name" are not currently distinguished. That is
    precisely the distinction the next pass has to decide about.
    """
    full = link(project, tmp_path, "full.fasta", FULL)
    edited = link(project, tmp_path, "edited.fasta", SUBSET_EDITED)
    set_id = focal_set(project)

    with pytest.raises(ProjectError) as raised:
        project.run_molecular_diagnosis(
            focal_set_id=set_id,
            fasta_file_ids=[full, edited],
            single_file=False,
            options={},
        )

    assert raised.value.code == "DUPLICATE_HEADER_ACROSS_FILES"


def test_each_shared_header_produces_its_own_problem_but_the_detail_is_capped(project, tmp_path):
    """
    One problem PER shared header, and the message says how many were elided.

    With a real alignment this is hundreds of near-identical sentences; the
    refusal used to concatenate all of them into the message the UI shows.
    """
    many = [(f"Leptacis_tipulae_{index}", "AACCGGTT") for index in range(20)]
    full = link(project, tmp_path, "full.fasta", many + [("Leptacis_other_1", "AGCCGGTT")])
    subset = link(project, tmp_path, "subset.fasta", many)

    scope = project.build_scope([full, subset])

    assert len(scope.problems) == 20
    assert {problem.code for problem in scope.problems} == {"DUPLICATE_HEADER_ACROSS_FILES"}

    set_id = focal_set(project, [header for header, _ in many])
    with pytest.raises(ProjectError) as raised:
        project.run_molecular_diagnosis(
            focal_set_id=set_id,
            fasta_file_ids=[full, subset],
            single_file=False,
            options={},
        )
    assert "(+12 more)" in (raised.value.detail or "")


# ---------------------------------------------------------------------------
# What the backend actually builds from several files
# ---------------------------------------------------------------------------


def test_files_are_pooled_into_one_header_keyed_dictionary(project, tmp_path):
    """
    Disjoint files are UNIONED in memory — no temporary combined FASTA, no
    re-parsing, and no record ids: the scope is a `header -> sequence` dict, in
    selection order.
    """
    full = link(project, tmp_path, "full.fasta", FULL)
    other = link(project, tmp_path, "outgroups.fasta", DISJOINT)

    scope = project.build_scope([full, other])

    assert scope.ok
    assert list(scope.sequences) == [
        "Leptacis_tipulae_1",
        "Leptacis_tipulae_2",
        "Leptacis_other_1",
        "Leptacis_other_2",
        "Outgroup_1",
        "Outgroup_2",
    ]
    assert scope.alignment_length == 8


def test_a_repeated_header_keeps_the_FIRST_file_s_record_in_the_pool(project, tmp_path):
    """
    Where the record came from, if the run were not refused.

    The combined dict keeps the first file's record and the second is skipped;
    the skip is what becomes the problem. Nothing is concatenated, renamed or
    merged. This is the mechanism the next pass will be choosing a policy for —
    the run is refused today precisely so that this silent first-wins outcome
    never reaches the science.
    """
    full = link(project, tmp_path, "full.fasta", FULL)
    edited = link(project, tmp_path, "edited.fasta", SUBSET_EDITED)

    scope = project.build_scope([full, edited])

    assert not scope.ok
    assert scope.sequences["Leptacis_tipulae_1"] == "AACCGGTT"  # from full.fasta
    # Selecting them the other way round keeps the other record — order is the
    # only thing deciding, which is exactly why it is refused.
    reversed_scope = project.build_scope([edited, full])
    assert reversed_scope.sequences["Leptacis_tipulae_1"] == "TTCCGGTT"


def test_every_file_is_validated_on_its_own_before_anything_is_pooled(project, tmp_path):
    """
    Per-file validity is checked FIRST, and a bad file is named.

    A file with repeated headers inside itself cannot take part, because
    `sequences` is a dict and one of its records would vanish silently.
    """
    full = link(project, tmp_path, "full.fasta", FULL)
    internal_duplicate = link(
        project,
        tmp_path,
        "dupes.fasta",
        [("Rep", "AACCGGTT"), ("Rep", "TTCCGGTT"), ("Unique", "AGCCGGTT")],
    )

    scope = project.build_scope([full, internal_duplicate])

    assert [problem.code for problem in scope.problems] == ["DUPLICATE_HEADER_IN_FILE"]
    assert "dupes.fasta" in (scope.problems[0].detail or "")
    # The offending file contributes NOTHING; the healthy one is still pooled.
    assert "Unique" not in scope.sequences
    assert "Leptacis_other_1" in scope.sequences


def test_files_of_different_alignment_lengths_are_refused_together(project, tmp_path):
    full = link(project, tmp_path, "full.fasta", FULL)
    longer = link(
        project,
        tmp_path,
        "longer.fasta",
        [("Other_1", "AACCGGTTAA"), ("Other_2", "AACCGGTTAC")],
    )

    scope = project.build_scope([full, longer])

    assert [problem.code for problem in scope.problems] == ["INCOMPATIBLE_ALIGNMENT_LENGTHS"]
    assert "full.fasta: 8" in (scope.problems[0].detail or "")


def test_a_multi_file_run_over_disjoint_files_still_works(project, tmp_path):
    """The All-files path itself is not broken — overlap is what refuses it."""
    full = link(project, tmp_path, "full.fasta", FULL)
    other = link(project, tmp_path, "outgroups.fasta", DISJOINT)
    set_id = focal_set(project)

    result = project.run_molecular_diagnosis(
        focal_set_id=set_id,
        fasta_file_ids=[full, other],
        single_file=False,
        options={},
    )

    assert result["sequenceCount"] == 6
    assert result["fastaFileIds"] == [full, other]


# ---------------------------------------------------------------------------
# Established All-files UI semantics, which must not have moved
# ---------------------------------------------------------------------------


def test_presence_under_all_files_is_green_when_the_header_is_in_any_file(project, tmp_path):
    """
    Green if present in ANY selected file, and no orange under All files.

    `selected_file_id=None` is what removes the "present, but elsewhere" case:
    with every file selected there is no elsewhere.
    """
    from molecular_diagnosis.project.locations import PresenceState

    full = link(project, tmp_path, "full.fasta", FULL)
    subset = link(project, tmp_path, "subset.fasta", SUBSET)
    set_id = focal_set(project, FOCAL_HEADERS + ["Leptacis_other_1"])

    presence = project.focal_presence(
        set_id, selected_file_id=None, fasta_file_ids=[full, subset]
    )

    states = {entry.header: entry.state for entry in presence}
    assert states["Leptacis_tipulae_1"] is PresenceState.PRESENT_CURRENT
    # Only in the full file, and still green: All files has no "elsewhere".
    assert states["Leptacis_other_1"] is PresenceState.PRESENT_CURRENT
    assert PresenceState.PRESENT_OTHER not in set(states.values())


def test_search_under_all_files_covers_the_union(project, tmp_path):
    full = link(project, tmp_path, "full.fasta", FULL)
    other = link(project, tmp_path, "outgroups.fasta", DISJOINT)

    hits = project.search_headers("_", fasta_file_ids=[full, other])["hits"]

    files = {hit["fastaFileId"] for hit in hits}
    assert files == {full, other}
