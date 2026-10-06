"""
Vetting a FASTA before it is linked, and the per-source lock.

Two invariants live here:

  * a file the user cannot use must never become a `fasta_file` row — not even
    briefly, and not even when the failure happens after the insert;
  * a locked source stays fully analysable but cannot be unlinked, and that
    refusal is the backend's, not a hidden button's.
"""

from __future__ import annotations

import pytest

from molecular_diagnosis.project.service import (
    ProjectError,
    ProjectService,
    validate_fasta_candidate,
)

ALIGNED = [
    ("focal_1|AU|Target", "AACC"),
    ("focal_2|GB|Target", "AACC"),
    ("other_1|FR|Contrast", "AGCT"),
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


# ---------------------------------------------------------------------------
# validate_fasta_candidate
# ---------------------------------------------------------------------------


def test_a_valid_aligned_fasta_reports_the_metadata_the_row_needs(tmp_path):
    path = write_fasta(tmp_path / "good.fasta", ALIGNED)

    candidate = validate_fasta_candidate(str(path))

    assert candidate["displayName"] == "good.fasta"
    assert candidate["sequenceCount"] == 3
    assert candidate["alignmentLength"] == 4
    assert candidate["path"] == str(path)


def test_validation_needs_no_project(tmp_path):
    """The new-project screen runs before any database exists."""
    path = write_fasta(tmp_path / "good.fasta", ALIGNED)
    # No ProjectService anywhere in this test.
    assert validate_fasta_candidate(str(path))["sequenceCount"] == 3


def test_a_missing_file_is_refused(tmp_path):
    with pytest.raises(ProjectError) as caught:
        validate_fasta_candidate(str(tmp_path / "nope.fasta"))
    assert caught.value.code == "FASTA_NOT_FOUND"


def test_a_directory_is_refused(tmp_path):
    with pytest.raises(ProjectError) as caught:
        validate_fasta_candidate(str(tmp_path))
    assert caught.value.code == "FASTA_NOT_FOUND"


def test_a_file_that_is_not_a_fasta_is_refused(tmp_path):
    """Judged by content, not by extension."""
    path = tmp_path / "notes.fasta"
    path.write_text("id,count\nA,1\nB,2\n", encoding="utf-8")

    with pytest.raises(ProjectError) as caught:
        validate_fasta_candidate(str(path))
    assert caught.value.code == "FASTA_EMPTY"


def test_an_empty_file_is_refused(tmp_path):
    path = tmp_path / "empty.fasta"
    path.write_text("", encoding="utf-8")

    with pytest.raises(ProjectError) as caught:
        validate_fasta_candidate(str(path))
    assert caught.value.code == "FASTA_EMPTY"


def test_an_unaligned_fasta_is_refused_with_the_phrasing_the_ui_shows(tmp_path):
    path = write_fasta(
        tmp_path / "ragged.fasta", [("a|X", "AACC"), ("b|X", "AACCTT"), ("c|X", "AG")]
    )

    with pytest.raises(ProjectError) as caught:
        validate_fasta_candidate(str(path))

    assert caught.value.code == "FASTA_NOT_ALIGNED"
    assert caught.value.message == "The FASTA file needs to be aligned."
    assert "ragged.fasta" in (caught.value.detail or "")


def test_a_valid_fasta_without_a_fasta_extension_is_accepted(tmp_path):
    """The converse of judging by extension: a good file must not be rejected."""
    path = write_fasta(tmp_path / "alignment.txt", ALIGNED)
    assert validate_fasta_candidate(str(path))["alignmentLength"] == 4


# ---------------------------------------------------------------------------
# A failed link leaves no source row
# ---------------------------------------------------------------------------


def test_linking_an_unaligned_fasta_leaves_no_source_row(project, tmp_path):
    path = write_fasta(tmp_path / "ragged.fasta", [("a|X", "AACC"), ("b|X", "AACCTT")])

    with pytest.raises(ProjectError) as caught:
        project.link_fasta(str(path))

    assert caught.value.code == "FASTA_NOT_ALIGNED"
    assert project.repository.list_fasta_files() == []


def test_linking_an_empty_fasta_leaves_no_source_row(project, tmp_path):
    path = tmp_path / "empty.fasta"
    path.write_text("", encoding="utf-8")

    with pytest.raises(ProjectError):
        project.link_fasta(str(path))

    assert project.repository.list_fasta_files() == []


def test_linking_a_missing_file_leaves_no_source_row(project, tmp_path):
    with pytest.raises(ProjectError):
        project.link_fasta(str(tmp_path / "gone.fasta"))

    assert project.repository.list_fasta_files() == []


def test_a_failure_after_the_insert_still_leaves_no_row(project, tmp_path, monkeypatch):
    """
    The race the rollback exists for.

    The candidate passes validation, then the file changes before it is
    indexed. The row was already inserted at that point, so it has to be
    removed again — otherwise the project keeps a link to something that was
    never successfully read.
    """
    path = write_fasta(tmp_path / "racy.fasta", ALIGNED)

    from molecular_diagnosis.project import service as service_module

    real_reindex = service_module.ProjectService.reindex

    def explode(self, file_id, **kwargs):
        raise ProjectError("FASTA_UNREADABLE", "The file vanished mid-index.")

    monkeypatch.setattr(service_module.ProjectService, "reindex", explode)

    with pytest.raises(ProjectError):
        project.link_fasta(str(path))

    assert project.repository.list_fasta_files() == []

    # And the same path links cleanly once the problem is gone.
    monkeypatch.setattr(service_module.ProjectService, "reindex", real_reindex)
    row, _status = project.link_fasta(str(path))
    assert [item.id for item in project.repository.list_fasta_files()] == [row.id]


def test_a_failed_relink_of_an_existing_row_leaves_that_row_alone(project, tmp_path, monkeypatch):
    """A row we did not create is not ours to delete."""
    path = write_fasta(tmp_path / "good.fasta", ALIGNED)
    row, _status = project.link_fasta(str(path))

    from molecular_diagnosis.project import service as service_module

    monkeypatch.setattr(
        service_module.ProjectService,
        "reindex",
        lambda self, file_id, **kwargs: (_ for _ in ()).throw(
            ProjectError("FASTA_UNREADABLE", "boom")
        ),
    )
    # Force the "already linked, needs reindexing" branch.
    project.repository.begin()
    project.repository.con.execute(
        "UPDATE fasta_file SET indexed_sha256 = NULL, indexed_size_bytes = NULL,"
        " indexed_mtime_ns = NULL, index_revision = 0 WHERE id = ?",
        (row.id,),
    )
    project.repository.commit()

    with pytest.raises(ProjectError):
        project.link_fasta(str(path))

    assert [item.id for item in project.repository.list_fasta_files()] == [row.id]


# ---------------------------------------------------------------------------
# The per-source lock
# ---------------------------------------------------------------------------


@pytest.fixture()
def linked(project, tmp_path):
    row, _status = project.link_fasta(str(write_fasta(tmp_path / "a.fasta", ALIGNED)))
    return row


def test_a_source_starts_unlocked(linked):
    assert linked.locked is False


def test_locking_persists_and_is_reported_in_the_status(project, linked):
    status = project.set_fasta_file_locked(linked.id, True)

    assert status.locked is True
    assert status.to_payload()["locked"] is True
    assert project.repository.get_fasta_file(linked.id).locked is True

    assert project.set_fasta_file_locked(linked.id, False).locked is False
    assert project.repository.get_fasta_file(linked.id).locked is False


def test_a_lock_survives_reopening_the_project(project, linked, tmp_path):
    project.set_fasta_file_locked(linked.id, True)
    project.close()

    reopened = ProjectService(tmp_path / "proj", create_if_missing=False)
    try:
        assert reopened.repository.get_fasta_file(linked.id).locked is True
    finally:
        reopened.close()


def test_a_locked_source_cannot_be_unlinked(project, linked):
    project.set_fasta_file_locked(linked.id, True)

    with pytest.raises(ProjectError) as caught:
        project.unlink_fasta(linked.id)

    assert caught.value.code == "FASTA_FILE_LOCKED"
    assert project.repository.get_fasta_file(linked.id) is not None

    # Unlocking restores the ability to remove it.
    project.set_fasta_file_locked(linked.id, False)
    project.unlink_fasta(linked.id)
    assert project.repository.get_fasta_file(linked.id) is None


def test_a_locked_source_is_still_analysable(project, linked):
    """Locking protects the LINK, not the data."""
    project.set_fasta_file_locked(linked.id, True)

    focal = project.save_focal_set(
        focal_set_id=None, title="S", headers=["focal_1|AU|Target", "focal_2|GB|Target"]
    )
    result = project.run_molecular_diagnosis(
        focal_set_id=focal["id"], fasta_file_ids=[linked.id], single_file=True, options={}
    )

    assert result["sequenceCount"] == 3


def test_locking_an_unknown_file_is_a_coded_refusal(project):
    with pytest.raises(ProjectError) as caught:
        project.set_fasta_file_locked("nope", True)
    assert caught.value.code == "UNKNOWN_FILE"
