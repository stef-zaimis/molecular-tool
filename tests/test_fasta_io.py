from pathlib import Path

import pytest

from molecular_diagnosis.fasta_io import (
    parse_fasta,
    split_focal_headers,
    validate_aligned_fasta,
)


def test_parse_fasta(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nACGT\n>seq2\nTGCA\n", encoding="utf-8")

    sequences = parse_fasta(fasta)

    assert sequences == {
        "seq1": "ACGT",
        "seq2": "TGCA",
    }


def test_parse_fasta_multiline_sequence(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nAC\nGT\n", encoding="utf-8")

    sequences = parse_fasta(fasta)

    assert sequences == {"seq1": "ACGT"}


def test_validate_aligned_fasta_returns_length() -> None:
    sequences = {
        "seq1": "ACGT",
        "seq2": "TGCA",
    }

    assert validate_aligned_fasta(sequences) == 4


def test_validate_aligned_fasta_rejects_empty() -> None:
    with pytest.raises(ValueError, match="No sequences"):
        validate_aligned_fasta({})


def test_validate_aligned_fasta_rejects_unequal_lengths() -> None:
    sequences = {
        "seq1": "ACGT",
        "seq2": "TGC",
    }

    with pytest.raises(ValueError, match="not all the same length"):
        validate_aligned_fasta(sequences)


def test_split_focal_headers() -> None:
    sequences = {
        "target_1": "ACGT",
        "target_2": "ACGT",
        "other": "TGCA",
    }

    focal, non_focal = split_focal_headers(sequences, "target")

    assert focal == ["target_1", "target_2"]
    assert non_focal == ["other"]


def test_split_focal_headers_rejects_no_focal() -> None:
    sequences = {
        "other": "ACGT",
    }

    with pytest.raises(ValueError, match="No sequences matched"):
        split_focal_headers(sequences, "target")


def test_split_focal_headers_rejects_no_non_focal() -> None:
    sequences = {
        "target_1": "ACGT",
        "target_2": "ACGT",
    }

    with pytest.raises(ValueError, match="All sequences match"):
        split_focal_headers(sequences, "target")