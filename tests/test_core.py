from molecular_diagnosis.core import (
    compute_match_score,
    compute_similarity,
    extract_sites,
    find_dmc_information,
    format_diag,
    score_state,
)


def test_extract_sites() -> None:
    assert extract_sites("ACGT", [0, 2, 3]) == ["A", "G", "T"]


def test_score_state_exact_match() -> None:
    assert score_state("A", "A") == 1.0


def test_score_state_ambiguous_match() -> None:
    assert score_state("A", "R") == 0.5


def test_score_state_no_match() -> None:
    assert score_state("C", "R") == 0.0


def test_compute_similarity() -> None:
    assert compute_similarity(["A", "C"], ["A", "N"]) == 0.625


def test_compute_match_score() -> None:
    assert compute_match_score(["A", "C"], ["A", "N"]) == 1.25


def test_find_dmc_information() -> None:
    sequences = {
        "focal_1": "AACC",
        "focal_2": "AACC",
        "other_1": "AGCT",
        "other_2": "CGTT",
    }

    result = find_dmc_information(sequences, "focal")

    assert result.fixed_count == 4
    assert result.candidate_count == 4
    assert 0 in result.single
    assert result.pairs_tested == 6


def test_format_diag() -> None:
    assert format_diag([0, 2], "ACGT") == "1:A, 3:G"


def test_format_diag_empty() -> None:
    assert format_diag([], "ACGT") == "None"