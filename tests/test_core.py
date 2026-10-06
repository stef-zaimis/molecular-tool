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
    assert result.single == [1, 3]


def test_find_dmc_information_stops_before_pairs_once_singles_are_found() -> None:
    """
    The search is a ladder, not an exhaustive sweep.

    Sites 1 and 3 are each diagnostic on their own, and the minimum combination
    length is 1, so the search stops at length 1 having satisfied its goal. It
    therefore never enumerates the 6 possible pairs, and `pairs_tested` is 0.

    An earlier version of this test asserted `pairs_tested == 6`, from before
    the early stop existed. The behaviour below is the verified current one:
    the assertion was stale, not the implementation. A smaller pair count is
    the point of the pruning, and "tested" means "actually evaluated".
    """
    sequences = {
        "focal_1": "AACC",
        "focal_2": "AACC",
        "other_1": "AGCT",
        "other_2": "CGTT",
    }

    result = find_dmc_information(sequences, "focal")

    assert result.stop_reason == "found_at_or_above_minimum_length"
    assert result.stopped_at_length == 1
    assert result.combinations_tested_by_length == {1: 4}
    assert result.pairs_tested == 0
    assert result.pairs == []

    # Forcing the search to start at length 2 does enumerate every pair, which
    # is what shows the zero above is early stopping rather than a broken count.
    from_pairs = find_dmc_information(sequences, "focal", start_combination_length=2)
    assert from_pairs.pairs_tested == 6


def test_format_diag() -> None:
    assert format_diag([0, 2], "ACGT") == "1:A, 3:G"


def test_format_diag_empty() -> None:
    assert format_diag([], "ACGT") == "None"