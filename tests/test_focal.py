"""
Focal-set selection tests.

The first block is CHARACTERIZATION: it pins the behaviour the single-string
pipeline had before multi-selector support, so the refactor of the four former
`target_string in header` call sites can be shown not to have changed it.
"""

import pytest

from molecular_diagnosis.core import find_dmc_information
from molecular_diagnosis.fasta_io import split_focal_headers
from molecular_diagnosis.focal import (
    focal_label,
    header_matches_focal,
    matching_focal_strings,
    normalise_focal_strings,
    partition_headers,
)

HEADERS = [
    "AACTA5253-20|AU|Leptacis",
    "BBHYQ1043-18|GB|Leptacis_tipulae",
    "BCJKR7756-20|NL|Leptacis_tipulae",
    "CGHIJ5502-19|ES|Leptacis_phantasmatica",
    "DJPQR8835-22|GR|Platygaster",
    "ELVWX1057-24|PL|Synopeas",
]

SEQUENCES = {
    "focal_1": "AACC",
    "focal_2": "AACC",
    "other_1": "AGCT",
    "other_2": "CGTT",
}


# --------------------------------------------------------------------------
# Characterization: single-string behaviour must not drift
# --------------------------------------------------------------------------


def test_split_focal_headers_single_string_is_substring_containment() -> None:
    sequences = {header: "ACGT" for header in HEADERS}
    focal, non_focal = split_focal_headers(sequences, "Leptacis_tipulae")

    assert focal == [
        "BBHYQ1043-18|GB|Leptacis_tipulae",
        "BCJKR7756-20|NL|Leptacis_tipulae",
    ]
    assert len(non_focal) == 4


def test_split_focal_headers_single_string_is_case_sensitive() -> None:
    sequences = {header: "ACGT" for header in HEADERS}

    with pytest.raises(ValueError, match="No sequences matched"):
        split_focal_headers(sequences, "leptacis_tipulae")


def test_find_dmc_information_accepts_a_plain_string() -> None:
    # The exact call shape used by the pre-existing test suite.
    result = find_dmc_information(SEQUENCES, "focal")

    assert result.fixed_count == 4
    assert result.candidate_count == 4
    assert result.single == [1, 3]


def test_find_dmc_information_string_and_one_element_list_agree() -> None:
    from_string = find_dmc_information(SEQUENCES, "focal")
    from_list = find_dmc_information(SEQUENCES, ["focal"])

    assert from_list.single == from_string.single
    assert from_list.states == from_string.states
    assert from_list.candidate_count == from_string.candidate_count
    assert from_list.stop_reason == from_string.stop_reason


# --------------------------------------------------------------------------
# normalise_focal_strings
# --------------------------------------------------------------------------


def test_normalise_accepts_a_single_string() -> None:
    assert normalise_focal_strings("Leptacis") == ["Leptacis"]


def test_normalise_accepts_a_list() -> None:
    assert normalise_focal_strings(["Leptacis", "Synopeas"]) == ["Leptacis", "Synopeas"]


def test_normalise_strips_surrounding_whitespace() -> None:
    # The original pipeline stripped its single target string; keep that.
    assert normalise_focal_strings(["  Leptacis  "]) == ["Leptacis"]


def test_normalise_deduplicates_preserving_first_seen_order() -> None:
    assert normalise_focal_strings(["b", "a", "b", "c", "a"]) == ["b", "a", "c"]


def test_normalise_deduplicates_after_stripping() -> None:
    assert normalise_focal_strings(["Leptacis", " Leptacis "]) == ["Leptacis"]


def test_normalise_rejects_an_empty_selector() -> None:
    with pytest.raises(ValueError, match="cannot be empty"):
        normalise_focal_strings(["Leptacis", "   "])


def test_normalise_rejects_an_empty_single_string() -> None:
    with pytest.raises(ValueError, match="cannot be empty"):
        normalise_focal_strings("")


def test_normalise_rejects_an_empty_list() -> None:
    with pytest.raises(ValueError, match="No identifier string entered."):
        normalise_focal_strings([])


def test_normalise_keeps_a_selector_containing_a_semicolon() -> None:
    # ';' is only a display separator in the UI; it is ordinary data here.
    assert normalise_focal_strings(["a;b"]) == ["a;b"]


# --------------------------------------------------------------------------
# header_matches_focal
# --------------------------------------------------------------------------


def test_one_selector_matches_by_substring() -> None:
    assert header_matches_focal("AACTA5253-20|AU|Leptacis", ["Leptacis"]) is True
    assert header_matches_focal("AACTA5253-20|AU|Leptacis", ["Lept"]) is True
    assert header_matches_focal("AACTA5253-20|AU|Leptacis", ["AU"]) is True


def test_one_selector_is_case_sensitive() -> None:
    assert header_matches_focal("AACTA5253-20|AU|Leptacis", ["leptacis"]) is False


def test_multiple_selectors_act_as_or() -> None:
    focal, non_focal = partition_headers(HEADERS, ["Platygaster", "Synopeas"])

    assert focal == [
        "DJPQR8835-22|GR|Platygaster",
        "ELVWX1057-24|PL|Synopeas",
    ]
    assert len(non_focal) == 4


def test_multiple_selectors_match_a_header_matched_by_either() -> None:
    header = "BBHYQ1043-18|GB|Leptacis_tipulae"
    assert header_matches_focal(header, ["Synopeas", "Leptacis_tipulae"]) is True
    # A header need only satisfy one selector, never all of them.
    assert header_matches_focal(header, ["Synopeas", "Platygaster"]) is False


def test_overlapping_selectors_do_not_double_count_headers() -> None:
    # "Leptacis" is a substring of "Leptacis_tipulae", so the sets overlap.
    focal, non_focal = partition_headers(HEADERS, ["Leptacis", "Leptacis_tipulae"])

    assert len(focal) == 4
    assert len(non_focal) == 2
    assert len(focal) == len(set(focal))


def test_selector_matching_no_header() -> None:
    focal, non_focal = partition_headers(HEADERS, ["Nothosaurus"])

    assert focal == []
    assert non_focal == HEADERS


def test_duplicate_selectors_behave_as_one() -> None:
    once = partition_headers(HEADERS, normalise_focal_strings(["Leptacis"]))
    twice = partition_headers(HEADERS, normalise_focal_strings(["Leptacis", "Leptacis"]))

    assert once == twice


def test_selector_containing_a_semicolon_matches_literally() -> None:
    headers = ["weird;header|AU|Leptacis", "AACTA5253-20|AU|Leptacis"]
    focal, non_focal = partition_headers(headers, ["weird;header"])

    assert focal == ["weird;header|AU|Leptacis"]
    assert non_focal == ["AACTA5253-20|AU|Leptacis"]


def test_partition_preserves_input_order() -> None:
    focal, non_focal = partition_headers(HEADERS, ["Leptacis_tipulae"])

    assert focal == [h for h in HEADERS if h in focal]
    assert non_focal == [h for h in HEADERS if h in non_focal]


def test_matching_focal_strings_reports_which_selectors_hit() -> None:
    header = "BBHYQ1043-18|GB|Leptacis_tipulae"
    assert matching_focal_strings(header, ["Leptacis", "Synopeas", "GB"]) == [
        "Leptacis",
        "GB",
    ]


# --------------------------------------------------------------------------
# focal_label
# --------------------------------------------------------------------------


def test_focal_label_of_a_single_selector_is_that_selector() -> None:
    # Guarantees single-selector runs keep byte-identical report/FASTA naming.
    assert focal_label(["Leptacis_tipulae"]) == "Leptacis_tipulae"


def test_focal_label_joins_multiple_selectors() -> None:
    assert focal_label(["Leptacis", "Synopeas"]) == "Leptacis+Synopeas"


# --------------------------------------------------------------------------
# Multi-selector behaviour through the public splitter
# --------------------------------------------------------------------------


def test_split_focal_headers_accepts_a_list() -> None:
    sequences = {header: "ACGT" for header in HEADERS}
    focal, non_focal = split_focal_headers(sequences, ["Platygaster", "Synopeas"])

    assert focal == [
        "DJPQR8835-22|GR|Platygaster",
        "ELVWX1057-24|PL|Synopeas",
    ]
    assert len(non_focal) == 4


def test_split_focal_headers_rejects_a_list_matching_nothing() -> None:
    sequences = {header: "ACGT" for header in HEADERS}

    with pytest.raises(ValueError, match="No sequences matched"):
        split_focal_headers(sequences, ["Nothosaurus", "Absent"])


def test_split_focal_headers_rejects_a_list_matching_everything() -> None:
    sequences = {header: "ACGT" for header in HEADERS}

    with pytest.raises(ValueError, match="All sequences match"):
        split_focal_headers(sequences, ["|"])


def test_find_dmc_information_with_multiple_selectors() -> None:
    sequences = {
        "alpha_1": "AACC",
        "beta_1": "AACC",
        "other_1": "AGCT",
        "other_2": "CGTT",
    }

    # Two selectors that together pick out exactly the same focal group as the
    # single-string case in the characterization test above.
    result = find_dmc_information(sequences, ["alpha", "beta"])

    assert result.fixed_count == 4
    assert result.candidate_count == 4
    assert result.single == [1, 3]
