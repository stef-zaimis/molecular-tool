"""
The exhaustive 5-site optimisation, pinned against a reference implementation.

`find_best_five_site_sets` is the slowest thing in the codebase and therefore
the thing most likely to be "improved". This file is the guard: it recomputes
every result the long way — through the untouched public `compute_metrics`,
one combination at a time, exactly as the original loop did — and demands
**exact** equality, floats included.

That matters more than it looks. The scores are sums of `1/len(possibilities)`
terms, so changing the ORDER of accumulation (or deduplicating sequences, or
reordering the search) can move the last bits of a score and silently change
which 5-site set is reported as best. Nothing here is written with `pytest.approx`
on purpose.

The behaviours pinned:

* every combination of `sites` is considered, in `itertools.combinations` order;
* the FIRST combination achieving a minimum wins (strict `<`), so ties are
  resolved by search order;
* strict bases, IUPAC ambiguity, gaps and `?` all score as they did;
* `ExactHeaders` selects by exact identity, plain selectors by substring;
* `diagnostic_states` supplied and not supplied are both honoured;
* the reference sequence and every focal sequence are excluded from scoring;
* fewer than five candidate sites returns the empty result without searching.
"""

from __future__ import annotations

import random
from itertools import combinations

import pytest

from molecular_diagnosis.core import compute_metrics, find_best_five_site_sets
from molecular_diagnosis.focal import ExactHeaders
from molecular_diagnosis.models import FiveSiteResult


def reference_five_site(
    sequences: dict[str, str],
    ref_id: str,
    sites: list[int],
    focal_strings,
    diagnostic_states: dict[int, str] | None = None,
) -> FiveSiteResult:
    """
    The original loop, verbatim, over the untouched `compute_metrics`.

    Deliberately naive: this is the definition of correct for this test file,
    so it must not share any of the implementation's shortcuts.
    """
    total_combinations_tested = 0
    best_gap_score = None
    best_gap_sites = None
    best_avg_score = None
    best_avg_sites = None

    if len(sites) < 5:
        return FiveSiteResult(
            total_combinations_tested=0,
            best_gap_score=None,
            best_gap_sites=None,
            best_avg_score=None,
            best_avg_sites=None,
        )

    for combo in combinations(sites, 5):
        total_combinations_tested += 1
        max_similarity, avg_similarity, _gap = compute_metrics(
            sequences=sequences,
            ref_id=ref_id,
            sites=combo,
            focal_strings=focal_strings,
            diagnostic_states=diagnostic_states,
        )
        if best_gap_score is None or max_similarity < best_gap_score:
            best_gap_score = max_similarity
            best_gap_sites = combo
        if best_avg_score is None or avg_similarity < best_avg_score:
            best_avg_score = avg_similarity
            best_avg_sites = combo

    return FiveSiteResult(
        total_combinations_tested=total_combinations_tested,
        best_gap_score=best_gap_score,
        best_gap_sites=best_gap_sites,
        best_avg_score=best_avg_score,
        best_avg_sites=best_avg_sites,
    )


def assert_identical(
    sequences, ref_id, sites, focal_strings, diagnostic_states=None, *, case: str = ""
) -> FiveSiteResult:
    """Run both implementations and demand bit-for-bit identical results."""
    expected = reference_five_site(
        sequences, ref_id, sites, focal_strings, diagnostic_states
    )
    actual = find_best_five_site_sets(
        sequences=sequences,
        ref_id=ref_id,
        sites=sites,
        focal_strings=focal_strings,
        diagnostic_states=diagnostic_states,
    )

    assert actual.total_combinations_tested == expected.total_combinations_tested, case
    # Exact float equality, and `repr` in the message so a last-bit difference
    # is readable rather than printing as the same number.
    assert actual.best_gap_score == expected.best_gap_score, (
        f"{case}: gap {actual.best_gap_score!r} != {expected.best_gap_score!r}"
    )
    assert actual.best_avg_score == expected.best_avg_score, (
        f"{case}: avg {actual.best_avg_score!r} != {expected.best_avg_score!r}"
    )
    # The winning combinations themselves, which is where a tie-break change
    # would show up even when the scores match.
    assert actual.best_gap_sites == expected.best_gap_sites, case
    assert actual.best_avg_sites == expected.best_avg_sites, case
    return actual


# ---------------------------------------------------------------------------
# Alphabets
# ---------------------------------------------------------------------------

STRICT = "ACGT"
#: Ambiguity codes, gaps and the unknown marker: everything `score_state`
#: distinguishes, including the states that score zero.
MIXED = "ACGTRYSWKMBDHVN-?"


def alignment(n_focal: int, n_other: int, length: int, alphabet: str, seed: int):
    rng = random.Random(seed)
    sequences: dict[str, str] = {}
    for index in range(n_focal):
        sequences[f"focal_{index}|AU|Target"] = "".join(
            rng.choice(alphabet) for _ in range(length)
        )
    for index in range(n_other):
        sequences[f"other_{index}|FR|Contrast"] = "".join(
            rng.choice(alphabet) for _ in range(length)
        )
    return sequences


# ---------------------------------------------------------------------------
# The states
# ---------------------------------------------------------------------------


def test_strict_bases_only():
    sequences = alignment(2, 6, 12, STRICT, seed=1)
    assert_identical(sequences, "focal_0|AU|Target", list(range(8)), ["Target"], case="strict")


def test_iupac_ambiguity_gaps_and_unknown():
    """
    The whole alphabet, including the states that score 0.

    `score_state` gives `1/len(possibilities)` for a match, so a `N` contributes
    0.25 and a gap contributes nothing. Those fractions are what make the sums
    float-sensitive.
    """
    sequences = alignment(2, 6, 14, MIXED, seed=2)
    assert_identical(sequences, "focal_0|AU|Target", list(range(9)), ["Target"], case="mixed")


def test_every_contrast_state_is_a_gap():
    """All-zero scores still have to produce the same first-wins answer."""
    sequences = {
        "focal_0|AU|Target": "ACGTACGTAC",
        "focal_1|GB|Target": "ACGTACGTAC",
        "other_0|FR|Contrast": "----------",
        "other_1|DE|Contrast": "----------",
    }
    result = assert_identical(
        sequences, "focal_0|AU|Target", list(range(7)), ["Target"], case="all gaps"
    )
    assert result.best_gap_score == 0.0
    assert result.best_avg_score == 0.0


# ---------------------------------------------------------------------------
# Focal selection
# ---------------------------------------------------------------------------


def test_substring_selector():
    sequences = alignment(3, 5, 12, MIXED, seed=3)
    assert_identical(sequences, "focal_0|AU|Target", list(range(8)), ["Target"], case="substring")


def test_several_substring_selectors_are_ored():
    sequences = alignment(2, 4, 12, MIXED, seed=4)
    sequences["other_9|XX|Second"] = "ACGTACGTACGT"
    assert_identical(
        sequences,
        "focal_0|AU|Target",
        list(range(8)),
        ["Target", "Second"],
        case="two selectors",
    )


def test_exact_headers_selector():
    """
    Exact identity, not substring.

    `other_0|FR|Contrast_extra` must stay in the comparison set even though the
    chosen header is a prefix of it — which is the entire reason `ExactHeaders`
    exists.
    """
    sequences = alignment(2, 4, 12, MIXED, seed=5)
    sequences["focal_0|AU|Target_extra"] = "ACGTACGTACGT"
    selector = ExactHeaders(["focal_0|AU|Target", "focal_1|AU|Target"])
    assert_identical(
        sequences, "focal_0|AU|Target", list(range(8)), selector, case="exact headers"
    )


def test_a_non_focal_reference_is_still_excluded_from_scoring():
    """`header == ref_id` is a separate exclusion from the focal test."""
    sequences = alignment(2, 5, 12, MIXED, seed=6)
    assert_identical(
        sequences, "other_0|FR|Contrast", list(range(8)), ["Target"], case="non-focal ref"
    )


# ---------------------------------------------------------------------------
# diagnostic_states
# ---------------------------------------------------------------------------


def test_without_diagnostic_states_the_reference_sequence_supplies_them():
    sequences = alignment(2, 5, 12, MIXED, seed=7)
    assert_identical(
        sequences, "focal_0|AU|Target", list(range(8)), ["Target"], None, case="no states"
    )


def test_with_diagnostic_states_those_states_are_used_instead():
    sequences = alignment(2, 5, 12, MIXED, seed=8)
    sites = list(range(8))
    # Deliberately NOT the reference's own bases, so using the wrong source
    # produces different scores.
    states = {site: "ACGT"[site % 4] for site in sites}
    assert_identical(
        sequences, "focal_0|AU|Target", sites, ["Target"], states, case="states given"
    )


def test_the_two_sources_of_reference_states_disagree_as_expected():
    """A guard on the guard: the fixture must actually distinguish the paths."""
    sequences = {
        "focal_0|AU|Target": "AAAAAAAA",
        "focal_1|GB|Target": "AAAAAAAA",
        "other_0|FR|Contrast": "ACACACAC",
        "other_1|DE|Contrast": "AGAGAGAG",
    }
    sites = list(range(8))
    from_ref = find_best_five_site_sets(
        sequences=sequences, ref_id="focal_0|AU|Target", sites=sites, focal_strings=["Target"]
    )
    from_states = find_best_five_site_sets(
        sequences=sequences,
        ref_id="focal_0|AU|Target",
        sites=sites,
        focal_strings=["Target"],
        diagnostic_states={site: "C" for site in sites},
    )
    assert from_ref.best_avg_score != from_states.best_avg_score


# ---------------------------------------------------------------------------
# Search shape: order, ties, size
# ---------------------------------------------------------------------------


def test_fewer_than_five_sites_searches_nothing():
    sequences = alignment(2, 3, 12, STRICT, seed=9)
    result = find_best_five_site_sets(
        sequences=sequences,
        ref_id="focal_0|AU|Target",
        sites=[0, 1, 2, 3],
        focal_strings=["Target"],
    )
    assert result == FiveSiteResult(
        total_combinations_tested=0,
        best_gap_score=None,
        best_gap_sites=None,
        best_avg_score=None,
        best_avg_sites=None,
    )


def test_exactly_five_sites_tests_exactly_one_combination():
    sequences = alignment(2, 3, 12, MIXED, seed=10)
    result = assert_identical(
        sequences, "focal_0|AU|Target", [1, 3, 5, 7, 9], ["Target"], case="exactly five"
    )
    assert result.total_combinations_tested == 1
    assert result.best_gap_sites == (1, 3, 5, 7, 9)


def test_every_combination_is_counted():
    sequences = alignment(2, 3, 20, MIXED, seed=11)
    sites = list(range(10))
    result = assert_identical(
        sequences, "focal_0|AU|Target", sites, ["Target"], case="counted"
    )
    assert result.total_combinations_tested == 252  # C(10,5)


def test_a_total_tie_is_won_by_the_first_combination():
    """
    Strict `<` means the earliest combination keeps the prize.

    Every site here scores identically for every contrast sequence, so all 56
    combinations tie and the answer is forced to be the first one
    `itertools.combinations` yields.
    """
    sequences = {
        "focal_0|AU|Target": "AAAAAAAA",
        "focal_1|GB|Target": "AAAAAAAA",
        "other_0|FR|Contrast": "CCCCCCCC",
        "other_1|DE|Contrast": "AAAAAAAA",
    }
    sites = list(range(8))
    result = assert_identical(sequences, "focal_0|AU|Target", sites, ["Target"], case="tie")
    assert result.total_combinations_tested == 56
    assert result.best_gap_sites == (0, 1, 2, 3, 4)
    assert result.best_avg_sites == (0, 1, 2, 3, 4)


def test_site_order_follows_the_supplied_list_not_sorted_order():
    """
    `combinations` walks the list as given.

    The search must not quietly sort its input: the reported tuples, and the
    tie-break, follow the caller's order.
    """
    sequences = alignment(2, 4, 20, MIXED, seed=12)
    sites = [9, 2, 7, 0, 15, 4]
    result = assert_identical(
        sequences, "focal_0|AU|Target", sites, ["Target"], case="unsorted sites"
    )
    assert result.best_gap_sites[0] in {9, 2}


# ---------------------------------------------------------------------------
# Differential sweep
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("seed", range(12))
def test_randomised_alignments_match_the_reference_exactly(seed):
    """
    The broad net.

    Varies the alphabet, the shape of the alignment and whether diagnostic
    states are supplied, and insists on exact equality every time.
    """
    rng = random.Random(1000 + seed)
    n_focal = rng.randint(1, 4)
    n_other = rng.randint(1, 7)
    length = rng.randint(8, 16)
    alphabet = MIXED if seed % 2 else STRICT
    sequences = alignment(n_focal, n_other, length, alphabet, seed=seed)

    sites = sorted(rng.sample(range(length), rng.randint(5, min(8, length))))
    states = {site: rng.choice("ACGT") for site in sites} if seed % 3 == 0 else None

    assert_identical(
        sequences,
        f"focal_0|AU|Target",
        sites,
        ["Target"],
        states,
        case=f"seed {seed}",
    )


@pytest.mark.parametrize("seed", range(6))
def test_randomised_alignments_with_exact_headers_match(seed):
    rng = random.Random(2000 + seed)
    sequences = alignment(3, rng.randint(2, 6), 14, MIXED, seed=seed)
    selector = ExactHeaders([header for header in sequences if "Target" in header][:2])
    sites = sorted(rng.sample(range(14), 7))
    assert_identical(
        sequences,
        "focal_0|AU|Target",
        sites,
        selector,
        case=f"exact headers seed {seed}",
    )


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------


class _Recorder:
    """Minimal observer: records events, observes nothing else."""

    def __init__(self):
        self.events = []

    def event(self, name, /, **fields):
        self.events.append((name, fields))

    def progress(self, stage, /, **kwargs):
        pass

    def stage(self, name, /, **fields):
        from contextlib import nullcontext

        return nullcontext()


def test_five_site_start_counts_the_sequences_actually_scored():
    """
    `comparison_sequences` is the size of the comparison set.

    It used to report `len(sequences)` — every record in the alignment,
    including the focal group and the reference, which are precisely the ones
    never scored. On a real alignment that overstated the work by the size of
    the focal set.
    """
    sequences = {
        "focal_0|AU|Target": "ACGTACGTAC",
        "focal_1|GB|Target": "ACGTACGTAC",
        "focal_2|FR|Target": "ACGTACGTAC",
        "other_0|DE|Contrast": "CGTACGTACG",
        "other_1|ES|Contrast": "TACGTACGTA",
    }
    recorder = _Recorder()

    find_best_five_site_sets(
        sequences=sequences,
        ref_id="focal_0|AU|Target",
        sites=list(range(6)),
        focal_strings=["Target"],
        observer=recorder,
    )

    start = dict(recorder.events)["five_site.start"]
    assert start["comparison_sequences"] == 2  # not 5
    assert start["candidate_sites"] == 6
    assert start["combinations"] == 6  # C(6,5)


def test_a_non_focal_reference_is_not_counted_as_a_comparison_sequence():
    sequences = {
        "focal_0|AU|Target": "ACGTACGTAC",
        "other_0|DE|Contrast": "CGTACGTACG",
        "other_1|ES|Contrast": "TACGTACGTA",
    }
    recorder = _Recorder()

    find_best_five_site_sets(
        sequences=sequences,
        ref_id="other_0|DE|Contrast",
        sites=list(range(6)),
        focal_strings=["Target"],
        observer=recorder,
    )

    # Three records, one focal, one of them the reference: one left to score.
    assert dict(recorder.events)["five_site.start"]["comparison_sequences"] == 1


def test_no_comparison_sequences_is_refused():
    """The same refusal `compute_metrics` raises, with the same wording."""
    sequences = {
        "focal_0|AU|Target": "ACGTACGTAC",
        "focal_1|GB|Target": "ACGTACGTAC",
    }

    with pytest.raises(ValueError, match="No non-focal sequences available"):
        find_best_five_site_sets(
            sequences=sequences,
            ref_id="focal_0|AU|Target",
            sites=list(range(6)),
            focal_strings=["Target"],
        )


# ---------------------------------------------------------------------------
# The float rule the implementation depends on
# ---------------------------------------------------------------------------


def test_builtin_sum_is_not_interchangeable_with_chained_addition():
    """
    Why the search still calls `sum()` for five numbers.

    CPython 3.12's `sum` compensates (Neumaier) when adding floats, so it does
    NOT agree with `a + b + c + d + e` — and the scores here are exactly the
    values `score_state` can return, so this is not a theoretical concern. If
    this test ever fails because the two agree, the fast path may be simplified;
    until then, simplifying it silently moves scores by one ulp and can hand
    `best_avg_sites` to a different combination.
    """
    from itertools import product

    # 0.0 plus 1/k for every IUPAC set size k in {1, 2, 3, 4}.
    values = [0.0, 1.0, 0.5, 1.0 / 3, 0.25]
    differing = [
        tuple(t)
        for t in product(values, repeat=5)
        if (t[0] + t[1] + t[2] + t[3] + t[4]) != sum(t)
    ]

    assert differing, "chained addition and sum() agree here; see the docstring"
    # The generator form `compute_similarity` uses is the same as the sequence
    # form the search uses, which is what makes the fast path exact.
    assert all(sum(x for x in t) == sum(t) for t in product(values, repeat=5))
