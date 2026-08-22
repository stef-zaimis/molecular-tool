from itertools import combinations
from math import comb

from molecular_diagnosis.constants import (
    BALANCING_EMPTY_LIMIT,
    IUPAC,
    STRICT_BASES,
)
from molecular_diagnosis.focal import (
    FocalSelector,
    header_matches_focal,
    normalise_focal_strings,
)
from molecular_diagnosis.models import DMCResult, FiveSiteResult
from molecular_diagnosis.progress import (
    NULL_OBSERVER,
    STAGE_DMC_SEARCH,
    STAGE_FIVE_SITE,
    ProgressTicker,
    RunObserver,
)


def state_possibilities(state: str) -> set[str]:
    return IUPAC.get(state.upper(), set())


def is_empty_state(state: str) -> bool:
    return not state_possibilities(state)


def is_ambiguous_state(state: str) -> bool:
    possibilities = state_possibilities(state)
    return bool(possibilities) and state.upper() not in STRICT_BASES


def column_empty_fraction(column: list[str]) -> float:
    if not column:
        return 0.0

    return sum(1 for state in column if is_empty_state(state)) / len(column)


def column_has_empty_state(column: list[str]) -> bool:
    return any(is_empty_state(state) for state in column)


def column_has_ambiguous_state(column: list[str]) -> bool:
    return any(is_ambiguous_state(state) for state in column)


def state_matches_dmc_base(
    state: str,
    base: str,
    *,
    include_ambiguous_dmc_bd: bool,
) -> bool:
    possibilities = state_possibilities(state)

    if not possibilities:
        return False

    if include_ambiguous_dmc_bd:
        return base in possibilities

    return state.upper() == base


def consensus_base_for_column(
    column: list[str],
    *,
    include_ambiguous_dmc_bd: bool,
    allow_gaps: bool,
) -> str | None:
    """
    Returns a single strict A/C/G/T consensus base for a focal column.

    Rules:
      - Empty states are allowed only when allow_gaps=True.
      - Ambiguous non-empty states are allowed only when
        include_ambiguous_dmc_bd=True.
      - With ambiguity enabled, all non-empty states must share exactly one
        possible strict base.
    """

    possibilities_by_state: list[set[str]] = []

    for state in column:
        state = state.upper()
        possibilities = state_possibilities(state)

        if not possibilities:
            if allow_gaps:
                continue
            return None

        if not include_ambiguous_dmc_bd and state not in STRICT_BASES:
            return None

        possibilities_by_state.append(possibilities)

    if not possibilities_by_state:
        return None

    if include_ambiguous_dmc_bd:
        common = set.intersection(*possibilities_by_state)

        if len(common) == 1:
            return sorted(common)[0]

        return None

    strict_bases = {next(iter(possibilities)) for possibilities in possibilities_by_state}

    if len(strict_bases) == 1:
        return next(iter(strict_bases))

    return None


def column_has_non_empty_variation_from_base(
    column: list[str],
    base: str,
    *,
    include_ambiguous_dmc_bd: bool,
) -> bool:
    for state in column:
        if is_empty_state(state):
            continue

        if not state_matches_dmc_base(
            state=state,
            base=base,
            include_ambiguous_dmc_bd=include_ambiguous_dmc_bd,
        ):
            return True

    return False


def extract_sites(seq: str, sites: list[int] | tuple[int, ...]) -> list[str]:
    return [seq[i] for i in sites]


def score_state(ref: str, query: str) -> float:
    possibilities = IUPAC.get(query, set())

    if not possibilities:
        return 0.0

    if ref in possibilities:
        return 1.0 / len(possibilities)

    return 0.0


def compute_similarity(ref_states: list[str], query_states: list[str]) -> float:
    if not ref_states:
        return 0.0

    total = sum(score_state(ref, query) for ref, query in zip(ref_states, query_states))
    return total / len(ref_states)


def compute_match_score(ref_states: list[str], query_states: list[str]) -> float:
    return sum(score_state(ref, query) for ref, query in zip(ref_states, query_states))


def compute_metrics(
    sequences: dict[str, str],
    ref_id: str,
    sites: list[int] | tuple[int, ...],
    focal_strings: FocalSelector,
    diagnostic_states: dict[int, str] | None = None,
) -> tuple[float, float, float]:
    selectors = normalise_focal_strings(focal_strings)

    if diagnostic_states is None:
        ref_states = extract_sites(sequences[ref_id], sites)
    else:
        ref_states = [diagnostic_states[site] for site in sites]

    similarities = []

    for header, seq in sequences.items():
        if header == ref_id or header_matches_focal(header, selectors):
            continue

        query_states = extract_sites(seq, sites)
        similarity = compute_similarity(ref_states, query_states)
        similarities.append(similarity)

    if not similarities:
        raise ValueError("No non-focal sequences available for similarity comparison.")

    max_similarity = max(similarities)
    avg_similarity = sum(similarities) / len(similarities)
    similarity_gap = 1 - max_similarity

    return max_similarity, avg_similarity, similarity_gap


def find_dmc_information(
    sequences: dict[str, str],
    focal_strings: FocalSelector,
    *,
    include_ambiguous_dmc_bd: bool = False,
    include_gappy_consensus_dmc_sites: bool = False,
    min_combination_length: int = 1,
    max_combination_length: int = 2,
    start_combination_length: int = 1,
    initial_diagnostic_combinations: list[tuple[int, ...]] | None = None,
    initial_combinations_tested_by_length: dict[int, int] | None = None,
    observer: RunObserver = NULL_OBSERVER,
) -> DMCResult:
    if min_combination_length < 1:
        raise ValueError("Minimum combination length must be at least 1.")

    if max_combination_length < 1:
        raise ValueError("Maximum combination length must be at least 1.")

    if start_combination_length < 1:
        raise ValueError("Start combination length must be at least 1.")

    if min_combination_length > max_combination_length:
        raise ValueError("Minimum combination length cannot exceed maximum combination length.")

    selectors = normalise_focal_strings(focal_strings)

    all_headers = list(sequences.keys())
    all_seqs = list(sequences.values())

    focal_headers = [
        header for header in all_headers if header_matches_focal(header, selectors)
    ]
    focal = [sequences[header] for header in focal_headers]
    non_focal = [
        seq
        for header, seq in sequences.items()
        if not header_matches_focal(header, selectors)
    ]

    if not focal:
        raise ValueError("No focal sequences found.")

    if not non_focal:
        raise ValueError("No non-focal sequences found.")

    length = len(focal[0])

    fixed: dict[int, str] = {}
    skipped_non_acgt = 0

    for site in range(length):
        focal_column = [seq[site] for seq in focal]
        full_column = [seq[site] for seq in all_seqs]

        if column_has_ambiguous_state(full_column) and not include_ambiguous_dmc_bd:
            skipped_non_acgt += 1
            continue

        focal_gap_dominated = column_empty_fraction(focal_column) >= BALANCING_EMPTY_LIMIT

        allow_gaps = (
            include_gappy_consensus_dmc_sites
            and not focal_gap_dominated
        )

        consensus_base = consensus_base_for_column(
            focal_column,
            include_ambiguous_dmc_bd=include_ambiguous_dmc_bd,
            allow_gaps=allow_gaps,
        )

        if consensus_base is None:
            skipped_non_acgt += 1
            continue

        fixed[site] = consensus_base

    candidates: dict[int, str] = {}
    globally_conserved_removed = 0
    gappy_consensus_flags: dict[int, bool] = {}

    for site, base in fixed.items():
        full_column = [seq[site] for seq in all_seqs]
        focal_column = [seq[site] for seq in focal]

        focal_gap_dominated = column_empty_fraction(focal_column) >= BALANCING_EMPTY_LIMIT

        has_non_empty_variation = column_has_non_empty_variation_from_base(
            full_column,
            base,
            include_ambiguous_dmc_bd=include_ambiguous_dmc_bd,
        )

        has_gaps = column_has_empty_state(full_column)

        gappy_consensus_signal = (
            include_gappy_consensus_dmc_sites
            and has_gaps
            and not focal_gap_dominated
            and not has_non_empty_variation
        )

        if has_non_empty_variation or gappy_consensus_signal:
            candidates[site] = base

            if gappy_consensus_signal:
                gappy_consensus_flags[site] = True
        else:
            globally_conserved_removed += 1

    candidate_sites = list(candidates.keys())
    candidate_site_set = set(candidate_sites)

    def sequence_matches_site(seq: str, site: int) -> bool:
        return state_matches_dmc_base(
            state=seq[site],
            base=candidates[site],
            include_ambiguous_dmc_bd=include_ambiguous_dmc_bd,
        )

    def sequence_matches_combination(seq: str, combo: tuple[int, ...]) -> bool:
        return all(sequence_matches_site(seq, site) for site in combo)

    def combination_is_diagnostic(combo: tuple[int, ...]) -> bool:
        return all(
            not sequence_matches_combination(seq, combo)
            for seq in non_focal
        )

    def normalise_combo(combo: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(sorted(set(combo)))

    diagnostic_combinations: list[tuple[int, ...]] = []
    seen_diagnostic_combinations: set[tuple[int, ...]] = set()

    for combo in initial_diagnostic_combinations or []:
        normalised = normalise_combo(combo)

        if not normalised:
            continue

        if not set(normalised).issubset(candidate_site_set):
            continue

        if normalised in seen_diagnostic_combinations:
            continue

        diagnostic_combinations.append(normalised)
        seen_diagnostic_combinations.add(normalised)

    diagnostic_combinations.sort(key=lambda combo: (len(combo), combo))

    combinations_tested_by_length: dict[int, int] = {}

    if initial_combinations_tested_by_length is not None:
        combinations_tested_by_length.update(initial_combinations_tested_by_length)

    def is_pruned_by_existing_dmc(combo: tuple[int, ...]) -> bool:
        combo_set = set(combo)

        return any(
            set(found_combo).issubset(combo_set)
            for found_combo in diagnostic_combinations
            if len(found_combo) < len(combo)
        )

    stop_reason = "reached_maximum_length"
    stopped_at_length = max_combination_length

    if not candidate_sites:
        stop_reason = "no_candidate_sites"
        stopped_at_length = 0

    elif start_combination_length > max_combination_length:
        stop_reason = "start_length_exceeds_maximum_length"
        stopped_at_length = start_combination_length - 1

    else:
        for combo_length in range(start_combination_length, max_combination_length + 1):
            combinations_tested_by_length.setdefault(combo_length, 0)
            found_this_length: list[tuple[int, ...]] = []

            # C(n, k) for this size: what the loop below is walking through.
            # Reported up front so a stage that will take an hour says so at
            # the start rather than after it.
            size_total = comb(len(candidate_sites), combo_length)
            observer.event(
                "dmc.size.start",
                size=combo_length,
                candidate_sites=len(candidate_sites),
                combinations=size_total,
            )
            ticker = ProgressTicker(
                observer,
                STAGE_DMC_SEARCH,
                total=size_total,
                detail=str(combo_length),
            )

            for combo in combinations(candidate_sites, combo_length):
                ticker.advance()

                if is_pruned_by_existing_dmc(combo):
                    continue

                combinations_tested_by_length[combo_length] += 1

                if combination_is_diagnostic(combo):
                    found_this_length.append(combo)

            ticker.finish()
            observer.event(
                "dmc.size.end",
                size=combo_length,
                tested=combinations_tested_by_length[combo_length],
                found=len(found_this_length),
            )

            for combo in found_this_length:
                if combo not in seen_diagnostic_combinations:
                    diagnostic_combinations.append(combo)
                    seen_diagnostic_combinations.add(combo)

            diagnostic_combinations.sort(key=lambda found: (len(found), found))

            stopped_at_length = combo_length

            if found_this_length and combo_length >= min_combination_length:
                stop_reason = "found_at_or_above_minimum_length"
                break

    diagnostic_combinations_by_length: dict[int, list[tuple[int, ...]]] = {}

    for combo in diagnostic_combinations:
        diagnostic_combinations_by_length.setdefault(len(combo), []).append(combo)

    singles = [
        combo[0]
        for combo in diagnostic_combinations_by_length.get(1, [])
    ]

    pairs = [
        combo
        for combo in diagnostic_combinations_by_length.get(2, [])
    ]

    unique = sorted({
        site
        for combo in diagnostic_combinations
        for site in combo
    })

    ambiguous_bd_sites_included = sorted(
        site
        for site in candidate_sites
        if column_has_ambiguous_state([seq[site] for seq in all_seqs])
    )

    gappy_consensus_sites_included = sorted(gappy_consensus_flags)

    total_combinations_tested = sum(combinations_tested_by_length.values())

    return DMCResult(
        single=singles,
        pairs=pairs,
        unique=unique,
        states=dict(candidates),
        diagnostic_combinations=diagnostic_combinations,
        diagnostic_combinations_by_length=diagnostic_combinations_by_length,
        min_combination_length=min_combination_length,
        max_combination_length=max_combination_length,
        start_combination_length=start_combination_length,
        stopped_at_length=stopped_at_length,
        stop_reason=stop_reason,
        combinations_tested_by_length=dict(combinations_tested_by_length),
        total_combinations_tested=total_combinations_tested,
        fixed_count=len(fixed),
        skipped_non_acgt=skipped_non_acgt,
        globally_conserved_removed=globally_conserved_removed,
        candidate_count=len(candidates),
        pairs_tested=combinations_tested_by_length.get(2, 0),
        ambiguous_bd_sites_included=ambiguous_bd_sites_included,
        gappy_consensus_sites_included=gappy_consensus_sites_included,
        include_ambiguous_dmc_bd=include_ambiguous_dmc_bd,
        include_gappy_consensus_dmc_sites=include_gappy_consensus_dmc_sites,
    )


def find_best_five_site_sets(
    sequences: dict[str, str],
    ref_id: str,
    sites: list[int],
    focal_strings: FocalSelector,
    diagnostic_states: dict[int, str] | None = None,
    observer: RunObserver = NULL_OBSERVER,
) -> FiveSiteResult:
    """
    Exhaustive best-of-`C(n, 5)` search.

    THE combinatorial cliff in this codebase, and the reason this function
    reports progress: the loop is `C(n, 5)` iterations and each one compares
    the reference against every non-focal sequence, so the cost is
    `C(n, 5) x non-focal count x 5`. Twenty candidate sites is 15,504
    combinations; forty is 658,008; sixty is 5,461,512. Against a couple of
    thousand contrast sequences the last of those is not a hang, but it is
    hours, and without progress it is indistinguishable from one.

    Nothing about the search is changed by observing it.
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

    total = comb(len(sites), 5)
    observer.event(
        "five_site.start",
        candidate_sites=len(sites),
        combinations=total,
        comparison_sequences=len(sequences),
    )
    ticker = ProgressTicker(observer, STAGE_FIVE_SITE, total=total, check_every=256)

    for combo in combinations(sites, 5):
        total_combinations_tested += 1
        ticker.advance()

        max_similarity, avg_similarity, similarity_gap = compute_metrics(
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

    ticker.finish()
    observer.event("five_site.end", tested=total_combinations_tested)

    return FiveSiteResult(
        total_combinations_tested=total_combinations_tested,
        best_gap_score=best_gap_score,
        best_gap_sites=best_gap_sites,
        best_avg_score=best_avg_score,
        best_avg_sites=best_avg_sites,
    )


def format_diag(sites: list[int] | tuple[int, ...], seq: str) -> str:
    if not sites:
        return "None"

    return ", ".join(f"{i + 1}:{seq[i]}" for i in sites)


def format_diag_from_states(
    sites: list[int] | tuple[int, ...],
    states: dict[int, str],
) -> str:
    if not sites:
        return "None"

    return ", ".join(f"{site + 1}:{states[site]}" for site in sites)