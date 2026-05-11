from itertools import combinations

from molecular_diagnosis.constants import (
    BALANCING_EMPTY_LIMIT,
    IUPAC,
    STRICT_BASES,
)
from molecular_diagnosis.models import DMCResult, FiveSiteResult


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
    target_string: str,
    diagnostic_states: dict[int, str] | None = None,
) -> tuple[float, float, float]:
    if diagnostic_states is None:
        ref_states = extract_sites(sequences[ref_id], sites)
    else:
        ref_states = [diagnostic_states[site] for site in sites]

    similarities = []

    for header, seq in sequences.items():
        if header == ref_id or target_string in header:
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
    target_string: str,
    *,
    include_ambiguous_dmc_bd: bool = False,
    include_gappy_consensus_dmc_sites: bool = False,
) -> DMCResult:
    all_headers = list(sequences.keys())
    all_seqs = list(sequences.values())

    focal_headers = [header for header in all_headers if target_string in header]
    focal = [sequences[header] for header in focal_headers]
    non_focal = [seq for header, seq in sequences.items() if target_string not in header]

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

    def sequence_matches_site(seq: str, site: int) -> bool:
        return state_matches_dmc_base(
            state=seq[site],
            base=candidates[site],
            include_ambiguous_dmc_bd=include_ambiguous_dmc_bd,
        )

    singles = [
        site
        for site in candidate_sites
        if all(not sequence_matches_site(seq, site) for seq in non_focal)
    ]

    pairs: list[tuple[int, int]] = []
    pairs_tested = 0

    for site_a_index in range(len(candidate_sites)):
        for site_b_index in range(site_a_index + 1, len(candidate_sites)):
            pairs_tested += 1

            site_a = candidate_sites[site_a_index]
            site_b = candidate_sites[site_b_index]

            is_diagnostic_pair = all(
                not (
                    sequence_matches_site(seq, site_a)
                    and sequence_matches_site(seq, site_b)
                )
                for seq in non_focal
            )

            if is_diagnostic_pair:
                pairs.append((site_a, site_b))

    unique = sorted({site for pair in pairs for site in pair})

    ambiguous_bd_sites_included = sorted(
        site
        for site in candidate_sites
        if column_has_ambiguous_state([seq[site] for seq in all_seqs])
    )

    gappy_consensus_sites_included = sorted(gappy_consensus_flags)

    return DMCResult(
        single=singles,
        pairs=pairs,
        unique=unique,
        states=dict(candidates),
        fixed_count=len(fixed),
        skipped_non_acgt=skipped_non_acgt,
        globally_conserved_removed=globally_conserved_removed,
        candidate_count=len(candidates),
        pairs_tested=pairs_tested,
        ambiguous_bd_sites_included=ambiguous_bd_sites_included,
        gappy_consensus_sites_included=gappy_consensus_sites_included,
        include_ambiguous_dmc_bd=include_ambiguous_dmc_bd,
        include_gappy_consensus_dmc_sites=include_gappy_consensus_dmc_sites,
    )


def find_best_five_site_sets(
    sequences: dict[str, str],
    ref_id: str,
    sites: list[int],
    target_string: str,
    diagnostic_states: dict[int, str] | None = None,
) -> FiveSiteResult:
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

        max_similarity, avg_similarity, similarity_gap = compute_metrics(
            sequences=sequences,
            ref_id=ref_id,
            sites=combo,
            target_string=target_string,
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