from itertools import combinations

from molecular_diagnosis.constants import IUPAC, STRICT_BASES
from molecular_diagnosis.models import DMCResult, FiveSiteResult


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
) -> tuple[float, float, float]:
    ref_states = extract_sites(sequences[ref_id], sites)

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
) -> DMCResult:
    all_seqs = list(sequences.values())
    focal = [seq for header, seq in sequences.items() if target_string in header]
    non_focal = [seq for header, seq in sequences.items() if target_string not in header]

    if not focal:
        raise ValueError("No focal sequences found.")

    if not non_focal:
        raise ValueError("No non-focal sequences found.")

    length = len(focal[0])

    fixed: dict[int, str] = {}
    skipped_non_acgt = 0

    for i in range(length):
        focal_column = [seq[i] for seq in focal]

        if any(base not in STRICT_BASES for base in focal_column):
            skipped_non_acgt += 1
            continue

        if len(set(focal_column)) == 1:
            fixed[i] = focal_column[0]

    candidates: dict[int, str] = {}
    globally_conserved_removed = 0

    for i, base in fixed.items():
        acgt_column = [seq[i] for seq in all_seqs if seq[i] in STRICT_BASES]

        if len(set(acgt_column)) > 1:
            candidates[i] = base
        else:
            globally_conserved_removed += 1

    candidate_sites = list(candidates.keys())

    singles = [
        i
        for i in candidate_sites
        if all(seq[i] != candidates[i] for seq in non_focal)
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
                    seq[site_a] == candidates[site_a]
                    and seq[site_b] == candidates[site_b]
                )
                for seq in non_focal
            )

            if is_diagnostic_pair:
                pairs.append((site_a, site_b))

    unique = sorted({site for pair in pairs for site in pair})

    return DMCResult(
        single=singles,
        pairs=pairs,
        unique=unique,
        fixed_count=len(fixed),
        skipped_non_acgt=skipped_non_acgt,
        globally_conserved_removed=globally_conserved_removed,
        candidate_count=len(candidates),
        pairs_tested=pairs_tested,
    )


def find_best_five_site_sets(
    sequences: dict[str, str],
    ref_id: str,
    sites: list[int],
    target_string: str,
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