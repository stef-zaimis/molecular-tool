from collections import defaultdict
from collections.abc import Mapping

from molecular_diagnosis.constants import (
    BALANCING_EMPTY_LIMIT,
    BALANCING_EMPTY_WEIGHT,
    BD_AMBIGUOUS_FRACTION_LIMIT,
    IUPAC,
    POLYMORPHISM_EMPTY_LIMIT,
    POLYMORPHISM_EMPTY_WEIGHT,
    PROLONGATION_WEIGHT,
    STRICT_BASES,
)
from molecular_diagnosis.models import PunishmentEvent, PunishmentResult


def state_possibilities(state: str) -> set[str]:
    return IUPAC.get(state.upper(), set())


def is_empty_state(state: str) -> bool:
    return not state_possibilities(state)


def empty_fraction(column: list[str]) -> float:
    if not column:
        return 0.0

    empty_count = sum(1 for state in column if is_empty_state(state))
    return empty_count / len(column)


def compatible_strict_base_for_bd(column: list[str]) -> str | None:
    """
    Benefit-of-doubt rule.

    BD applies only if:
      1. exactly one strict A/C/G/T base is present among non-empty states;
      2. every ambiguous non-empty state includes that strict base; and
      3. ambiguous non-empty states account for no more than the configured
         BD ambiguity fraction among non-empty states.

    Empty states are ignored for the ambiguity-fraction calculation because
    gaps/missing states are handled separately by empty-site punishment.
    """

    non_empty_states = [
        state.upper()
        for state in column
        if not is_empty_state(state)
    ]

    if not non_empty_states:
        return None

    strict_states = {
        state
        for state in non_empty_states
        if state in STRICT_BASES
    }

    if len(strict_states) != 1:
        return None

    strict_base = next(iter(strict_states))

    ambiguous_states = [
        state
        for state in non_empty_states
        if state not in STRICT_BASES
    ]

    ambiguous_fraction = len(ambiguous_states) / len(non_empty_states)

    if ambiguous_fraction > BD_AMBIGUOUS_FRACTION_LIMIT:
        return None

    for state in ambiguous_states:
        possibilities = state_possibilities(state)

        if strict_base not in possibilities:
            return None

    return strict_base

def weighted_counts_for_column(
    column: list[str],
    *,
    focal_index: int | None = None,
    assumed_base: str | None = None,
) -> dict[str, float]:
    """
    Counts bases fractionally.

    Example:
        R contributes 0.5 A + 0.5 G.
        N contributes 0.25 A + 0.25 C + 0.25 G + 0.25 T.

    When scoring an ambiguous sequence itself, focal_index and assumed_base are
    used to replace that sequence's ambiguous state with one full assumed base.
    """

    counts = {base: 0.0 for base in STRICT_BASES}

    for index, state in enumerate(column):
        possibilities = state_possibilities(state)

        if not possibilities:
            continue

        if focal_index is not None and index == focal_index:
            if assumed_base is None:
                raise ValueError("assumed_base is required when focal_index is used.")

            counts[assumed_base] += 1.0
            continue

        contribution = 1.0 / len(possibilities)

        for base in possibilities:
            counts[base] += contribution

    return counts


def score_sequence_as_base(
    column: list[str],
    sequence_index: int,
    assumed_base: str,
) -> float:
    """
    Scores one sequence at one column as if its state were assumed_base.

    Bdif / Bsame, where:
      Bsame includes the sequence being scored.
      Empty states in other sequences are excluded.
      Ambiguous states in other sequences are fractionally distributed.
    """

    counts = weighted_counts_for_column(
        column,
        focal_index=sequence_index,
        assumed_base=assumed_base,
    )

    total_non_empty = sum(counts.values())
    b_same = counts[assumed_base]

    if b_same == 0:
        return 0.0

    b_dif = total_non_empty - b_same
    return b_dif / b_same


def score_non_empty_state(
    column: list[str],
    sequence_index: int,
    state: str,
) -> float:
    """
    Scores an unambiguous or ambiguous non-empty state.

    Ambiguous focal states are scored as the average of their possible strict
    base interpretations.
    """

    possibilities = state_possibilities(state)

    if not possibilities:
        return 0.0

    scores = [
        score_sequence_as_base(
            column=column,
            sequence_index=sequence_index,
            assumed_base=base,
        )
        for base in sorted(possibilities)
    ]

    return sum(scores) / len(scores)


def add_score_event(
    *,
    events: list[PunishmentEvent],
    headers: list[str],
    scores: dict[str, float],
    category_scores: dict[str, float],
    sequence_index: int,
    site: int,
    state: str,
    category: str,
    score: float,
    empty_frac: float,
    direction: str | None = None,
    note: str = "",
) -> None:
    sequence_id = headers[sequence_index]

    scores[sequence_id] += score
    category_scores[sequence_id] += score

    events.append(
        PunishmentEvent(
            sequence_id=sequence_id,
            site=site,
            state=state,
            category=category,
            score=score,
            empty_fraction=empty_frac,
            direction=direction,
            note=note,
        )
    )


def process_polymorphism_or_balancing_column(
    *,
    site: int,
    column: list[str],
    headers: list[str],
    category: str,
    empty_weight: float,
    total_scores: dict[str, float],
    category_scores: dict[str, float],
    bd_counts: dict[str, int],
    events: list[PunishmentEvent],
) -> None:
    empty_frac = empty_fraction(column)

    b_empty = sum(1 for state in column if is_empty_state(state))
    b_any = len(column) - b_empty

    empty_score = 0.0

    if b_empty > 0:
        empty_score = (b_any / b_empty) * empty_weight

    bd_base = compatible_strict_base_for_bd(column)

    if bd_base is not None:
        for sequence_index, state in enumerate(column):
            possibilities = state_possibilities(state)

            if not possibilities:
                if empty_score > 0:
                    add_score_event(
                        events=events,
                        headers=headers,
                        scores=total_scores,
                        category_scores=category_scores,
                        sequence_index=sequence_index,
                        site=site,
                        state=state,
                        category=category,
                        score=empty_score,
                        empty_frac=empty_frac,
                        note="Empty state in otherwise benefit-of-doubt-compatible column.",
                    )

                continue

            if state.upper() not in STRICT_BASES:
                sequence_id = headers[sequence_index]
                bd_counts[sequence_id] += 1

                events.append(
                    PunishmentEvent(
                        sequence_id=sequence_id,
                        site=site,
                        state=state,
                        category="BD",
                        score=0.0,
                        empty_fraction=empty_frac,
                        note=f"Ambiguous state treated as compatible with {bd_base}.",
                    )
                )

        return

    for sequence_index, state in enumerate(column):
        if is_empty_state(state):
            if empty_score > 0:
                add_score_event(
                    events=events,
                    headers=headers,
                    scores=total_scores,
                    category_scores=category_scores,
                    sequence_index=sequence_index,
                    site=site,
                    state=state,
                    category=category,
                    score=empty_score,
                    empty_frac=empty_frac,
                    note="Empty state.",
                )

            continue

        score = score_non_empty_state(
            column=column,
            sequence_index=sequence_index,
            state=state,
        )

        if score > 0:
            add_score_event(
                events=events,
                headers=headers,
                scores=total_scores,
                category_scores=category_scores,
                sequence_index=sequence_index,
                site=site,
                state=state,
                category=category,
                score=score,
                empty_frac=empty_frac,
                note="Base-state disagreement within focal group.",
            )


def process_gap_dominated_column(
    *,
    site: int,
    column: list[str],
    headers: list[str],
    category: str,
    direction: str,
    total_scores: dict[str, float],
    prolongation_scores: dict[str, float],
    insertion_scores: dict[str, float],
    prl_counts: dict[str, int],
    ins_counts: dict[str, int],
    events: list[PunishmentEvent],
) -> None:
    empty_frac = empty_fraction(column)

    b_empty = sum(1 for state in column if is_empty_state(state))
    b_any = len(column) - b_empty

    if b_any == 0:
        return

    if category == "PRL":
        score = (b_empty / b_any) * PROLONGATION_WEIGHT
        category_scores = prolongation_scores
        point_counts = prl_counts
        note = "Non-empty state in terminal mostly-empty region."
    elif category == "INS":
        score = float(b_any)
        category_scores = insertion_scores
        point_counts = ins_counts
        note = "Non-empty state in internal mostly-empty region."
    else:
        raise ValueError(f"Unsupported gap-dominated category: {category}")

    for sequence_index, state in enumerate(column):
        if is_empty_state(state):
            continue

        sequence_id = headers[sequence_index]
        point_counts[sequence_id] += 1

        add_score_event(
            events=events,
            headers=headers,
            scores=total_scores,
            category_scores=category_scores,
            sequence_index=sequence_index,
            site=site,
            state=state,
            category=category,
            score=score,
            empty_frac=empty_frac,
            direction=direction,
            note=note,
        )


def find_focal_punishments(
    sequences: Mapping[str, str],
    focal_headers: list[str],
) -> PunishmentResult:
    """
    Computes focal punishment/anomaly scores.

    This only uses the focal group. It does not use non-focal sequences and does
    not alter the DMC analysis.
    """

    if not focal_headers:
        raise ValueError("No focal headers supplied for punishment analysis.")

    headers = list(focal_headers)
    focal_sequences = [sequences[header] for header in headers]

    alignment_length = len(focal_sequences[0])

    total_scores: dict[str, float] = defaultdict(float)
    polymorphism_scores: dict[str, float] = defaultdict(float)
    balancing_scores: dict[str, float] = defaultdict(float)
    prolongation_scores: dict[str, float] = defaultdict(float)
    insertion_scores: dict[str, float] = defaultdict(float)

    bd_counts: dict[str, int] = defaultdict(int)
    prl_counts: dict[str, int] = defaultdict(int)
    ins_counts: dict[str, int] = defaultdict(int)

    events: list[PunishmentEvent] = []

    for header in headers:
        total_scores[header] = 0.0
        polymorphism_scores[header] = 0.0
        balancing_scores[header] = 0.0
        prolongation_scores[header] = 0.0
        insertion_scores[header] = 0.0
        bd_counts[header] = 0
        prl_counts[header] = 0
        ins_counts[header] = 0

    def column_at(site: int) -> list[str]:
        return [seq[site] for seq in focal_sequences]

    def process_scan_positions(positions: range, direction: str) -> None:
        terminal_gap_run_open = True

        for site in positions:
            column = column_at(site)
            empty_frac = empty_fraction(column)

            if empty_frac >= BALANCING_EMPTY_LIMIT:
                category = "PRL" if terminal_gap_run_open else "INS"

                process_gap_dominated_column(
                    site=site,
                    column=column,
                    headers=headers,
                    category=category,
                    direction=direction,
                    total_scores=total_scores,
                    prolongation_scores=prolongation_scores,
                    insertion_scores=insertion_scores,
                    prl_counts=prl_counts,
                    ins_counts=ins_counts,
                    events=events,
                )

                continue

            terminal_gap_run_open = False

            if empty_frac < POLYMORPHISM_EMPTY_LIMIT:
                process_polymorphism_or_balancing_column(
                    site=site,
                    column=column,
                    headers=headers,
                    category="POLY",
                    empty_weight=POLYMORPHISM_EMPTY_WEIGHT,
                    total_scores=total_scores,
                    category_scores=polymorphism_scores,
                    bd_counts=bd_counts,
                    events=events,
                )
            else:
                process_polymorphism_or_balancing_column(
                    site=site,
                    column=column,
                    headers=headers,
                    category="BAL",
                    empty_weight=BALANCING_EMPTY_WEIGHT,
                    total_scores=total_scores,
                    category_scores=balancing_scores,
                    bd_counts=bd_counts,
                    events=events,
                )

    left_stop = (alignment_length + 1) // 2

    left_positions = range(0, left_stop)
    right_positions = range(alignment_length - 1, left_stop - 1, -1)

    process_scan_positions(left_positions, direction="left-to-right")
    process_scan_positions(right_positions, direction="right-to-left")

    return PunishmentResult(
        total_scores=dict(total_scores),
        polymorphism_scores=dict(polymorphism_scores),
        balancing_scores=dict(balancing_scores),
        prolongation_scores=dict(prolongation_scores),
        insertion_scores=dict(insertion_scores),
        bd_counts=dict(bd_counts),
        prl_counts=dict(prl_counts),
        ins_counts=dict(ins_counts),
        events=events,
    )