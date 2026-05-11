from dataclasses import dataclass
from pathlib import Path
from collections.abc import Sequence

from molecular_diagnosis.constants import (
    BALANCING_EMPTY_LIMIT,
    IUPAC,
    IUPAC_FROM_BASES,
    STRICT_BASES,
)


@dataclass(frozen=True)
class ConsensusResult:
    untrimmed_sequence: str
    trimmed_sequence: str
    kept_indices: tuple[int, ...]
    removed_prl_indices: tuple[int, ...]
    removed_ins_indices: tuple[int, ...]
    untrimmed_n_count: int
    trimmed_n_count: int


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(
        sequence[index:index + width]
        for index in range(0, len(sequence), width)
    )


def is_empty_state(state: str) -> bool:
    return not IUPAC.get(state.upper(), set())


def empty_fraction(column: Sequence[str]) -> float:
    if not column:
        return 0.0

    return sum(1 for state in column if is_empty_state(state)) / len(column)


def assign_consensus_state(non_gap_column: Sequence[str]) -> str:
    """
    Consensus rule copied from the older script.

    - Ignore gaps/missing states before calling this function.
    - Count only strict A/C/G/T states.
    - If one base dominates the second-most common base by > 0.5, use it.
    - Otherwise use the exact IUPAC code for the observed strict base set.
    - If no strict A/C/G/T bases are available, return N.
    """
    counts: dict[str, int] = {}

    for state in non_gap_column:
        state = state.upper()

        if state in STRICT_BASES:
            counts[state] = counts.get(state, 0) + 1

    if not counts:
        return "N"

    total = sum(counts.values())
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))

    if len(ranked) == 1:
        return ranked[0][0]

    base_1, count_1 = ranked[0]
    _base_2, count_2 = ranked[1]

    frequency_1 = count_1 / total
    frequency_2 = count_2 / total

    if frequency_1 - frequency_2 > 0.5:
        return base_1

    base_set = frozenset(counts.keys())
    return IUPAC_FROM_BASES.get(base_set, "N")


def build_untrimmed_consensus(
    focal_sequences: Sequence[str],
) -> tuple[str, int]:
    if not focal_sequences:
        raise ValueError("Cannot build consensus without focal sequences.")

    alignment_length = len(focal_sequences[0])
    consensus: list[str] = []
    n_count = 0

    for site in range(alignment_length):
        column = [sequence[site] for sequence in focal_sequences]
        non_gap_column = [
            state
            for state in column
            if not is_empty_state(state)
        ]

        state = assign_consensus_state(non_gap_column)

        if state == "N":
            n_count += 1

        consensus.append(state)

    return "".join(consensus), n_count


def classify_prl_ins_sites(
    focal_sequences: Sequence[str],
) -> tuple[set[int], set[int]]:
    """
    Classify focal columns using the same scan idea as punishment scoring.

    A column is PRL/INS-like if:
    empty_fraction >= BALANCING_EMPTY_LIMIT

    Terminal gap-dominated runs from each end are PRL.
    Gap-dominated columns after the terminal run has closed are INS.
    """
    if not focal_sequences:
        raise ValueError("Cannot classify PRL/INS sites without focal sequences.")

    alignment_length = len(focal_sequences[0])

    removed_prl: set[int] = set()
    removed_ins: set[int] = set()

    def column_at(site: int) -> list[str]:
        return [sequence[site] for sequence in focal_sequences]

    def process_scan_positions(positions: range) -> None:
        terminal_gap_run_open = True

        for site in positions:
            column = column_at(site)
            empty_frac = empty_fraction(column)

            if empty_frac >= BALANCING_EMPTY_LIMIT:
                if terminal_gap_run_open:
                    removed_prl.add(site)
                else:
                    removed_ins.add(site)

                continue

            terminal_gap_run_open = False

    left_stop = (alignment_length + 1) // 2
    left_positions = range(0, left_stop)
    right_positions = range(alignment_length - 1, left_stop - 1, -1)

    process_scan_positions(left_positions)
    process_scan_positions(right_positions)

    return removed_prl, removed_ins


def build_trimmed_consensus(
    focal_sequences: Sequence[str],
    removed_prl: set[int],
    removed_ins: set[int],
) -> tuple[str, tuple[int, ...], int]:
    if not focal_sequences:
        raise ValueError("Cannot build consensus without focal sequences.")

    alignment_length = len(focal_sequences[0])
    removed = removed_prl | removed_ins

    consensus: list[str] = []
    kept_indices: list[int] = []
    n_count = 0

    for site in range(alignment_length):
        if site in removed:
            continue

        column = [sequence[site] for sequence in focal_sequences]
        non_gap_column = [
            state
            for state in column
            if not is_empty_state(state)
        ]

        state = assign_consensus_state(non_gap_column)

        if state == "N":
            n_count += 1

        consensus.append(state)
        kept_indices.append(site)

    return "".join(consensus), tuple(kept_indices), n_count


def build_focal_consensus_result(
    focal_sequences: Sequence[str],
) -> ConsensusResult:
    if not focal_sequences:
        raise ValueError("Cannot build focal consensus without focal sequences.")

    lengths = {len(sequence) for sequence in focal_sequences}

    if len(lengths) != 1:
        raise ValueError("Focal sequences are not all the same aligned length.")

    untrimmed_sequence, untrimmed_n_count = build_untrimmed_consensus(
        focal_sequences=focal_sequences,
    )

    removed_prl, removed_ins = classify_prl_ins_sites(
        focal_sequences=focal_sequences,
    )

    trimmed_sequence, kept_indices, trimmed_n_count = build_trimmed_consensus(
        focal_sequences=focal_sequences,
        removed_prl=removed_prl,
        removed_ins=removed_ins,
    )

    return ConsensusResult(
        untrimmed_sequence=untrimmed_sequence,
        trimmed_sequence=trimmed_sequence,
        kept_indices=kept_indices,
        removed_prl_indices=tuple(sorted(removed_prl)),
        removed_ins_indices=tuple(sorted(removed_ins)),
        untrimmed_n_count=untrimmed_n_count,
        trimmed_n_count=trimmed_n_count,
    )


def write_consensus_text_report(
    output_path: str | Path,
    *,
    target_string: str,
    focal_headers: Sequence[str],
    alignment_length: int,
    consensus_result: ConsensusResult,
    dmc_sites: Sequence[int] | None = None,
) -> None:
    output_path = Path(output_path)

    kept_index_to_trimmed_position = {
        full_index: trimmed_index + 1
        for trimmed_index, full_index in enumerate(consensus_result.kept_indices)
    }

    removed_prl = set(consensus_result.removed_prl_indices)
    removed_ins = set(consensus_result.removed_ins_indices)

    with open(output_path, "w", encoding="utf-8") as file:
        file.write("Focal consensus report\n\n")

        file.write("Input settings:\n")
        file.write(f"Focal identifier string: {target_string}\n")
        file.write(f"Input focal sequences: {len(focal_headers)}\n")
        file.write(f"Original alignment length: {alignment_length}\n\n")

        file.write("Consensus rule:\n")
        file.write("Gaps/missing states are ignored when assigning consensus states.\n")
        file.write("Only strict A/C/G/T states are counted for consensus assignment.\n")
        file.write(
            "A strict base is used only when its frequency exceeds the "
            "second-most frequent base by > 0.5.\n"
        )
        file.write(
            "Otherwise, the exact IUPAC ambiguity code for the observed strict "
            "base set is used.\n"
        )
        file.write("If no strict A/C/G/T state exists at a site, N is used.\n\n")

        file.write("Trimming rule:\n")
        file.write(
            "Trimmed consensus removes focal columns classified as PRL or INS "
            "by the punishment gap-dominated-column logic.\n"
        )
        file.write(
            f"A column is PRL/INS-like when empty_fraction >= "
            f"{BALANCING_EMPTY_LIMIT:.6f}.\n\n"
        )

        file.write("Summary:\n")
        file.write(f"Untrimmed length: {len(consensus_result.untrimmed_sequence)}\n")
        file.write(f"Untrimmed N sites: {consensus_result.untrimmed_n_count}\n")
        file.write(f"Trimmed length: {len(consensus_result.trimmed_sequence)}\n")
        file.write(f"Trimmed N sites: {consensus_result.trimmed_n_count}\n")
        file.write(f"Removed PRL sites: {len(consensus_result.removed_prl_indices)}\n")
        file.write(f"Removed INS sites: {len(consensus_result.removed_ins_indices)}\n\n")

        file.write("Removed PRL sites, 1-based:\n")
        if consensus_result.removed_prl_indices:
            file.write(
                ", ".join(
                    str(site + 1)
                    for site in consensus_result.removed_prl_indices
                )
                + "\n\n"
            )
        else:
            file.write("None\n\n")

        file.write("Removed INS sites, 1-based:\n")
        if consensus_result.removed_ins_indices:
            file.write(
                ", ".join(
                    str(site + 1)
                    for site in consensus_result.removed_ins_indices
                )
                + "\n\n"
            )
        else:
            file.write("None\n\n")

        file.write("Untrimmed focal consensus:\n")
        file.write(f">{target_string}_focal_consensus_untrimmed\n")
        file.write(wrap_sequence(consensus_result.untrimmed_sequence))
        file.write("\n\n")

        file.write("Trimmed focal consensus:\n")
        file.write(f">{target_string}_focal_consensus_trimmed\n")
        file.write(wrap_sequence(consensus_result.trimmed_sequence))
        file.write("\n\n")

        if dmc_sites is not None:
            file.write("DMC site mapping into trimmed consensus:\n")

            if not dmc_sites:
                file.write("No DMC sites supplied.\n")
            else:
                for site in sorted(dmc_sites):
                    site_1based = site + 1

                    if site in kept_index_to_trimmed_position:
                        trimmed_position = kept_index_to_trimmed_position[site]
                        base = consensus_result.trimmed_sequence[trimmed_position - 1]
                        file.write(
                            f"{site_1based} -> trimmed {trimmed_position} -> {base}\n"
                        )
                    elif site in removed_prl:
                        file.write(f"{site_1based} -> DROPPED as PRL\n")
                    elif site in removed_ins:
                        file.write(f"{site_1based} -> DROPPED as INS\n")
                    else:
                        file.write(f"{site_1based} -> DROPPED\n")