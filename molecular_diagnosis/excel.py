from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.worksheet.worksheet import Worksheet

from molecular_diagnosis.constants import COLORS
from molecular_diagnosis.core import (
    compute_match_score,
    compute_similarity,
    extract_sites,
)


def autosize_columns(ws: Worksheet) -> None:
    for column in ws.columns:
        max_len = 0
        column_letter = column[0].column_letter

        for cell in column:
            if cell.value is not None:
                max_len = max(max_len, len(str(cell.value)))

        ws.column_dimensions[column_letter].width = max_len + 2


def build_sheet(
    workbook: Workbook,
    name: str,
    sequences: dict[str, str],
    ref_id: str,
    sites: list[int] | tuple[int, ...],
    target_string: str,
) -> None:
    ws = workbook.create_sheet(name)

    header = ["ID"] + [site + 1 for site in sites] + ["matches", "% similarity", "avg similarity"]
    ws.append(header)

    for cell in ws[1]:
        cell.font = Font(bold=True)

    ws.freeze_panes = "B3"

    ref_states = extract_sites(sequences[ref_id], sites)

    rows = []

    for sequence_id, sequence in sequences.items():
        if sequence_id != ref_id and target_string in sequence_id:
            continue

        states = extract_sites(sequence, sites)
        match_score = compute_match_score(ref_states, states)
        similarity = compute_similarity(ref_states, states)

        rows.append((sequence_id, states, match_score, similarity))

    ref_row = next(row for row in rows if row[0] == ref_id)
    non_ref_rows = [row for row in rows if row[0] != ref_id]

    avg_similarity = sum(row[3] for row in non_ref_rows) / len(non_ref_rows)

    non_ref_rows = sorted(non_ref_rows, key=lambda row: (-row[2], -row[3], row[0]))
    final_rows = [ref_row] + non_ref_rows

    for row_index, (sequence_id, states, match_score, similarity) in enumerate(final_rows, start=2):
        avg_value = round(avg_similarity * 100, 2) if row_index == 2 else ""

        ws.append(
            [
                sequence_id,
                *states,
                round(match_score, 2),
                round(similarity * 100, 2),
                avg_value,
            ]
        )

        for column_index, base in enumerate(states, start=2):
            color = COLORS.get(base, "FFFFFF")

            ws.cell(row_index, column_index).fill = PatternFill(
                start_color=color,
                end_color=color,
                fill_type="solid",
            )

    autosize_columns(ws)


def write_excel_report(
    output_path: str | Path,
    sequences: dict[str, str],
    ref_id: str,
    full_sites: list[int],
    target_string: str,
    best_gap_sites: tuple[int, ...] | None,
    best_avg_sites: tuple[int, ...] | None,
) -> None:
    workbook = Workbook()
    workbook.remove(workbook.active)

    build_sheet(
        workbook=workbook,
        name="Full",
        sequences=sequences,
        ref_id=ref_id,
        sites=full_sites,
        target_string=target_string,
    )

    if best_gap_sites is not None:
        build_sheet(
            workbook=workbook,
            name="Gap5",
            sequences=sequences,
            ref_id=ref_id,
            sites=best_gap_sites,
            target_string=target_string,
        )

    if best_avg_sites is not None:
        build_sheet(
            workbook=workbook,
            name="Avg5",
            sequences=sequences,
            ref_id=ref_id,
            sites=best_avg_sites,
            target_string=target_string,
        )

    workbook.save(output_path)