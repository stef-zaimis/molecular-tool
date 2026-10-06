"""
Parity driver: dumps a complete, JSON-comparable picture of one Molecular
Diagnosis run.

This file is executed by `tests/test_parity.py` in TWO trees — the current
working tree and a read-only git worktree of the pre-integration commit — and
the two dumps are compared.

It therefore has to run unchanged against BOTH versions of the package, which
constrains it in two ways:

  * every call uses POSITIONAL arguments, because the focal parameter was
    renamed (`target_string` -> `focal_strings`) during the integration work;
  * it imports nothing that only exists on one side (no `focal`, no `service`).

Run it as:  PYTHONPATH=<tree> python tests/parity_driver.py <spec.json> <out.json>
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def jsonify(value):
    """Make results comparable: tuples -> lists, int keys -> str keys."""
    if isinstance(value, dict):
        return {str(key): jsonify(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        items = [jsonify(item) for item in value]
        return sorted(items, key=repr) if isinstance(value, (set, frozenset)) else items
    return value


def dump_dmc(dmc) -> dict:
    return {
        "single": jsonify(dmc.single),
        "pairs": jsonify(dmc.pairs),
        "unique": jsonify(dmc.unique),
        "states": jsonify(dmc.states),
        "diagnostic_combinations": jsonify(dmc.diagnostic_combinations),
        "diagnostic_combinations_by_length": jsonify(dmc.diagnostic_combinations_by_length),
        "min_combination_length": dmc.min_combination_length,
        "max_combination_length": dmc.max_combination_length,
        "start_combination_length": dmc.start_combination_length,
        "stopped_at_length": dmc.stopped_at_length,
        "stop_reason": dmc.stop_reason,
        "combinations_tested_by_length": jsonify(dmc.combinations_tested_by_length),
        "total_combinations_tested": dmc.total_combinations_tested,
        "fixed_count": dmc.fixed_count,
        "skipped_non_acgt": dmc.skipped_non_acgt,
        "globally_conserved_removed": dmc.globally_conserved_removed,
        "candidate_count": dmc.candidate_count,
        "pairs_tested": dmc.pairs_tested,
        "ambiguous_bd_sites_included": jsonify(dmc.ambiguous_bd_sites_included),
        "gappy_consensus_sites_included": jsonify(dmc.gappy_consensus_sites_included),
        "include_ambiguous_dmc_bd": dmc.include_ambiguous_dmc_bd,
        "include_gappy_consensus_dmc_sites": dmc.include_gappy_consensus_dmc_sites,
    }


def dump_five_site(result) -> dict:
    return {
        "total_combinations_tested": result.total_combinations_tested,
        "best_gap_score": result.best_gap_score,
        "best_gap_sites": jsonify(result.best_gap_sites),
        "best_avg_score": result.best_avg_score,
        "best_avg_sites": jsonify(result.best_avg_sites),
    }


def dump_consensus(result) -> dict:
    return {
        "untrimmed_sequence": result.untrimmed_sequence,
        "trimmed_sequence": result.trimmed_sequence,
        "kept_indices": jsonify(result.kept_indices),
        "removed_prl_indices": jsonify(result.removed_prl_indices),
        "removed_ins_indices": jsonify(result.removed_ins_indices),
        "untrimmed_n_count": result.untrimmed_n_count,
        "trimmed_n_count": result.trimmed_n_count,
    }


def dump_workbook(path: Path) -> dict:
    """
    Logical contents of an xlsx: sheets, cell values and the styling the
    reports actually set. Deliberately NOT a byte comparison — the zip carries
    timestamps and openpyxl metadata that differ between runs.
    """
    from openpyxl import load_workbook

    workbook = load_workbook(path)
    sheets = {}

    for worksheet in workbook.worksheets:
        cells = []
        for row in worksheet.iter_rows():
            for cell in row:
                if cell.value is None and (cell.fill is None or cell.fill.fill_type is None):
                    continue
                fill = cell.fill
                colour = None
                if fill is not None and fill.fill_type == "solid":
                    colour = getattr(fill.start_color, "rgb", None)
                    if not isinstance(colour, str):
                        colour = None
                cells.append(
                    {
                        "ref": cell.coordinate,
                        "value": cell.value,
                        "number_format": cell.number_format,
                        "bold": bool(cell.font.bold) if cell.font else False,
                        "fill": colour,
                        "alignment": cell.alignment.horizontal if cell.alignment else None,
                    }
                )

        sheets[worksheet.title] = {
            "max_row": worksheet.max_row,
            "max_column": worksheet.max_column,
            "freeze_panes": worksheet.freeze_panes,
            "auto_filter": worksheet.auto_filter.ref,
            "column_widths": {
                letter: round(dim.width, 4) if dim.width is not None else None
                for letter, dim in sorted(worksheet.column_dimensions.items())
            },
            "cells": cells,
        }

    return {"sheet_names": workbook.sheetnames, "sheets": sheets}


def main() -> int:
    spec = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    out_path = Path(sys.argv[2])

    from molecular_diagnosis.consensus import build_focal_consensus_result
    from molecular_diagnosis.core import find_best_five_site_sets, find_dmc_information
    from molecular_diagnosis.fasta_io import (
        parse_fasta,
        split_focal_headers,
        validate_aligned_fasta,
    )
    from molecular_diagnosis.pipeline import run_pipeline_core

    fasta_path = spec["fasta_path"]
    focal = spec["focal"]
    output_dir = Path(spec["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    result: dict = {}

    # --- selection -------------------------------------------------------
    sequences = parse_fasta(fasta_path)
    result["alignment_length"] = validate_aligned_fasta(sequences)
    result["all_headers"] = list(sequences)

    focal_headers, non_focal_headers = split_focal_headers(sequences, focal)
    result["focal_headers"] = focal_headers
    result["non_focal_headers"] = non_focal_headers
    result["ref_id"] = focal_headers[0]

    # --- consensus -------------------------------------------------------
    focal_sequences = [sequences[header] for header in focal_headers]
    result["consensus"] = dump_consensus(build_focal_consensus_result(focal_sequences))

    # --- DMC search ------------------------------------------------------
    resume = spec.get("resume")
    dmc_kwargs = {
        "include_ambiguous_dmc_bd": spec["include_ambiguous_dmc_bd"],
        "include_gappy_consensus_dmc_sites": spec["include_gappy_consensus_dmc_sites"],
        "min_combination_length": spec["min_combination_length"],
        "max_combination_length": spec["max_combination_length"],
    }
    if resume:
        dmc_kwargs["start_combination_length"] = resume["start_combination_length"]
        dmc_kwargs["initial_diagnostic_combinations"] = [
            tuple(combo) for combo in resume["initial_diagnostic_combinations"]
        ]
        dmc_kwargs["initial_combinations_tested_by_length"] = {
            int(length): count
            for length, count in resume["initial_combinations_tested_by_length"].items()
        }

    dmc = find_dmc_information(sequences, focal, **dmc_kwargs)
    result["dmc"] = dump_dmc(dmc)

    # --- five-site optimisation -----------------------------------------
    five_site = find_best_five_site_sets(
        sequences, focal_headers[0], dmc.unique, focal, dmc.states
    )
    result["five_site"] = dump_five_site(five_site)

    # --- full pipeline, including every written file ---------------------
    pipeline_kwargs = dict(dmc_kwargs)
    pipeline = run_pipeline_core(fasta_path, focal, str(output_dir), **pipeline_kwargs)

    report_path = Path(pipeline.txt_output_path)
    workbook_path = Path(pipeline.xlsx_output_path)
    consensus_path = (
        Path(pipeline.consensus_txt_output_path)
        if pipeline.consensus_txt_output_path is not None
        else None
    )

    result["outputs"] = {
        # Names only: the absolute paths differ because the trees differ.
        "report_name": report_path.name,
        "workbook_name": workbook_path.name,
        "consensus_name": consensus_path.name if consensus_path else None,
    }
    # The DMC report echoes the output directory back, and the two sides
    # deliberately write to different directories so their files cannot
    # collide. Redact that one path so the comparison is about content.
    def redact(text: str) -> str:
        return text.replace(str(output_dir), "<OUTPUT_DIR>").replace(
            str(output_dir).replace("\\", "/"), "<OUTPUT_DIR>"
        )

    result["report_text"] = redact(report_path.read_text(encoding="utf-8"))
    result["consensus_text"] = (
        redact(consensus_path.read_text(encoding="utf-8")) if consensus_path else None
    )
    result["workbook"] = dump_workbook(workbook_path)
    result["pipeline_dmc"] = dump_dmc(pipeline.dmc)

    out_path.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
