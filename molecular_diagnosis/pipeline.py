from pathlib import Path

from molecular_diagnosis.constants import (
    CONSENSUS_TXT_OUTPUT_BASENAME,
    PUNISHMENT_XLSX_OUTPUT_BASENAME,
    SEQUENCE_SUBSETS_XLSX_OUTPUT_BASENAME,
    TXT_OUTPUT_BASENAME,
    XLSX_OUTPUT_BASENAME,
)
from molecular_diagnosis.consensus import (
    build_focal_consensus_result,
    write_consensus_text_report,
)
from molecular_diagnosis.sequence_subsets import (
    group_sequence_subsets,
    write_sequence_subset_excel_report,
)
from molecular_diagnosis.core import find_best_five_site_sets, find_dmc_information
from molecular_diagnosis.focal import FocalSelector, normalise_focal_strings
from molecular_diagnosis.excel import write_excel_report, write_punishment_excel_report
from molecular_diagnosis.fasta_io import (
    parse_fasta,
    split_focal_headers,
    validate_aligned_fasta,
)
from molecular_diagnosis.models import PipelineResult, PunishmentPipelineResult
from molecular_diagnosis.punishments import find_focal_punishments
from molecular_diagnosis.reports import write_text_report
from molecular_diagnosis.utils import next_available_filename


def load_inputs(
    fasta_path: str | Path,
    focal_strings: FocalSelector,
    output_dir: str | Path,
) -> tuple[
    Path,
    list[str],
    Path,
    dict[str, str],
    int,
    list[str],
    list[str],
]:
    fasta_path_text = str(fasta_path).strip()
    output_dir_text = str(output_dir).strip()

    if not fasta_path_text:
        raise ValueError("No FASTA file selected.")

    # Keeps the original wording for the "nothing supplied" case, and keeps
    # the original order in which these three errors fire.
    selectors = normalise_focal_strings(focal_strings)

    if not output_dir_text:
        raise ValueError("No output directory selected.")

    fasta_path = Path(fasta_path_text)
    output_dir = Path(output_dir_text)

    if not fasta_path.exists():
        raise ValueError("FASTA file does not exist.")

    if not fasta_path.is_file():
        raise ValueError("FASTA path is not a file.")

    if not output_dir.exists():
        raise ValueError("Output directory does not exist.")

    if not output_dir.is_dir():
        raise ValueError("Output path is not a directory.")

    sequences = parse_fasta(fasta_path)
    alignment_length = validate_aligned_fasta(sequences)

    focal_headers, non_focal_headers = split_focal_headers(
        sequences=sequences,
        focal_strings=selectors,
    )

    return (
        fasta_path,
        selectors,
        output_dir,
        sequences,
        alignment_length,
        focal_headers,
        non_focal_headers,
    )


def run_pipeline_core(
    fasta_path: str | Path,
    focal_strings: FocalSelector,
    output_dir: str | Path,
    *,
    include_ambiguous_dmc_bd: bool = False,
    include_gappy_consensus_dmc_sites: bool = False,
    min_combination_length: int = 1,
    max_combination_length: int = 2,
    start_combination_length: int = 1,
    initial_diagnostic_combinations: list[tuple[int, ...]] | None = None,
    initial_combinations_tested_by_length: dict[int, int] | None = None,
) -> PipelineResult:
    (
        fasta_path,
        selectors,
        output_dir,
        sequences,
        alignment_length,
        focal_headers,
        non_focal_headers,
    ) = load_inputs(
        fasta_path=fasta_path,
        focal_strings=focal_strings,
        output_dir=output_dir,
    )

    return run_pipeline_on_sequences(
        sequences=sequences,
        selectors=selectors,
        focal_headers=focal_headers,
        non_focal_headers=non_focal_headers,
        alignment_length=alignment_length,
        source_label=str(fasta_path),
        output_dir=output_dir,
        include_ambiguous_dmc_bd=include_ambiguous_dmc_bd,
        include_gappy_consensus_dmc_sites=include_gappy_consensus_dmc_sites,
        min_combination_length=min_combination_length,
        max_combination_length=max_combination_length,
        start_combination_length=start_combination_length,
        initial_diagnostic_combinations=initial_diagnostic_combinations,
        initial_combinations_tested_by_length=initial_combinations_tested_by_length,
    )


def run_pipeline_on_sequences(
    *,
    sequences: dict[str, str],
    selectors: FocalSelector,
    focal_headers: list[str],
    non_focal_headers: list[str],
    alignment_length: int,
    source_label: str,
    output_dir: str | Path,
    include_ambiguous_dmc_bd: bool = False,
    include_gappy_consensus_dmc_sites: bool = False,
    min_combination_length: int = 1,
    max_combination_length: int = 2,
    start_combination_length: int = 1,
    initial_diagnostic_combinations: list[tuple[int, ...]] | None = None,
    initial_combinations_tested_by_length: dict[int, int] | None = None,
) -> PipelineResult:
    """
    The pipeline body, over sequences that are ALREADY parsed and verified.

    Extracted from `run_pipeline_core` so a multi-file run can combine several
    verified alignments in memory instead of writing and re-parsing a temporary
    combined FASTA. `run_pipeline_core` calls straight through to it, so the
    single-file path is unchanged.

    `source_label` is what the report prints as the input; it is a label only.
    """
    output_dir = Path(output_dir)
    ref_id = focal_headers[0]

    focal_sequences = [
        sequences[header]
        for header in focal_headers
    ]

    consensus_result = build_focal_consensus_result(
        focal_sequences=focal_sequences,
    )

    dmc = find_dmc_information(
        sequences=sequences,
        focal_strings=selectors,
        include_ambiguous_dmc_bd=include_ambiguous_dmc_bd,
        include_gappy_consensus_dmc_sites=include_gappy_consensus_dmc_sites,
        min_combination_length=min_combination_length,
        max_combination_length=max_combination_length,
        start_combination_length=start_combination_length,
        initial_diagnostic_combinations=initial_diagnostic_combinations,
        initial_combinations_tested_by_length=initial_combinations_tested_by_length,
    )

    five_site_result = find_best_five_site_sets(
        sequences=sequences,
        ref_id=ref_id,
        sites=dmc.unique,
        focal_strings=selectors,
        diagnostic_states=dmc.states,
    )

    txt_output_path = next_available_filename(output_dir / TXT_OUTPUT_BASENAME)
    xlsx_output_path = next_available_filename(output_dir / XLSX_OUTPUT_BASENAME)
    consensus_txt_output_path = next_available_filename(
        output_dir / CONSENSUS_TXT_OUTPUT_BASENAME
    )

    write_text_report(
        output_path=txt_output_path,
        fasta_path=source_label,
        output_dir=output_dir,
        focal_strings=selectors,
        sequences=sequences,
        alignment_length=alignment_length,
        focal_headers=focal_headers,
        non_focal_headers=non_focal_headers,
        ref_id=ref_id,
        dmc=dmc,
        five_site_result=five_site_result,
        punishment_result=None,
    )

    write_consensus_text_report(
        output_path=consensus_txt_output_path,
        focal_strings=selectors,
        focal_headers=focal_headers,
        alignment_length=alignment_length,
        consensus_result=consensus_result,
        dmc_sites=dmc.unique,
    )

    write_excel_report(
        output_path=xlsx_output_path,
        sequences=sequences,
        ref_id=ref_id,
        full_sites=dmc.unique,
        focal_strings=selectors,
        best_gap_sites=five_site_result.best_gap_sites,
        best_avg_sites=five_site_result.best_avg_sites,
        diagnostic_states=dmc.states,
    )

    return PipelineResult(
        txt_output_path=txt_output_path,
        xlsx_output_path=xlsx_output_path,
        dmc=dmc,
        consensus_txt_output_path=consensus_txt_output_path,
    )


def run_punishment_core(
    fasta_path: str | Path,
    focal_strings: FocalSelector,
    output_dir: str | Path,
) -> PunishmentPipelineResult:
    (
        _fasta_path,
        _selectors,
        output_dir,
        sequences,
        alignment_length,
        focal_headers,
        _non_focal_headers,
    ) = load_inputs(
        fasta_path=fasta_path,
        focal_strings=focal_strings,
        output_dir=output_dir,
    )

    punishment_result = find_focal_punishments(
        sequences=sequences,
        focal_headers=focal_headers,
    )

    focal_sequences = {
        header: sequences[header]
        for header in focal_headers
    }

    sequence_subset_groups = group_sequence_subsets(
        sequences=focal_sequences,
    )

    xlsx_output_path = next_available_filename(
        output_dir / PUNISHMENT_XLSX_OUTPUT_BASENAME
    )

    sequence_subsets_xlsx_output_path = next_available_filename(
        output_dir / SEQUENCE_SUBSETS_XLSX_OUTPUT_BASENAME
    )

    write_sequence_subset_excel_report(
        output_path=sequence_subsets_xlsx_output_path,
        sequence_subset_groups=sequence_subset_groups,
    )

    write_punishment_excel_report(
        output_path=xlsx_output_path,
        focal_headers=focal_headers,
        alignment_length=alignment_length,
        punishment_result=punishment_result,
    )

    return PunishmentPipelineResult(
        xlsx_output_path=xlsx_output_path,
        sequence_subsets_xlsx_output_path=sequence_subsets_xlsx_output_path,
    )