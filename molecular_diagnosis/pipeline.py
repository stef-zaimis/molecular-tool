from pathlib import Path

from molecular_diagnosis.constants import (
    PUNISHMENT_XLSX_OUTPUT_BASENAME,
    TXT_OUTPUT_BASENAME,
    XLSX_OUTPUT_BASENAME,
)
from molecular_diagnosis.core import find_best_five_site_sets, find_dmc_information
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
    target_string: str,
    output_dir: str | Path,
) -> tuple[
    Path,
    str,
    Path,
    dict[str, str],
    int,
    list[str],
    list[str],
]:
    fasta_path_text = str(fasta_path).strip()
    target_string = str(target_string).strip()
    output_dir_text = str(output_dir).strip()

    if not fasta_path_text:
        raise ValueError("No FASTA file selected.")

    if not target_string:
        raise ValueError("No identifier string entered.")

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
        target_string=target_string,
    )

    return (
        fasta_path,
        target_string,
        output_dir,
        sequences,
        alignment_length,
        focal_headers,
        non_focal_headers,
    )


def run_pipeline_core(
    fasta_path: str | Path,
    target_string: str,
    output_dir: str | Path,
) -> PipelineResult:
    (
        fasta_path,
        target_string,
        output_dir,
        sequences,
        alignment_length,
        focal_headers,
        non_focal_headers,
    ) = load_inputs(
        fasta_path=fasta_path,
        target_string=target_string,
        output_dir=output_dir,
    )

    ref_id = focal_headers[0]

    dmc = find_dmc_information(
        sequences=sequences,
        target_string=target_string,
    )

    five_site_result = find_best_five_site_sets(
        sequences=sequences,
        ref_id=ref_id,
        sites=dmc.unique,
        target_string=target_string,
    )

    txt_output_path = next_available_filename(output_dir / TXT_OUTPUT_BASENAME)
    xlsx_output_path = next_available_filename(output_dir / XLSX_OUTPUT_BASENAME)

    write_text_report(
        output_path=txt_output_path,
        fasta_path=fasta_path,
        output_dir=output_dir,
        target_string=target_string,
        sequences=sequences,
        alignment_length=alignment_length,
        focal_headers=focal_headers,
        non_focal_headers=non_focal_headers,
        ref_id=ref_id,
        dmc=dmc,
        five_site_result=five_site_result,
        punishment_result=None,
    )

    write_excel_report(
        output_path=xlsx_output_path,
        sequences=sequences,
        ref_id=ref_id,
        full_sites=dmc.unique,
        target_string=target_string,
        best_gap_sites=five_site_result.best_gap_sites,
        best_avg_sites=five_site_result.best_avg_sites,
    )

    return PipelineResult(
        txt_output_path=txt_output_path,
        xlsx_output_path=xlsx_output_path,
    )


def run_punishment_core(
    fasta_path: str | Path,
    target_string: str,
    output_dir: str | Path,
) -> PunishmentPipelineResult:
    (
        _fasta_path,
        _target_string,
        output_dir,
        sequences,
        alignment_length,
        focal_headers,
        _non_focal_headers,
    ) = load_inputs(
        fasta_path=fasta_path,
        target_string=target_string,
        output_dir=output_dir,
    )

    punishment_result = find_focal_punishments(
        sequences=sequences,
        focal_headers=focal_headers,
    )

    xlsx_output_path = next_available_filename(
        output_dir / PUNISHMENT_XLSX_OUTPUT_BASENAME
    )

    write_punishment_excel_report(
        output_path=xlsx_output_path,
        focal_headers=focal_headers,
        alignment_length=alignment_length,
        punishment_result=punishment_result,
    )

    return PunishmentPipelineResult(
        xlsx_output_path=xlsx_output_path,
    )