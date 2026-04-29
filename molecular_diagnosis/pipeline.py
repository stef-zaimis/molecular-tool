from pathlib import Path

from molecular_diagnosis.constants import TXT_OUTPUT_BASENAME, XLSX_OUTPUT_BASENAME
from molecular_diagnosis.core import find_best_five_site_sets, find_dmc_information
from molecular_diagnosis.excel import write_excel_report
from molecular_diagnosis.fasta_io import (
    parse_fasta,
    split_focal_headers,
    validate_aligned_fasta,
)
from molecular_diagnosis.models import PipelineResult
from molecular_diagnosis.reports import write_text_report
from molecular_diagnosis.utils import next_available_filename


def run_pipeline_core(
    fasta_path: str | Path,
    target_string: str,
    output_dir: str | Path,
) -> PipelineResult:
    fasta_path = Path(str(fasta_path).strip())
    target_string = str(target_string).strip()
    output_dir = Path(str(output_dir).strip())

    if not str(fasta_path):
        raise ValueError("No FASTA file selected.")

    if not target_string:
        raise ValueError("No identifier string entered.")

    if not str(output_dir):
        raise ValueError("No output directory selected.")

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