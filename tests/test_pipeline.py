from pathlib import Path

from molecular_diagnosis.pipeline import run_pipeline_core


def test_run_pipeline_core_creates_outputs(tmp_path: Path) -> None:
    fasta = tmp_path / "input.fasta"

    fasta.write_text(
        ">focal_1\nAACC\n"
        ">focal_2\nAACC\n"
        ">other_1\nAGCT\n"
        ">other_2\nCGTT\n",
        encoding="utf-8",
    )

    result = run_pipeline_core(
        fasta_path=fasta,
        target_string="focal",
        output_dir=tmp_path,
    )

    assert result.txt_output_path.exists()
    assert result.xlsx_output_path.exists()

    assert result.txt_output_path.name == "DMCs_output.txt"
    assert result.xlsx_output_path.name == "comparison_output.xlsx"


def test_run_pipeline_core_uses_next_available_filename(tmp_path: Path) -> None:
    fasta = tmp_path / "input.fasta"

    fasta.write_text(
        ">focal_1\nAACC\n"
        ">focal_2\nAACC\n"
        ">other_1\nAGCT\n"
        ">other_2\nCGTT\n",
        encoding="utf-8",
    )

    (tmp_path / "DMCs_output.txt").write_text("existing", encoding="utf-8")
    (tmp_path / "comparison_output.xlsx").write_text("existing", encoding="utf-8")

    result = run_pipeline_core(
        fasta_path=fasta,
        target_string="focal",
        output_dir=tmp_path,
    )

    assert result.txt_output_path.name == "DMCs_output(2).txt"
    assert result.xlsx_output_path.name == "comparison_output(2).xlsx"