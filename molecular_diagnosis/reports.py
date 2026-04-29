from pathlib import Path

from molecular_diagnosis.core import format_diag
from molecular_diagnosis.models import DMCResult, FiveSiteResult


def write_text_report(
    output_path: str | Path,
    fasta_path: str | Path,
    output_dir: str | Path,
    target_string: str,
    sequences: dict[str, str],
    alignment_length: int,
    focal_headers: list[str],
    non_focal_headers: list[str],
    ref_id: str,
    dmc: DMCResult,
    five_site_result: FiveSiteResult,
) -> None:
    output_path = Path(output_path)
    ref_seq = sequences[ref_id]
    sites = dmc.unique

    with open(output_path, "w", encoding="utf-8") as file:
        file.write("Molecular diagnostic report\n\n")

        file.write("Input settings:\n")
        file.write(f"FASTA file: {fasta_path}\n")
        file.write(f"Output directory: {output_dir}\n")
        file.write(f"Focal identifier string: {target_string}\n\n")

        file.write("Sequence summary:\n")
        file.write(f"Total sequences read: {len(sequences)}\n")
        file.write(f"Alignment length: {alignment_length}\n")
        file.write(f"Focal sequences: {len(focal_headers)}\n")
        file.write(f"Non-focal sequences: {len(non_focal_headers)}\n")
        file.write(f"Reference sequence selected: {ref_id}\n\n")

        file.write("Run diagnostics:\n")
        file.write(f"Sites skipped because focal column contains non-ACGT: {dmc.skipped_non_acgt}\n")
        file.write(f"Sites fixed in focal set (strict A/C/G/T only): {dmc.fixed_count}\n")
        file.write(f"Fixed sites removed as globally conserved: {dmc.globally_conserved_removed}\n")
        file.write(f"Candidate diagnostic sites retained: {dmc.candidate_count}\n")
        file.write(f"1-site diagnoses found: {len(dmc.single)}\n")
        file.write(f"2-site combinations tested: {dmc.pairs_tested}\n")
        file.write(f"Valid 2-site diagnostic combinations recovered: {len(dmc.pairs)}\n")
        file.write(f"Unique sites participating in valid 2-site combinations: {len(sites)}\n")
        file.write(f"5-site combinations tested: {five_site_result.total_combinations_tested}\n\n")

        file.write("1-site diagnosis:\n")
        if not dmc.single:
            file.write("None found\n\n")
        else:
            file.write(format_diag(dmc.single, ref_seq) + "\n\n")

        file.write("2-site diagnostic combinations:\n")
        if not dmc.pairs:
            file.write("None found\n")
        else:
            for site_a, site_b in dmc.pairs:
                file.write(f"[{site_a + 1}:{ref_seq[site_a]}, {site_b + 1}:{ref_seq[site_b]}]\n")

        file.write(f"\nTotal: {len(dmc.pairs)}\n\n")

        file.write("Unique diagnostic sites:\n")
        file.write(format_diag(sites, ref_seq) + "\n\n")

        file.write(f"Full diagnosis ({len(sites)} sites):\n")
        file.write(format_diag(sites, ref_seq) + "\n\n")

        file.write("5-site diagnosis (for similarity gap optimisation):\n")
        if five_site_result.best_gap_sites is None:
            file.write("Not computed: fewer than 5 unique diagnostic sites.\n\n")
        else:
            file.write(format_diag(five_site_result.best_gap_sites, ref_seq) + "\n\n")

        file.write("5-site diagnosis (for average similarity optimisation):\n")
        if five_site_result.best_avg_sites is None:
            file.write("Not computed: fewer than 5 unique diagnostic sites.\n\n")
        else:
            file.write(format_diag(five_site_result.best_avg_sites, ref_seq) + "\n\n")