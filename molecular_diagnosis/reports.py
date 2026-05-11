from pathlib import Path

from molecular_diagnosis.core import format_diag_from_states
from molecular_diagnosis.models import DMCResult, FiveSiteResult, PunishmentResult


def format_combo(
    combo: tuple[int, ...],
    states: dict[int, str],
) -> str:
    return "[" + ", ".join(f"{site + 1}:{states[site]}" for site in combo) + "]"


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
    punishment_result: PunishmentResult | None = None,
) -> None:
    output_path = Path(output_path)
    sites = dmc.unique

    with open(output_path, "w", encoding="utf-8") as file:
        file.write("Molecular diagnostic report\n\n")

        file.write("Input settings:\n")
        file.write(f"FASTA file: {fasta_path}\n")
        file.write(f"Output directory: {output_dir}\n")
        file.write(f"Focal identifier string: {target_string}\n")
        file.write(f"DMC ambiguous-site benefit of doubt: {dmc.include_ambiguous_dmc_bd}\n")
        file.write(f"DMC gappy consensus sites included: {dmc.include_gappy_consensus_dmc_sites}\n")
        file.write(f"Minimum combination length: {dmc.min_combination_length}\n")
        file.write(f"Maximum combination length: {dmc.max_combination_length}\n")
        file.write(f"Search started at combination length: {dmc.start_combination_length}\n\n")

        file.write("Sequence summary:\n")
        file.write(f"Total sequences read: {len(sequences)}\n")
        file.write(f"Alignment length: {alignment_length}\n")
        file.write(f"Focal sequences: {len(focal_headers)}\n")
        file.write(f"Non-focal sequences: {len(non_focal_headers)}\n")
        file.write(f"Reference sequence selected: {ref_id}\n\n")

        file.write("Run diagnostics:\n")
        file.write(
            "Sites skipped because focal column lacked an allowed consensus state "
            f"or contained disabled ambiguity/gap states: {dmc.skipped_non_acgt}\n"
        )
        file.write(f"Focal consensus sites retained before global filtering: {dmc.fixed_count}\n")
        file.write(f"Consensus sites removed as globally conserved: {dmc.globally_conserved_removed}\n")
        file.write(f"Candidate diagnostic sites retained: {dmc.candidate_count}\n")
        file.write(f"Candidate sites included using ambiguous-site BD: {len(dmc.ambiguous_bd_sites_included)}\n")
        file.write(f"Gappy consensus candidate sites included: {len(dmc.gappy_consensus_sites_included)}\n")
        file.write(f"Search stopped at combination length: {dmc.stopped_at_length}\n")
        file.write(f"Stop reason: {dmc.stop_reason}\n")

        file.write("Combinations tested by length:\n")
        if not dmc.combinations_tested_by_length:
            file.write("None\n")
        else:
            for combo_length in sorted(dmc.combinations_tested_by_length):
                file.write(
                    f"{combo_length}-site combinations tested after pruning: "
                    f"{dmc.combinations_tested_by_length[combo_length]}\n"
                )

        file.write(f"Total combinations tested after pruning: {dmc.total_combinations_tested}\n")
        file.write(f"Diagnostic combinations recovered: {len(dmc.diagnostic_combinations)}\n")
        file.write(f"1-site diagnoses found: {len(dmc.single)}\n")
        file.write(f"2-site diagnostic combinations recovered: {len(dmc.pairs)}\n")
        file.write(f"Unique sites participating in recovered diagnostic combinations: {len(sites)}\n")
        file.write(f"5-site combinations tested: {five_site_result.total_combinations_tested}\n\n")

        file.write("Diagnostic combinations by length:\n")

        if not dmc.diagnostic_combinations_by_length:
            file.write("None found\n\n")
        else:
            for combo_length in sorted(dmc.diagnostic_combinations_by_length):
                combos = dmc.diagnostic_combinations_by_length[combo_length]

                file.write(f"{combo_length}-site diagnostic combinations:\n")

                for combo in combos:
                    file.write(format_combo(combo, dmc.states) + "\n")

                file.write(f"Total {combo_length}-site combinations: {len(combos)}\n\n")

        file.write("Unique diagnostic sites:\n")
        file.write(format_diag_from_states(sites, dmc.states) + "\n\n")

        file.write(f"Full diagnosis ({len(sites)} sites):\n")
        file.write(format_diag_from_states(sites, dmc.states) + "\n\n")

        file.write("5-site diagnosis (for similarity gap optimisation):\n")
        if five_site_result.best_gap_sites is None:
            file.write("Not computed: fewer than 5 unique diagnostic sites.\n\n")
        else:
            file.write(
                format_diag_from_states(
                    five_site_result.best_gap_sites,
                    dmc.states,
                )
                + "\n\n"
            )

        file.write("5-site diagnosis (for average similarity optimisation):\n")
        if five_site_result.best_avg_sites is None:
            file.write("Not computed: fewer than 5 unique diagnostic sites.\n\n")
        else:
            file.write(
                format_diag_from_states(
                    five_site_result.best_avg_sites,
                    dmc.states,
                )
                + "\n\n"
            )

        if punishment_result is not None:
            file.write("Punishment summary:\n")
            file.write(
                "Punishment results were supplied to the report writer, "
                "but detailed punishment text reporting is not currently implemented.\n"
            )