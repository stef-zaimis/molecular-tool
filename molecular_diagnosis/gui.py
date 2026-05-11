from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog

from molecular_diagnosis.fasta_io import parse_fasta
from molecular_diagnosis.pipeline import run_pipeline_core, run_punishment_core
from molecular_diagnosis.viewer import open_fasta_viewer


def launch_gui() -> None:
    root = tk.Tk()
    root.title("Molecular Diagnosis Tool")
    root.geometry("650x470")

    fasta_var = tk.StringVar()
    target_var = tk.StringVar()
    output_dir_var = tk.StringVar()
    
    include_ambiguous_dmc_bd_var = tk.BooleanVar(value=False)
    include_gappy_consensus_dmc_sites_var = tk.BooleanVar(value=False)

    min_combination_length_var = tk.IntVar(value=1)
    max_combination_length_var = tk.IntVar(value=2)

    def browse_fasta() -> None:
        path = filedialog.askopenfilename(
            initialdir=str(Path.home()),
            title="Select aligned FASTA file",
            filetypes=[
                ("FASTA files", "*.fasta *.fa *.fna"),
                ("All files", "*.*"),
            ],
        )

        if path:
            fasta_var.set(path)

    def browse_output() -> None:
        folder = filedialog.askdirectory(
            initialdir=str(Path.home()),
            title="Select output folder",
        )

        if folder:
            output_dir_var.set(folder)

    def use_fasta_location() -> None:
        fasta_path = fasta_var.get().strip()

        if not fasta_path:
            messagebox.showerror("Error", "Select a FASTA file first.")
            return

        output_dir_var.set(str(Path(fasta_path).parent))

    def view_fasta() -> None:
        fasta_path = fasta_var.get().strip()

        if not fasta_path:
            messagebox.showerror("Error", "Please select a FASTA file first.")
            return

        try:
            sequences = parse_fasta(fasta_path)
        except Exception as error:
            messagebox.showerror(
                "Error",
                f"Unable to read the file.\n\n"
                f"Error type: {type(error).__name__}\n"
                f"Details: {error}",
            )
            return

        if not sequences:
            messagebox.showerror(
                "Error",
                "Unable to read the file.\n\n"
                "Error type: EmptyFASTA\n"
                "Details: No sequences were read from the FASTA file.",
            )
            return

        open_fasta_viewer(root, sequences)
    
    def dmc_summary_message(result) -> str:
        dmc = result.dmc

        found_parts = []

        for combo_length in sorted(dmc.diagnostic_combinations_by_length):
            count = len(dmc.diagnostic_combinations_by_length[combo_length])
            found_parts.append(f"{combo_length}-site: {count}")

        if found_parts:
            found_text = ", ".join(found_parts)
        else:
            found_text = "none"

        return (
            f"Analysis completed.\n\n"
            f"Text output:\n{result.txt_output_path}\n\n"
            f"Excel output:\n{result.xlsx_output_path}\n\n"
            f"DMC search summary:\n"
            f"Found combinations: {found_text}\n"
            f"Stopped at length: {dmc.stopped_at_length}\n"
            f"Stop reason: {dmc.stop_reason}"
        )

    def run() -> None:
        fasta_path = fasta_var.get().strip()
        target_string = target_var.get().strip()
        output_dir = output_dir_var.get().strip()

        if not fasta_path:
            messagebox.showerror("Error", "Please select a FASTA file.")
            return

        if not target_string:
            messagebox.showerror("Error", "Please enter an identifier string.")
            return

        if not output_dir:
            messagebox.showerror("Error", "Please select an output directory.")
            return

        try:
            min_combination_length = int(min_combination_length_var.get())
            max_combination_length = int(max_combination_length_var.get())

            if min_combination_length < 1:
                messagebox.showerror("Error", "Minimum combination length must be at least 1.")
                return

            if max_combination_length < 1:
                messagebox.showerror("Error", "Maximum combination length must be at least 1.")
                return

            if min_combination_length > max_combination_length:
                messagebox.showerror(
                    "Error",
                    "Minimum combination length cannot exceed maximum combination length.",
                )
                return

            start_combination_length = 1
            current_max_combination_length = max_combination_length
            initial_diagnostic_combinations = None
            initial_combinations_tested_by_length = None

            while True:
                result = run_pipeline_core(
                    fasta_path=fasta_path,
                    target_string=target_string,
                    output_dir=output_dir,
                    include_ambiguous_dmc_bd=include_ambiguous_dmc_bd_var.get(),
                    include_gappy_consensus_dmc_sites=include_gappy_consensus_dmc_sites_var.get(),
                    min_combination_length=min_combination_length,
                    max_combination_length=current_max_combination_length,
                    start_combination_length=start_combination_length,
                    initial_diagnostic_combinations=initial_diagnostic_combinations,
                    initial_combinations_tested_by_length=initial_combinations_tested_by_length,
                )

                message = dmc_summary_message(result)

                if result.dmc.stop_reason != "reached_maximum_length":
                    messagebox.showinfo("Done", message)
                    break

                continue_search = messagebox.askyesno(
                    "Continue DMC search?",
                    message
                    + "\n\n"
                    + "The configured maximum combination length was reached. "
                    + "Continue from the next length with a new maximum?",
                )

                if not continue_search:
                    break

                next_start = result.dmc.stopped_at_length + 1

                new_max = simpledialog.askinteger(
                    "New maximum combination length",
                    "Enter the new maximum combination length:",
                    initialvalue=next_start,
                    minvalue=next_start,
                    maxvalue=100,
                )

                if new_max is None:
                    break

                start_combination_length = next_start
                current_max_combination_length = new_max
                initial_diagnostic_combinations = result.dmc.diagnostic_combinations
                initial_combinations_tested_by_length = result.dmc.combinations_tested_by_length

        except Exception as error:
            messagebox.showerror("Error", str(error))

    def run_punishments() -> None:
        fasta_path = fasta_var.get().strip()
        target_string = target_var.get().strip()
        output_dir = output_dir_var.get().strip()

        if not fasta_path:
            messagebox.showerror("Error", "Please select a FASTA file.")
            return

        if not target_string:
            messagebox.showerror("Error", "Please enter an identifier string.")
            return

        if not output_dir:
            messagebox.showerror("Error", "Please select an output directory.")
            return

        try:
            result = run_punishment_core(
                fasta_path=fasta_path,
                target_string=target_string,
                output_dir=output_dir,
            )

            messagebox.showinfo(
                "Done",
                f"Punishment analysis completed.\n\n"
                f"Excel output:\n{result.xlsx_output_path}",
            )

        except Exception as error:
            messagebox.showerror("Error", str(error))

    tk.Label(root, text="Aligned FASTA file").pack(pady=(10, 0))
    tk.Entry(root, textvariable=fasta_var, width=80).pack(padx=10)

    fasta_button_frame = tk.Frame(root)
    fasta_button_frame.pack(pady=4)

    tk.Button(
        fasta_button_frame,
        text="Browse FASTA",
        command=browse_fasta,
    ).grid(row=0, column=0, padx=5)

    tk.Button(
        fasta_button_frame,
        text="View FASTA",
        command=view_fasta,
    ).grid(row=0, column=1, padx=5)

    tk.Label(root, text="Focal identifier string").pack(pady=(10, 0))
    tk.Entry(root, textvariable=target_var, width=80).pack(padx=10)

    options_frame = tk.LabelFrame(root, text="DMC rule options")
    options_frame.pack(fill="x", padx=10, pady=(10, 0))

    tk.Checkbutton(
        options_frame,
        text="Give benefit of doubt to ambiguous bases",
        variable=include_ambiguous_dmc_bd_var,
    ).pack(anchor="w", padx=10, pady=(4, 0))

    tk.Checkbutton(
        options_frame,
        text="Ignore gaps",
        variable=include_gappy_consensus_dmc_sites_var,
    ).pack(anchor="w", padx=10, pady=(0, 4))

    combo_frame = tk.Frame(options_frame)
    combo_frame.pack(anchor="w", padx=10, pady=(6, 4))

    tk.Label(combo_frame, text="Minimum combination length").grid(row=0, column=0, padx=(0, 5))

    tk.Spinbox(
        combo_frame,
        from_=1,
        to=20,
        width=5,
        textvariable=min_combination_length_var,
    ).grid(row=0, column=1, padx=(0, 15))

    tk.Label(combo_frame, text="Maximum combination length").grid(row=0, column=2, padx=(0, 5))

    tk.Spinbox(
        combo_frame,
        from_=1,
        to=20,
        width=5,
        textvariable=max_combination_length_var,
    ).grid(row=0, column=3)

    tk.Label(root, text="Output directory").pack(pady=(10, 0))
    tk.Entry(root, textvariable=output_dir_var, width=80).pack(padx=10)

    button_frame = tk.Frame(root)
    button_frame.pack(pady=4)

    tk.Button(
        button_frame,
        text="Browse Output Folder",
        command=browse_output,
    ).grid(row=0, column=0, padx=5)

    tk.Label(button_frame, text="or").grid(row=0, column=1, padx=5)

    tk.Button(
        button_frame,
        text="Use FASTA Location",
        command=use_fasta_location,
    ).grid(row=0, column=2, padx=5)

    action_frame = tk.Frame(root)
    action_frame.pack(pady=15)

    tk.Button(
        action_frame,
        text="Run DMC Analysis",
        command=run,
    ).grid(row=0, column=0, padx=5)

    tk.Button(
        action_frame,
        text="Run Punishment Analysis",
        command=run_punishments,
    ).grid(row=0, column=1, padx=5)

    root.mainloop()