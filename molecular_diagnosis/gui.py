from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox

from molecular_diagnosis.fasta_io import parse_fasta
from molecular_diagnosis.pipeline import run_pipeline_core
from molecular_diagnosis.viewer import open_fasta_viewer


def launch_gui() -> None:
    root = tk.Tk()
    root.title("Molecular Diagnosis Tool")
    root.geometry("650x360")

    fasta_var = tk.StringVar()
    target_var = tk.StringVar()
    output_dir_var = tk.StringVar()

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
            result = run_pipeline_core(
                fasta_path=fasta_path,
                target_string=target_string,
                output_dir=output_dir,
            )

            messagebox.showinfo(
                "Done",
                f"Analysis completed.\n\n"
                f"Text output:\n{result.txt_output_path}\n\n"
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

    tk.Button(root, text="Run", command=run).pack(pady=15)

    root.mainloop()