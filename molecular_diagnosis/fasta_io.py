from pathlib import Path


def parse_fasta(path: str | Path) -> dict[str, str]:
    sequences: dict[str, str] = {}

    with open(path, encoding="utf-8") as file:
        header: str | None = None
        seq_chunks: list[str] = []

        for raw_line in file:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    sequences[header] = "".join(seq_chunks).upper()

                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if header is not None:
            sequences[header] = "".join(seq_chunks).upper()

    return sequences


def validate_aligned_fasta(sequences: dict[str, str]) -> int:
    if not sequences:
        raise ValueError("No sequences were read from the FASTA file.")

    lengths = {len(seq) for seq in sequences.values()}

    if len(lengths) != 1:
        raise ValueError("Sequences are not all the same length. Input must be an aligned FASTA.")

    return next(iter(lengths))


def split_focal_headers(
    sequences: dict[str, str],
    target_string: str,
) -> tuple[list[str], list[str]]:
    focal_headers = [header for header in sequences if target_string in header]
    non_focal_headers = [header for header in sequences if target_string not in header]

    if not focal_headers:
        raise ValueError("No sequences matched the identifier.")

    if not non_focal_headers:
        raise ValueError("All sequences match the identifier; need contrast set.")

    return focal_headers, non_focal_headers