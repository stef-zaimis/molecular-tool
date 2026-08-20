"""Path normalisation for `fasta_file.source_path_key`."""

from __future__ import annotations

import os
import unicodedata
from pathlib import Path

__all__ = ["normalise_path_key", "display_name_for"]


def normalise_path_key(path: str | os.PathLike[str]) -> str:
    """
    A stable key for "is this the same file we already linked?".

    Absolute, symlink-resolved where possible, forward-slashed, NFC-normalised,
    and case-folded only on case-insensitive filesystems. Case folding is
    deliberately platform-dependent: two paths differing only in case are the
    same file on Windows and macOS but different files on Linux, and folding
    unconditionally would make the UNIQUE constraint reject legitimate files.
    """
    candidate = Path(path)
    try:
        resolved = candidate.resolve(strict=False)
    except OSError:
        resolved = candidate.absolute()

    text = unicodedata.normalize("NFC", str(resolved)).replace("\\", "/")

    # os.path.normcase folds case on Windows and is a no-op on POSIX.
    return os.path.normcase(text)


def display_name_for(path: str | os.PathLike[str]) -> str:
    name = Path(path).name
    return name or str(path)
