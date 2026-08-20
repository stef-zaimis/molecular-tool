"""
Byte-stable FASTA scanning.

One open file, one read: the SHA-256, the record table (with byte offsets) and
the parsed sequences all come out of the same pass, so a scientific run never
pays for a second traversal just to fingerprint the file.

Parsing semantics are deliberately identical to `fasta_io.parse_fasta`:
a header line starts with '>', the header is the rest of that line stripped,
sequence lines are concatenated and upper-cased, and blank lines are ignored.
The difference is that this module keeps EVERY record, including repeated
headers, because `fasta_record` is a positional index of the file rather than a
dictionary of it.
"""

from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "FastaScan",
    "ScannedRecord",
    "SourceChangedDuringRead",
    "read_header_at",
    "scan_fasta",
    "stat_signature",
]

CHUNK = 1 << 20


class SourceChangedDuringRead(RuntimeError):
    """The file's metadata moved between the start and end of the read."""


@dataclass(frozen=True)
class ScannedRecord:
    ordinal: int
    header: str
    header_casefold: str
    sequence_length: int
    record_start_byte: int
    record_end_byte: int


@dataclass
class FastaScan:
    """Everything one complete, stable pass over a FASTA yields."""

    path: str
    size_bytes: int
    mtime_ns: int
    sha256: str
    records: list[ScannedRecord] = field(default_factory=list)
    #: header -> sequence, collapsing duplicates exactly as parse_fasta does.
    sequences: dict[str, str] = field(default_factory=dict)

    @property
    def sequence_count(self) -> int:
        """Distinct headers — what the analysis will actually see."""
        return len(self.sequences)

    @property
    def record_count(self) -> int:
        return len(self.records)

    @property
    def duplicate_header_count(self) -> int:
        return len(self.records) - len(self.sequences)

    @property
    def alignment_length(self) -> int | None:
        lengths = {len(value) for value in self.sequences.values()}
        return next(iter(lengths)) if len(lengths) == 1 else None

    @property
    def is_aligned(self) -> bool:
        return bool(self.sequences) and self.alignment_length is not None


def stat_signature(stat: os.stat_result) -> tuple[int, int, int, int]:
    """
    The cheap change detector: size, mtime, and identity when available.

    st_ino/st_dev are 0 on some Windows filesystems; comparing them anyway is
    harmless because both sides come from the same platform.
    """
    return (stat.st_size, stat.st_mtime_ns, getattr(stat, "st_ino", 0), getattr(stat, "st_dev", 0))


def _decode(raw: bytes) -> str:
    # Mirrors text-mode reading with encoding="utf-8"; surrogateescape keeps a
    # malformed byte round-trippable instead of exploding mid-index.
    return raw.decode("utf-8", errors="surrogateescape")


def scan_fasta(path: str | Path) -> FastaScan:
    """
    Read a FASTA once, producing its fingerprint, record index and sequences.

    Raises `SourceChangedDuringRead` if the file's stat differs before and
    after the read, so a snapshot assembled from a file that was being written
    is never committed. Callers retry or fail; they must not use the result.
    """
    file_path = Path(path)

    with open(file_path, "rb") as handle:
        before = os.fstat(handle.fileno())

        digest = hashlib.sha256()
        records: list[ScannedRecord] = []
        sequences: dict[str, str] = {}

        pending_header: str | None = None
        pending_start = 0
        pending_chunks: list[str] = []
        last_content_end = 0

        offset = 0
        carry = b""

        def close_record(end_byte: int) -> None:
            nonlocal pending_header, pending_chunks
            if pending_header is None:
                return
            sequence = "".join(pending_chunks).upper()
            header_casefold = pending_header.casefold()
            records.append(
                ScannedRecord(
                    ordinal=len(records),
                    header=pending_header,
                    header_casefold=header_casefold,
                    sequence_length=len(sequence),
                    record_start_byte=pending_start,
                    record_end_byte=max(pending_start, end_byte),
                )
            )
            # Last occurrence wins, matching parse_fasta's dict assignment.
            sequences[pending_header] = sequence
            pending_header = None
            pending_chunks = []

        while True:
            chunk = handle.read(CHUNK)
            if not chunk:
                break
            digest.update(chunk)

            buffer = carry + chunk
            carry = b""
            start = 0
            while True:
                newline = buffer.find(b"\n", start)
                if newline == -1:
                    carry = buffer[start:]
                    break
                line = buffer[start:newline]
                line_start = offset + start
                start = newline + 1

                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith(b">"):
                    close_record(last_content_end)
                    pending_header = _decode(stripped[1:]).strip()
                    pending_start = line_start
                else:
                    pending_chunks.append(_decode(stripped))
                last_content_end = line_start + len(line.rstrip(b"\r\n")) - 1
                if last_content_end < line_start:
                    last_content_end = line_start

            offset += len(buffer) - len(carry)

        if carry:
            stripped = carry.strip()
            line_start = offset
            if stripped:
                if stripped.startswith(b">"):
                    close_record(last_content_end)
                    pending_header = _decode(stripped[1:]).strip()
                    pending_start = line_start
                else:
                    pending_chunks.append(_decode(stripped))
                last_content_end = line_start + len(carry.rstrip(b"\r\n")) - 1

        close_record(last_content_end)

        after = os.fstat(handle.fileno())

    if stat_signature(before) != stat_signature(after):
        raise SourceChangedDuringRead(
            f"{file_path} changed while it was being read; the index was discarded."
        )

    return FastaScan(
        path=str(file_path),
        size_bytes=after.st_size,
        mtime_ns=after.st_mtime_ns,
        sha256=digest.hexdigest(),
        records=records,
        sequences=sequences,
    )


def read_header_at(path: str | Path, offset: int) -> str | None:
    """
    Read the FASTA header that starts at `offset`, or None if there is not one.

    Used to confirm a cached byte offset still points at the record it claims
    to, before anything is displayed or used from it.
    """
    if offset < 0:
        return None
    try:
        with open(path, "rb") as handle:
            handle.seek(offset)
            line = handle.readline()
    except OSError:
        return None

    stripped = line.strip()
    if not stripped.startswith(b">"):
        return None
    return _decode(stripped[1:]).strip()
