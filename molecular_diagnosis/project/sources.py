"""
Source-file verification lifecycle.

The linked FASTA is always authoritative. `fasta_record`, byte offsets,
ordinals and `focal_entry_location` are derived caches describing exactly the
bytes fingerprinted by `fasta_file.indexed_*`.

Two tiers, deliberately:

* CHEAP  — `os.stat`: size + mtime_ns. Fast enough to run before every
  operation. It detects change; it does not prove identity.
* STRONG — SHA-256 of the raw bytes. The real identity, and the only thing
  that can tell "restored to the original size and mtime but edited" from
  "genuinely untouched".

Every file is strongly verified at least once per service session before its
indexed contents are relied on. The result is cached in memory against the
stat it was observed with, so typing in the focal box does not re-hash a
multi-megabyte alignment on every keystroke. Any stat change invalidates that
session verification.

Nothing here persists a current/missing/changed status: availability is
runtime state, recomputed from the filesystem.
"""

from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from molecular_diagnosis.project.indexing import (
    FastaScan,
    SourceChangedDuringRead,
    scan_fasta,
    stat_signature,
)
from molecular_diagnosis.project.repository import FastaFileRow, Repository

__all__ = ["SourceState", "SourceStatus", "SourceVerifier"]

CHUNK = 1 << 20


class SourceState(str, Enum):
    """Runtime state of a linked source. Never persisted."""

    #: The path does not exist. The row stays; cached data must not be used.
    MISSING = "missing"
    #: Exists but is not a readable regular file.
    UNREADABLE = "unreadable"
    #: Linked but never indexed.
    NEVER_INDEXED = "never_indexed"
    #: Stat matches the indexed fingerprint, but not strongly verified yet
    #: this session. Show it as provisional, not confidently current.
    UNVERIFIED = "unverified"
    #: Strongly verified this session and unchanged since. Index is usable.
    CURRENT = "current"
    #: The current bytes differ from the indexed snapshot. The complete header
    #: index is stale and must be rebuilt before it can be trusted.
    STALE = "stale"


@dataclass(frozen=True)
class SourceStatus:
    fasta_file_id: str
    source_path: str
    display_name: str
    state: SourceState
    exists: bool
    current_size_bytes: int | None
    current_mtime_ns: int | None
    indexed_size_bytes: int | None
    indexed_mtime_ns: int | None
    index_revision: int
    sequence_count: int | None
    alignment_length: int | None
    duplicate_header_count: int
    message: str | None = None

    @property
    def available(self) -> bool:
        """Readable right now. Says nothing about whether the index matches."""
        return self.state not in (SourceState.MISSING, SourceState.UNREADABLE)

    @property
    def index_usable(self) -> bool:
        """May the complete header index be trusted for arbitrary searches?"""
        return self.state is SourceState.CURRENT

    def to_payload(self) -> dict[str, object]:
        return {
            "fastaFileId": self.fasta_file_id,
            "sourcePath": self.source_path,
            "displayName": self.display_name,
            "state": self.state.value,
            "available": self.available,
            "indexUsable": self.index_usable,
            "exists": self.exists,
            "currentSizeBytes": self.current_size_bytes,
            "currentMtimeNs": self.current_mtime_ns,
            "indexedSizeBytes": self.indexed_size_bytes,
            "indexedMtimeNs": self.indexed_mtime_ns,
            "indexRevision": self.index_revision,
            "sequenceCount": self.sequence_count,
            "alignmentLength": self.alignment_length,
            "duplicateHeaderCount": self.duplicate_header_count,
            "message": self.message,
        }


@dataclass
class _SessionVerification:
    """A strong verification, valid only while the stat still matches."""

    size_bytes: int
    mtime_ns: int
    sha256: str
    matched_index: bool


def hash_file(path: str | Path) -> tuple[str, int, int]:
    """SHA-256 plus the stat observed around the read. Raises on a racing write."""
    with open(path, "rb") as handle:
        before = os.fstat(handle.fileno())
        digest = hashlib.sha256()
        while True:
            chunk = handle.read(CHUNK)
            if not chunk:
                break
            digest.update(chunk)
        after = os.fstat(handle.fileno())

    if stat_signature(before) != stat_signature(after):
        raise SourceChangedDuringRead(f"{path} changed while it was being hashed.")
    return digest.hexdigest(), after.st_size, after.st_mtime_ns


class SourceVerifier:
    """
    The one place that decides whether a linked FASTA is currently usable.

    Every consumer — search, focal presence, analysis — goes through
    `ensure_current`, so the rules cannot drift between call sites.
    """

    def __init__(self, repository: Repository) -> None:
        self.repository = repository
        self._verified: dict[str, _SessionVerification] = {}

    # ------------------------------------------------------------------

    def forget(self, file_id: str) -> None:
        """Drop the session verification, e.g. after a reindex or relink."""
        self._verified.pop(file_id, None)

    def forget_all(self) -> None:
        self._verified.clear()

    def session_verified(self, file_id: str) -> bool:
        return file_id in self._verified

    # ------------------------------------------------------------------

    @staticmethod
    def _stat(path: str) -> os.stat_result | None:
        try:
            stat = os.stat(path)
        except OSError:
            return None
        return stat

    def cheap_status(self, row: FastaFileRow) -> SourceStatus:
        """
        Stat-only status. Never hashes, so it is safe on every keystroke and at
        every lifecycle point.
        """
        stat = self._stat(row.source_path)

        if stat is None:
            self.forget(row.id)
            return self._status(row, SourceState.MISSING, None, None,
                                message="The linked FASTA file was not found.")

        if not os.path.isfile(row.source_path) or not os.access(row.source_path, os.R_OK):
            self.forget(row.id)
            return self._status(row, SourceState.UNREADABLE, stat.st_size, stat.st_mtime_ns,
                                message="The linked path is not a readable file.")

        size, mtime = stat.st_size, stat.st_mtime_ns

        if not row.has_index:
            return self._status(row, SourceState.NEVER_INDEXED, size, mtime)

        cached = self._verified.get(row.id)
        if cached is not None:
            if (cached.size_bytes, cached.mtime_ns) == (size, mtime):
                # Verified this session and untouched since.
                return self._status(
                    row,
                    SourceState.CURRENT if cached.matched_index else SourceState.STALE,
                    size, mtime,
                )
            # Stat moved: the session verification no longer describes the file.
            self.forget(row.id)

        if (row.indexed_size_bytes, row.indexed_mtime_ns) == (size, mtime):
            # Looks unchanged, but cheap evidence only.
            return self._status(row, SourceState.UNVERIFIED, size, mtime)

        return self._status(row, SourceState.STALE, size, mtime,
                            message="The file has changed since it was indexed.")

    def ensure_current(self, row: FastaFileRow, *, strong: bool = True) -> SourceStatus:
        """
        The single verification entry point.

        With `strong=True` (the default before any indexed data is relied on),
        an UNVERIFIED file is hashed once and the result cached for the session.
        With `strong=False` the caller accepts UNVERIFIED, which is what the
        UI uses for cheap refreshes.
        """
        status = self.cheap_status(row)

        if not strong:
            return status
        if self.session_verified(row.id):
            # Already proven this session, either way.
            return status
        if status.state not in (SourceState.UNVERIFIED, SourceState.STALE):
            return status

        # A differing stat is only a SUSPICION of staleness. Hash before
        # believing it: a file that was touched, copied, or restored keeps the
        # same bytes, and reindexing it would be pure waste.

        try:
            digest, size, mtime = hash_file(row.source_path)
        except SourceChangedDuringRead as error:
            return self._status(row, SourceState.STALE, status.current_size_bytes,
                                status.current_mtime_ns, message=str(error))
        except OSError as error:
            return self._status(row, SourceState.UNREADABLE, None, None, message=str(error))

        matched = digest == row.indexed_sha256
        self._verified[row.id] = _SessionVerification(size, mtime, digest, matched)

        if matched:
            # Same bytes, new stat: refresh the cheap detector only. The index
            # still describes the file, so index_revision must not move.
            if (row.indexed_size_bytes, row.indexed_mtime_ns) != (size, mtime):
                self.repository.touch_indexed_stat(row.id, size, mtime)
            return self._status(row, SourceState.CURRENT, size, mtime)

        return self._status(row, SourceState.STALE, size, mtime,
                            message="The file's contents differ from the indexed version.")

    def note_scan(self, file_id: str, scan: FastaScan, *, matched_index: bool = True) -> None:
        """Record the verification a full scan already produced, for free."""
        self._verified[file_id] = _SessionVerification(
            scan.size_bytes, scan.mtime_ns, scan.sha256, matched_index
        )

    def scan_current(self, row: FastaFileRow) -> FastaScan:
        """
        Read the file completely, computing its hash in the same pass.

        Used when the bytes are needed anyway (reindex, analysis), so
        verification costs nothing extra.
        """
        return scan_fasta(row.source_path)

    # ------------------------------------------------------------------

    def _status(
        self,
        row: FastaFileRow,
        state: SourceState,
        size: int | None,
        mtime: int | None,
        *,
        message: str | None = None,
    ) -> SourceStatus:
        return SourceStatus(
            fasta_file_id=row.id,
            source_path=row.source_path,
            display_name=row.display_name,
            state=state,
            exists=size is not None,
            current_size_bytes=size,
            current_mtime_ns=mtime,
            indexed_size_bytes=row.indexed_size_bytes,
            indexed_mtime_ns=row.indexed_mtime_ns,
            index_revision=row.index_revision,
            sequence_count=row.sequence_count,
            alignment_length=row.alignment_length,
            duplicate_header_count=row.duplicate_header_count,
            message=message,
        )
