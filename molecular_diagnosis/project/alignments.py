"""
In-memory parsed-alignment cache.

SQLite holds headers and record positions, not sequences. Molecular Diagnosis
needs every residue of every sequence, so the alignment still has to be parsed;
what this cache avoids is parsing the SAME unchanged file again for the next
operation in the same session.

The key is (fasta file id, index_revision, sha256): a content reindex changes
the revision, and any edit changes the hash, so a stale entry cannot be
selected. Before an entry is handed out the caller must still have established
that the source is current — the cache never decides that for itself.

Deliberately NOT persistent. A binary/memory-mapped on-disk sequence cache is
a plausible later optimisation, but only once profiling shows parsing is
actually the bottleneck.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass

from molecular_diagnosis.project.indexing import FastaScan

__all__ = ["AlignmentCache", "CachedAlignment"]

#: Alignments are large; keep only a handful of the most recently used.
DEFAULT_CAPACITY = 8


@dataclass(frozen=True)
class CachedAlignment:
    fasta_file_id: str
    index_revision: int
    sha256: str
    sequences: dict[str, str]
    alignment_length: int | None
    duplicate_header_count: int

    @property
    def headers(self) -> list[str]:
        return list(self.sequences)


class AlignmentCache:
    def __init__(self, capacity: int = DEFAULT_CAPACITY) -> None:
        self._entries: OrderedDict[tuple[str, int, str], CachedAlignment] = OrderedDict()
        self._capacity = capacity
        self.hits = 0
        self.misses = 0

    @staticmethod
    def _key(file_id: str, revision: int, sha256: str) -> tuple[str, int, str]:
        return (file_id, revision, sha256)

    def get(self, file_id: str, revision: int, sha256: str | None) -> CachedAlignment | None:
        if sha256 is None:
            self.misses += 1
            return None
        key = self._key(file_id, revision, sha256)
        entry = self._entries.get(key)
        if entry is None:
            self.misses += 1
            return None
        self._entries.move_to_end(key)
        self.hits += 1
        return entry

    def put(self, file_id: str, revision: int, scan: FastaScan) -> CachedAlignment:
        entry = CachedAlignment(
            fasta_file_id=file_id,
            index_revision=revision,
            sha256=scan.sha256,
            sequences=dict(scan.sequences),
            alignment_length=scan.alignment_length,
            duplicate_header_count=scan.duplicate_header_count,
        )
        key = self._key(file_id, revision, scan.sha256)
        self._entries[key] = entry
        self._entries.move_to_end(key)
        while len(self._entries) > self._capacity:
            self._entries.popitem(last=False)
        return entry

    def invalidate(self, file_id: str) -> None:
        """Drop every entry for a file, whatever revision or hash it had."""
        for key in [key for key in self._entries if key[0] == file_id]:
            del self._entries[key]

    def clear(self) -> None:
        self._entries.clear()

    def __len__(self) -> int:
        return len(self._entries)
