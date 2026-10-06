"""
Focal-entry location verification and repair.

`focal_entry_location` is a persistent cache that outlives both application
sessions and header-index rebuilds. It is NOT discarded merely because the
FASTA changed.

Verification is local first, and only escalates when it has to:

1. A saved location carries a byte offset -> seek there, read that record's
   header, compare it exactly with the focal entry's stored header. A match
   keeps the location, stamps `last_verified_at_ms`, and reads nothing else.
   The whole-file hash is irrelevant to this: an edit elsewhere in the file
   does not move this record.
2. Anything that fails direct verification is repaired from ONE scan of the
   current file, shared across every entry that needs it. Found -> the
   location, ordinal and offsets are updated permanently. Not found -> that
   entry/file relationship is removed.

Duplicate identical headers in one file are kept as separate locations.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from molecular_diagnosis.project.indexing import FastaScan, read_header_at
from molecular_diagnosis.project.repository import FastaFileRow, FocalEntryRow, Repository

__all__ = ["EntryPresence", "PresenceState", "verify_entry_locations"]


class PresenceState(str, Enum):
    """Per §10. The UI maps these to colours; the colours are never stored."""

    PRESENT_CURRENT = "present_current"
    PRESENT_OTHER = "present_other"
    MISSING = "missing"
    #: The relevant source is unavailable, so presence cannot be asserted.
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class EntryPresence:
    entry_id: str
    header: str
    state: PresenceState
    #: file id -> number of occurrences currently verified in that file
    occurrences: dict[str, int]

    def to_payload(self) -> dict[str, object]:
        return {
            "entryId": self.entry_id,
            "header": self.header,
            "state": self.state.value,
            "occurrences": dict(self.occurrences),
        }


def verify_entry_locations(
    repository: Repository,
    file_row: FastaFileRow,
    entries: list[FocalEntryRow],
    *,
    scan: FastaScan | None = None,
) -> dict[str, list[int]]:
    """
    Verify and repair every entry's locations within ONE file.

    Returns entry_id -> ordinals currently holding that exact header.

    `scan` lets a caller that already read the file (a reindex, an analysis)
    donate its result instead of paying for another read.
    """
    if not entries:
        return {}

    saved: dict[str, list] = {}
    for entry in entries:
        saved[entry.id] = [
            location
            for location in repository.locations_for_entry(entry.id)
            if location.fasta_file_id == file_row.id
        ]

    resolved: dict[str, list[int]] = {}
    needs_search: list[FocalEntryRow] = []

    if scan is None:
        # Try each saved location directly against the current file first.
        for entry in entries:
            locations = saved[entry.id]
            if not locations:
                needs_search.append(entry)
                continue

            verified: list[int] = []
            all_ok = True
            for location in locations:
                if location.record_start_byte is None:
                    all_ok = False
                    break
                found = read_header_at(file_row.source_path, location.record_start_byte)
                if found == entry.header:
                    verified.append(location.ordinal)
                else:
                    all_ok = False
                    break

            if all_ok and verified:
                for ordinal in verified:
                    repository.mark_location_verified(entry.id, file_row.id, ordinal)
                resolved[entry.id] = sorted(verified)
            else:
                needs_search.append(entry)
    else:
        needs_search = list(entries)

    if not needs_search:
        return resolved

    # One shared read repairs everything that failed direct verification.
    if scan is None:
        try:
            from molecular_diagnosis.project.indexing import scan_fasta

            scan = scan_fasta(file_row.source_path)
        except OSError:
            # Source unavailable: leave the persisted rows alone. They are not
            # currently usable, but they are not proven wrong either.
            return resolved

    placements: dict[str, list[tuple[int, int | None, int | None]]] = {}
    for record in scan.records:
        placements.setdefault(record.header, []).append(
            (record.ordinal, record.record_start_byte, record.record_end_byte)
        )

    for entry in needs_search:
        found = placements.get(entry.header, [])
        repository.replace_entry_locations_in_file(entry.id, file_row.id, found)
        if found:
            resolved[entry.id] = sorted(ordinal for ordinal, _s, _e in found)
        else:
            resolved.pop(entry.id, None)

    return resolved


def presence_for_entries(
    entries: list[FocalEntryRow],
    per_file: dict[str, dict[str, list[int]]],
    *,
    selected_file_id: str | None,
    unavailable_files: set[str],
) -> list[EntryPresence]:
    """
    Fold per-file verification results into the §10 presence states.

    `selected_file_id=None` means the "All files" scope, which has no orange
    state: an entry either exists somewhere or it does not.
    """
    results: list[EntryPresence] = []

    for entry in entries:
        occurrences = {
            file_id: len(resolved[entry.id])
            for file_id, resolved in per_file.items()
            if entry.id in resolved and resolved[entry.id]
        }

        if selected_file_id is None:
            if occurrences:
                state = PresenceState.PRESENT_CURRENT
            elif unavailable_files:
                state = PresenceState.UNKNOWN
            else:
                state = PresenceState.MISSING
        elif selected_file_id in unavailable_files:
            state = PresenceState.UNKNOWN
        elif occurrences.get(selected_file_id):
            state = PresenceState.PRESENT_CURRENT
        elif occurrences:
            state = PresenceState.PRESENT_OTHER
        elif unavailable_files:
            state = PresenceState.UNKNOWN
        else:
            state = PresenceState.MISSING

        results.append(EntryPresence(entry.id, entry.header, state, occurrences))

    return results
