"""
Every SQL statement in the application lives here.

Nothing else in the codebase — and certainly nothing in Electron — executes
SQL. Callers pass and receive plain dataclasses.
"""

from __future__ import annotations

import sqlite3
import time
import uuid
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

from molecular_diagnosis.project.indexing import FastaScan
from molecular_diagnosis.project.paths import display_name_for, normalise_path_key

__all__ = [
    "FastaFileRow",
    "FocalEntryRow",
    "FocalLocationRow",
    "FocalSetRow",
    "Repository",
    "now_ms",
]


def now_ms() -> int:
    return int(time.time() * 1000)


def new_id() -> str:
    return str(uuid.uuid4())


@dataclass(frozen=True)
class FastaFileRow:
    id: str
    source_path: str
    source_path_key: str
    display_name: str
    sort_order: int
    indexed_size_bytes: int | None
    indexed_mtime_ns: int | None
    indexed_sha256: str | None
    sequence_count: int | None
    alignment_length: int | None
    duplicate_header_count: int
    index_revision: int
    indexed_at_ms: int | None

    @property
    def has_index(self) -> bool:
        return self.indexed_sha256 is not None and self.index_revision > 0


@dataclass(frozen=True)
class FocalSetRow:
    id: str
    title: str
    locked: bool
    sort_order: int


@dataclass(frozen=True)
class FocalEntryRow:
    id: str
    focal_set_id: str
    header: str
    sort_order: int


@dataclass(frozen=True)
class FocalLocationRow:
    focal_entry_id: str
    fasta_file_id: str
    ordinal: int
    record_start_byte: int | None
    record_end_byte: int | None
    last_verified_at_ms: int | None


def _file_row(row: sqlite3.Row) -> FastaFileRow:
    return FastaFileRow(
        id=row["id"],
        source_path=row["source_path"],
        source_path_key=row["source_path_key"],
        display_name=row["display_name"],
        sort_order=row["sort_order"],
        indexed_size_bytes=row["indexed_size_bytes"],
        indexed_mtime_ns=row["indexed_mtime_ns"],
        indexed_sha256=row["indexed_sha256"],
        sequence_count=row["sequence_count"],
        alignment_length=row["alignment_length"],
        duplicate_header_count=row["duplicate_header_count"],
        index_revision=row["index_revision"],
        indexed_at_ms=row["indexed_at_ms"],
    )


class Repository:
    def __init__(self, con: sqlite3.Connection, *, fts_enabled: bool) -> None:
        self.con = con
        self.fts_enabled = fts_enabled

    # ------------------------------------------------------------------
    # Transactions
    # ------------------------------------------------------------------

    def begin(self) -> None:
        """IMMEDIATE so a writer takes the lock up front rather than mid-work."""
        self.con.execute("BEGIN IMMEDIATE")

    def commit(self) -> None:
        self.con.commit()

    def rollback(self) -> None:
        self.con.rollback()

    # ------------------------------------------------------------------
    # Project metadata
    # ------------------------------------------------------------------

    def ensure_project(self, title: str) -> dict[str, Any]:
        row = self.con.execute("SELECT * FROM project_metadata WHERE singleton = 1").fetchone()
        if row is not None:
            return dict(row)

        stamp = now_ms()
        self.con.execute(
            "INSERT INTO project_metadata(singleton, project_uuid, title, created_at_ms,"
            " updated_at_ms) VALUES (1, ?, ?, ?, ?)",
            (new_id(), title, stamp, stamp),
        )
        return dict(self.con.execute("SELECT * FROM project_metadata WHERE singleton = 1").fetchone())

    def set_project_title(self, title: str) -> None:
        self.con.execute(
            "UPDATE project_metadata SET title = ?, updated_at_ms = ? WHERE singleton = 1",
            (title, now_ms()),
        )

    # ------------------------------------------------------------------
    # FASTA files
    # ------------------------------------------------------------------

    def list_fasta_files(self) -> list[FastaFileRow]:
        rows = self.con.execute(
            "SELECT * FROM fasta_file ORDER BY sort_order, id"
        ).fetchall()
        return [_file_row(row) for row in rows]

    def get_fasta_file(self, file_id: str) -> FastaFileRow | None:
        row = self.con.execute("SELECT * FROM fasta_file WHERE id = ?", (file_id,)).fetchone()
        return _file_row(row) if row else None

    def find_by_path(self, path: str) -> FastaFileRow | None:
        row = self.con.execute(
            "SELECT * FROM fasta_file WHERE source_path_key = ?", (normalise_path_key(path),)
        ).fetchone()
        return _file_row(row) if row else None

    def add_fasta_file(self, path: str) -> FastaFileRow:
        """Register a link. Indexing is a separate, verified step."""
        existing = self.find_by_path(path)
        if existing is not None:
            return existing

        stamp = now_ms()
        next_order = self.con.execute(
            "SELECT COALESCE(MAX(sort_order) + 1, 0) FROM fasta_file"
        ).fetchone()[0]
        file_id = new_id()
        self.con.execute(
            "INSERT INTO fasta_file(id, source_path, source_path_key, display_name, sort_order,"
            " duplicate_header_count, index_revision, created_at_ms, updated_at_ms)"
            " VALUES (?, ?, ?, ?, ?, 0, 0, ?, ?)",
            (file_id, str(path), normalise_path_key(path), display_name_for(path),
             next_order, stamp, stamp),
        )
        row = self.get_fasta_file(file_id)
        assert row is not None
        return row

    def remove_fasta_file(self, file_id: str) -> None:
        """Cascades to fasta_record and to that file's focal locations."""
        self.con.execute("DELETE FROM fasta_file WHERE id = ?", (file_id,))

    def update_source_path(self, file_id: str, path: str) -> None:
        self.con.execute(
            "UPDATE fasta_file SET source_path = ?, source_path_key = ?, display_name = ?,"
            " updated_at_ms = ? WHERE id = ?",
            (str(path), normalise_path_key(path), display_name_for(path), now_ms(), file_id),
        )

    def touch_indexed_stat(self, file_id: str, size_bytes: int, mtime_ns: int) -> None:
        """
        Refresh only the cheap change-detector fields.

        Used when the strong hash proved the bytes are unchanged: the index
        still describes the file, so index_revision must NOT move.
        """
        self.con.execute(
            "UPDATE fasta_file SET indexed_size_bytes = ?, indexed_mtime_ns = ?,"
            " updated_at_ms = ? WHERE id = ?",
            (size_bytes, mtime_ns, now_ms(), file_id),
        )

    # ------------------------------------------------------------------
    # Header index
    # ------------------------------------------------------------------

    def replace_index(self, file_id: str, scan: FastaScan, *, bump_revision: bool) -> int:
        """
        Swap in a complete header snapshot for one file, atomically.

        Caller owns the transaction. The whole index is replaced; it is never
        edited record-by-record, because a partial index cannot be reasoned
        about. Focal locations are untouched here by design — they are a
        separate cache with their own verification.
        """
        if self.fts_enabled:
            self.con.execute(
                "DELETE FROM fasta_record_search WHERE rowid IN"
                " (SELECT id FROM fasta_record WHERE fasta_file_id = ?)",
                (file_id,),
            )
        self.con.execute("DELETE FROM fasta_record WHERE fasta_file_id = ?", (file_id,))

        self.con.executemany(
            "INSERT INTO fasta_record(fasta_file_id, ordinal, header, header_casefold,"
            " sequence_length, record_start_byte, record_end_byte)"
            " VALUES (?, ?, ?, ?, ?, ?, ?)",
            [
                (file_id, r.ordinal, r.header, r.header_casefold, r.sequence_length,
                 r.record_start_byte, r.record_end_byte)
                for r in scan.records
            ],
        )

        if self.fts_enabled:
            self.con.execute(
                "INSERT INTO fasta_record_search(rowid, header_casefold)"
                " SELECT id, header_casefold FROM fasta_record WHERE fasta_file_id = ?",
                (file_id,),
            )

        current = self.get_fasta_file(file_id)
        revision = (current.index_revision if current else 0) + (1 if bump_revision else 0)
        if revision == 0:
            revision = 1

        self.con.execute(
            "UPDATE fasta_file SET indexed_size_bytes = ?, indexed_mtime_ns = ?,"
            " indexed_sha256 = ?, sequence_count = ?, alignment_length = ?,"
            " duplicate_header_count = ?, index_revision = ?, indexed_at_ms = ?,"
            " updated_at_ms = ? WHERE id = ?",
            (scan.size_bytes, scan.mtime_ns, scan.sha256, scan.sequence_count,
             scan.alignment_length, scan.duplicate_header_count, revision,
             now_ms(), now_ms(), file_id),
        )
        return revision

    def record_count(self, file_id: str) -> int:
        return self.con.execute(
            "SELECT COUNT(*) FROM fasta_record WHERE fasta_file_id = ?", (file_id,)
        ).fetchone()[0]

    def headers_for_file(self, file_id: str) -> list[tuple[int, str]]:
        return [
            (row["ordinal"], row["header"])
            for row in self.con.execute(
                "SELECT ordinal, header FROM fasta_record WHERE fasta_file_id = ?"
                " ORDER BY ordinal", (file_id,)
            )
        ]

    def locate_header(self, file_id: str, header: str) -> list[sqlite3.Row]:
        """Every occurrence of an exact header in one file."""
        return self.con.execute(
            "SELECT ordinal, record_start_byte, record_end_byte FROM fasta_record"
            " WHERE header = ? AND fasta_file_id = ? ORDER BY ordinal",
            (header, file_id),
        ).fetchall()

    def files_containing_header(self, header: str) -> list[str]:
        return [
            row["fasta_file_id"]
            for row in self.con.execute(
                "SELECT DISTINCT fasta_file_id FROM fasta_record WHERE header = ?", (header,)
            )
        ]

    # ------------------------------------------------------------------
    # Focal sets
    # ------------------------------------------------------------------

    def list_focal_sets(self) -> list[FocalSetRow]:
        return [
            FocalSetRow(row["id"], row["title"], bool(row["locked"]), row["sort_order"])
            for row in self.con.execute("SELECT * FROM focal_set ORDER BY sort_order, id")
        ]

    def get_focal_set(self, set_id: str) -> FocalSetRow | None:
        row = self.con.execute("SELECT * FROM focal_set WHERE id = ?", (set_id,)).fetchone()
        if row is None:
            return None
        return FocalSetRow(row["id"], row["title"], bool(row["locked"]), row["sort_order"])

    def create_focal_set(self, title: str) -> FocalSetRow:
        stamp = now_ms()
        next_order = self.con.execute(
            "SELECT COALESCE(MAX(sort_order) + 1, 0) FROM focal_set"
        ).fetchone()[0]
        set_id = new_id()
        self.con.execute(
            "INSERT INTO focal_set(id, title, locked, sort_order, created_at_ms, updated_at_ms)"
            " VALUES (?, ?, 0, ?, ?, ?)",
            (set_id, title, next_order, stamp, stamp),
        )
        return FocalSetRow(set_id, title, False, next_order)

    def rename_focal_set(self, set_id: str, title: str) -> None:
        self.con.execute(
            "UPDATE focal_set SET title = ?, updated_at_ms = ? WHERE id = ?",
            (title, now_ms(), set_id),
        )

    def set_focal_set_locked(self, set_id: str, locked: bool) -> None:
        self.con.execute(
            "UPDATE focal_set SET locked = ?, updated_at_ms = ? WHERE id = ?",
            (1 if locked else 0, now_ms(), set_id),
        )

    def delete_focal_set(self, set_id: str) -> None:
        self.con.execute("DELETE FROM focal_set WHERE id = ?", (set_id,))

    def list_focal_entries(self, set_id: str) -> list[FocalEntryRow]:
        return [
            FocalEntryRow(row["id"], row["focal_set_id"], row["header"], row["sort_order"])
            for row in self.con.execute(
                "SELECT id, focal_set_id, header, sort_order FROM focal_set_entry"
                " WHERE focal_set_id = ? ORDER BY sort_order, id", (set_id,)
            )
        ]

    def add_focal_entries(self, set_id: str, headers: Sequence[str]) -> list[str]:
        """
        Insert exact headers, ignoring ones already present.

        Returns the headers actually added, in order.
        """
        stamp = now_ms()
        next_order = self.con.execute(
            "SELECT COALESCE(MAX(sort_order) + 1, 0) FROM focal_set_entry WHERE focal_set_id = ?",
            (set_id,),
        ).fetchone()[0]

        existing = {row.header for row in self.list_focal_entries(set_id)}
        added: list[str] = []
        for header in headers:
            if header in existing:
                continue
            existing.add(header)
            self.con.execute(
                "INSERT INTO focal_set_entry(id, focal_set_id, header, sort_order,"
                " created_at_ms, updated_at_ms) VALUES (?, ?, ?, ?, ?, ?)",
                (new_id(), set_id, header, next_order, stamp, stamp),
            )
            next_order += 1
            added.append(header)
        return added

    def replace_focal_entries(self, set_id: str, headers: Sequence[str]) -> dict[str, list[str]]:
        """
        Make the set's explicit contents exactly `headers`, as a DIFF.

        Caller owns the transaction and has already trimmed and deduplicated.
        Headers that are staying keep their existing `focal_set_entry.id`, and
        therefore keep the `focal_entry_location` rows that hang off it — a
        delete-and-recreate would throw away a verified location cache on every
        debounced keystroke. Only the genuinely new and the genuinely gone are
        written; everything else moves at most its `sort_order`, which is set to
        the caller's order so the textbox and the stored set agree.
        """
        existing = {row.header: row for row in self.list_focal_entries(set_id)}
        wanted = list(headers)
        wanted_set = set(wanted)

        removed = [header for header in existing if header not in wanted_set]
        for header in removed:
            self.con.execute(
                "DELETE FROM focal_set_entry WHERE focal_set_id = ? AND header = ?",
                (set_id, header),
            )

        stamp = now_ms()
        added: list[str] = []
        for position, header in enumerate(wanted):
            row = existing.get(header)
            if row is None:
                self.con.execute(
                    "INSERT INTO focal_set_entry(id, focal_set_id, header, sort_order,"
                    " created_at_ms, updated_at_ms) VALUES (?, ?, ?, ?, ?, ?)",
                    (new_id(), set_id, header, position, stamp, stamp),
                )
                added.append(header)
            elif row.sort_order != position:
                self.con.execute(
                    "UPDATE focal_set_entry SET sort_order = ?, updated_at_ms = ? WHERE id = ?",
                    (position, stamp, row.id),
                )

        return {
            "added": added,
            "removed": removed,
            "kept": [header for header in wanted if header in existing],
        }

    def remove_focal_entries(self, set_id: str, headers: Iterable[str]) -> list[str]:
        removed: list[str] = []
        for header in headers:
            cursor = self.con.execute(
                "DELETE FROM focal_set_entry WHERE focal_set_id = ? AND header = ?",
                (set_id, header),
            )
            if cursor.rowcount:
                removed.append(header)
        return removed

    # ------------------------------------------------------------------
    # Focal entry locations (persistent, independent cache)
    # ------------------------------------------------------------------

    def locations_for_set(self, set_id: str) -> list[FocalLocationRow]:
        return [
            FocalLocationRow(
                row["focal_entry_id"], row["fasta_file_id"], row["ordinal"],
                row["record_start_byte"], row["record_end_byte"], row["last_verified_at_ms"],
            )
            for row in self.con.execute(
                "SELECT l.* FROM focal_entry_location l"
                " JOIN focal_set_entry e ON e.id = l.focal_entry_id"
                " WHERE e.focal_set_id = ?", (set_id,)
            )
        ]

    def locations_for_entry(self, entry_id: str) -> list[FocalLocationRow]:
        return [
            FocalLocationRow(
                row["focal_entry_id"], row["fasta_file_id"], row["ordinal"],
                row["record_start_byte"], row["record_end_byte"], row["last_verified_at_ms"],
            )
            for row in self.con.execute(
                "SELECT * FROM focal_entry_location WHERE focal_entry_id = ?", (entry_id,)
            )
        ]

    def replace_entry_locations_in_file(
        self,
        entry_id: str,
        file_id: str,
        placements: Sequence[tuple[int, int | None, int | None]],
    ) -> None:
        """Repair one entry's locations within one file. Empty clears them."""
        self.con.execute(
            "DELETE FROM focal_entry_location WHERE focal_entry_id = ? AND fasta_file_id = ?",
            (entry_id, file_id),
        )
        stamp = now_ms()
        self.con.executemany(
            "INSERT INTO focal_entry_location(focal_entry_id, fasta_file_id, ordinal,"
            " record_start_byte, record_end_byte, last_verified_at_ms) VALUES (?, ?, ?, ?, ?, ?)",
            [(entry_id, file_id, ordinal, start, end, stamp) for ordinal, start, end in placements],
        )

    def mark_location_verified(self, entry_id: str, file_id: str, ordinal: int) -> None:
        self.con.execute(
            "UPDATE focal_entry_location SET last_verified_at_ms = ?"
            " WHERE focal_entry_id = ? AND fasta_file_id = ? AND ordinal = ?",
            (now_ms(), entry_id, file_id, ordinal),
        )

    def delete_locations_in_file(self, file_id: str) -> None:
        self.con.execute("DELETE FROM focal_entry_location WHERE fasta_file_id = ?", (file_id,))
