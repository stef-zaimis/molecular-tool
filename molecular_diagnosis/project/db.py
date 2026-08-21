"""
Connection setup, capability detection and `PRAGMA user_version` migrations.

Connection-level policy lives here and NOT in the migration SQL, because those
settings are per-connection: a migration file that set them would only affect
the connection that happened to run the migration.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "Capabilities",
    "MIGRATIONS_DIR",
    "SCHEMA_VERSION",
    "connect",
    "detect_capabilities",
    "install_fts",
    "migrate",
    "open_project_db",
    "optimize",
]

MIGRATIONS_DIR = Path(__file__).resolve().parents[2] / "db"

#: Highest migration this build knows how to apply.
SCHEMA_VERSION = 2

#: Ordered migrations. The index is the version the migration produces.
MIGRATIONS: tuple[tuple[int, str], ...] = (
    (1, "001_initial.sql"),
    (2, "002_fasta_file_locked.sql"),
)

FTS_SCRIPT = "optional_fts5_trigram.sql"

#: Long enough to ride out another connection's write, short enough that a
#: genuinely stuck lock surfaces as an error instead of hanging the UI.
BUSY_TIMEOUT_MS = 5_000


@dataclass(frozen=True)
class Capabilities:
    """What the SQLite actually linked into this Python can do."""

    sqlite_version: str
    fts5: bool
    trigram: bool

    @property
    def accelerated_search(self) -> bool:
        """FTS is usable only if both the module and the tokenizer are present."""
        return self.fts5 and self.trigram


def connect(path: str | Path) -> sqlite3.Connection:
    """
    Open a project database with the connection policy this app requires.

    * foreign_keys ON — the schema relies on ON DELETE CASCADE, and SQLite
      leaves enforcement off by default.
    * trusted_schema OFF — the database is a user-supplied file; this stops
      functions/virtual tables named inside schema structures from being run
      implicitly.
    * busy_timeout — the Electron main process and a test may both touch the
      file; without it a concurrent write raises immediately.

    WAL is deliberately NOT forced here: it is a persistent database property,
    it changes the on-disk file set, and nothing in this pass needs it.
    """
    con = sqlite3.connect(str(path), isolation_level=None)
    con.row_factory = sqlite3.Row

    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    # Must be outside a transaction to take effect; isolation_level=None keeps
    # us in autocommit until an explicit BEGIN.
    con.execute("PRAGMA foreign_keys = ON")
    try:
        con.execute("PRAGMA trusted_schema = OFF")
    except sqlite3.Error:
        # Older SQLite without the pragma: the rest of the policy still holds.
        pass

    return con


def detect_capabilities(con: sqlite3.Connection | None = None) -> Capabilities:
    """
    Probe FTS5 and the trigram tokenizer by actually creating a table.

    Compile options are not sufficient: a build can report ENABLE_FTS5 while
    lacking the tokenizer, and the tokenizer is what the accelerator needs.
    """
    probe = sqlite3.connect(":memory:")
    try:
        fts5 = False
        trigram = False
        try:
            probe.execute("CREATE VIRTUAL TABLE __probe_fts USING fts5(x)")
            fts5 = True
        except sqlite3.Error:
            return Capabilities(sqlite3.sqlite_version, False, False)

        try:
            probe.execute(
                "CREATE VIRTUAL TABLE __probe_tri USING fts5(x, tokenize='trigram case_sensitive 0')"
            )
            probe.execute("INSERT INTO __probe_tri(rowid, x) VALUES (1, 'abcdef')")
            hit = probe.execute('SELECT rowid FROM __probe_tri WHERE x MATCH \'"cde"\'').fetchall()
            trigram = hit == [(1,)]
        except sqlite3.Error:
            trigram = False

        return Capabilities(sqlite3.sqlite_version, fts5, trigram)
    finally:
        probe.close()


def _read_script(name: str) -> str:
    return (MIGRATIONS_DIR / name).read_text(encoding="utf-8")


def migrate(con: sqlite3.Connection) -> int:
    """
    Bring the database up to `SCHEMA_VERSION`, applying only what is missing.

    The monolithic schema is never re-run on an already-migrated database; that
    is the entire point of tracking `PRAGMA user_version`.

    Returns the resulting version. Raises if the database is NEWER than this
    build understands, because silently opening a future schema is how data
    gets corrupted.
    """
    version = int(con.execute("PRAGMA user_version").fetchone()[0])

    if version > SCHEMA_VERSION:
        raise RuntimeError(
            f"This project database is version {version}, but this build only "
            f"understands up to {SCHEMA_VERSION}. Update the application."
        )

    for target, script in MIGRATIONS:
        if version >= target:
            continue
        sql = _read_script(script)
        try:
            # executescript() commits any open transaction first; each migration
            # file carries its own BEGIN IMMEDIATE/COMMIT so it is atomic.
            con.executescript(sql)
        except sqlite3.Error:
            if con.in_transaction:
                con.rollback()
            raise
        version = int(con.execute("PRAGMA user_version").fetchone()[0])
        if version != target:
            raise RuntimeError(
                f"Migration {script} finished with user_version={version}, expected {target}."
            )

    return version


def install_fts(con: sqlite3.Connection, capabilities: Capabilities) -> bool:
    """
    Attach the optional FTS accelerator. No-op when unsupported.

    Returns whether the FTS table is present afterwards.
    """
    if not capabilities.accelerated_search:
        return False
    try:
        con.executescript(_read_script(FTS_SCRIPT))
    except sqlite3.Error:
        return False
    return fts_available(con)


def fts_available(con: sqlite3.Connection) -> bool:
    row = con.execute(
        "SELECT 1 FROM sqlite_schema WHERE type='table' AND name='fasta_record_search'"
    ).fetchone()
    return row is not None


def optimize(con: sqlite3.Connection) -> None:
    """
    SQLite's own post-change tuning hook.

    Cheap when there is nothing to do, so it is safe to call after schema or
    bulk index changes and before closing.
    """
    try:
        con.execute("PRAGMA optimize")
    except sqlite3.Error:
        pass


def open_project_db(path: str | Path) -> tuple[sqlite3.Connection, Capabilities, bool]:
    """Open (creating if needed), migrate, and attach FTS when supported."""
    con = connect(path)
    capabilities = detect_capabilities(con)
    migrate(con)
    fts = install_fts(con, capabilities)
    optimize(con)
    return con, capabilities, fts
