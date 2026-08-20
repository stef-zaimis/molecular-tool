"""
Schema, migration and query-plan tests for the project database.

These pin the design decisions the audit produced, so a future edit to the SQL
that reintroduces a coupling or loses an index fails here rather than in the
field.
"""

from __future__ import annotations

import sqlite3
import time
import uuid

import pytest

from molecular_diagnosis.project.db import (
    SCHEMA_VERSION,
    Capabilities,
    connect,
    detect_capabilities,
    fts_available,
    install_fts,
    migrate,
    open_project_db,
)


@pytest.fixture()
def db_path(tmp_path):
    return tmp_path / "project.sqlite"


def stamp() -> int:
    return int(time.time() * 1000)


def seed(con: sqlite3.Connection) -> tuple[str, str, str]:
    """A file with a duplicated header, a focal set, and one entry."""
    now = stamp()
    con.execute(
        "INSERT INTO project_metadata(singleton, project_uuid, title, created_at_ms, updated_at_ms)"
        " VALUES (1, ?, 'T', ?, ?)", (str(uuid.uuid4()), now, now))

    file_id = str(uuid.uuid4())
    con.execute(
        "INSERT INTO fasta_file(id, source_path, source_path_key, display_name, sort_order,"
        " indexed_size_bytes, indexed_mtime_ns, indexed_sha256, sequence_count, alignment_length,"
        " duplicate_header_count, index_revision, indexed_at_ms, created_at_ms, updated_at_ms)"
        " VALUES (?,?,?,?,0,10,1,?,2,4,1,1,?,?,?)",
        (file_id, "C:/x/a.fasta", "c:/x/a.fasta", "a.fasta", "a" * 64, now, now, now))

    for ordinal, header in enumerate(["dup|H", "dup|H", "solo|H"]):
        con.execute(
            "INSERT INTO fasta_record(fasta_file_id, ordinal, header, header_casefold,"
            " sequence_length, record_start_byte, record_end_byte) VALUES (?,?,?,?,4,?,?)",
            (file_id, ordinal, header, header.casefold(), ordinal * 10, ordinal * 10 + 9))

    set_id = str(uuid.uuid4())
    con.execute(
        "INSERT INTO focal_set(id,title,locked,sort_order,created_at_ms,updated_at_ms)"
        " VALUES (?,'S',0,0,?,?)", (set_id, now, now))
    entry_id = str(uuid.uuid4())
    con.execute(
        "INSERT INTO focal_set_entry(id,focal_set_id,header,sort_order,created_at_ms,updated_at_ms)"
        " VALUES (?,?,'dup|H',0,?,?)", (entry_id, set_id, now, now))
    for ordinal in (0, 1):
        con.execute(
            "INSERT INTO focal_entry_location(focal_entry_id, fasta_file_id, ordinal,"
            " record_start_byte, record_end_byte, last_verified_at_ms) VALUES (?,?,?,?,?,?)",
            (entry_id, file_id, ordinal, ordinal * 10, ordinal * 10 + 9, now))
    con.commit()
    return file_id, set_id, entry_id


# ---------------------------------------------------------------------------
# Connection policy
# ---------------------------------------------------------------------------


def test_connection_enables_foreign_keys_and_locks_down_schema(db_path):
    con = connect(db_path)
    assert con.execute("PRAGMA foreign_keys").fetchone()[0] == 1
    assert con.execute("PRAGMA trusted_schema").fetchone()[0] == 0
    assert con.execute("PRAGMA busy_timeout").fetchone()[0] >= 1000
    con.close()


def test_wal_is_not_forced(db_path):
    # WAL is a persistent database property; this pass deliberately leaves it.
    con, _caps, _fts = open_project_db(db_path)
    assert con.execute("PRAGMA journal_mode").fetchone()[0].lower() != "wal"
    con.close()


# ---------------------------------------------------------------------------
# Capabilities
# ---------------------------------------------------------------------------


def test_capabilities_are_probed_not_assumed():
    capabilities = detect_capabilities()
    assert capabilities.sqlite_version == sqlite3.sqlite_version
    # Whatever the answer, the two flags must be coherent.
    if capabilities.trigram:
        assert capabilities.fts5
    assert capabilities.accelerated_search == (capabilities.fts5 and capabilities.trigram)


def test_fts_is_not_installed_when_unsupported(db_path):
    con = connect(db_path)
    migrate(con)
    assert install_fts(con, Capabilities("0.0.0", False, False)) is False
    assert fts_available(con) is False
    con.close()


# ---------------------------------------------------------------------------
# Migration
# ---------------------------------------------------------------------------


def test_initial_migration_creates_the_schema(db_path):
    con = connect(db_path)
    assert con.execute("PRAGMA user_version").fetchone()[0] == 0
    assert migrate(con) == SCHEMA_VERSION
    tables = {
        row[0]
        for row in con.execute("SELECT name FROM sqlite_schema WHERE type='table'")
    }
    assert {
        "project_metadata", "fasta_file", "fasta_record",
        "focal_set", "focal_set_entry", "focal_entry_location",
    } <= tables
    con.close()


def test_reopening_does_not_rerun_the_migration(db_path):
    con, _caps, _fts = open_project_db(db_path)
    file_id, _set_id, _entry_id = seed(con)
    con.close()

    con2, _caps, _fts = open_project_db(db_path)
    assert con2.execute("PRAGMA user_version").fetchone()[0] == SCHEMA_VERSION
    # Re-running the monolithic schema would have thrown or wiped the data.
    assert con2.execute("SELECT COUNT(*) FROM fasta_record").fetchone()[0] == 3
    assert con2.execute("SELECT id FROM fasta_file").fetchone()[0] == file_id
    con2.close()


def test_migrate_is_idempotent(db_path):
    con = connect(db_path)
    assert migrate(con) == SCHEMA_VERSION
    assert migrate(con) == SCHEMA_VERSION
    con.close()


def test_a_newer_database_is_refused(db_path):
    con = connect(db_path)
    migrate(con)
    con.execute(f"PRAGMA user_version = {SCHEMA_VERSION + 5}")
    with pytest.raises(RuntimeError, match="only"):
        migrate(con)
    con.close()


# ---------------------------------------------------------------------------
# Integrity
# ---------------------------------------------------------------------------


def test_foreign_key_and_integrity_checks_are_clean(db_path):
    con, _caps, _fts = open_project_db(db_path)
    seed(con)
    assert con.execute("PRAGMA foreign_key_check").fetchall() == []
    assert con.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
    con.close()


def test_foreign_keys_are_actually_enforced(db_path):
    con, _caps, _fts = open_project_db(db_path)
    seed(con)
    with pytest.raises(sqlite3.IntegrityError):
        con.execute(
            "INSERT INTO fasta_record(fasta_file_id, ordinal, header, header_casefold,"
            " sequence_length) VALUES ('no-such-file', 0, 'h', 'h', 1)")
    con.close()


def test_fingerprint_columns_are_all_or_nothing(db_path):
    con, _caps, _fts = open_project_db(db_path)
    now = stamp()
    with pytest.raises(sqlite3.IntegrityError):
        con.execute(
            "INSERT INTO fasta_file(id, source_path, source_path_key, display_name, sort_order,"
            " indexed_size_bytes, created_at_ms, updated_at_ms) VALUES ('x','p','k','d',0,10,?,?)",
            (now, now))
    con.close()


# ---------------------------------------------------------------------------
# The decoupling the audit corrected
# ---------------------------------------------------------------------------


def test_focal_locations_survive_a_full_header_index_rebuild(db_path):
    """
    The whole point of the correction: rebuilding fasta_record must not take
    focal membership locations with it.
    """
    con, _caps, _fts = open_project_db(db_path)
    file_id, _set_id, _entry_id = seed(con)

    before = con.execute("SELECT COUNT(*) FROM focal_entry_location").fetchone()[0]
    con.execute("DELETE FROM fasta_record WHERE fasta_file_id = ?", (file_id,))
    con.commit()
    after = con.execute("SELECT COUNT(*) FROM focal_entry_location").fetchone()[0]

    assert before == 2
    assert after == 2
    con.close()


def test_a_location_carries_everything_needed_to_verify_itself(db_path):
    con, _caps, _fts = open_project_db(db_path)
    columns = {row[1] for row in con.execute("PRAGMA table_info(focal_entry_location)")}
    assert {"fasta_file_id", "ordinal", "record_start_byte", "last_verified_at_ms"} <= columns
    # It must NOT depend on the header index.
    assert "fasta_record_id" not in columns
    con.close()


def test_duplicate_headers_get_separate_locations(db_path):
    con, _caps, _fts = open_project_db(db_path)
    _file_id, _set_id, entry_id = seed(con)
    ordinals = [
        row[0]
        for row in con.execute(
            "SELECT ordinal FROM focal_entry_location WHERE focal_entry_id = ? ORDER BY ordinal",
            (entry_id,))
    ]
    assert ordinals == [0, 1]
    con.close()


def test_removing_a_file_removes_its_locations(db_path):
    con, _caps, _fts = open_project_db(db_path)
    file_id, _set_id, _entry_id = seed(con)
    con.execute("DELETE FROM fasta_file WHERE id = ?", (file_id,))
    con.commit()
    assert con.execute("SELECT COUNT(*) FROM focal_entry_location").fetchone()[0] == 0
    assert con.execute("SELECT COUNT(*) FROM fasta_record").fetchone()[0] == 0
    con.close()


def test_deleting_a_focal_entry_removes_its_locations(db_path):
    con, _caps, _fts = open_project_db(db_path)
    _file_id, _set_id, entry_id = seed(con)
    con.execute("DELETE FROM focal_set_entry WHERE id = ?", (entry_id,))
    con.commit()
    assert con.execute("SELECT COUNT(*) FROM focal_entry_location").fetchone()[0] == 0
    con.close()


# ---------------------------------------------------------------------------
# Query plans
# ---------------------------------------------------------------------------


def plan(con: sqlite3.Connection, sql: str, params=()) -> str:
    rows = con.execute("EXPLAIN QUERY PLAN " + sql, params).fetchall()
    return " | ".join(row[-1] for row in rows)


def test_exact_header_lookups_use_an_index(db_path):
    con, _caps, _fts = open_project_db(db_path)
    file_id, _set_id, _entry_id = seed(con)

    single = plan(
        con,
        "SELECT id, ordinal, record_start_byte FROM fasta_record"
        " WHERE fasta_file_id = ? AND header = ?", (file_id, "dup|H"))
    assert "SCAN" not in single
    assert "idx_fasta_record_header_file" in single

    project_wide = plan(
        con, "SELECT fasta_file_id, ordinal FROM fasta_record WHERE header = ?", ("dup|H",))
    assert "SCAN" not in project_wide
    assert "idx_fasta_record_header_file" in project_wide
    con.close()


def test_ordinal_lookup_uses_the_unique_index(db_path):
    con, _caps, _fts = open_project_db(db_path)
    file_id, _set_id, _entry_id = seed(con)
    text = plan(
        con, "SELECT header FROM fasta_record WHERE fasta_file_id = ? AND ordinal = ?",
        (file_id, 0))
    assert "SCAN" not in text
    con.close()


def test_location_lookups_use_indexes(db_path):
    con, _caps, _fts = open_project_db(db_path)
    file_id, _set_id, entry_id = seed(con)

    by_entry = plan(
        con, "SELECT fasta_file_id, ordinal FROM focal_entry_location WHERE focal_entry_id = ?",
        (entry_id,))
    assert "SCAN" not in by_entry

    by_file = plan(
        con, "SELECT focal_entry_id, ordinal FROM focal_entry_location WHERE fasta_file_id = ?",
        (file_id,))
    assert "SCAN" not in by_file
    assert "idx_focal_entry_location_file" in by_file
    con.close()


def test_substring_search_is_a_scan_by_design(db_path):
    """
    A B-tree cannot accelerate %substring%, so the fallback scans. This is
    pinned so nobody 'optimises' it with an index that cannot help.
    """
    con, _caps, _fts = open_project_db(db_path)
    seed(con)
    text = plan(con, "SELECT id FROM fasta_record WHERE instr(header_casefold, ?) > 0", ("up",))
    assert "SCAN" in text
    con.close()


def test_no_redundant_index_on_fasta_record(db_path):
    con, _caps, _fts = open_project_db(db_path)
    names = {
        row[0]
        for row in con.execute(
            "SELECT name FROM sqlite_schema WHERE type='index' AND tbl_name='fasta_record'"
            " AND sql IS NOT NULL")
    }
    # One explicit index serves both exact-header lookups; the second one the
    # audit found was pure write cost.
    assert names == {"idx_fasta_record_header_file"}
    con.close()
