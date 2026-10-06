-- 001_initial.sql
-- Molecular Diagnosis project database schema v1.
-- One SQLite database == one project.
--
-- IMPORTANT:
--   * Connection-level PRAGMAs (foreign_keys, trusted_schema, busy_timeout,
--     journal_mode, synchronous) belong in Python connection setup, not here.
--   * PRAGMA user_version is the application migration version.
--   * FTS5 trigram search is optional and installed separately only when the
--     bundled/runtime SQLite reports FTS5 support.
--
-- REVIEW NOTES (see REPO_MAP/UI_NOTES for the full audit). Four corrections
-- were made to the originally proposed v1 after testing it against the
-- requirements on the real runtime (SQLite 3.43.1):
--
--   1. focal_entry_location no longer references fasta_record. It referenced
--      it ON DELETE CASCADE, which meant a full header-index rebuild
--      (DELETE FROM fasta_record ...) silently destroyed every focal-entry
--      location. That was demonstrated, not assumed. Focal locations are a
--      persistent cache that must outlive header-index rebuilds.
--   2. focal_entry_location now carries the coordinates needed to verify a
--      location directly against the current source file (fasta_file_id,
--      ordinal, byte offsets) plus last_verified_at_ms. Previously a location
--      could only be resolved by joining fasta_record, i.e. only while the
--      header index was itself trustworthy.
--   3. Its primary key includes ordinal, so several identical headers in one
--      file are preserved as separate locations.
--   4. The fingerprint columns are named indexed_* and are nullable. They
--      describe the source version the index was built from, NOT the file's
--      current stat; the old names (size_bytes/mtime_ns/sha256) invited
--      exactly that confusion, and NOT NULL forced a placeholder value for a
--      file that has never been indexed.
--
-- Also dropped: idx_fasta_record_file_header. EXPLAIN QUERY PLAN shows the
-- planner satisfies both the single-file and project-wide exact-header lookups
-- from idx_fasta_record_header_file, so the second index only added write cost
-- to every reindex.

BEGIN IMMEDIATE;

CREATE TABLE project_metadata (
    singleton       INTEGER PRIMARY KEY CHECK (singleton = 1),
    project_uuid    TEXT    NOT NULL UNIQUE,
    title           TEXT    NOT NULL CHECK (length(trim(title)) > 0),
    created_at_ms   INTEGER NOT NULL CHECK (created_at_ms >= 0),
    updated_at_ms   INTEGER NOT NULL CHECK (updated_at_ms >= created_at_ms)
);

CREATE TABLE fasta_file (
    id                      TEXT    PRIMARY KEY,
    source_path             TEXT    NOT NULL,
    source_path_key         TEXT    NOT NULL UNIQUE,
    display_name            TEXT    NOT NULL CHECK (length(display_name) > 0),
    sort_order              INTEGER NOT NULL CHECK (sort_order >= 0),

    -- Fingerprint of the source version that fasta_record/FTS were built from.
    -- NULL across all three means "never indexed". These are NOT the file's
    -- current stat; current stat is runtime state and is never persisted.
    indexed_size_bytes      INTEGER CHECK (indexed_size_bytes IS NULL OR indexed_size_bytes >= 0),
    indexed_mtime_ns        INTEGER CHECK (indexed_mtime_ns IS NULL OR indexed_mtime_ns >= 0),
    indexed_sha256          TEXT
        CHECK (indexed_sha256 IS NULL OR
               (length(indexed_sha256) = 64 AND indexed_sha256 NOT GLOB '*[^0-9a-fA-F]*')),

    sequence_count          INTEGER CHECK (sequence_count IS NULL OR sequence_count >= 0),
    alignment_length        INTEGER CHECK (alignment_length IS NULL OR alignment_length >= 0),
    duplicate_header_count  INTEGER NOT NULL DEFAULT 0 CHECK (duplicate_header_count >= 0),

    -- Incremented only on a true content reindex, never for a metadata-only
    -- refresh or a byte-identical relink.
    index_revision          INTEGER NOT NULL DEFAULT 0 CHECK (index_revision >= 0),
    indexed_at_ms           INTEGER CHECK (indexed_at_ms IS NULL OR indexed_at_ms >= 0),

    created_at_ms           INTEGER NOT NULL CHECK (created_at_ms >= 0),
    updated_at_ms           INTEGER NOT NULL CHECK (updated_at_ms >= created_at_ms),

    -- The fingerprint is all-or-nothing: a half-written one would let a stale
    -- index look current on one axis and unknown on another.
    CHECK (
        (indexed_size_bytes IS NULL AND indexed_mtime_ns IS NULL AND indexed_sha256 IS NULL)
        OR
        (indexed_size_bytes IS NOT NULL AND indexed_mtime_ns IS NOT NULL
         AND indexed_sha256 IS NOT NULL)
    )
);

CREATE INDEX idx_fasta_file_sort
    ON fasta_file(sort_order, id);

-- Complete header snapshot of ONE indexed version of a FASTA. Rebuilt whole,
-- never mutated record-by-record. Trustworthy only while fasta_file.indexed_*
-- still describes the current source.
CREATE TABLE fasta_record (
    id                  INTEGER PRIMARY KEY,
    fasta_file_id       TEXT    NOT NULL
        REFERENCES fasta_file(id) ON DELETE CASCADE,
    ordinal             INTEGER NOT NULL CHECK (ordinal >= 0),

    header              TEXT    NOT NULL CHECK (length(header) > 0),
    -- Produced with Python str.casefold(); the canonical search form.
    header_casefold     TEXT    NOT NULL CHECK (length(header_casefold) > 0),
    sequence_length     INTEGER NOT NULL CHECK (sequence_length >= 0),

    -- Byte offsets are derived from the exact indexed file version.
    -- They allow future direct record retrieval / viewer jumps without
    -- rescanning from the beginning of the FASTA. NULL is allowed for
    -- importers that cannot provide offsets yet.
    record_start_byte   INTEGER CHECK (record_start_byte IS NULL OR record_start_byte >= 0),
    record_end_byte     INTEGER CHECK (
        record_end_byte IS NULL OR
        (record_end_byte >= 0 AND
         (record_start_byte IS NULL OR record_end_byte >= record_start_byte))
    ),

    UNIQUE (fasta_file_id, ordinal)
);

-- Exact-header resolution, both project-wide (header = ?) and for one file
-- (header = ? AND fasta_file_id = ?). Verified with EXPLAIN QUERY PLAN.
CREATE INDEX idx_fasta_record_header_file
    ON fasta_record(header, fasta_file_id, id);

CREATE TABLE focal_set (
    id              TEXT    PRIMARY KEY,
    title           TEXT    NOT NULL CHECK (length(trim(title)) > 0),
    locked          INTEGER NOT NULL DEFAULT 0 CHECK (locked IN (0, 1)),
    sort_order      INTEGER NOT NULL CHECK (sort_order >= 0),
    created_at_ms   INTEGER NOT NULL CHECK (created_at_ms >= 0),
    updated_at_ms   INTEGER NOT NULL CHECK (updated_at_ms >= created_at_ms)
);

CREATE INDEX idx_focal_set_sort
    ON focal_set(sort_order, id);

CREATE TABLE focal_set_entry (
    id              TEXT    PRIMARY KEY,
    focal_set_id    TEXT    NOT NULL
        REFERENCES focal_set(id) ON DELETE CASCADE,

    -- This is a resolved exact FASTA header, never a substring query.
    header          TEXT    NOT NULL CHECK (length(header) > 0),
    sort_order      INTEGER NOT NULL CHECK (sort_order >= 0),

    created_at_ms   INTEGER NOT NULL CHECK (created_at_ms >= 0),
    updated_at_ms   INTEGER NOT NULL CHECK (updated_at_ms >= created_at_ms),

    UNIQUE (focal_set_id, header)
);

CREATE INDEX idx_focal_set_entry_order
    ON focal_set_entry(focal_set_id, sort_order, id);

-- Persistent cache of where a focal entry's exact header occurs.
--
-- Deliberately independent of fasta_record: it survives application sessions
-- AND full header-index rebuilds, and each row carries everything needed to
-- re-verify itself directly against the current source file (seek to
-- record_start_byte, read the header, compare it to focal_set_entry.header).
-- The source FASTA stays authoritative; a saved location is an optimisation.
--
-- ordinal is part of the key so several identical headers in one file are kept
-- as separate locations.
CREATE TABLE focal_entry_location (
    focal_entry_id      TEXT    NOT NULL
        REFERENCES focal_set_entry(id) ON DELETE CASCADE,
    fasta_file_id       TEXT    NOT NULL
        REFERENCES fasta_file(id) ON DELETE CASCADE,
    ordinal             INTEGER NOT NULL CHECK (ordinal >= 0),

    record_start_byte   INTEGER CHECK (record_start_byte IS NULL OR record_start_byte >= 0),
    record_end_byte     INTEGER CHECK (
        record_end_byte IS NULL OR
        (record_end_byte >= 0 AND
         (record_start_byte IS NULL OR record_end_byte >= record_start_byte))
    ),

    -- When this individual location was last confirmed against the source.
    last_verified_at_ms INTEGER CHECK (last_verified_at_ms IS NULL OR last_verified_at_ms >= 0),

    PRIMARY KEY (focal_entry_id, fasta_file_id, ordinal)
) WITHOUT ROWID;

-- Reverse lookup: every location in a file, for revalidation after that file
-- changes, and for the per-file presence status.
CREATE INDEX idx_focal_entry_location_file
    ON focal_entry_location(fasta_file_id, focal_entry_id);

PRAGMA user_version = 1;

COMMIT;
