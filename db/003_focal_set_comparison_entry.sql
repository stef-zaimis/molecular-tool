-- 003_focal_set_comparison_entry.sql
-- An optional explicit Comparison Set, stored WITH the focal set it belongs to.
--
-- A saved focal set may name exact headers to compare against. Zero rows for
-- a set means "no explicit comparison": the run compares against every
-- non-focal sequence in the selected scope, exactly as before this table
-- existed. That is why existing projects need no backfill — every set they
-- hold already means what it meant.
--
-- Deliberately NOT a second library: there is no comparison_set table and no
-- title. A comparison list has no identity apart from its focal set, and it is
-- saved, locked and deleted with it (ON DELETE CASCADE).
--
-- There is no location cache for these entries. The run checks them against
-- the verified sequences it has already loaded, and editing feedback comes
-- from the header index, so a cache would only be one more thing to keep
-- consistent.
--
-- Additive, so an existing project opens unchanged.

BEGIN IMMEDIATE;

CREATE TABLE focal_set_comparison_entry (
    id              TEXT    PRIMARY KEY,
    focal_set_id    TEXT    NOT NULL
        REFERENCES focal_set(id) ON DELETE CASCADE,

    -- A resolved exact FASTA header, never a substring query.
    header          TEXT    NOT NULL CHECK (length(header) > 0),
    sort_order      INTEGER NOT NULL CHECK (sort_order >= 0),

    created_at_ms   INTEGER NOT NULL CHECK (created_at_ms >= 0),
    updated_at_ms   INTEGER NOT NULL CHECK (updated_at_ms >= created_at_ms),

    UNIQUE (focal_set_id, header)
);

CREATE INDEX idx_focal_set_comparison_entry_order
    ON focal_set_comparison_entry(focal_set_id, sort_order, id);

PRAGMA user_version = 3;

COMMIT;
