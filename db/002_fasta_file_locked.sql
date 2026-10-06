-- 002_fasta_file_locked.sql
-- Per-source lock, so a FASTA the user has settled on cannot be unlinked or
-- renamed by accident.
--
-- A locked source stays fully usable for analysis: locking says "stop editing
-- this link", not "stop reading this file". That is why nothing here touches
-- the index, the fingerprint or the record table.
--
-- Additive and defaulted, so an existing project opens unchanged.

BEGIN IMMEDIATE;

ALTER TABLE fasta_file
    ADD COLUMN locked INTEGER NOT NULL DEFAULT 0 CHECK (locked IN (0, 1));

PRAGMA user_version = 2;

COMMIT;
