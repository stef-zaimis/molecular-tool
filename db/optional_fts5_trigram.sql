-- optional_fts5_trigram.sql
-- Apply only after Python verifies that the runtime SQLite supports FTS5 and
-- the trigram tokenizer. This is a performance accelerator, not authoritative
-- project data; the application must be fully correct without it.
--
-- Keep rowid equal to fasta_record.id. The repository layer is responsible
-- for inserting/deleting these rows in the same transaction as fasta_record.
--
-- The indexed column holds header_casefold (Python str.casefold()), because
-- that is the canonical case-insensitive search representation. Queries must
-- casefold the user's text the same way, wrap it as a single quoted FTS
-- string so it cannot be read as FTS query syntax, and then post-filter the
-- candidates with the real Python rule
--     query.casefold() in header_casefold
-- so the accelerated and fallback paths return identical results.
--
-- The trigram tokenizer cannot serve queries shorter than 3 characters: it
-- returns no rows rather than an error. Those queries MUST take the fallback
-- path. Verified on SQLite 3.43.1.

CREATE VIRTUAL TABLE IF NOT EXISTS fasta_record_search
USING fts5(
    header_casefold,
    tokenize = 'trigram case_sensitive 0'
);
