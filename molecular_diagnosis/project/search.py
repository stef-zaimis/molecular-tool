"""
Header search: plain case-insensitive substring matching.

The rule the whole application obeys is exactly this one line of Python:

    query.casefold() in header_casefold

FTS5 with the trigram tokenizer is an accelerator for it, never a different
search. Two things make the accelerated path safe:

* the user's text is wrapped as a single quoted FTS string, so `a OR b`,
  `lept*`, `a-b` and an unbalanced quote are matched literally instead of
  being parsed as FTS query syntax (all four were verified to either mismatch
  or raise when passed raw);
* every candidate FTS returns is re-checked in Python with the rule above.

The trigram tokenizer cannot answer queries shorter than three characters — it
silently returns nothing — so those always take the fallback.
"""

from __future__ import annotations

from dataclasses import dataclass

from molecular_diagnosis.project.repository import Repository

__all__ = ["HeaderHit", "MIN_TRIGRAM_LENGTH", "quote_fts_query", "search_headers"]

#: Below this the trigram index has no term to look up.
MIN_TRIGRAM_LENGTH = 3

#: Results come back in the project's own file order, then in file order within
#: a file. Ordering by `fasta_file_id` instead — as this once did — orders by a
#: random UUID, which made `+` expand to the same headers in a DIFFERENT order
#: from one run to the next, and that order is what a focal set stores.
ORDER_BY = " ORDER BY f.sort_order, f.id, r.ordinal"


@dataclass(frozen=True)
class HeaderHit:
    fasta_file_id: str
    ordinal: int
    header: str
    record_start_byte: int | None

    def to_payload(self) -> dict[str, object]:
        return {
            "fastaFileId": self.fasta_file_id,
            "ordinal": self.ordinal,
            "header": self.header,
            "recordStartByte": self.record_start_byte,
        }


def quote_fts_query(text: str) -> str:
    """Wrap user text as one literal FTS string; `"` is escaped by doubling."""
    return '"' + text.replace('"', '""') + '"'


def _rows_to_hits(rows) -> list[HeaderHit]:
    return [
        HeaderHit(row["fasta_file_id"], row["ordinal"], row["header"], row["record_start_byte"])
        for row in rows
    ]


def search_headers(
    repository: Repository,
    query: str,
    *,
    fasta_file_ids: list[str] | None = None,
    limit: int | None = None,
    force_fallback: bool = False,
) -> list[HeaderHit]:
    """
    Every record whose header contains `query`, case-insensitively.

    `fasta_file_ids=None` searches the whole project. Results are ordered by
    file then ordinal so the caller sees them in file order.
    """
    needle = query.casefold()
    if not needle:
        return []

    con = repository.con
    scope_sql = ""
    scope_params: list[object] = []
    if fasta_file_ids is not None:
        if not fasta_file_ids:
            return []
        placeholders = ",".join("?" for _ in fasta_file_ids)
        scope_sql = f" AND r.fasta_file_id IN ({placeholders})"
        scope_params = list(fasta_file_ids)

    use_fts = (
        repository.fts_enabled
        and not force_fallback
        and len(needle) >= MIN_TRIGRAM_LENGTH
    )

    if use_fts:
        sql = (
            "SELECT r.fasta_file_id, r.ordinal, r.header, r.header_casefold, r.record_start_byte"
            " FROM fasta_record_search s"
            " JOIN fasta_record r ON r.id = s.rowid"
            " JOIN fasta_file f ON f.id = r.fasta_file_id"
            " WHERE fasta_record_search MATCH ?" + scope_sql +
            ORDER_BY
        )
        try:
            rows = con.execute(sql, [quote_fts_query(needle), *scope_params]).fetchall()
        except Exception:
            # Any FTS surprise degrades to the fallback rather than to no results.
            rows = None
        if rows is not None:
            # Post-filter with the real rule: trigram matching is an
            # approximation of substring containment, never the definition.
            hits = _rows_to_hits(row for row in rows if needle in row["header_casefold"])
            return hits[:limit] if limit else hits

    sql = (
        "SELECT r.fasta_file_id, r.ordinal, r.header, r.header_casefold, r.record_start_byte"
        " FROM fasta_record r"
        " JOIN fasta_file f ON f.id = r.fasta_file_id"
        " WHERE instr(r.header_casefold, ?) > 0" + scope_sql +
        ORDER_BY
    )
    rows = con.execute(sql, [needle, *scope_params]).fetchall()
    hits = _rows_to_hits(row for row in rows if needle in row["header_casefold"])
    return hits[:limit] if limit else hits


def resolve_headers(
    repository: Repository,
    query: str,
    *,
    fasta_file_ids: list[str] | None = None,
) -> list[str]:
    """
    The complete, deduplicated exact headers a query selects.

    This is the `+` expansion: a transient query becomes a set of exact
    headers, and only those exact headers ever enter a focal set.
    """
    seen: set[str] = set()
    ordered: list[str] = []
    for hit in search_headers(repository, query, fasta_file_ids=fasta_file_ids):
        if hit.header in seen:
            continue
        seen.add(hit.header)
        ordered.append(hit.header)
    return ordered
