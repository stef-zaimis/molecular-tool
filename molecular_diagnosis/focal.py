"""
Focal-set selection.

THE single place that decides which FASTA headers belong to the focal set.

Before this module existed the rule `target_string in header` was written out
independently in four places (`fasta_io.split_focal_headers`,
`core.find_dmc_information`, `core.compute_metrics`, `excel.build_sheet`), which
meant any change to focal semantics had to be made four times and kept in step
by hand. Every one of those call sites now routes through `header_matches_focal`.

Semantics, unchanged from the original single-string behaviour except for the
addition of multiple selectors:

- A header belongs to the focal set when ANY focal string occurs literally
  within it (logical OR across the selectors).
- Matching is CASE-SENSITIVE literal substring containment. It is not exact-id
  equality, not case-insensitive, and not a pattern match.
- Selectors are stripped of surrounding whitespace, as the original pipeline
  did, and are then deduplicated while preserving first-seen order.
- An empty (or whitespace-only) selector is invalid: under substring semantics
  it would match every header, which is never what a caller means.

A single `str` is still accepted everywhere a selector list is, so existing
callers and tests keep working; it is treated as a one-element list.
"""

from collections.abc import Sequence

__all__ = [
    "FocalSelector",
    "focal_label",
    "header_matches_focal",
    "matching_focal_strings",
    "normalise_focal_strings",
    "partition_headers",
]

#: A single selector or a list of them. Both are accepted everywhere.
FocalSelector = str | Sequence[str]


def normalise_focal_strings(
    value: FocalSelector,
    *,
    empty_message: str = "No identifier string entered.",
) -> list[str]:
    """
    Coerce a selector argument into a clean, ordered, duplicate-free list.

    `empty_message` lets callers keep their existing user-facing wording for the
    "nothing supplied" case rather than inventing a new one.

    Raises:
        ValueError: if any individual selector is empty/whitespace-only, or if
            no selectors remain.
    """
    candidates: list[str]

    if isinstance(value, str):
        candidates = [value]
    else:
        candidates = [str(item) for item in value]

    if not candidates:
        raise ValueError(empty_message)

    seen: set[str] = set()
    cleaned: list[str] = []

    for raw in candidates:
        selector = raw.strip()

        if not selector:
            raise ValueError("Focal strings cannot be empty.")

        if selector in seen:
            continue

        seen.add(selector)
        cleaned.append(selector)

    if not cleaned:
        raise ValueError(empty_message)

    return cleaned


def header_matches_focal(header: str, focal_strings: Sequence[str]) -> bool:
    """True when any focal string occurs literally inside `header`."""
    return any(selector in header for selector in focal_strings)


def matching_focal_strings(header: str, focal_strings: Sequence[str]) -> list[str]:
    """Which selectors matched this header. Useful for reporting, not for filtering."""
    return [selector for selector in focal_strings if selector in header]


def partition_headers(
    headers: Sequence[str],
    focal_strings: Sequence[str],
) -> tuple[list[str], list[str]]:
    """Split headers into (focal, non-focal), preserving input order in both."""
    focal: list[str] = []
    non_focal: list[str] = []

    for header in headers:
        if header_matches_focal(header, focal_strings):
            focal.append(header)
        else:
            non_focal.append(header)

    return focal, non_focal


def focal_label(focal_strings: Sequence[str]) -> str:
    """
    Human-readable name for the focal set, used in report text and in the
    consensus FASTA record names.

    A single selector formats exactly as it did before multi-selector support,
    so single-string runs produce byte-identical output.
    """
    return "+".join(focal_strings)
