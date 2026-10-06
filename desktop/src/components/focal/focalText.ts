/**
 * The focal set as TEXT, and the text as a focal set.
 *
 * The editor is a semicolon-separated document the user types into directly,
 * so the array and the text have to round-trip predictably. Two rules make
 * that possible:
 *
 *  - `;` separates entries and is never part of one. A FASTA header may
 *    legitimately contain spaces, so whitespace INSIDE an entry is preserved;
 *    only whitespace immediately around a separator is dropped.
 *  - The text is the working document, the array is the working membership.
 *    Parsing is lossy on purpose (blanks vanish, duplicates collapse), which
 *    is why the text is kept as its own state rather than re-serialised from
 *    the array on every keystroke — doing that would yank the caret around
 *    while someone is still typing an entry.
 *
 * This mirrors `parse_focal_text` in `project/service.py`. Both ends need the
 * answer; the rule is small enough that they cannot drift, and Save returns
 * the canonical result anyway.
 */

export const ENTRY_SEPARATOR = ';';

export interface FocalTextSpan {
  /** Entry text with surrounding whitespace trimmed off. */
  readonly text: string;
  /** Offsets of the TRIMMED text within the document. */
  readonly from: number;
  readonly to: number;
  /** False for blank segments and for repeats of an earlier entry. */
  readonly isMember: boolean;
}

/**
 * Split a document into entry spans, keeping document offsets.
 *
 * Offsets are what lets the editor colour each entry in place, so they are
 * computed here rather than re-derived by a decorator that would have to
 * re-implement the same trimming.
 */
export function scanFocalText(text: string): FocalTextSpan[] {
  const spans: FocalTextSpan[] = [];
  const seen = new Set<string>();
  let cursor = 0;

  for (const segment of text.split(ENTRY_SEPARATOR)) {
    const start = cursor;
    cursor += segment.length + ENTRY_SEPARATOR.length;

    const leading = segment.length - segment.trimStart().length;
    const trimmed = segment.trim();
    const from = start + leading;

    // A duplicate is still shown — it is text the user typed — but it is not a
    // second member, and colouring it as one would imply the set has two.
    const isMember = trimmed.length > 0 && !seen.has(trimmed);
    if (isMember) seen.add(trimmed);

    spans.push({ text: trimmed, from, to: from + trimmed.length, isMember });
  }

  return spans;
}

/** The membership a document expresses: trimmed, non-blank, first-wins. */
export function parseFocalText(text: string): string[] {
  return scanFocalText(text)
    .filter((span) => span.isMember)
    .map((span) => span.text);
}

/** Render a membership array as a document, in the design's `a; b; c` form. */
export function serialiseFocalText(headers: readonly string[]): string {
  return headers.join(`${ENTRY_SEPARATOR} `);
}

/**
 * Do a document and an array express the same membership?
 *
 * Used to decide whether an incoming array (from `+`, `−`, undo, or a save
 * response) actually needs to replace what is in the editor. If it does not,
 * the document is left exactly as typed — including trailing separators and
 * half-finished entries the user is still working on.
 */
export function textMatchesHeaders(text: string, headers: readonly string[]): boolean {
  const parsed = parseFocalText(text);
  return parsed.length === headers.length && parsed.every((entry, i) => entry === headers[i]);
}
