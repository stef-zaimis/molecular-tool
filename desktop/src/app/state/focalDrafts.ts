/**
 * Focal sets as SAVED DATA plus a WORKING COPY.
 *
 * The distinction is the point of this file. A focal set in the project
 * database is a scientific commitment — it names the sequences an analysis
 * will treat as focal. Typing in the workspace is not that commitment, so
 * editing writes nothing: every change lands in a draft here, and only
 * `project.saveFocalSet` moves a draft's contents into the database.
 *
 * That is what makes "unsaved changes" a real, checkable state rather than a
 * label stuck on data that was already committed, and it is why Run refuses a
 * dirty draft instead of quietly saving it first. A run reports which focal set
 * it used; if Run could auto-save, that name would describe something the user
 * never chose.
 */

import type { FocalSetPayload } from '../../backendContract';

/**
 * Undo/redo over the working header list.
 *
 * Snapshots, not inverse operations: focal sets are tens of entries, so a
 * complete copy per step is simpler and cannot drift out of sync with the
 * edit that produced it. Title edits are deliberately NOT in the history —
 * the undo/redo controls sit on the header box in the design, and folding an
 * unrelated rename into a header undo would surprise.
 */
export interface FocalHistory {
  readonly past: readonly (readonly string[])[];
  readonly future: readonly (readonly string[])[];
}

export const EMPTY_HISTORY: FocalHistory = { past: [], future: [] };

/**
 * Which of a draft's two header lists an edit addresses.
 *
 * `focal` is the focal set proper. `comparison` is its optional explicit
 * Comparison Set: saved, locked and deleted WITH the focal set, never on its
 * own, which is why it lives on the same draft rather than in a second
 * library. Each list keeps its own undo stack, because each editor has its own
 * undo/redo controls.
 */
export type HeaderListKind = 'focal' | 'comparison';

export interface FocalDraft {
  /**
   * Stable local identity, for React keys and for addressing a draft that has
   * never been saved. NOT the database id — see `persistedId`.
   */
  readonly key: string;
  /** The row this draft edits, or null while it exists only in this session. */
  readonly persistedId: string | null;
  /** Last values the backend confirmed. The baseline `dirty` compares against. */
  readonly savedTitle: string;
  readonly savedHeaders: readonly string[];
  /** What the user is editing right now. */
  readonly title: string;
  readonly headers: readonly string[];
  /** A locked set is readable, selectable and runnable, but never editable. */
  readonly locked: boolean;
  readonly history: FocalHistory;
  /**
   * True while a continuous manual typing burst is open.
   *
   * Typing must not make every keystroke its own undo step — undoing a pasted
   * header one character at a time is useless. While a burst is open the
   * working headers are replaced WITHOUT pushing history, so the whole burst
   * collapses to the one step recorded when it began. `+` and `−` never
   * coalesce: each is one deliberate action and one step.
   */
  readonly burstOpen: boolean;
  /**
   * The explicit Comparison Set: saved, working, history and burst, exactly
   * parallel to the focal fields above. Empty means "compare against every
   * non-focal sequence in the FASTA pool".
   */
  readonly savedComparison: readonly string[];
  readonly comparison: readonly string[];
  readonly comparisonHistory: FocalHistory;
  readonly comparisonBurstOpen: boolean;
}

/** One list's editable state, so every edit operation is written once. */
interface HeaderListView {
  readonly headers: readonly string[];
  readonly history: FocalHistory;
  readonly burstOpen: boolean;
}

function viewOf(draft: FocalDraft, list: HeaderListKind): HeaderListView {
  return list === 'focal'
    ? { headers: draft.headers, history: draft.history, burstOpen: draft.burstOpen }
    : {
        headers: draft.comparison,
        history: draft.comparisonHistory,
        burstOpen: draft.comparisonBurstOpen,
      };
}

function withView(draft: FocalDraft, list: HeaderListKind, view: HeaderListView): FocalDraft {
  return list === 'focal'
    ? { ...draft, headers: view.headers, history: view.history, burstOpen: view.burstOpen }
    : {
        ...draft,
        comparison: view.headers,
        comparisonHistory: view.history,
        comparisonBurstOpen: view.burstOpen,
      };
}

/**
 * Trim, drop blanks, and collapse duplicates onto their first occurrence.
 *
 * Mirrors `dedupe_headers` in `project/service.py`. Both exist because both
 * ends need the answer: Python to store it, the renderer to compare against
 * what was stored without a round trip. The rule is plain enough that the two
 * cannot drift, and `saveFocalSet` returns the canonical result anyway.
 */
export function normaliseHeaders(headers: readonly string[]): readonly string[] {
  const seen = new Set<string>();
  const ordered: string[] = [];
  for (const raw of headers) {
    const header = raw.trim();
    if (!header || seen.has(header)) continue;
    seen.add(header);
    ordered.push(header);
  }
  return ordered;
}

function sameHeaders(left: readonly string[], right: readonly string[]): boolean {
  // Order matters: `sort_order` is persisted, so a reorder is a real change.
  return left.length === right.length && left.every((value, index) => value === right[index]);
}

let draftCounter = 0;

/** A local key. Never sent anywhere: the backend mints the real id on save. */
export function nextDraftKey(): string {
  draftCounter += 1;
  return `draft-${draftCounter}`;
}

/** An empty, never-saved draft. No database row exists until Save. */
export function blankDraft(title = ''): FocalDraft {
  return {
    key: nextDraftKey(),
    persistedId: null,
    savedTitle: '',
    savedHeaders: [],
    title,
    headers: [],
    locked: false,
    history: EMPTY_HISTORY,
    burstOpen: false,
    savedComparison: [],
    comparison: [],
    comparisonHistory: EMPTY_HISTORY,
    comparisonBurstOpen: false,
  };
}

/**
 * A draft whose working copy starts out identical to what is stored.
 *
 * A set saved before comparison sets existed carries no comparison list; it
 * loads as empty, which is exactly the behaviour it always had.
 */
export function draftFromPayload(payload: FocalSetPayload): FocalDraft {
  const headers = payload.entries.map((entry) => entry.header);
  const comparison = payload.comparisonHeaders ?? [];
  return {
    key: `set-${payload.id}`,
    persistedId: payload.id,
    savedTitle: payload.title,
    savedHeaders: headers,
    title: payload.title,
    headers,
    locked: payload.locked,
    history: EMPTY_HISTORY,
    burstOpen: false,
    savedComparison: comparison,
    comparison,
    comparisonHistory: EMPTY_HISTORY,
    comparisonBurstOpen: false,
  };
}

/**
 * Adopt a backend response as the new saved AND working state.
 *
 * The response is canonical — it carries the trimming, deduplication and
 * ordering the database actually applied — so the working copy is replaced by
 * it rather than left as whatever the user typed. History is cleared: the steps
 * that led here describe a draft that no longer exists.
 */
export function draftSaved(draft: FocalDraft, payload: FocalSetPayload): FocalDraft {
  return { ...draftFromPayload(payload), key: draft.key };
}

/** True until the draft has been saved even once. */
export function isNewDraft(draft: FocalDraft): boolean {
  return draft.persistedId === null;
}

/** Derived, never stored: comparing normalised working values against saved ones. */
export function isDirty(draft: FocalDraft): boolean {
  if (isNewDraft(draft)) return true;
  if (draft.title.trim() !== draft.savedTitle.trim()) return true;
  if (!sameHeaders(normaliseHeaders(draft.comparison), normaliseHeaders(draft.savedComparison))) {
    return true;
  }
  return !sameHeaders(normaliseHeaders(draft.headers), normaliseHeaders(draft.savedHeaders));
}

/**
 * Exact headers that are in BOTH the focal and the comparison list.
 *
 * Complete-header equality, the rule the backend refuses a run on. Neither
 * side is silently trimmed of them: both editors mark them and Run waits
 * until the user decides which set each one belongs to.
 */
export function overlappingHeaders(draft: FocalDraft): readonly string[] {
  const focal = new Set(normaliseHeaders(draft.headers));
  return normaliseHeaders(draft.comparison).filter((header) => focal.has(header));
}

/**
 * Replace the working headers, pushing one undo step.
 *
 * `coalesce` marks a keystroke inside a continuous manual burst: the first one
 * records a step and opens the burst, and the rest fold into it. Anything
 * else — `+`, `−`, a programmatic replacement — closes the burst and records
 * its own step, so the toolbar's undo always walks back one whole action.
 */
export function withHeaders(
  draft: FocalDraft,
  headers: readonly string[],
  { coalesce = false, list = 'focal' }: { coalesce?: boolean; list?: HeaderListKind } = {},
): FocalDraft {
  const current = viewOf(draft, list);
  const next = normaliseHeaders(headers);
  if (sameHeaders(next, current.headers)) {
    return current.burstOpen === coalesce
      ? draft
      : withView(draft, list, { ...current, burstOpen: coalesce });
  }

  if (coalesce && current.burstOpen) {
    // Inside an open burst: the step recorded when it began still describes
    // the state to go back to.
    return withView(draft, list, { ...current, headers: next });
  }

  return withView(draft, list, {
    headers: next,
    burstOpen: coalesce,
    // Any fresh edit invalidates the redo branch.
    history: { past: [...current.history.past, current.headers], future: [] },
  });
}

/** Close an open typing burst, so the next keystroke starts a new undo step. */
export function endBurst(draft: FocalDraft, list: HeaderListKind = 'focal'): FocalDraft {
  const current = viewOf(draft, list);
  return current.burstOpen ? withView(draft, list, { ...current, burstOpen: false }) : draft;
}

/** Close the bursts of BOTH lists — what leaving a draft does. */
export function endAllBursts(draft: FocalDraft): FocalDraft {
  return endBurst(endBurst(draft, 'focal'), 'comparison');
}

/** The working headers of one list. */
export function headersOf(draft: FocalDraft, list: HeaderListKind): readonly string[] {
  return viewOf(draft, list).headers;
}

/** `+`: append exact headers already resolved by the backend, in hit order. */
export function appendHeaders(
  draft: FocalDraft,
  incoming: readonly string[],
  list: HeaderListKind = 'focal',
): FocalDraft {
  return withHeaders(draft, [...headersOf(draft, list), ...incoming], { list });
}

/** `-`: drop the headers the backend matched against this working copy. */
export function removeHeaders(
  draft: FocalDraft,
  doomed: readonly string[],
  list: HeaderListKind = 'focal',
): FocalDraft {
  const remove = new Set(doomed);
  const headers = headersOf(draft, list);
  if (!headers.some((header) => remove.has(header))) return draft;
  return withHeaders(
    draft,
    headers.filter((header) => !remove.has(header)),
    { list },
  );
}

export function canUndo(draft: FocalDraft, list: HeaderListKind = 'focal'): boolean {
  return viewOf(draft, list).history.past.length > 0;
}

export function canRedo(draft: FocalDraft, list: HeaderListKind = 'focal'): boolean {
  return viewOf(draft, list).history.future.length > 0;
}

export function undoDraft(draft: FocalDraft, list: HeaderListKind = 'focal'): FocalDraft {
  const current = viewOf(draft, list);
  const { past, future } = current.history;
  if (past.length === 0) return draft;
  return withView(draft, list, {
    headers: past[past.length - 1],
    // Undoing ends the burst: the next keystroke is a new step, not a
    // continuation of the one just reverted.
    burstOpen: false,
    history: { past: past.slice(0, -1), future: [current.headers, ...future] },
  });
}

export function redoDraft(draft: FocalDraft, list: HeaderListKind = 'focal'): FocalDraft {
  const current = viewOf(draft, list);
  const { past, future } = current.history;
  if (future.length === 0) return draft;
  const [next, ...rest] = future;
  return withView(draft, list, {
    headers: next,
    burstOpen: false,
    history: { past: [...past, current.headers], future: rest },
  });
}

/**
 * The payload for `project.saveFocalSet`, from a draft.
 *
 * Focal and comparison travel together so the backend writes both in one
 * transaction. `comparisonHeaders` is always sent — `[]` included — so saving
 * a cleared comparison list really clears it.
 */
export function saveRequestFor(draft: FocalDraft): {
  focalSetId: string | null;
  title: string;
  headers: readonly string[];
  comparisonHeaders: readonly string[];
} {
  return {
    focalSetId: draft.persistedId,
    title: draft.title.trim(),
    headers: normaliseHeaders(draft.headers),
    comparisonHeaders: normaliseHeaders(draft.comparison),
  };
}

/**
 * Why locking is unavailable, or null when it may be locked.
 *
 * Locking marks a set as final, so there must BE a final version: locking a
 * dirty draft would either lock the stale saved copy (silently discarding the
 * edits) or auto-save first (committing a version the user never approved).
 * Neither is acceptable, so it is refused with an explanation instead.
 */
export function lockBlockedReason(draft: FocalDraft): string | null {
  if (isNewDraft(draft)) return 'Save this focal set before locking it.';
  if (isDirty(draft)) return 'Save the changes to this focal set before locking it.';
  return null;
}

/** Why Save is unavailable, or null when it is available. */
export function saveBlockedReason(draft: FocalDraft): string | null {
  if (draft.locked) return 'This focal set is locked. Unlock it before editing.';
  if (!draft.title.trim()) return 'Give the focal set a title before saving it.';
  if (!isDirty(draft)) return 'No unsaved changes.';
  return null;
}
