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
  };
}

/** A draft whose working copy starts out identical to what is stored. */
export function draftFromPayload(payload: FocalSetPayload): FocalDraft {
  const headers = payload.entries.map((entry) => entry.header);
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
  return !sameHeaders(normaliseHeaders(draft.headers), normaliseHeaders(draft.savedHeaders));
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
  { coalesce = false }: { coalesce?: boolean } = {},
): FocalDraft {
  const next = normaliseHeaders(headers);
  if (sameHeaders(next, draft.headers)) {
    return draft.burstOpen === coalesce ? draft : { ...draft, burstOpen: coalesce };
  }

  if (coalesce && draft.burstOpen) {
    // Inside an open burst: the step recorded when it began still describes
    // the state to go back to.
    return { ...draft, headers: next };
  }

  return {
    ...draft,
    headers: next,
    burstOpen: coalesce,
    // Any fresh edit invalidates the redo branch.
    history: { past: [...draft.history.past, draft.headers], future: [] },
  };
}

/** Close an open typing burst, so the next keystroke starts a new undo step. */
export function endBurst(draft: FocalDraft): FocalDraft {
  return draft.burstOpen ? { ...draft, burstOpen: false } : draft;
}

/** `+`: append exact headers already resolved by the backend, in hit order. */
export function appendHeaders(draft: FocalDraft, incoming: readonly string[]): FocalDraft {
  return withHeaders(draft, [...draft.headers, ...incoming]);
}

/** `-`: drop the headers the backend matched against this working copy. */
export function removeHeaders(draft: FocalDraft, doomed: readonly string[]): FocalDraft {
  const remove = new Set(doomed);
  if (!draft.headers.some((header) => remove.has(header))) return draft;
  return withHeaders(
    draft,
    draft.headers.filter((header) => !remove.has(header)),
  );
}

export function canUndo(draft: FocalDraft): boolean {
  return draft.history.past.length > 0;
}

export function canRedo(draft: FocalDraft): boolean {
  return draft.history.future.length > 0;
}

export function undoDraft(draft: FocalDraft): FocalDraft {
  const { past, future } = draft.history;
  if (past.length === 0) return draft;
  return {
    ...draft,
    headers: past[past.length - 1],
    // Undoing ends the burst: the next keystroke is a new step, not a
    // continuation of the one just reverted.
    burstOpen: false,
    history: { past: past.slice(0, -1), future: [draft.headers, ...future] },
  };
}

export function redoDraft(draft: FocalDraft): FocalDraft {
  const { past, future } = draft.history;
  if (future.length === 0) return draft;
  const [next, ...rest] = future;
  return {
    ...draft,
    headers: next,
    burstOpen: false,
    history: { past: [...past, draft.headers], future: rest },
  };
}

/** The payload for `project.saveFocalSet`, from a draft. */
export function saveRequestFor(draft: FocalDraft): {
  focalSetId: string | null;
  title: string;
  headers: readonly string[];
} {
  return {
    focalSetId: draft.persistedId,
    title: draft.title.trim(),
    headers: normaliseHeaders(draft.headers),
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
