import { describe, expect, it } from 'vitest';
import { activeDraft, appReducer, createInitialState } from './projectState';
import type { AppAction, AppState } from './projectState';
import {
  canRedo,
  canUndo,
  draftFromPayload,
  isDirty,
  overlappingHeaders,
  saveBlockedReason,
  saveRequestFor,
} from './focalDrafts';
import type { FocalDraft } from './focalDrafts';
import { evaluateRunGate } from './runGate';
import type { FocalPresenceState } from './projectState';
import type {
  FocalPresenceState as PresenceVerdict,
  FocalSetPayload,
  SourceStatusPayload,
} from '../../backendContract';

/**
 * The optional Comparison Set, as part of a focal-set draft.
 *
 * It is saved, dirtied, locked and loaded WITH its focal set, has its own undo
 * stack, and an empty list must behave exactly like a set saved before
 * comparison sets existed.
 */

const OLD_SET: FocalSetPayload = {
  id: 'set-1',
  title: 'Targets',
  locked: false,
  // No `comparisonHeaders`: a set saved before the feature existed.
  entries: [
    { id: 'e1', header: 'focal_1' },
    { id: 'e2', header: 'focal_2' },
  ],
};

const WITH_COMPARISON: FocalSetPayload = {
  ...OLD_SET,
  comparisonHeaders: ['other_1', 'other_2'],
};

function run(state: AppState, actions: readonly AppAction[]): AppState {
  return actions.reduce(appReducer, state);
}

function loaded(payload: FocalSetPayload = OLD_SET): AppState {
  return appReducer(createInitialState(), { type: 'focalSetsLoaded', focalSets: [payload] });
}

describe('loading saved sets', () => {
  it('loads a set saved before comparison sets as an empty, clean comparison list', () => {
    const draft = activeDraft(loaded(OLD_SET));
    expect(draft.comparison).toEqual([]);
    expect(draft.savedComparison).toEqual([]);
    expect(isDirty(draft)).toBe(false);
    expect(overlappingHeaders(draft)).toEqual([]);
  });

  it('restores both lists of a saved set', () => {
    const draft = activeDraft(loaded(WITH_COMPARISON));
    expect(draft.headers).toEqual(['focal_1', 'focal_2']);
    expect(draft.comparison).toEqual(['other_1', 'other_2']);
    expect(isDirty(draft)).toBe(false);
  });

  it('restores the right comparison list when switching between saved sets', () => {
    const second: FocalSetPayload = {
      id: 'set-2',
      title: 'Second',
      locked: false,
      entries: [{ id: 'x', header: 'focal_9' }],
      comparisonHeaders: ['other_9'],
    };
    const state = appReducer(createInitialState(), {
      type: 'focalSetsLoaded',
      focalSets: [WITH_COMPARISON, second],
    });
    const switched = appReducer(state, { type: 'selectFocalDraft', key: 'set-set-2' });
    expect(activeDraft(switched).comparison).toEqual(['other_9']);
    const back = appReducer(switched, { type: 'selectFocalDraft', key: 'set-set-1' });
    expect(activeDraft(back).comparison).toEqual(['other_1', 'other_2']);
  });
});

describe('editing the comparison list', () => {
  it('edits only the comparison list and makes the draft dirty', () => {
    const state = loaded();
    const key = activeDraft(state).key;
    const edited = run(state, [
      { type: 'setDraftHeaders', key, headers: ['other_1', ' other_1 ', ''], list: 'comparison' },
    ]);
    const draft = activeDraft(edited);
    expect(draft.comparison).toEqual(['other_1']);
    expect(draft.headers).toEqual(['focal_1', 'focal_2']);
    expect(isDirty(draft)).toBe(true);
    expect(saveBlockedReason(draft)).toBeNull();
  });

  it('is clean again when the comparison edit is undone', () => {
    const state = loaded(WITH_COMPARISON);
    const key = activeDraft(state).key;
    const edited = run(state, [
      { type: 'setDraftHeaders', key, headers: ['other_2'], list: 'comparison' },
      { type: 'undoDraftEdit', key, list: 'comparison' },
    ]);
    expect(activeDraft(edited).comparison).toEqual(['other_1', 'other_2']);
    expect(isDirty(activeDraft(edited))).toBe(false);
  });

  it('keeps separate undo stacks for the two lists', () => {
    const state = loaded();
    const key = activeDraft(state).key;
    const edited = run(state, [
      { type: 'setDraftHeaders', key, headers: ['focal_1'] },
      { type: 'setDraftHeaders', key, headers: ['other_1'], list: 'comparison' },
    ]);
    let draft = activeDraft(edited);
    expect(canUndo(draft, 'focal')).toBe(true);
    expect(canUndo(draft, 'comparison')).toBe(true);

    // Undoing the focal list leaves the comparison list alone.
    const focalUndone = appReducer(edited, { type: 'undoDraftEdit', key });
    draft = activeDraft(focalUndone);
    expect(draft.headers).toEqual(['focal_1', 'focal_2']);
    expect(draft.comparison).toEqual(['other_1']);
    expect(canRedo(draft, 'focal')).toBe(true);
    expect(canRedo(draft, 'comparison')).toBe(false);
  });

  it('coalesces a comparison typing burst into one undo step', () => {
    const state = loaded();
    const key = activeDraft(state).key;
    const typed = run(state, [
      { type: 'setDraftHeaders', key, headers: ['o'], coalesce: true, list: 'comparison' },
      { type: 'setDraftHeaders', key, headers: ['ot'], coalesce: true, list: 'comparison' },
      { type: 'setDraftHeaders', key, headers: ['oth'], coalesce: true, list: 'comparison' },
      { type: 'endDraftBurst', key, list: 'comparison' },
    ]);
    expect(activeDraft(typed).comparisonHistory.past).toHaveLength(1);
    const undone = appReducer(typed, { type: 'undoDraftEdit', key, list: 'comparison' });
    expect(activeDraft(undone).comparison).toEqual([]);
  });

  it('refuses comparison edits on a locked set', () => {
    const state = loaded({ ...WITH_COMPARISON, locked: true });
    const key = activeDraft(state).key;
    const edited = run(state, [
      { type: 'setDraftHeaders', key, headers: [], list: 'comparison' },
    ]);
    expect(activeDraft(edited).comparison).toEqual(['other_1', 'other_2']);
  });

  it('adopts the saved comparison list from the save response', () => {
    const state = loaded();
    const key = activeDraft(state).key;
    const edited = run(state, [
      { type: 'setDraftHeaders', key, headers: ['other_1'], list: 'comparison' },
      { type: 'focalDraftSaved', key, focalSet: { ...OLD_SET, comparisonHeaders: ['other_1'] } },
    ]);
    const draft = activeDraft(edited);
    expect(draft.savedComparison).toEqual(['other_1']);
    expect(isDirty(draft)).toBe(false);
  });

  it('sends both lists in one save request', () => {
    const state = loaded(WITH_COMPARISON);
    expect(saveRequestFor(activeDraft(state))).toEqual({
      focalSetId: 'set-1',
      title: 'Targets',
      headers: ['focal_1', 'focal_2'],
      comparisonHeaders: ['other_1', 'other_2'],
    });
  });
});

describe('overlap', () => {
  const draftWith = (focal: string[], comparison: string[]): FocalDraft => ({
    ...draftFromPayload(OLD_SET),
    headers: focal,
    comparison,
  });

  it('is exact complete-header equality', () => {
    expect(overlappingHeaders(draftWith(['a', 'b'], ['b', 'c']))).toEqual(['b']);
    // Substrings are not overlap.
    expect(overlappingHeaders(draftWith(['abc'], ['ab', 'abcd']))).toEqual([]);
  });

  it('ignores surrounding whitespace, as storage does', () => {
    expect(overlappingHeaders(draftWith([' a '], ['a']))).toEqual(['a']);
  });

  it('never removes the header from either list', () => {
    const state = loaded();
    const key = activeDraft(state).key;
    const edited = run(state, [
      { type: 'setDraftHeaders', key, headers: ['focal_1'], list: 'comparison' },
    ]);
    const draft = activeDraft(edited);
    expect(draft.headers).toContain('focal_1');
    expect(draft.comparison).toEqual(['focal_1']);
    expect(overlappingHeaders(draft)).toEqual(['focal_1']);
  });
});

/* ------------------------------------------------------------------ */
/* Run gate                                                            */
/* ------------------------------------------------------------------ */

const SOURCE = {
  fastaFileId: 'f1',
  displayName: 'a.fasta',
  available: true,
} as SourceStatusPayload;

function presence(states: Record<string, PresenceVerdict>): FocalPresenceState {
  const byHeader = Object.fromEntries(
    Object.entries(states).map(([header, state]) => [
      header,
      {
        header,
        state,
        occurrences: (state === 'present_current' ? { f1: 1 } : {}) as Record<string, number>,
      },
    ]),
  );
  return { status: 'loaded', byHeader, generation: 1, error: null };
}

const GREEN_FOCAL = { focal_1: 'present_current', focal_2: 'present_current' } as const;

function gate(draft: FocalDraft, states: Record<string, PresenceVerdict>, singleFile = true) {
  return evaluateRunGate({
    draft,
    presence: presence(states),
    singleFileScope: singleFile,
    scopeSources: [SOURCE],
    running: false,
  });
}

describe('the run gate with a comparison set', () => {
  it('runs with a blank comparison exactly as before', () => {
    expect(gate(draftFromPayload(OLD_SET), GREEN_FOCAL)).toEqual({ canRun: true, reason: null });
  });

  it('runs when every comparison entry is green', () => {
    const verdict = gate(draftFromPayload(WITH_COMPARISON), {
      ...GREEN_FOCAL,
      other_1: 'present_current',
      other_2: 'present_current',
    });
    expect(verdict.canRun).toBe(true);
  });

  it('blocks on focal/comparison overlap, before anything else about presence', () => {
    const overlapping = draftFromPayload({ ...OLD_SET, comparisonHeaders: ['focal_2'] });
    const verdict = gate(overlapping, GREEN_FOCAL);
    expect(verdict.canRun).toBe(false);
    expect(verdict.reason).toMatch(/focal_2 is in both the focal and the comparison set/);
  });

  it('counts several overlaps', () => {
    const overlapping = draftFromPayload({
      ...OLD_SET,
      comparisonHeaders: ['focal_1', 'focal_2'],
    });
    expect(gate(overlapping, GREEN_FOCAL).reason).toMatch(/^2 headers are in both/);
  });

  it('blocks on a comparison entry missing from the selected file', () => {
    const verdict = gate(draftFromPayload(WITH_COMPARISON), {
      ...GREEN_FOCAL,
      other_1: 'present_current',
      other_2: 'present_other',
    });
    expect(verdict.canRun).toBe(false);
    expect(verdict.reason).toMatch(/Comparison entry other_2 is not in the selected FASTA/);
  });

  it('blocks on a comparison entry missing from every selected file', () => {
    const verdict = gate(
      draftFromPayload(WITH_COMPARISON),
      { ...GREEN_FOCAL, other_1: 'present_current', other_2: 'missing' },
      false,
    );
    expect(verdict.reason).toMatch(/not in any of the selected FASTA files/);
  });

  it('blocks when a comparison entry cannot be located', () => {
    const verdict = gate(draftFromPayload(WITH_COMPARISON), {
      ...GREEN_FOCAL,
      other_1: 'unknown',
      other_2: 'present_current',
    });
    expect(verdict.reason).toMatch(/Comparison entry other_1 cannot be located/);
  });

  it('waits while comparison presence has not been answered', () => {
    const verdict = gate(draftFromPayload(WITH_COMPARISON), GREEN_FOCAL);
    expect(verdict.canRun).toBe(false);
    expect(verdict.reason).toMatch(/Checking/);
  });

  it('still reports focal problems first', () => {
    const verdict = gate(draftFromPayload(WITH_COMPARISON), {
      focal_1: 'missing',
      focal_2: 'present_current',
      other_1: 'missing',
      other_2: 'missing',
    });
    expect(verdict.reason).toMatch(/^focal_1 is not in the selected FASTA/);
  });

  it('blocks a dirty comparison edit like any other unsaved change', () => {
    const draft: FocalDraft = { ...draftFromPayload(OLD_SET), comparison: ['other_1'] };
    const verdict = gate(draft, { ...GREEN_FOCAL, other_1: 'present_current' });
    expect(verdict.reason).toMatch(/unsaved changes/);
  });
});
