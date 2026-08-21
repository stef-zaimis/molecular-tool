import { describe, expect, it } from 'vitest';
import {
  DEFAULT_VIEWER_SNAP,
  VIEWER_SNAPS,
  activeDraft,
  appReducer,
  createInitialState,
  scopeFileIds,
  selectedAnalyses,
} from './projectState';
import type { AppAction, AppState } from './projectState';
import {
  blankDraft,
  draftFromPayload,
  isDirty,
  isNewDraft,
  lockBlockedReason,
  normaliseHeaders,
  saveBlockedReason,
  saveRequestFor,
  withHeaders as withHeadersFor,
} from './focalDrafts';
import type { FocalSetPayload, SourceStatusPayload } from '../../backendContract';

/**
 * Focal editing as SAVED DATA plus a WORKING COPY.
 *
 * The distinction is what these tests are about: editing must never write, and
 * the difference between "what is on screen" and "what a run would use" must be
 * derivable at any moment.
 */

function run(state: AppState, actions: readonly AppAction[]): AppState {
  return actions.reduce(appReducer, state);
}

function headers(state: AppState): readonly string[] {
  return activeDraft(state).headers;
}

const SAVED_SET: FocalSetPayload = {
  id: 'set-1',
  title: 'Targets',
  locked: false,
  entries: [
    { id: 'e1', header: 'focal_1|AU|Target' },
    { id: 'e2', header: 'focal_2|GB|Target' },
  ],
};

function source(id: string, displayName: string): SourceStatusPayload {
  return {
    fastaFileId: id,
    sourcePath: `/data/${displayName}`,
    displayName,
    state: 'current',
    available: true,
    indexUsable: true,
    exists: true,
    currentSizeBytes: 10,
    currentMtimeNs: 1,
    indexedSizeBytes: 10,
    indexedMtimeNs: 1,
    indexRevision: 1,
    sequenceCount: 4,
    alignmentLength: 8,
    duplicateHeaderCount: 0,
    locked: false,
    message: null,
  };
}

/** A state with one saved focal set loaded, as after opening a project. */
function withSavedSet(): AppState {
  return appReducer(createInitialState(), { type: 'focalSetsLoaded', focalSets: [SAVED_SET] });
}

describe('normaliseHeaders', () => {
  it('trims, drops blanks and collapses duplicates onto the first occurrence', () => {
    expect(normaliseHeaders([' b ', 'a', 'b', '', '   '])).toEqual(['b', 'a']);
  });

  it('keeps whitespace inside a header, which FASTA headers legitimately have', () => {
    expect(normaliseHeaders(['gi|1 Homo sapiens'])).toEqual(['gi|1 Homo sapiens']);
  });
});

describe('a project with no focal sets', () => {
  it('gets one blank local draft rather than a database row', () => {
    const state = appReducer(createInitialState(), { type: 'focalSetsLoaded', focalSets: [] });
    const draft = activeDraft(state);

    expect(state.focalDrafts).toHaveLength(1);
    expect(draft.persistedId).toBeNull();
    expect(isNewDraft(draft)).toBe(true);
    expect(draft.headers).toEqual([]);
  });
});

describe('dirty is derived, not stored', () => {
  it('a freshly loaded set is clean', () => {
    expect(isDirty(activeDraft(withSavedSet()))).toBe(false);
  });

  it('a never-saved draft is always dirty, even when empty', () => {
    expect(isDirty(blankDraft())).toBe(true);
  });

  it('changing the title makes it dirty', () => {
    const state = withSavedSet();
    const next = appReducer(state, {
      type: 'setDraftTitle',
      key: activeDraft(state).key,
      title: 'Renamed',
    });
    expect(isDirty(activeDraft(next))).toBe(true);
  });

  it('ignores whitespace-only title differences', () => {
    const state = withSavedSet();
    const next = appReducer(state, {
      type: 'setDraftTitle',
      key: activeDraft(state).key,
      title: '  Targets  ',
    });
    expect(isDirty(activeDraft(next))).toBe(false);
  });

  it('reordering entries counts as a change, because the order is persisted', () => {
    const state = withSavedSet();
    const next = appReducer(state, {
      type: 'setDraftHeaders',
      key: activeDraft(state).key,
      headers: ['focal_2|GB|Target', 'focal_1|AU|Target'],
    });
    expect(isDirty(activeDraft(next))).toBe(true);
  });

  it('is clean again once the backend response is adopted', () => {
    const state = withSavedSet();
    const key = activeDraft(state).key;
    const edited = appReducer(state, {
      type: 'setDraftHeaders',
      key,
      headers: ['focal_1|AU|Target'],
    });
    expect(isDirty(activeDraft(edited))).toBe(true);

    const saved = appReducer(edited, {
      type: 'focalDraftSaved',
      key,
      focalSet: { ...SAVED_SET, entries: [{ id: 'e1', header: 'focal_1|AU|Target' }] },
    });

    expect(isDirty(activeDraft(saved))).toBe(false);
    expect(headers(saved)).toEqual(['focal_1|AU|Target']);
  });

  it('adopts the backend ordering rather than what the user typed', () => {
    const state = withSavedSet();
    const key = activeDraft(state).key;

    // The backend deduplicated and reordered; the working copy must follow it,
    // otherwise the screen would keep claiming unsaved changes forever.
    const saved = appReducer(state, {
      type: 'focalDraftSaved',
      key,
      focalSet: { ...SAVED_SET, entries: [{ id: 'e2', header: 'focal_2|GB|Target' }] },
    });

    expect(headers(saved)).toEqual(['focal_2|GB|Target']);
    expect(isDirty(activeDraft(saved))).toBe(false);
  });
});

describe('editing the working copy', () => {
  it('appends headers and deduplicates against what is already there', () => {
    const state = withSavedSet();
    const key = activeDraft(state).key;
    const next = appReducer(state, {
      type: 'setDraftHeaders',
      key,
      headers: [...headers(state), 'focal_1|AU|Target', 'other_1|FR|Contrast'],
    });

    expect(headers(next)).toEqual([
      'focal_1|AU|Target',
      'focal_2|GB|Target',
      'other_1|FR|Contrast',
    ]);
  });

  it('undoes and redoes one edit at a time', () => {
    const state = withSavedSet();
    const key = activeDraft(state).key;

    const edited = run(state, [
      { type: 'setDraftHeaders', key, headers: ['a'] },
      { type: 'setDraftHeaders', key, headers: ['a', 'b'] },
    ]);
    expect(headers(edited)).toEqual(['a', 'b']);

    const undone = appReducer(edited, { type: 'undoDraftEdit', key });
    expect(headers(undone)).toEqual(['a']);

    const redone = appReducer(undone, { type: 'redoDraftEdit', key });
    expect(headers(redone)).toEqual(['a', 'b']);
  });

  it('a fresh edit after undo invalidates the redo branch', () => {
    const state = withSavedSet();
    const key = activeDraft(state).key;

    const next = run(state, [
      { type: 'setDraftHeaders', key, headers: ['a'] },
      { type: 'undoDraftEdit', key },
      { type: 'setDraftHeaders', key, headers: ['z'] },
      { type: 'redoDraftEdit', key },
    ]);

    expect(headers(next)).toEqual(['z']);
  });

  it('a no-op edit records no history step', () => {
    const state = withSavedSet();
    const key = activeDraft(state).key;
    const same = appReducer(state, { type: 'setDraftHeaders', key, headers: [...headers(state)] });

    expect(activeDraft(same).history.past).toHaveLength(0);
  });
});

describe('locked focal sets', () => {
  const locked = draftFromPayload({ ...SAVED_SET, locked: true });

  function lockedState(): AppState {
    return appReducer(createInitialState(), {
      type: 'focalSetsLoaded',
      focalSets: [{ ...SAVED_SET, locked: true }],
    });
  }

  it('refuse edits in the reducer, not only in the backend', () => {
    const state = lockedState();
    const key = activeDraft(state).key;

    const next = run(state, [
      { type: 'setDraftHeaders', key, headers: ['anything'] },
      { type: 'setDraftTitle', key, title: 'Renamed' },
    ]);

    expect(headers(next)).toEqual(['focal_1|AU|Target', 'focal_2|GB|Target']);
    expect(activeDraft(next).title).toBe('Targets');
  });

  it('cannot be saved', () => {
    expect(saveBlockedReason(locked)).toMatch(/locked/i);
  });

  it('remain selectable and runnable', () => {
    // Nothing about being locked makes a set unusable — only uneditable.
    expect(isDirty(locked)).toBe(false);
    expect(locked.headers).toHaveLength(2);
  });
});

describe('switching drafts', () => {
  it('does not save the one being left', () => {
    const base = appReducer(createInitialState(), {
      type: 'focalSetsLoaded',
      focalSets: [SAVED_SET, { ...SAVED_SET, id: 'set-2', title: 'Second', entries: [] }],
    });
    const firstKey = base.focalDrafts[0].key;

    const next = run(base, [
      { type: 'setDraftHeaders', key: firstKey, headers: ['edited'] },
      { type: 'selectFocalDraft', key: base.focalDrafts[1].key },
    ]);

    const left = next.focalDrafts.find((draft) => draft.key === firstKey);
    expect(left?.headers).toEqual(['edited']);
    // Still dirty: leaving a draft is not a decision to commit it.
    expect(isDirty(left!)).toBe(true);
    expect(left?.savedHeaders).toEqual(['focal_1|AU|Target', 'focal_2|GB|Target']);
  });
});

describe('save requests', () => {
  it('carry the persisted id, a trimmed title and normalised headers', () => {
    const state = withSavedSet();
    const key = activeDraft(state).key;
    const edited = run(state, [
      { type: 'setDraftTitle', key, title: '  Renamed  ' },
      { type: 'setDraftHeaders', key, headers: [' a ', 'b', 'a'] },
    ]);

    expect(saveRequestFor(activeDraft(edited))).toEqual({
      focalSetId: 'set-1',
      title: 'Renamed',
      headers: ['a', 'b'],
    });
  });

  it('carry a null id for a set that has never been saved', () => {
    expect(saveRequestFor(blankDraft('New')).focalSetId).toBeNull();
  });

  it('are blocked without a title', () => {
    expect(saveBlockedReason(blankDraft(''))).toMatch(/title/i);
  });
});

describe('FASTA scope', () => {
  const sources = [source('f1', 'a.fasta'), source('f2', 'b.fasta')];

  it('defaults to the first linked file when a project opens', () => {
    const state = appReducer(createInitialState(), {
      type: 'projectOpened',
      result: {
        projectDir: '/p',
        outputsDir: '/p/outputs',
        metadata: { projectUuid: 'u', title: 'P' },
        capabilities: { sqliteVersion: '3', fts5: true, trigram: true, acceleratedSearch: true },
        sources,
      },
    });

    expect(state.fastaScope).toEqual({ kind: 'file', fastaFileId: 'f1' });
    expect(scopeFileIds(state)).toEqual(['f1']);
  });

  it('All files resolves to every linked file, in project order', () => {
    const state = run(createInitialState(), [
      { type: 'sourcesRefreshed', sources },
      { type: 'setFastaScope', scope: { kind: 'all' } },
    ]);

    expect(scopeFileIds(state)).toEqual(['f1', 'f2']);
  });

  it('falls back to the first file when the selected one is unlinked', () => {
    const state = run(createInitialState(), [
      { type: 'sourcesRefreshed', sources },
      { type: 'setFastaScope', scope: { kind: 'file', fastaFileId: 'f2' } },
      { type: 'sourceRemoved', fastaFileId: 'f2' },
    ]);

    // Never silently promoted to All files, which would widen the analysis.
    expect(state.fastaScope).toEqual({ kind: 'file', fastaFileId: 'f1' });
  });

  it('changing scope drops presence answers, which were scope-specific', () => {
    const state = run(createInitialState(), [
      { type: 'sourcesRefreshed', sources },
      { type: 'presenceRequested', generation: 1 },
      {
        type: 'presenceLoaded',
        generation: 1,
        entries: [{ header: 'a', state: 'present_current', occurrences: { f1: 1 } }],
        runEntries: [{ header: 'a', state: 'present_current', occurrences: { f1: 1 } }],
      },
      { type: 'setFastaScope', scope: { kind: 'file', fastaFileId: 'f2' } },
    ]);

    expect(state.focalPresence.byHeader).toEqual({});
    expect(state.focalPresence.status).toBe('idle');
  });

  it('re-selecting the SAME scope changes nothing at all', () => {
    /*
     * Clearing presence here would drop the colours while leaving the presence
     * effect's inputs identical — nothing would ask again, and the set would
     * sit neutral until something unrelated moved.
     */
    const base = run(createInitialState(), [
      { type: 'sourcesRefreshed', sources },
      { type: 'setFastaScope', scope: { kind: 'file', fastaFileId: 'f1' } },
      { type: 'presenceRequested', generation: 1 },
      {
        type: 'presenceLoaded',
        generation: 1,
        entries: [{ header: 'a', state: 'present_current', occurrences: { f1: 1 } }],
        runEntries: [{ header: 'a', state: 'present_current', occurrences: { f1: 1 } }],
      },
    ]);

    const again = appReducer(base, {
      type: 'setFastaScope',
      scope: { kind: 'file', fastaFileId: 'f1' },
    });

    expect(again).toBe(base);
    expect(again.focalPresence.byHeader.a.state).toBe('present_current');
  });
});

describe('presence responses', () => {
  it('are keyed by header so reordering cannot mis-colour entries', () => {
    const state = run(createInitialState(), [
      { type: 'presenceRequested', generation: 1 },
      {
        type: 'presenceLoaded',
        generation: 1,
        entries: [
          { header: 'a', state: 'present_current', occurrences: { f1: 1 } },
          { header: 'b', state: 'missing', occurrences: {} },
        ],
        runEntries: [
          { header: 'a', state: 'present_current', occurrences: { f1: 1 } },
          { header: 'b', state: 'missing', occurrences: {} },
        ],
      },
    ]);

    expect(state.focalPresence.byHeader.a.state).toBe('present_current');
    expect(state.focalPresence.byHeader.b.state).toBe('missing');
  });

  it('ignores a stale response from an older generation', () => {
    const state = run(createInitialState(), [
      { type: 'presenceRequested', generation: 1 },
      { type: 'presenceRequested', generation: 2 },
      {
        type: 'presenceLoaded',
        generation: 1,
        entries: [{ header: 'stale', state: 'present_current', occurrences: {} }],
        runEntries: [{ header: 'stale', state: 'present_current', occurrences: {} }],
      },
    ]);

    expect(state.focalPresence.byHeader).toEqual({});
    expect(state.focalPresence.status).toBe('checking');
  });
});

describe('analysis selection drives the tab list', () => {
  it('lists only selected analyses, in canonical order', () => {
    const state = run(createInitialState(), [
      { type: 'toggleAnalysis', analysis: 'consensusSequenceGeneration' },
    ]);

    expect(selectedAnalyses(state.analyses)).toEqual([
      'molecularDiagnosis',
      'consensusSequenceGeneration',
    ]);
  });

  it('refuses to activate an analysis that is not selected', () => {
    const state = appReducer(createInitialState(), {
      type: 'setActiveAnalysis',
      analysis: 'sequencePunishmentTest',
    });
    expect(state.activeAnalysis).toBe('molecularDiagnosis');
  });

  it('moves the active tab off an analysis that gets deselected', () => {
    const state = run(createInitialState(), [
      { type: 'toggleAnalysis', analysis: 'sequencePunishmentTest' },
      { type: 'setActiveAnalysis', analysis: 'sequencePunishmentTest' },
      { type: 'toggleAnalysis', analysis: 'sequencePunishmentTest' },
    ]);

    expect(state.activeAnalysis).toBe('molecularDiagnosis');
  });
});

describe('the project-to-be', () => {
  /** Candidates, not bare paths: each was vetted by the backend first. */
  const candidate = (path: string) => ({
    path,
    displayName: path.replace('/', ''),
    sequenceCount: 3,
    alignmentLength: 4,
    duplicateHeaderCount: 0,
  });

  it('collects vetted candidates without duplicates', () => {
    const state = run(createInitialState(), [
      { type: 'addPendingCandidates', candidates: [candidate('/a.fasta'), candidate('/b.fasta')] },
      { type: 'addPendingCandidates', candidates: [candidate('/b.fasta'), candidate('/c.fasta')] },
    ]);

    expect(state.pending.candidates.map((item) => item.path)).toEqual([
      '/a.fasta',
      '/b.fasta',
      '/c.fasta',
    ]);
  });

  it('drops one candidate without disturbing the others', () => {
    const state = run(createInitialState(), [
      { type: 'addPendingCandidates', candidates: [candidate('/a.fasta'), candidate('/b.fasta')] },
      { type: 'removePendingCandidate', path: '/a.fasta' },
    ]);

    expect(state.pending.candidates.map((item) => item.path)).toEqual(['/b.fasta']);
  });

  it('carries the metadata the shared FASTA table shows', () => {
    const state = appReducer(createInitialState(), {
      type: 'addPendingCandidates',
      candidates: [candidate('/a.fasta')],
    });

    expect(state.pending.candidates[0]).toMatchObject({ sequenceCount: 3, alignmentLength: 4 });
  });
});


describe('typing bursts collapse into one undo step', () => {
  function typed(state: AppState, key: string, values: readonly (readonly string[])[]): AppState {
    // Each entry is one keystroke's worth of membership, as the editor reports it.
    return values.reduce(
      (acc, headers) => appReducer(acc, { type: 'setDraftHeaders', key, headers, coalesce: true }),
      state,
    );
  }

  it('records one step for a continuous burst, not one per keystroke', () => {
    const state = withSavedSet();
    const key = activeDraft(state).key;

    const next = typed(state, key, [['a'], ['ab'], ['abc']]);

    expect(headers(next)).toEqual(['abc']);
    expect(activeDraft(next).history.past).toHaveLength(1);

    // One undo returns the whole burst to where it started.
    const undone = appReducer(next, { type: 'undoDraftEdit', key });
    expect(headers(undone)).toEqual(['focal_1|AU|Target', 'focal_2|GB|Target']);
  });

  it('starts a new step after the burst is closed', () => {
    const state = withSavedSet();
    const key = activeDraft(state).key;

    const next = run(state, [
      { type: 'setDraftHeaders', key, headers: ['a'], coalesce: true },
      { type: 'setDraftHeaders', key, headers: ['ab'], coalesce: true },
      { type: 'endDraftBurst', key },
      { type: 'setDraftHeaders', key, headers: ['abc'], coalesce: true },
    ]);

    expect(activeDraft(next).history.past).toHaveLength(2);
    expect(headers(appReducer(next, { type: 'undoDraftEdit', key }))).toEqual(['ab']);
  });

  it('keeps `+` and `−` as one step each, never coalesced into typing', () => {
    const state = withSavedSet();
    const key = activeDraft(state).key;

    const next = run(state, [
      { type: 'setDraftHeaders', key, headers: ['a'], coalesce: true },
      // A `+` arrives without `coalesce`, so it closes the burst and records
      // its own step.
      { type: 'setDraftHeaders', key, headers: ['a', 'b'] },
      { type: 'setDraftHeaders', key, headers: ['a', 'b', 'c'], coalesce: true },
    ]);

    expect(activeDraft(next).history.past).toHaveLength(3);
  });

  it('switching drafts closes an open burst', () => {
    const base = appReducer(createInitialState(), {
      type: 'focalSetsLoaded',
      focalSets: [SAVED_SET, { ...SAVED_SET, id: 'set-2', title: 'Second', entries: [] }],
    });
    const first = base.focalDrafts[0].key;

    const next = run(base, [
      { type: 'setDraftHeaders', key: first, headers: ['a'], coalesce: true },
      { type: 'selectFocalDraft', key: base.focalDrafts[1].key },
    ]);

    expect(next.focalDrafts.find((d) => d.key === first)?.burstOpen).toBe(false);
  });
});

describe('the focal-set library', () => {
  function twoSets(): AppState {
    return appReducer(createInitialState(), {
      type: 'focalSetsLoaded',
      focalSets: [
        SAVED_SET,
        { ...SAVED_SET, id: 'set-2', title: 'Second', entries: [{ id: 'x', header: 'h' }] },
      ],
    });
  }

  it('adopts only the lock from a lock response, never the membership', () => {
    const state = twoSets();
    const key = state.focalDrafts[0].key;

    const edited = appReducer(state, { type: 'setDraftHeaders', key, headers: ['edited'] });
    const locked = appReducer(edited, {
      type: 'focalDraftLockChanged',
      key,
      focalSet: { ...SAVED_SET, locked: true },
    });

    const draft = locked.focalDrafts.find((d) => d.key === key);
    expect(draft?.locked).toBe(true);
    // The unsaved edit is still there: locking is not a save.
    expect(draft?.headers).toEqual(['edited']);
  });

  it('removing a draft selects a neighbour', () => {
    const state = twoSets();
    const [first, second] = state.focalDrafts;

    const next = appReducer(state, { type: 'removeFocalDraft', key: first.key });

    expect(next.focalDrafts).toHaveLength(1);
    expect(next.activeFocalKey).toBe(second.key);
  });

  it('removing the last draft leaves one blank LOCAL draft, not a database row', () => {
    const state = appReducer(createInitialState(), {
      type: 'focalSetsLoaded',
      focalSets: [SAVED_SET],
    });

    const next = appReducer(state, { type: 'removeFocalDraft', key: state.focalDrafts[0].key });

    expect(next.focalDrafts).toHaveLength(1);
    expect(next.focalDrafts[0].persistedId).toBeNull();
    expect(next.focalDrafts[0].headers).toEqual([]);
    expect(next.activeFocalKey).toBe(next.focalDrafts[0].key);
  });

  it('removing an inactive draft leaves the selection alone', () => {
    const state = twoSets();
    const [first, second] = state.focalDrafts;

    const next = appReducer(state, { type: 'removeFocalDraft', key: second.key });

    expect(next.activeFocalKey).toBe(first.key);
    expect(next.focalPresence).toBe(state.focalPresence);
  });

  it('the pencil marks a draft for renaming, and saving clears it', () => {
    const state = twoSets();
    const key = state.focalDrafts[0].key;

    const renaming = appReducer(state, { type: 'setTitleEditing', key });
    expect(renaming.titleEditingKey).toBe(key);

    const saved = appReducer(renaming, { type: 'focalDraftSaved', key, focalSet: SAVED_SET });
    expect(saved.titleEditingKey).toBeNull();
  });

  it('a new draft opens with its title ready to type', () => {
    const next = appReducer(twoSets(), { type: 'addFocalDraft' });
    expect(next.titleEditingKey).toBe(next.activeFocalKey);
  });
});

describe('locking is refused for anything unsaved', () => {
  it('explains that a new draft must be saved first', () => {
    expect(lockBlockedReason(blankDraft('New'))).toMatch(/save/i);
  });

  it('explains that a dirty draft must be saved first', () => {
    const dirtyDraft = withHeadersFor(draftFromPayload(SAVED_SET), ['changed']);
    expect(lockBlockedReason(dirtyDraft)).toMatch(/save/i);
  });

  it('allows a saved, clean draft', () => {
    expect(lockBlockedReason(draftFromPayload(SAVED_SET))).toBeNull();
  });
});

describe('the sequence visualizer position', () => {
  it('starts at the design\'s default snap', () => {
    expect(createInitialState().viewerSnap).toBe(DEFAULT_VIEWER_SNAP);
  });

  it('clamps to the available snaps, so the handle stays reachable', () => {
    const low = appReducer(createInitialState(), { type: 'setViewerSnap', snap: -5 });
    expect(low.viewerSnap).toBe(0);

    const high = appReducer(createInitialState(), { type: 'setViewerSnap', snap: 99 });
    expect(high.viewerSnap).toBe(VIEWER_SNAPS.length - 1);
  });

  it('toggles the library without touching the viewer position', () => {
    const shown = appReducer(createInitialState(), { type: 'setLibraryVisible', visible: true });
    expect(shown.libraryVisible).toBe(true);
    expect(shown.viewerSnap).toBe(DEFAULT_VIEWER_SNAP);
  });
});


describe('re-selecting the active focal set', () => {
  it('is a strict no-op', () => {
    const base = run(
      appReducer(createInitialState(), { type: 'focalSetsLoaded', focalSets: [SAVED_SET] }),
      [
        { type: 'presenceRequested', generation: 1 },
        {
          type: 'presenceLoaded',
          generation: 1,
          entries: [{ header: 'focal_1|AU|Target', state: 'present_current', occurrences: {} }],
          runEntries: [
            { header: 'focal_1|AU|Target', state: 'present_current', occurrences: {} },
          ],
        },
        { type: 'setTitleEditing', key: activeDraft(createInitialState()).key },
      ],
    );

    const withEditing = appReducer(base, {
      type: 'setTitleEditing',
      key: base.activeFocalKey,
    });

    const again = appReducer(withEditing, {
      type: 'selectFocalDraft',
      key: withEditing.activeFocalKey,
    });

    // Same object: nothing cleared, nothing recomputed, no colour flicker.
    expect(again).toBe(withEditing);
    expect(again.focalPresence.byHeader['focal_1|AU|Target'].state).toBe('present_current');
    expect(again.titleEditingKey).toBe(withEditing.activeFocalKey);
  });

  it('still switches to a different draft', () => {
    const base = appReducer(createInitialState(), {
      type: 'focalSetsLoaded',
      focalSets: [SAVED_SET, { ...SAVED_SET, id: 'set-2', title: 'Second', entries: [] }],
    });

    const next = appReducer(base, {
      type: 'selectFocalDraft',
      key: base.focalDrafts[1].key,
    });

    expect(next.activeFocalKey).toBe(base.focalDrafts[1].key);
  });
});

describe('the sequence visualizer and the library are independent', () => {
  it('can be shown and hidden separately', () => {
    const both = run(createInitialState(), [
      { type: 'setLibraryVisible', visible: true },
      { type: 'setViewerVisible', visible: true },
    ]);
    expect(both.libraryVisible).toBe(true);
    expect(both.viewerVisible).toBe(true);

    const libraryOnly = appReducer(both, { type: 'setViewerVisible', visible: false });
    expect(libraryOnly.libraryVisible).toBe(true);
    expect(libraryOnly.viewerVisible).toBe(false);

    const neither = appReducer(libraryOnly, { type: 'setLibraryVisible', visible: false });
    expect(neither.libraryVisible).toBe(false);
    expect(neither.viewerVisible).toBe(false);
  });

  it('hiding the viewer preserves its width for when it returns', () => {
    const widened = run(createInitialState(), [
      { type: 'setViewerSnap', snap: 0 },
      { type: 'setViewerVisible', visible: false },
    ]);
    expect(widened.viewerSnap).toBe(0);

    const back = appReducer(widened, { type: 'setViewerVisible', visible: true });
    expect(back.viewerSnap).toBe(0);
  });

  it('offers a full-expansion snap reaching the workspace edge', () => {
    expect(VIEWER_SNAPS[0]).toBe(0);
    // And never travels outside the window at the other end.
    expect(VIEWER_SNAPS[VIEWER_SNAPS.length - 1]).toBeLessThan(1);
  });
});
