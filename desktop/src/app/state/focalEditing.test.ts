import { describe, expect, it } from 'vitest';
import {
  appReducer,
  canEnterWorkspace,
  createInitialState,
  focalHistoryFor,
  selectedAnalyses,
} from './projectState';
import type { AppAction, AppState } from './projectState';

function setup() {
  const state = createInitialState();
  return { state, focalSetId: state.activeFocalSetId };
}

function run(state: AppState, actions: readonly AppAction[]): AppState {
  return actions.reduce(appReducer, state);
}

const strings = (state: AppState) => state.focalSets[0].strings;

describe('focal set — add and remove', () => {
  it('adds a string in add mode', () => {
    const { state, focalSetId } = setup();
    const next = appReducer(state, { type: 'addFocalString', focalSetId, value: 'Leptacis' });
    expect(strings(next)).toEqual(['Leptacis']);
  });

  it('trims whitespace around an added string', () => {
    const { state, focalSetId } = setup();
    const next = appReducer(state, { type: 'addFocalString', focalSetId, value: '  Leptacis  ' });
    expect(strings(next)).toEqual(['Leptacis']);
  });

  it('keeps entries as an array, preserving insertion order', () => {
    const { state, focalSetId } = setup();
    const next = run(state, [
      { type: 'addFocalString', focalSetId, value: 'b' },
      { type: 'addFocalString', focalSetId, value: 'a' },
      { type: 'addFocalString', focalSetId, value: 'c' },
    ]);
    expect(strings(next)).toEqual(['b', 'a', 'c']);
  });

  it('removes an exact entry in remove mode', () => {
    const { state, focalSetId } = setup();
    const next = run(state, [
      { type: 'addFocalString', focalSetId, value: 'Leptacis' },
      { type: 'addFocalString', focalSetId, value: 'Synopeas' },
      { type: 'removeFocalString', focalSetId, value: 'Leptacis' },
    ]);
    expect(strings(next)).toEqual(['Synopeas']);
  });

  it('removes only on an exact match, not a substring of an entry', () => {
    const { state, focalSetId } = setup();
    const added = appReducer(state, {
      type: 'addFocalString',
      focalSetId,
      value: 'Leptacis_tipulae',
    });
    const next = appReducer(added, { type: 'removeFocalString', focalSetId, value: 'Leptacis' });
    expect(strings(next)).toEqual(['Leptacis_tipulae']);
  });
});

describe('focal set — no-ops record no history', () => {
  it('ignores an empty add', () => {
    const { state, focalSetId } = setup();
    const next = appReducer(state, { type: 'addFocalString', focalSetId, value: '   ' });
    expect(strings(next)).toEqual([]);
    expect(focalHistoryFor(next, focalSetId).past).toHaveLength(0);
  });

  it('ignores a duplicate add and records no history entry', () => {
    const { state, focalSetId } = setup();
    const once = appReducer(state, { type: 'addFocalString', focalSetId, value: 'Leptacis' });
    const twice = appReducer(once, { type: 'addFocalString', focalSetId, value: 'Leptacis' });
    expect(strings(twice)).toEqual(['Leptacis']);
    expect(focalHistoryFor(twice, focalSetId).past).toHaveLength(1);
  });

  it('a removal that matches nothing is not destructive and adds no history', () => {
    const { state, focalSetId } = setup();
    const added = appReducer(state, { type: 'addFocalString', focalSetId, value: 'Leptacis' });
    const missed = appReducer(added, { type: 'removeFocalString', focalSetId, value: 'Nope' });

    expect(strings(missed)).toEqual(['Leptacis']);
    expect(focalHistoryFor(missed, focalSetId).past).toHaveLength(1);
    // The state object is returned unchanged, so nothing downstream re-renders.
    expect(missed).toBe(added);
  });
});

describe('focal set — undo and redo', () => {
  it('undoes one add at a time', () => {
    const { state, focalSetId } = setup();
    const next = run(state, [
      { type: 'addFocalString', focalSetId, value: 'a' },
      { type: 'addFocalString', focalSetId, value: 'b' },
      { type: 'undoFocalEdit', focalSetId },
    ]);
    expect(strings(next)).toEqual(['a']);
  });

  it('undoes a removal, restoring the entry', () => {
    const { state, focalSetId } = setup();
    const next = run(state, [
      { type: 'addFocalString', focalSetId, value: 'a' },
      { type: 'removeFocalString', focalSetId, value: 'a' },
      { type: 'undoFocalEdit', focalSetId },
    ]);
    expect(strings(next)).toEqual(['a']);
  });

  it('redoes an undone edit', () => {
    const { state, focalSetId } = setup();
    const next = run(state, [
      { type: 'addFocalString', focalSetId, value: 'a' },
      { type: 'addFocalString', focalSetId, value: 'b' },
      { type: 'undoFocalEdit', focalSetId },
      { type: 'redoFocalEdit', focalSetId },
    ]);
    expect(strings(next)).toEqual(['a', 'b']);
  });

  it('walks a whole sequence back and forward', () => {
    const { state, focalSetId } = setup();
    const built = run(state, [
      { type: 'addFocalString', focalSetId, value: 'a' },
      { type: 'addFocalString', focalSetId, value: 'b' },
      { type: 'removeFocalString', focalSetId, value: 'a' },
    ]);
    expect(strings(built)).toEqual(['b']);

    const back = run(built, [
      { type: 'undoFocalEdit', focalSetId },
      { type: 'undoFocalEdit', focalSetId },
      { type: 'undoFocalEdit', focalSetId },
    ]);
    expect(strings(back)).toEqual([]);

    const forward = run(back, [
      { type: 'redoFocalEdit', focalSetId },
      { type: 'redoFocalEdit', focalSetId },
      { type: 'redoFocalEdit', focalSetId },
    ]);
    expect(strings(forward)).toEqual(['b']);
  });

  it('a new edit after undo invalidates the redo branch', () => {
    const { state, focalSetId } = setup();
    const undone = run(state, [
      { type: 'addFocalString', focalSetId, value: 'a' },
      { type: 'addFocalString', focalSetId, value: 'b' },
      { type: 'undoFocalEdit', focalSetId },
    ]);
    expect(focalHistoryFor(undone, focalSetId).future).toHaveLength(1);

    const branched = appReducer(undone, { type: 'addFocalString', focalSetId, value: 'c' });
    expect(strings(branched)).toEqual(['a', 'c']);
    expect(focalHistoryFor(branched, focalSetId).future).toHaveLength(0);

    // Redo must now do nothing at all.
    expect(appReducer(branched, { type: 'redoFocalEdit', focalSetId })).toBe(branched);
  });

  it('undo and redo are inert when their stacks are empty', () => {
    const { state, focalSetId } = setup();
    expect(appReducer(state, { type: 'undoFocalEdit', focalSetId })).toBe(state);
    expect(appReducer(state, { type: 'redoFocalEdit', focalSetId })).toBe(state);
  });
});

describe('focal edit mode', () => {
  it('starts in add mode and is mutually exclusive', () => {
    const { state } = setup();
    expect(state.focalMode).toBe('add');

    const removing = appReducer(state, { type: 'setFocalMode', mode: 'remove' });
    expect(removing.focalMode).toBe('remove');

    const adding = appReducer(removing, { type: 'setFocalMode', mode: 'add' });
    expect(adding.focalMode).toBe('add');
  });
});

describe('navigation gates', () => {
  it('blocks the workspace until a FASTA is selected', () => {
    const { state } = setup();
    expect(state.draft.analyses.molecularDiagnosis).toBe(true);
    expect(canEnterWorkspace(state)).toBe(false);

    const withFasta = appReducer(state, { type: 'setFastaPath', path: '/tmp/x.fasta' });
    expect(canEnterWorkspace(withFasta)).toBe(true);
  });

  it('blocks the workspace when Molecular Diagnosis is not selected', () => {
    const { state } = setup();
    const next = run(state, [
      { type: 'setFastaPath', path: '/tmp/x.fasta' },
      { type: 'toggleAnalysis', analysis: 'molecularDiagnosis' },
    ]);
    expect(canEnterWorkspace(next)).toBe(false);
  });

  it('clearing the FASTA drops the loaded headers', () => {
    const { state } = setup();
    const loaded = run(state, [
      { type: 'setFastaPath', path: '/tmp/x.fasta' },
      {
        type: 'alignmentLoaded',
        data: { path: '/tmp/x.fasta', headers: ['a'], duplicateCount: 0 },
      },
    ]);
    expect(loaded.alignment.status).toBe('loaded');

    const cleared = appReducer(loaded, { type: 'setFastaPath', path: null });
    expect(cleared.alignment.status).toBe('idle');
  });
});

describe('analysis selection drives the tab list', () => {
  it('lists only selected analyses, in canonical order', () => {
    const { state } = setup();
    expect(selectedAnalyses(state.draft.analyses)).toEqual(['molecularDiagnosis']);

    const withConsensus = appReducer(state, {
      type: 'toggleAnalysis',
      analysis: 'consensusSequenceGeneration',
    });
    expect(selectedAnalyses(withConsensus.draft.analyses)).toEqual([
      'molecularDiagnosis',
      'consensusSequenceGeneration',
    ]);
  });

  it('is empty when nothing is selected', () => {
    const { state } = setup();
    const none = appReducer(state, { type: 'toggleAnalysis', analysis: 'molecularDiagnosis' });
    expect(selectedAnalyses(none.draft.analyses)).toEqual([]);
  });

  it('refuses to activate an analysis that is not selected', () => {
    const { state } = setup();
    const next = appReducer(state, {
      type: 'setActiveAnalysis',
      analysis: 'sequencePunishmentTest',
    });
    expect(next.activeAnalysis).toBe('molecularDiagnosis');
  });

  it('moves the active tab off an analysis that gets deselected', () => {
    const { state } = setup();
    const next = run(state, [
      { type: 'toggleAnalysis', analysis: 'sequencePunishmentTest' },
      { type: 'setActiveAnalysis', analysis: 'sequencePunishmentTest' },
      { type: 'toggleAnalysis', analysis: 'sequencePunishmentTest' },
    ]);
    expect(next.activeAnalysis).toBe('molecularDiagnosis');
  });
});

describe('project draft survives navigation', () => {
  it('keeps name, FASTA, analyses and focal strings across screens', () => {
    const { state, focalSetId } = setup();
    const next = run(state, [
      { type: 'navigate', screen: 'projectCreation' },
      { type: 'setProjectName', name: 'European Leptacis' },
      { type: 'setFastaPath', path: '/tmp/x.fasta' },
      { type: 'addFocalString', focalSetId, value: 'Leptacis_tipulae' },
      { type: 'navigate', screen: 'workspace' },
      { type: 'navigate', screen: 'projectCreation' },
      { type: 'navigate', screen: 'workspace' },
    ]);

    expect(next.draft.name).toBe('European Leptacis');
    expect(next.draft.fastaPath).toBe('/tmp/x.fasta');
    expect(strings(next)).toEqual(['Leptacis_tipulae']);
  });
});
