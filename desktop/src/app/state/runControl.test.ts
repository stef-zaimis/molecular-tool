import { describe, expect, it } from 'vitest';
import { appReducer, createInitialState, runInputKey } from './projectState';
import type { AppAction, AppState } from './projectState';
import type {
  BackendError,
  FocalSetPayload,
  OpenProjectResult,
  ProjectDiagnosisResult,
  SourceStatusPayload,
} from '../../backendContract';

/**
 * Run lifecycle rules in the reducer: Stop, the cancelled outcome, results
 * applying only to their own run, and what a continuation is bound to.
 */

function run(state: AppState, actions: readonly AppAction[]): AppState {
  return actions.reduce(appReducer, state);
}

const RESULT = {
  focalSetId: 'set-1',
  focalHeaders: ['f'],
  fastaFileIds: ['f1'],
  sequenceCount: 4,
  alignmentLength: 8,
  outputs: { reportTxt: '/r.txt', workbookXlsx: '/w.xlsx', consensusTxt: null },
  dmc: {} as ProjectDiagnosisResult['dmc'],
  canContinue: true,
  resume: {
    startCombinationLength: 3,
    diagnosticCombinations: [],
    combinationsTestedByLength: {},
    inputsFingerprint: 'abc',
  },
} as ProjectDiagnosisResult;

const CANCELLED: BackendError = {
  code: 'RUN_CANCELLED',
  message: 'Analysis stopped.',

};

const started = (runToken: string): AppAction => ({
  type: 'diagnosisStarted',
  continuing: false,
  runToken,
  startedAt: 0,
  inputKey: 'key',
});

describe('Stop', () => {
  it('marks only the run in flight as stopping, once', () => {
    const state = run(createInitialState(), [
      started('A'),
      { type: 'diagnosisStopRequested', runToken: 'A' },
    ]);
    expect(state.diagnosisRun).toMatchObject({ status: 'running', stopping: true });
    // A repeat is a no-op (same object back).
    expect(appReducer(state, { type: 'diagnosisStopRequested', runToken: 'A' })).toBe(state);
    // A stop for some other run changes nothing.
    const other = run(createInitialState(), [
      started('A'),
      { type: 'diagnosisStopRequested', runToken: 'B' },
    ]);
    expect(other.diagnosisRun).toMatchObject({ stopping: false });
  });

  it('becomes running again if the stop request could not be delivered', () => {
    const state = run(createInitialState(), [
      started('A'),
      { type: 'diagnosisStopRequested', runToken: 'A' },
      { type: 'diagnosisStopFailed', runToken: 'A' },
    ]);
    expect(state.diagnosisRun).toMatchObject({ status: 'running', stopping: false });
  });

  it('ends in its own cancelled state, not as a failure and not as a result', () => {
    const state = run(createInitialState(), [
      started('A'),
      { type: 'diagnosisStopRequested', runToken: 'A' },
      { type: 'diagnosisFailed', runToken: 'A', error: CANCELLED },
    ]);
    expect(state.diagnosisRun).toEqual({ status: 'cancelled' });
  });

  it('still reports a genuine failure as a failure', () => {
    const state = run(createInitialState(), [
      started('A'),
      {
        type: 'diagnosisFailed',
        runToken: 'A',
        error: { code: 'BACKEND_UNAVAILABLE', message: 'gone' },
      },
    ]);
    expect(state.diagnosisRun).toMatchObject({ status: 'failed' });
  });

  it('allows a new run after a cancelled one', () => {
    const state = run(createInitialState(), [
      started('A'),
      { type: 'diagnosisFailed', runToken: 'A', error: CANCELLED },
      started('B'),
      { type: 'diagnosisSucceeded', runToken: 'B', result: RESULT },
    ]);
    expect(state.diagnosisRun).toMatchObject({ status: 'succeeded' });
  });
});

describe('a run ending applies only to that run', () => {
  it('ignores a late result or failure from an earlier run', () => {
    const base = run(createInitialState(), [started('A'), started('B')]);
    expect(appReducer(base, { type: 'diagnosisSucceeded', runToken: 'A', result: RESULT })).toBe(
      base,
    );
    expect(appReducer(base, { type: 'diagnosisFailed', runToken: 'A', error: CANCELLED })).toBe(
      base,
    );
  });

  it('ignores an ending when nothing is running', () => {
    const done = run(createInitialState(), [
      started('A'),
      { type: 'diagnosisFailed', runToken: 'A', error: CANCELLED },
    ]);
    expect(appReducer(done, { type: 'diagnosisSucceeded', runToken: 'A', result: RESULT })).toBe(
      done,
    );
  });

  it('ignores progress from a cancelled run', () => {
    const done = run(createInitialState(), [
      started('A'),
      { type: 'diagnosisFailed', runToken: 'A', error: CANCELLED },
    ]);
    const next = appReducer(done, {
      type: 'diagnosisProgress',
      progress: {
        runId: 'x',
        runToken: 'A',
        stage: 'five_site',
        current: 1,
        total: 2,
        detail: null,
        elapsedMs: 1,
      },
    });
    expect(next).toBe(done);
  });

  it('keeps the input key of the run that produced the continuation', () => {
    const state = run(createInitialState(), [
      { ...started('A'), inputKey: 'the-inputs' } as AppAction,
      { type: 'diagnosisSucceeded', runToken: 'A', result: RESULT },
    ]);
    expect(state.diagnosisRun).toMatchObject({ inputKey: 'the-inputs' });
  });
});

/* ------------------------------------------------------------------ */
/* What a continuation is bound to                                      */
/* ------------------------------------------------------------------ */

function source(id: string): SourceStatusPayload {
  return {
    fastaFileId: id,
    sourcePath: `/data/${id}.fasta`,
    displayName: `${id}.fasta`,
    state: 'current',
    available: true,
    indexUsable: true,
    exists: true,
    currentSizeBytes: 1,
    currentMtimeNs: 1,
    indexedSizeBytes: 1,
    indexedMtimeNs: 1,
    indexRevision: 1,
    sequenceCount: 4,
    alignmentLength: 8,
    duplicateHeaderCount: 0,
    locked: false,
    message: null,
  };
}

const OPEN: OpenProjectResult = {
  projectDir: '/p',
  outputsDir: '/p/o',
  metadata: { projectUuid: 'u', title: 'p' },
  capabilities: { sqliteVersion: '3', fts5: true, trigram: true, acceleratedSearch: true },
  sources: [source('f1'), source('f2')],
};

const SET_A: FocalSetPayload = {
  id: 'set-a',
  title: 'A',
  locked: false,
  entries: [{ id: 'e', header: 'focal_1' }],
  comparisonHeaders: ['other_1'],
};
const SET_B: FocalSetPayload = { ...SET_A, id: 'set-b', entries: [{ id: 'e', header: 'focal_2' }] };

function opened(): AppState {
  return run(createInitialState(), [
    { type: 'projectOpened', result: OPEN },
    { type: 'focalSetsLoaded', focalSets: [SET_A, SET_B] },
  ]);
}

describe('runInputKey', () => {
  const key = opened();
  const original = runInputKey(key);

  it('ignores the maximum candidate size, which Continue is meant to raise', () => {
    const raised = appReducer(key, {
      type: 'updateMolecularDiagnosis',
      patch: { maxCandidateSize: 9 },
    });
    expect(runInputKey(raised)).toBe(original);
  });

  it.each<[string, AppAction]>([
    ['another focal set', { type: 'selectFocalDraft', key: 'set-set-b' }],
    ['another FASTA scope', { type: 'setFastaScope', scope: { kind: 'all' } }],
    ['Ignore Gaps', { type: 'updateMolecularDiagnosis', patch: { ignoreGaps: true } }],
    [
      'benefit of the doubt',
      { type: 'updateMolecularDiagnosis', patch: { giveBenefitOfDoubtToAmbiguousBases: true } },
    ],
    ['the minimum size', { type: 'updateMolecularDiagnosis', patch: { minCandidateSize: 2 } }],
    [
      'saved comparison membership',
      { type: 'focalDraftSaved', key: 'set-set-a', focalSet: { ...SET_A, comparisonHeaders: [] } },
    ],
    [
      'saved focal membership',
      {
        type: 'focalDraftSaved',
        key: 'set-set-a',
        focalSet: { ...SET_A, entries: [{ id: 'x', header: 'focal_9' }] },
      },
    ],
  ])('changes with %s', (_label, action) => {
    expect(runInputKey(appReducer(key, action))).not.toBe(original);
  });

  it('does not change with unsaved typing (the saved set is what runs)', () => {
    const typed = appReducer(key, {
      type: 'setDraftHeaders',
      key: 'set-set-a',
      headers: ['focal_1', 'more'],
    });
    expect(runInputKey(typed)).toBe(original);
  });
});
