import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act, cleanup, render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { App } from './App';
import { ProjectProvider } from './state/ProjectContext';
import type {
  BackendResult,
  DiagnosisProgress,
  DiagnosisResumeState,
  FocalSetPayload,
  OpenProjectResult,
  ProjectDiagnosisResult,
  SourceStatusPayload,
} from '../backendContract';

/**
 * Stop and Continue in the mounted workspace, against a stubbed backend.
 *
 * The stub plays the real protocol's part: `cancelDiagnosis` only records the
 * request, and the RUN's own promise settles with RUN_CANCELLED later — the
 * run is not over until then. (The real out-of-band path is covered on the
 * Python side by tests/test_run_cancellation.py, which drives the process.)
 */

const SOURCE: SourceStatusPayload = {
  fastaFileId: 'f1',
  sourcePath: '/data/a.fasta',
  displayName: 'a.fasta',
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

const OPEN_RESULT: OpenProjectResult = {
  projectDir: '/p',
  outputsDir: '/p/outputs',
  metadata: { projectUuid: 'u', title: 'Stop project' },
  capabilities: { sqliteVersion: '3', fts5: true, trigram: true, acceleratedSearch: true },
  sources: [SOURCE],
};

const SAVED_SET: FocalSetPayload = {
  id: 'set-1',
  title: 'Targets',
  locked: false,
  entries: [{ id: 'e1', header: 'focal_1|AU|Target' }],
};

const RESUME: DiagnosisResumeState = {
  startCombinationLength: 3,
  diagnosticCombinations: [],
  combinationsTestedByLength: { '1': 4, '2': 6 },
  inputsFingerprint: 'fp-1',
};

function result(overrides: Partial<ProjectDiagnosisResult> = {}): ProjectDiagnosisResult {
  return {
    focalSetId: 'set-1',
    focalHeaders: ['focal_1|AU|Target'],
    fastaFileIds: ['f1'],
    sequenceCount: 4,
    alignmentLength: 8,
    outputs: {
      reportTxt: '/p/outputs/DMCs_output.txt',
      workbookXlsx: '/p/outputs/comparison_output.xlsx',
      consensusTxt: null,
    },
    dmc: {
      combinationsByLength: {},
      singleSites: [],
      pairs: [],
      uniqueSites: [],
      states: {},
      stopReason: 'reached_maximum_length',
      stoppedAtLength: 2,
      minCombinationLength: 1,
      maxCombinationLength: 2,
      startCombinationLength: 1,
      diagnostics: {
        fixedCount: 4,
        skippedSites: 0,
        globallyConservedRemoved: 0,
        candidateCount: 4,
        pairsTested: 6,
        totalCombinationsTested: 10,
        combinationsTestedByLength: { '1': 4, '2': 6 },
        ambiguousBdSitesIncluded: [],
        gappyConsensusSitesIncluded: [],
      },
    },
    canContinue: true,
    resume: RESUME,
    ...overrides,
  };
}

const ok = <T,>(value: T): BackendResult<T> => ({ ok: true, result: value });
const cancelled: BackendResult<ProjectDiagnosisResult> = {
  ok: false,
  error: { code: 'RUN_CANCELLED', message: 'Analysis stopped.' },
};

function makeBackend() {
  /** One settle function per run request, in order. */
  const settlers: Array<(value: BackendResult<ProjectDiagnosisResult>) => void> = [];
  const runRequests: Array<{ runToken?: string; resume?: DiagnosisResumeState | null }> = [];
  const cancelTokens: string[] = [];
  let cancelResponse: BackendResult<{ runToken: string; cancelRequested: boolean }> | null = null;

  const api = {
    window: {
      minimize: vi.fn(),
      toggleMaximize: vi.fn(),
      close: vi.fn(),
      isMaximized: () => Promise.resolve(false),
      onMaximizedChanged: () => () => undefined,
    },
    shell: {
      enterWorkspaceLayout: vi.fn(),
      enterLauncherLayout: vi.fn(),
      showItemInFolder: () => Promise.resolve(true),
    },
    dialog: {
      selectFastaFile: () => Promise.resolve(null),
      selectFastaFiles: () => Promise.resolve([]),
      exportFocalSet: () => Promise.resolve({ ok: true as const, path: '/x' }),
    },
    analysis: {
      ping: vi.fn(),
      loadFasta: vi.fn(),
      validateFocalStrings: vi.fn(),
      runMolecularDiagnosis: vi.fn(),
    },
    project: {
      create: () => Promise.resolve(ok(OPEN_RESULT)),
      open: () => Promise.resolve(ok(OPEN_RESULT)),
      close: () => Promise.resolve(ok({ closed: true })),
      setTitle: vi.fn(),
      refreshSources: () => Promise.resolve(ok({ sources: [SOURCE] })),
      validateFastaCandidate: vi.fn(),
      linkFasta: vi.fn(),
      setFastaFileLocked: vi.fn(),
      unlinkFasta: vi.fn(),
      relinkFasta: vi.fn(),
      reindexFasta: vi.fn(),
      searchHeaders: vi.fn(),
      resolveFocalAddQuery: (query: string) => Promise.resolve(ok({ query, headers: [] })),
      listFocalSets: () => Promise.resolve(ok({ focalSets: [SAVED_SET] })),
      getFocalSet: vi.fn(),
      createFocalSet: vi.fn(),
      renameFocalSet: vi.fn(),
      setFocalSetLocked: vi.fn(),
      deleteFocalSet: vi.fn(),
      replaceFocalEntries: vi.fn(),
      replaceFocalEntryHeaders: vi.fn(),
      saveFocalSet: vi.fn(),
      headerPresence: (headers: readonly string[]) =>
        Promise.resolve(
          ok({
            entries: headers.map((header) => ({
              header,
              state: 'present_current' as const,
              occurrences: { f1: 1 },
            })),
          }),
        ),
      matchFocalHeaders: vi.fn(),
      runMolecularDiagnosis: (request: {
        runToken?: string;
        resume?: DiagnosisResumeState | null;
      }) => {
        runRequests.push(request);
        return new Promise<BackendResult<ProjectDiagnosisResult>>((resolve) => {
          settlers.push(resolve);
        });
      },
      cancelDiagnosis: (runToken: string) => {
        cancelTokens.push(runToken);
        return Promise.resolve(cancelResponse ?? ok({ runToken, cancelRequested: true }));
      },
      onDiagnosisProgress: (_listener: (progress: DiagnosisProgress) => void) => () => undefined,
    },
    projectDialog: { selectDirectory: () => Promise.resolve('/p') },
  };

  return {
    api,
    runRequests,
    cancelTokens,
    failCancel: (message: string) => {
      cancelResponse = { ok: false, error: { code: 'BACKEND_UNAVAILABLE', message } };
    },
    finishRun: async (index: number, value: BackendResult<ProjectDiagnosisResult>) => {
      await act(async () => {
        settlers[index](value);
        await Promise.resolve();
      });
    },
  };
}

let backend: ReturnType<typeof makeBackend>;

beforeEach(() => {
  vi.useFakeTimers({ shouldAdvanceTime: true });
  backend = makeBackend();
  (window as unknown as { desktop: unknown }).desktop = backend.api;
});

afterEach(() => {
  cleanup();
  vi.useRealTimers();
  delete (window as unknown as { desktop?: unknown }).desktop;
});

async function settle() {
  await act(async () => {
    vi.advanceTimersByTime(400);
    await Promise.resolve();
  });
}

const runButton = () =>
  screen.getByRole('button', { name: /run molecular diagnosis|running/i }) as HTMLButtonElement;

async function openWorkspace() {
  const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
  render(
    <ProjectProvider>
      <App />
    </ProjectProvider>,
  );
  await user.click(screen.getByRole('button', { name: /open existing/i }));
  await waitFor(() => expect(screen.getByText('Stop project')).toBeDefined());
  await user.click(screen.getByRole('tab', { name: 'Molecular Diagnosis' }));
  await settle();
  await waitFor(() => expect(runButton().disabled).toBe(false));
  return user;
}

async function startRun(user: ReturnType<typeof userEvent.setup>) {
  await user.click(runButton());
  await waitFor(() => expect(screen.getByText('Running')).toBeDefined());
}

describe('Stop analysis', () => {
  it('is offered only while a run is going', async () => {
    const user = await openWorkspace();
    expect(screen.queryByRole('button', { name: 'Stop analysis' })).toBeNull();
    await startRun(user);
    expect(screen.getByRole('button', { name: 'Stop analysis' })).toBeDefined();
  });

  it('sends one stop for this run, then waits in "Stopping…" for the backend', async () => {
    const user = await openWorkspace();
    await startRun(user);
    const token = backend.runRequests[0].runToken;

    await user.click(screen.getByRole('button', { name: 'Stop analysis' }));
    await waitFor(() => expect(screen.getByRole('button', { name: 'Stopping…' })).toBeDefined());
    const stopping = screen.getByRole('button', { name: 'Stopping…' }) as HTMLButtonElement;
    expect(stopping.disabled).toBe(true);
    await user.click(stopping);

    expect(backend.cancelTokens).toEqual([token]);
    // Still a run in flight: Run stays unavailable until the backend answers.
    expect(runButton().disabled).toBe(true);
    expect(screen.getByText('Stopping')).toBeDefined();
  });

  it('ends as "Analysis stopped.", not as an error or a result, and allows another run', async () => {
    const user = await openWorkspace();
    await startRun(user);
    await user.click(screen.getByRole('button', { name: 'Stop analysis' }));
    await backend.finishRun(0, cancelled);

    await waitFor(() => expect(screen.getByText('Analysis stopped.')).toBeDefined());
    expect(screen.queryByRole('alert')).toBeNull();
    expect(screen.queryByText(/unique diagnostic sites/)).toBeNull();
    expect(screen.queryByRole('button', { name: /stop analysis|stopping/i })).toBeNull();

    // A fresh run, with a fresh token, works normally.
    await waitFor(() => expect(runButton().disabled).toBe(false));
    await startRun(user);
    expect(backend.runRequests).toHaveLength(2);
    expect(backend.runRequests[1].runToken).not.toBe(backend.runRequests[0].runToken);
    await backend.finishRun(1, ok(result({ canContinue: false, resume: null })));
    await waitFor(() => expect(screen.getByText(/unique diagnostic sites/)).toBeDefined());
    expect(screen.queryByText('Analysis stopped.')).toBeNull();
  });

  it('returns to running if the stop could not be delivered', async () => {
    const user = await openWorkspace();
    await startRun(user);
    backend.failCancel('The analysis backend is not running.');

    await user.click(screen.getByRole('button', { name: 'Stop analysis' }));

    await waitFor(() =>
      expect(screen.getByText(/could not be stopped: The analysis backend is not running/)).toBeDefined(),
    );
    expect(screen.getByRole('button', { name: 'Stop analysis' })).toBeDefined();
  });

  it('shows a result if the run finished before the stop reached it', async () => {
    const user = await openWorkspace();
    await startRun(user);
    await user.click(screen.getByRole('button', { name: 'Stop analysis' }));
    await backend.finishRun(0, ok(result({ canContinue: false, resume: null })));

    await waitFor(() => expect(screen.getByText(/unique diagnostic sites/)).toBeDefined());
    expect(screen.queryByText('Analysis stopped.')).toBeNull();
  });
});

describe('Continue search', () => {
  async function stoppedAtMax() {
    const user = await openWorkspace();
    await startRun(user);
    await backend.finishRun(0, ok(result()));
    await waitFor(() =>
      expect(screen.getByRole('button', { name: 'Continue from size 3' })).toBeDefined(),
    );
    return user;
  }

  it('continues the same run, sending its resume state back unchanged', async () => {
    const user = await stoppedAtMax();
    // Raising the maximum is part of the workflow and keeps Continue available.
    await user.click(screen.getByRole('button', { name: 'Increase Maximum candidate DNC size' }));
    await user.click(screen.getByRole('button', { name: 'Increase Maximum candidate DNC size' }));
    await user.click(screen.getByRole('button', { name: 'Continue from size 3' }));
    expect(backend.runRequests).toHaveLength(2);
    expect(backend.runRequests[1].resume).toEqual(RESUME);
    expect((backend.runRequests[1] as { options?: { maxCandidateSize?: number } }).options)
      .toMatchObject({ maxCandidateSize: 4 });
  });

  it('is withdrawn when a search option changes, and offered again when it is restored', async () => {
    const user = await stoppedAtMax();
    await user.click(screen.getByRole('checkbox', { name: 'Ignore gaps' }));

    await waitFor(() =>
      expect(screen.queryByRole('button', { name: 'Continue from size 3' })).toBeNull(),
    );
    expect(screen.getByText(/can no longer be continued/)).toBeDefined();

    await user.click(screen.getByRole('checkbox', { name: 'Ignore gaps' }));
    await waitFor(() =>
      expect(screen.getByRole('button', { name: 'Continue from size 3' })).toBeDefined(),
    );
  });

  it('is withdrawn when the FASTA pool changes', async () => {
    const user = await stoppedAtMax();
    await user.selectOptions(screen.getByLabelText('FASTA POOL'), 'all');
    await waitFor(() =>
      expect(screen.queryByRole('button', { name: 'Continue from size 3' })).toBeNull(),
    );
    expect(backend.runRequests).toHaveLength(1);
  });
});
