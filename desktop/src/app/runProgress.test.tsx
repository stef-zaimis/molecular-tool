import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act, cleanup, render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { App } from './App';
import { ProjectProvider } from './state/ProjectContext';
import type {
  BackendError,
  BackendResult,
  DiagnosisProgress,
  FocalSetPayload,
  OpenProjectResult,
  ProjectDiagnosisResult,
  SourceStatusPayload,
} from '../backendContract';

/**
 * What the user sees while a run is running, and when one dies.
 *
 * The rules under test are the ones that decide whether the window can lie:
 * progress must belong to the run in flight, a notification from a replaced or
 * finished run must be ignored, and a backend that stops answering must not
 * leave "Running" on screen forever.
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
  metadata: { projectUuid: 'u', title: 'Progress project' },
  capabilities: { sqliteVersion: '3', fts5: true, trigram: true, acceleratedSearch: true },
  sources: [SOURCE],
};

/** Saved and clean, so the run gate is open. */
const SAVED_SET: FocalSetPayload = {
  id: 'set-1',
  title: 'Targets',
  locked: false,
  entries: [{ id: 'e1', header: 'focal_1|AU|Target' }],
};

const RUN_RESULT: ProjectDiagnosisResult = {
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
    combinationsByLength: { '1': [[1]] },
    singleSites: [1],
    pairs: [],
    uniqueSites: [1],
    states: { '1': 'A' },
    stopReason: 'found_at_or_above_minimum_length',
    stoppedAtLength: 1,
    minCombinationLength: 1,
    maxCombinationLength: 2,
    startCombinationLength: 1,
    diagnostics: {
      fixedCount: 1,
      skippedSites: 0,
      globallyConservedRemoved: 0,
      candidateCount: 1,
      pairsTested: 0,
      totalCombinationsTested: 1,
      combinationsTestedByLength: { '1': 1 },
      ambiguousBdSitesIncluded: [],
      gappyConsensusSitesIncluded: [],
    },
  },
  canContinue: false,
  resume: null,
};

const ok = <T,>(result: T): BackendResult<T> => ({ ok: true, result });

function makeBackend() {
  const progressListeners: Array<(progress: DiagnosisProgress) => void> = [];
  /** Resolves the pending run, so a test decides when it finishes. */
  let settleRun: (value: BackendResult<ProjectDiagnosisResult>) => void = () => undefined;
  const runTokens: (string | undefined)[] = [];

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
      runMolecularDiagnosis: (request: { runToken?: string }) => {
        runTokens.push(request.runToken);
        return new Promise<BackendResult<ProjectDiagnosisResult>>((resolve) => {
          settleRun = resolve;
        });
      },
      onDiagnosisProgress: (listener: (progress: DiagnosisProgress) => void) => {
        progressListeners.push(listener);
        return () => {
          const index = progressListeners.indexOf(listener);
          if (index >= 0) progressListeners.splice(index, 1);
        };
      },
    },
    projectDialog: { selectDirectory: () => Promise.resolve('/p') },
  };

  return {
    api,
    runTokens,
    emitProgress: (progress: DiagnosisProgress) => {
      for (const listener of [...progressListeners]) listener(progress);
    },
    finishRun: (value: BackendResult<ProjectDiagnosisResult>) => settleRun(value),
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

function progress(overrides: Partial<DiagnosisProgress> = {}): DiagnosisProgress {
  return {
    runId: 'abc12345',
    runToken: backend.runTokens[backend.runTokens.length - 1] ?? 'unknown',
    stage: 'five_site',
    current: 428_000,
    total: 2_118_760,
    detail: null,
    elapsedMs: 41_000,
    ...overrides,
  };
}

async function settle() {
  await act(async () => {
    vi.advanceTimersByTime(400);
    await Promise.resolve();
  });
}

async function startRun() {
  const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
  render(
    <ProjectProvider>
      <App />
    </ProjectProvider>,
  );
  await user.click(screen.getByRole('button', { name: /open existing/i }));
  await waitFor(() => expect(screen.getByText('Progress project')).toBeDefined());
  await user.click(screen.getByRole('tab', { name: 'Molecular Diagnosis' }));
  await settle();

  const run = screen.getByRole('button', { name: /run molecular diagnosis/i });
  await waitFor(() => expect((run as HTMLButtonElement).disabled).toBe(false));
  await user.click(run);
  await waitFor(() => expect(screen.getByText('Running')).toBeDefined());
  return user;
}

const statusText = () =>
  (document.querySelector('.run-panel__status--running') as HTMLElement | null)?.textContent ?? '';

describe('live progress', () => {
  it('sends a run token and shows the stage the backend reports', async () => {
    await startRun();
    expect(backend.runTokens).toHaveLength(1);
    expect(backend.runTokens[0]).toBeTruthy();

    act(() => backend.emitProgress(progress()));

    await waitFor(() => expect(statusText()).toMatch(/5-site optimisation/));
    // The counts are the backend's numbers, formatted here.
    expect(statusText()).toMatch(/428,000\/2,118,760/);
  });

  it('names the DMC search size, which is the expensive dimension', async () => {
    await startRun();

    act(() =>
      backend.emitProgress(
        progress({ stage: 'dmc_search', detail: '3', current: 1_000, total: 266_916 }),
      ),
    );

    await waitFor(() => expect(statusText()).toMatch(/DMC search, size 3/));
    expect(statusText()).toMatch(/1,000 \/ 266,916/);
  });

  it('shows a stage with nothing to count as just its name', async () => {
    await startRun();

    act(() =>
      backend.emitProgress(
        progress({ stage: 'writing_workbook', current: null, total: null }),
      ),
    );

    await waitFor(() => expect(statusText()).toMatch(/writing comparison workbook/));
  });

  it('keeps a clock that advances without any progress at all', async () => {
    await startRun();
    expect(statusText()).toMatch(/00:0\d/);

    await act(async () => {
      vi.advanceTimersByTime(65_000);
      await Promise.resolve();
    });

    await waitFor(() => expect(statusText()).toMatch(/01:0\d/));
  });
});

describe('progress belongs to exactly one run', () => {
  it('ignores a notification carrying a different run token', async () => {
    await startRun();
    act(() => backend.emitProgress(progress({ stage: 'consensus' })));
    await waitFor(() => expect(statusText()).toMatch(/focal consensus/));

    // A late notification from a run that has been replaced.
    act(() =>
      backend.emitProgress(progress({ runToken: 'a-different-run', stage: 'writing_report' })),
    );

    await settle();
    expect(statusText()).toMatch(/focal consensus/);
    expect(statusText()).not.toMatch(/writing report/);
  });

  it('ignores progress that arrives after the run has finished', async () => {
    await startRun();

    await act(async () => {
      backend.finishRun(ok(RUN_RESULT));
      await Promise.resolve();
    });
    await waitFor(() => expect(document.querySelector('.run-panel__status--running')).toBeNull());

    act(() => backend.emitProgress(progress({ stage: 'five_site' })));

    await settle();
    // Still finished: no running line came back.
    expect(document.querySelector('.run-panel__status--running')).toBeNull();
    expect(screen.getByText(/unique diagnostic sites/)).toBeDefined();
  });
});

describe('a run that ends badly', () => {
  it('clears Running and reports the error when the backend dies', async () => {
    await startRun();
    act(() => backend.emitProgress(progress({ stage: 'dmc_search', detail: '2' })));
    await waitFor(() => expect(statusText()).toMatch(/DMC search/));

    /*
     * What the bridge produces when the child process exits: the pending
     * request is resolved with a failure rather than left hanging, which is
     * what stops the UI sitting on "Running" forever.
     */
    const error: BackendError = {
      code: 'BACKEND_UNAVAILABLE',
      message: 'The analysis backend stopped unexpectedly.',
      detail: 'backend exited (code 1, signal none) (during project.runMolecularDiagnosis)',
    };
    await act(async () => {
      backend.finishRun({ ok: false, error });
      await Promise.resolve();
    });

    await waitFor(() =>
      expect(screen.getByRole('alert').textContent).toMatch(/stopped unexpectedly/),
    );
    expect(document.querySelector('.run-panel__status--running')).toBeNull();
    // And the button is offered again, rather than stuck disabled.
    await waitFor(() =>
      expect(
        (screen.getByRole('button', { name: /run molecular diagnosis/i }) as HTMLButtonElement)
          .disabled,
      ).toBe(false),
    );
  });
});
