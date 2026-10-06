import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act, cleanup, render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { App } from './App';
import { ProjectProvider } from './state/ProjectContext';
import type {
  BackendResult,
  FocalSetPayload,
  HeaderPresencePayload,
  OpenProjectResult,
  ProjectDiagnosisResult,
  SourceStatusPayload,
} from '../backendContract';

/**
 * The three scopes, kept apart, and the run gate enforced where it matters.
 *
 * These are regressions for three specific mistakes:
 *
 *  A. expanding `+` through the CAPPED preview, which silently truncated a
 *     query matching more headers than the cap;
 *  B. comparing presence against the RUN scope, which painted a header red
 *     when the project plainly contained it in another file;
 *  C. relying on the Run button's disabled state, which is a hint rather than
 *     a guarantee.
 */

const FILE_A = 'fileA';
const FILE_B = 'fileB';

/** fileA and fileB, with one header each that the other does not have. */
const CONTENTS: Record<string, readonly string[]> = {
  [FILE_A]: ['a_focal|AU|Target', 'a_other|FR|Contrast'],
  [FILE_B]: ['b_only|GB|Target', 'b_other|DE|Contrast'],
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
    sequenceCount: 2,
    alignmentLength: 8,
    duplicateHeaderCount: 0,
    locked: false,
    message: null,
  };
}

const SOURCES = [source(FILE_A, 'a.fasta'), source(FILE_B, 'b.fasta')];

const OPEN_RESULT: OpenProjectResult = {
  projectDir: '/projects/p',
  outputsDir: '/projects/p/outputs',
  metadata: { projectUuid: 'uuid', title: 'Scopes' },
  capabilities: { sqliteVersion: '3', fts5: true, trigram: true, acceleratedSearch: true },
  sources: SOURCES,
};

const ok = <T,>(result: T): BackendResult<T> => ({ ok: true, result });

/** A minimal but complete run result, so the result panel can render it. */
const RUN_RESULT: ProjectDiagnosisResult = {
  focalSetId: 'set-1',
  focalHeaders: ['a_focal|AU|Target'],
  fastaFileIds: [FILE_A],
  sequenceCount: 2,
  alignmentLength: 8,
  outputs: { reportTxt: '/p/outputs/report.txt', workbookXlsx: '/p/outputs/w.xlsx', consensusTxt: null },
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

interface PresenceCall {
  readonly headers: readonly string[];
  readonly selected: string | null | undefined;
  readonly scope: readonly string[] | undefined;
}

interface AddCall {
  readonly query: string;
  readonly scope: readonly string[] | undefined;
}

function makeBackend(
  saved: FocalSetPayload[] = [],
  { unavailable = [] as readonly string[] } = {},
) {
  const presenceCalls: PresenceCall[] = [];
  const addCalls: AddCall[] = [];
  const runCalls: unknown[] = [];
  let focalSets = saved;

  /**
   * The real `headerPresence` semantics, so the renderer is tested against
   * what Python actually answers rather than a convenient stub.
   */
  const presenceFor = (
    headers: readonly string[],
    selected: string | null | undefined,
    scope: readonly string[] | undefined,
  ): HeaderPresencePayload[] => {
    const requested = scope ?? [FILE_A, FILE_B];
    // A file that cannot be consulted contributes nothing and makes an
    // otherwise-absent answer UNPROVEN, exactly as Python does.
    const answerable = requested.filter((file) => !unavailable.includes(file));
    const blind = requested.length !== answerable.length;

    return headers.map((header) => {
      const occurrences: Record<string, number> = {};
      for (const file of answerable) {
        if (CONTENTS[file]?.includes(header)) occurrences[file] = 1;
      }
      const found = Object.keys(occurrences);

      if (!selected) {
        if (found.length > 0) return { header, state: 'present_current', occurrences };
        return { header, state: blind ? 'unknown' : 'missing', occurrences };
      }
      if (unavailable.includes(selected)) return { header, state: 'unknown', occurrences };
      if (occurrences[selected]) return { header, state: 'present_current', occurrences };
      if (found.length > 0) return { header, state: 'present_other', occurrences };
      return { header, state: blind ? 'unknown' : 'missing', occurrences };
    });
  };

  /** Listeners registered through `project.onDiagnosisProgress`. */
  const progressListeners: Array<(progress: unknown) => void> = [];
  /** Push one progress notification at the renderer, as the backend would. */
  const emitProgress = (progress: unknown) => {
    for (const listener of [...progressListeners]) listener(progress);
  };

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
      exportFocalSet: () => Promise.resolve({ ok: true as const, path: '/tmp/x.txt' }),
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
      refreshSources: () => Promise.resolve(ok({ sources: SOURCES })),
      linkFasta: vi.fn(),
      unlinkFasta: vi.fn(),
      relinkFasta: vi.fn(),
      reindexFasta: vi.fn(),
      searchHeaders: vi.fn(),
      resolveFocalAddQuery: (query: string, fastaFileIds?: readonly string[]) => {
        addCalls.push({ query, scope: fastaFileIds });
        const files = fastaFileIds ?? [FILE_A, FILE_B];
        const headers: string[] = [];
        for (const file of files) {
          for (const header of CONTENTS[file] ?? []) {
            if (header.toLowerCase().includes(query.toLowerCase()) && !headers.includes(header)) {
              headers.push(header);
            }
          }
        }
        return Promise.resolve(ok({ query, headers }));
      },
      listFocalSets: () => Promise.resolve(ok({ focalSets })),
      getFocalSet: vi.fn(),
      createFocalSet: vi.fn(),
      renameFocalSet: vi.fn(),
      setFocalSetLocked: vi.fn(),
      deleteFocalSet: vi.fn(),
      replaceFocalEntries: vi.fn(),
      replaceFocalEntryHeaders: vi.fn(),
      saveFocalSet: (request: {
        focalSetId?: string | null;
        title: string;
        headers: readonly string[];
      }) => {
        const set: FocalSetPayload = {
          id: request.focalSetId ?? 'set-1',
          title: request.title,
          locked: false,
          entries: request.headers.map((header, index) => ({ id: `e${index}`, header })),
        };
        focalSets = [set];
        return Promise.resolve(ok({ focalSet: set }));
      },
      headerPresence: (
        headers: readonly string[],
        selected?: string | null,
        fastaFileIds?: readonly string[],
      ) => {
        presenceCalls.push({ headers, selected, scope: fastaFileIds });
        return Promise.resolve(ok({ entries: presenceFor(headers, selected, fastaFileIds) }));
      },
      matchFocalHeaders: (query: string, headers: readonly string[]) =>
        Promise.resolve(
          ok({
            query,
            matched: headers.filter((h) => h.toLowerCase().includes(query.toLowerCase())),
          }),
        ),
      runMolecularDiagnosis: (request: unknown) => {
        runCalls.push(request);
        return Promise.resolve(ok(RUN_RESULT));
      },
      /*
       * The progress channel. Tests that care drive it through
       * `backend.emitProgress(...)`; the rest simply need it to exist, because
       * the provider subscribes to it on mount.
       */
      onDiagnosisProgress: (listener: (progress: unknown) => void) => {
        progressListeners.push(listener);
        return () => {
          const index = progressListeners.indexOf(listener);
          if (index >= 0) progressListeners.splice(index, 1);
        };
      },
    },
    projectDialog: { selectDirectory: () => Promise.resolve('/projects/p') },
  };

  return { api, presenceCalls, addCalls, runCalls, emitProgress };
}

let backend: ReturnType<typeof makeBackend>;

function install(
  saved: FocalSetPayload[] = [],
  options: { unavailable?: readonly string[] } = {},
) {
  backend = makeBackend(saved, options);
  (window as unknown as { desktop: unknown }).desktop = backend.api;
}

beforeEach(() => {
  vi.useFakeTimers({ shouldAdvanceTime: true });
  install();
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

async function openWorkspace(user: ReturnType<typeof userEvent.setup>) {
  render(
    <ProjectProvider>
      <App />
    </ProjectProvider>,
  );
  await user.click(screen.getByRole('button', { name: /open existing/i }));
  await waitFor(() => expect(screen.getByText('Scopes')).toBeDefined());
  await user.click(screen.getByRole('tab', { name: 'Molecular Diagnosis' }));
}

function setup() {
  return userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
}

/** Type a header straight into the editable focal field. */
async function typeIntoEditor(user: ReturnType<typeof userEvent.setup>, text: string) {
  const content = document.querySelector('.cm-content') as HTMLElement;
  await user.click(content);
  await user.type(content, text);
}

/** The presence tone the editor is painting a header with. */
function toneOf(header: string): string | undefined {
  const spans = Array.from(document.querySelectorAll('.cm-focal'));
  const span = spans.find((element) => element.textContent === header);
  return Array.from(span?.classList ?? []).find((name) => name.startsWith('cm-focal--'));
}

async function selectPool(user: ReturnType<typeof userEvent.setup>, value: string) {
  await user.selectOptions(screen.getByLabelText('FASTA POOL'), value);
  await settle();
}

describe('presence compares against the whole project, not the run scope', () => {
  it('paints a header absent from the selected file but present elsewhere ORANGE', async () => {
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);

    // Typed by hand: it is in fileB, and fileA is selected.
    await typeIntoEditor(user, 'b_only|GB|Target');
    await settle();

    await waitFor(() => expect(toneOf('b_only|GB|Target')).toBe('cm-focal--elsewhere'));

    /*
     * Two calls per check: the DISPLAY one against every linked file (which is
     * what makes orange sayable) and the RUN one against the analysis scope.
     */
    const [display, run] = backend.presenceCalls.slice(-2);
    expect(display.selected).toBe(FILE_A);
    expect([...(display.scope ?? [])].sort()).toEqual([FILE_A, FILE_B].sort());
    expect(run.selected).toBe(FILE_A);
    expect(run.scope).toEqual([FILE_A]);
  });

  it('blocks a run against fileA while an entry is only in fileB', async () => {
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);

    await typeIntoEditor(user, 'b_only|GB|Target');
    await settle();

    await waitFor(() => expect(toneOf('b_only|GB|Target')).toBe('cm-focal--elsewhere'));
    const run = screen.getByRole('button', {
      name: /run molecular diagnosis/i,
    }) as HTMLButtonElement;
    expect(run.disabled).toBe(true);
  });

  it('turns green when the selection moves to the file that has it', async () => {
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);

    await typeIntoEditor(user, 'b_only|GB|Target');
    await settle();
    await waitFor(() => expect(toneOf('b_only|GB|Target')).toBe('cm-focal--elsewhere'));

    await selectPool(user, FILE_B);
    await waitFor(() => expect(toneOf('b_only|GB|Target')).toBe('cm-focal--match'));
  });

  it('turns green under All files, where there is no "elsewhere"', async () => {
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);

    await typeIntoEditor(user, 'b_only|GB|Target');
    await settle();
    await waitFor(() => expect(toneOf('b_only|GB|Target')).toBe('cm-focal--elsewhere'));

    await selectPool(user, 'all');
    await waitFor(() => expect(toneOf('b_only|GB|Target')).toBe('cm-focal--match'));

    // All files means no selected file, so `present_other` cannot arise.
    const last = backend.presenceCalls[backend.presenceCalls.length - 1];
    expect(last.selected).toBeNull();
  });

  it('paints a header in NO linked file red', async () => {
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);

    await typeIntoEditor(user, 'nowhere|ZZ');
    await settle();

    await waitFor(() => expect(toneOf('nowhere|ZZ')).toBe('cm-focal--missing'));
  });
});

describe('`+` searches only the selected file', () => {
  it('does not reach into other project files just to obtain green', async () => {
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);

    await user.type(screen.getByLabelText(/add every matching header/i), 'Target');
    await user.click(screen.getByRole('button', { name: 'Apply add' }));

    // fileA's Target header only — b_only is a Target too, and stays out.
    await waitFor(() => expect(toneOf('a_focal|AU|Target')).toBeDefined());
    expect(document.querySelector('.cm-content')?.textContent).toBe('a_focal|AU|Target');

    const add = backend.addCalls[backend.addCalls.length - 1];
    expect(add.scope).toEqual([FILE_A]);
  });

  it('searches every file under All files', async () => {
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, 'all');

    await user.type(screen.getByLabelText(/add every matching header/i), 'Target');
    await user.click(screen.getByRole('button', { name: 'Apply add' }));

    await waitFor(() =>
      expect(document.querySelector('.cm-content')?.textContent).toBe(
        'a_focal|AU|Target; b_only|GB|Target',
      ),
    );
    expect([...(backend.addCalls[backend.addCalls.length - 1].scope ?? [])].sort()).toEqual(
      [FILE_A, FILE_B].sort(),
    );
  });
});

describe('the run gate is enforced inside runDiagnosis', () => {
  it('never calls the backend for a dirty saved draft, even when invoked directly', async () => {
    // A saved, clean set to start from, so the ONLY thing blocking the run is
    // the edit made below.
    install([
      {
        id: 'set-1',
        title: 'Targets',
        locked: false,
        entries: [{ id: 'e0', header: 'a_focal|AU|Target' }],
      },
    ]);

    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);
    await settle();

    const run = () =>
      screen.getByRole('button', { name: /run molecular diagnosis/i }) as HTMLButtonElement;
    await waitFor(() => expect(run().disabled).toBe(false));

    // Make it dirty by typing.
    await typeIntoEditor(user, '; a_other|FR|Contrast');
    await settle();
    await waitFor(() => expect(screen.getByText('Unsaved changes')).toBeDefined());

    /*
     * Invoke Run directly rather than through the disabled button. A disabled
     * attribute is a hint; the guarantee has to live in the handler, so this
     * fires the click listener the way a stale render or a keyboard path
     * could.
     */
    await act(async () => {
      run().removeAttribute('disabled');
      run().click();
      await Promise.resolve();
    });

    expect(backend.runCalls).toHaveLength(0);
    expect(screen.getByText(/unsaved changes\. save it before running/i)).toBeDefined();
  });

  it('refuses when an entry is not in the run scope, naming it', async () => {
    install([
      {
        id: 'set-1',
        title: 'Targets',
        locked: false,
        entries: [{ id: 'e0', header: 'b_only|GB|Target' }],
      },
    ]);

    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);
    await settle();

    const run = screen.getByRole('button', {
      name: /run molecular diagnosis/i,
    }) as HTMLButtonElement;

    await act(async () => {
      run.removeAttribute('disabled');
      run.click();
      await Promise.resolve();
    });

    expect(backend.runCalls).toHaveLength(0);
  });

  it('does run once the set is saved, clean and fully present', async () => {
    install([
      {
        id: 'set-1',
        title: 'Targets',
        locked: false,
        entries: [{ id: 'e0', header: 'a_focal|AU|Target' }],
      },
    ]);

    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);
    await settle();

    const run = screen.getByRole('button', {
      name: /run molecular diagnosis/i,
    }) as HTMLButtonElement;
    await waitFor(() => expect(run.disabled).toBe(false));

    await user.click(run);
    await waitFor(() => expect(backend.runCalls).toHaveLength(1));

    // The RUN scope is the selected file, not the presence comparison scope.
    expect(backend.runCalls[0]).toMatchObject({
      focalSetId: 'set-1',
      fastaFileIds: [FILE_A],
      singleFile: true,
    });
  });
});


/**
 * The bug this whole split exists for.
 *
 * Display colouring and run gating answer different questions, and gating on
 * the display answer let an UNRELATED unavailable FASTA claim that "a selected
 * FASTA is unavailable" while the selected one was perfectly healthy.
 */
describe('display presence and RUN presence are separate answers', () => {
  const SAVED = (header: string): FocalSetPayload => ({
    id: 'set-1',
    title: 'Targets',
    locked: false,
    entries: [{ id: 'e0', header }],
  });

  function runReason(): string | null {
    const button = screen.getByRole('button', {
      name: /run molecular diagnosis/i,
    }) as HTMLButtonElement;
    return button.getAttribute('title');
  }

  it('1. absent from the selected file AND from every file reads red, blamed on the scope', async () => {
    install([SAVED('nowhere|ZZ')]);
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);
    await settle();

    await waitFor(() => expect(toneOf('nowhere|ZZ')).toBe('cm-focal--missing'));
    expect(runReason()).toMatch(/not in the selected fasta/i);
    expect(runReason()).not.toMatch(/unavailable/i);
  });

  it('2. absent from A but present in B reads orange, and a single-file run is blocked', async () => {
    install([SAVED('b_only|GB|Target')]);
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);
    await settle();

    await waitFor(() => expect(toneOf('b_only|GB|Target')).toBe('cm-focal--elsewhere'));

    const run = screen.getByRole('button', {
      name: /run molecular diagnosis/i,
    }) as HTMLButtonElement;
    expect(run.disabled).toBe(true);
    expect(runReason()).toMatch(/not in the selected fasta/i);

    await act(async () => {
      run.removeAttribute('disabled');
      run.click();
      await Promise.resolve();
    });
    expect(backend.runCalls).toHaveLength(0);
  });

  it('3. an unrelated unavailable B never blames the healthy selected A', async () => {
    install([SAVED('a_focal|AU|Target')], { unavailable: [FILE_B] });
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);
    await settle();

    // The header IS in A, which the run reads, so the run is allowed.
    await waitFor(() => expect(toneOf('a_focal|AU|Target')).toBe('cm-focal--match'));
    const run = screen.getByRole('button', {
      name: /run molecular diagnosis/i,
    }) as HTMLButtonElement;
    await waitFor(() => expect(run.disabled).toBe(false));

    await user.click(run);
    await waitFor(() => expect(backend.runCalls).toHaveLength(1));
    expect(backend.runCalls[0]).toMatchObject({ fastaFileIds: [FILE_A], singleFile: true });
  });

  it('3b. a header absent everywhere while an unrelated B is blind stays a scope refusal', async () => {
    install([SAVED('nowhere|ZZ')], { unavailable: [FILE_B] });
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);
    await settle();

    // Display cannot prove absence (B was not readable), so it is grey...
    await waitFor(() => expect(toneOf('nowhere|ZZ')).toBe('cm-focal--unknown'));
    // ...but the RUN scope is A alone, which answered: it is simply not there.
    expect(runReason()).toMatch(/not in the selected fasta/i);
    expect(runReason()).not.toMatch(/unavailable/i);
  });

  it('4. the SELECTED file being unavailable is reported as exactly that', async () => {
    install([SAVED('a_focal|AU|Target')], { unavailable: [FILE_A] });
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, FILE_A);
    await settle();

    await waitFor(() => expect(runReason()).toMatch(/unavailable/i));
    expect(runReason()).toMatch(/selected fasta is unavailable/i);
  });

  it('5. All files: present anywhere is green, absent everywhere is red, never orange', async () => {
    install([
      {
        id: 'set-1',
        title: 'Targets',
        locked: false,
        entries: [
          { id: 'e0', header: 'b_only|GB|Target' },
          { id: 'e1', header: 'nowhere|ZZ' },
        ],
      },
    ]);
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, 'all');
    await settle();

    await waitFor(() => expect(toneOf('b_only|GB|Target')).toBe('cm-focal--match'));
    expect(toneOf('nowhere|ZZ')).toBe('cm-focal--missing');
    expect(document.querySelector('.cm-focal--elsewhere')).toBeNull();
    expect(runReason()).toMatch(/not in any of the selected fasta files/i);

    /*
     * Under All files the display and run scopes coincide, so ONE lookup serves
     * both — asking twice for the same answer would be pure waste on a path
     * that runs on every debounced keystroke.
     */
    const last = backend.presenceCalls[backend.presenceCalls.length - 1];
    expect(last.selected).toBeNull();
    expect([...(last.scope ?? [])].sort()).toEqual([FILE_A, FILE_B].sort());
  });

  it('5b. All files: genuinely uncheckable reads grey, not red', async () => {
    install([SAVED('nowhere|ZZ')], { unavailable: [FILE_B] });
    const user = setup();
    await openWorkspace(user);
    await selectPool(user, 'all');
    await settle();

    await waitFor(() => expect(toneOf('nowhere|ZZ')).toBe('cm-focal--unknown'));
    expect(runReason()).toMatch(/unavailable/i);
  });
});
