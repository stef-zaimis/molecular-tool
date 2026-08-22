import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act, cleanup, render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { App } from './App';
import { ProjectProvider } from './state/ProjectContext';
import { TYPING_BURST_IDLE_MS } from './state/typingBurst';
import type {
  BackendResult,
  FocalSetPayload,
  HeaderPresencePayload,
  OpenProjectResult,
  SourceStatusPayload,
} from '../backendContract';

/**
 * The project flow, mounted for real against a stubbed backend.
 *
 * These are the checks a typecheck cannot make: that the screens render, that
 * the launcher actually reaches a project, and — most importantly — that
 * editing a focal set writes NOTHING until Save, and that Run stays refused
 * until it has been saved.
 */

const SOURCE: SourceStatusPayload = {
  fastaFileId: 'f1',
  sourcePath: '/data/a.fasta',
  displayName: 'a.fasta',
  state: 'current',
  available: true,
  indexUsable: true,
  exists: true,
  currentSizeBytes: 100,
  currentMtimeNs: 1,
  indexedSizeBytes: 100,
  indexedMtimeNs: 1,
  indexRevision: 1,
  sequenceCount: 375,
  alignmentLength: 712,
  duplicateHeaderCount: 0,
  locked: false,
  message: null,
};

const OPEN_RESULT: OpenProjectResult = {
  projectDir: '/projects/leptacis',
  outputsDir: '/projects/leptacis/outputs',
  metadata: { projectUuid: 'uuid-1', title: 'European Leptacis' },
  capabilities: { sqliteVersion: '3.43.1', fts5: true, trigram: true, acceleratedSearch: true },
  sources: [SOURCE],
};

const ok = <T,>(result: T): BackendResult<T> => ({ ok: true, result });

/** What the fixture project's single FASTA contains. */
const ALL_HEADERS = ['focal_1|AU|Target', 'focal_2|GB|Target', 'other_1|FR|Contrast'];

/** Records every backend call so a test can assert what was NOT called. */
function makeBackend() {
  const calls: string[] = [];
  let focalSets: FocalSetPayload[] = [];

  const record =
    <A extends unknown[], R>(name: string, impl: (...args: A) => R) =>
    (...args: A): R => {
      calls.push(name);
      return impl(...args);
    };

  /** Every header the fixture FASTA would match, uncapped. */
  const matchesFor = (query: string) =>
    ALL_HEADERS.filter((header) => header.toLowerCase().includes(query.toLowerCase()));

  const presenceFor = (headers: readonly string[]) =>
    headers.map<HeaderPresencePayload>((header) => {
      const present = header.startsWith('focal');
      const occurrences: Record<string, number> = present ? { f1: 1 } : {};
      return { header, state: present ? 'present_current' : 'missing', occurrences };
    });

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
      selectFastaFile: () => Promise.resolve('/data/b.fasta'),
      selectFastaFiles: () => Promise.resolve(['/data/b.fasta']),
      exportFocalSet: () => Promise.resolve({ ok: true as const, path: '/tmp/x.txt' }),
    },
    analysis: {
      ping: vi.fn(),
      loadFasta: vi.fn(),
      validateFocalStrings: vi.fn(),
      runMolecularDiagnosis: vi.fn(),
    },
    project: {
      create: record('create', () => Promise.resolve(ok(OPEN_RESULT))),
      open: record('open', () => Promise.resolve(ok(OPEN_RESULT))),
      close: () => Promise.resolve(ok({ closed: true })),
      setTitle: record('setTitle', (title: string) =>
        Promise.resolve(ok({ metadata: { projectUuid: 'uuid-1', title } })),
      ),
      refreshSources: record('refreshSources', () => Promise.resolve(ok({ sources: [SOURCE] }))),
      validateFastaCandidate: record('validateFastaCandidate', (path: string) =>
        Promise.resolve(
          ok({
            candidate: {
              path,
              displayName: path.split('/').pop() ?? path,
              sequenceCount: 4,
              alignmentLength: 8,
              duplicateHeaderCount: 0,
            },
          }),
        ),
      ),
      linkFasta: record('linkFasta', () =>
        Promise.resolve(ok({ fastaFileId: 'f1', source: SOURCE })),
      ),
      setFastaFileLocked: record('setFastaFileLocked', () =>
        Promise.resolve(ok({ source: SOURCE })),
      ),
      unlinkFasta: record('unlinkFasta', () => Promise.resolve(ok({ removed: 'f1' }))),
      relinkFasta: vi.fn(),
      reindexFasta: vi.fn(),
      searchHeaders: record('searchHeaders', (query: string) =>
        Promise.resolve(ok({ query, hits: [], unavailable: [] })),
      ),
      resolveFocalAddQuery: record('resolveFocalAddQuery', (query: string) =>
        Promise.resolve(ok({ query, headers: matchesFor(query) })),
      ),
      listFocalSets: record('listFocalSets', () => Promise.resolve(ok({ focalSets }))),
      getFocalSet: vi.fn(),
      createFocalSet: record('createFocalSet', vi.fn()),
      renameFocalSet: record('renameFocalSet', vi.fn()),
      setFocalSetLocked: vi.fn(),
      deleteFocalSet: vi.fn(),
      replaceFocalEntries: record('replaceFocalEntries', vi.fn()),
      replaceFocalEntryHeaders: record('replaceFocalEntryHeaders', vi.fn()),
      saveFocalSet: record(
        'saveFocalSet',
        (request: { focalSetId?: string | null; title: string; headers: readonly string[] }) => {
          const saved: FocalSetPayload = {
            id: request.focalSetId ?? 'set-1',
            title: request.title,
            locked: false,
            entries: request.headers.map((header, index) => ({ id: `e${index}`, header })),
          };
          focalSets = [saved];
          return Promise.resolve(ok({ focalSet: saved }));
        },
      ),
      headerPresence: record('headerPresence', (headers: readonly string[]) =>
        Promise.resolve(ok({ entries: presenceFor(headers) })),
      ),
      matchFocalHeaders: record('matchFocalHeaders', (query: string, headers: readonly string[]) =>
        Promise.resolve(
          ok({
            query,
            matched: headers.filter((header) =>
              header.toLowerCase().includes(query.toLowerCase()),
            ),
          }),
        ),
      ),
      runMolecularDiagnosis: record('runMolecularDiagnosis', vi.fn()),
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
    projectDialog: {
      selectDirectory: () => Promise.resolve('/projects/leptacis'),
    },
  };

  return { api, calls, emitProgress };
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

function renderApp() {
  return render(
    <ProjectProvider>
      <App />
    </ProjectProvider>,
  );
}

/** Let the debounced presence request fire and settle. */
async function settlePresence() {
  await act(async () => {
    vi.advanceTimersByTime(400);
    await Promise.resolve();
  });
}

async function openProjectPage(user: ReturnType<typeof userEvent.setup>) {
  renderApp();
  await user.click(screen.getByRole('button', { name: /open existing/i }));
  await waitFor(() => expect(screen.getByText('European Leptacis')).toBeDefined());
}

describe('the launcher', () => {
  it('offers Create and Open as separate actions', () => {
    renderApp();
    expect(screen.getByRole('button', { name: /create new/i })).toBeDefined();
    expect(screen.getByRole('button', { name: /open existing/i })).toBeDefined();
  });

  it('Create New goes to the new-project form without touching the backend', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    renderApp();

    await user.click(screen.getByRole('button', { name: /create new/i }));

    expect(screen.getByLabelText('Project title')).toBeDefined();
    // No project exists until Create is pressed.
    expect(backend.calls).not.toContain('create');
  });

  it('Open Existing opens a project and lands on the project page', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openProjectPage(user);

    expect(backend.calls).toContain('open');
    expect(backend.calls).not.toContain('create');
    // Focal sets are loaded straight after opening.
    expect(backend.calls).toContain('listFocalSets');
  });
});

describe('the new-project form', () => {
  it('requires a title before Create is available', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    renderApp();
    await user.click(screen.getByRole('button', { name: /create new/i }));

    const create = screen.getByRole('button', { name: 'Create project' }) as HTMLButtonElement;
    expect(create.disabled).toBe(true);

    await user.type(screen.getByLabelText('Project title'), 'European Leptacis');
    expect(
      (screen.getByRole('button', { name: 'Create project' }) as HTMLButtonElement).disabled,
    ).toBe(false);
  });

  it('creates the project first, then links the chosen FASTAs', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    renderApp();
    await user.click(screen.getByRole('button', { name: /create new/i }));
    await user.type(screen.getByLabelText('Project title'), 'European Leptacis');
    await user.click(screen.getByRole('button', { name: 'Browse' }));

    await waitFor(() => expect(screen.getByText('b.fasta')).toBeDefined());

    await user.click(screen.getByRole('button', { name: 'Create project' }));

    await waitFor(() => expect(backend.calls).toContain('linkFasta'));
    // Order matters: there is nothing to link into until the project exists.
    expect(backend.calls.indexOf('create')).toBeLessThan(backend.calls.indexOf('linkFasta'));
  });
});

describe('the project page', () => {
  it('shows the persistent title and the linked files with their counts', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openProjectPage(user);

    expect(screen.getByText('European Leptacis')).toBeDefined();
    expect(screen.getByText('a.fasta')).toBeDefined();
    expect(screen.getByText('375')).toBeDefined();
    expect(screen.getByText('712')).toBeDefined();
    // PIS is not computed in this build, so it is a dash rather than a number.
    expect(screen.getByText('—')).toBeDefined();
  });

  it('renames the project through the backend, not in local state', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openProjectPage(user);

    await user.click(screen.getByRole('button', { name: /rename this project/i }));
    const field = screen.getByLabelText('Project title');
    await user.clear(field);
    await user.type(field, 'Renamed');
    await user.click(screen.getByRole('button', { name: 'Save project title' }));

    await waitFor(() => expect(backend.calls).toContain('setTitle'));
    await waitFor(() => expect(screen.getByText('Renamed')).toBeDefined());
  });

  it('asks before unlinking a source', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openProjectPage(user);

    await user.click(screen.getByRole('button', { name: /remove a\.fasta/i }));
    expect(screen.getByRole('alertdialog', { name: /confirm removal/i })).toBeDefined();
    expect(backend.calls).not.toContain('unlinkFasta');

    await user.click(screen.getByRole('button', { name: 'Remove' }));
    await waitFor(() => expect(backend.calls).toContain('unlinkFasta'));
  });

  it('reaches the workspace through the analysis tab, with no Continue button', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openProjectPage(user);

    expect(screen.queryByRole('button', { name: 'Continue' })).toBeNull();

    await user.click(screen.getByRole('tab', { name: 'Molecular Diagnosis' }));
    expect(screen.getByLabelText('Focal set title')).toBeDefined();
  });
});

describe('the Molecular Diagnosis workspace', () => {
  async function openWorkspace(user: ReturnType<typeof userEvent.setup>) {
    await openProjectPage(user);
    await user.click(screen.getByRole('tab', { name: 'Molecular Diagnosis' }));
  }

  it('adds exact headers with `+`, never the query itself', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);

    await user.type(screen.getByLabelText(/add every matching header/i), 'Target');
    await user.click(screen.getByRole('button', { name: 'Apply add' }));

    await waitFor(() => expect(screen.getByText('focal_1|AU|Target')).toBeDefined());
    expect(screen.getByText('focal_2|GB|Target')).toBeDefined();
    // The search string is not a member.
    expect(screen.queryByText('Target')).toBeNull();
  });

  it('writes nothing to the database while editing', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);

    await user.type(screen.getByLabelText('Focal set title'), 'Targets');
    await user.type(screen.getByLabelText(/add every matching header/i), 'Target');
    await user.click(screen.getByRole('button', { name: 'Apply add' }));
    await waitFor(() => expect(screen.getByText('focal_1|AU|Target')).toBeDefined());

    // Editing is a working copy: no write of any kind has happened.
    expect(backend.calls).not.toContain('saveFocalSet');
    expect(backend.calls).not.toContain('createFocalSet');
    expect(backend.calls).not.toContain('replaceFocalEntries');
    // A never-saved set shows the screen-4 title form, not the saved heading.
    expect(screen.getByLabelText('Focal set title')).toBeDefined();
    expect(screen.queryByRole('button', { name: 'SAVE' })).toBeNull();
  });

  it('refuses to run an unsaved focal set, then allows it once saved', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);

    await user.type(screen.getByLabelText('Focal set title'), 'Targets');
    await user.type(screen.getByLabelText(/add every matching header/i), 'Target');
    await user.click(screen.getByRole('button', { name: 'Apply add' }));
    await waitFor(() => expect(screen.getByText('focal_1|AU|Target')).toBeDefined());
    await settlePresence();

    const runButton = () =>
      screen.getByRole('button', { name: /run molecular diagnosis/i }) as HTMLButtonElement;
    expect(runButton().disabled).toBe(true);
    expect(screen.getByText(/save this focal set before running/i)).toBeDefined();

    await user.click(screen.getByRole('button', { name: 'Save focal set' }));
    await waitFor(() => expect(backend.calls).toContain('saveFocalSet'));
    await settlePresence();

    // Saved: the title becomes the screen-6 heading with SAVE beside it.
    await waitFor(() => expect(screen.getByText('Saved')).toBeDefined());
    await waitFor(() => expect(runButton().disabled).toBe(false));
    // Nothing ran on the way there: Run never auto-saves.
    expect(backend.calls).not.toContain('runMolecularDiagnosis');
  });

  it('goes dirty again after a further edit, and blocks the run', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);

    await user.type(screen.getByLabelText('Focal set title'), 'Targets');
    await user.type(screen.getByLabelText(/add every matching header/i), 'Target');
    await user.click(screen.getByRole('button', { name: 'Apply add' }));
    await waitFor(() => expect(screen.getByText('focal_1|AU|Target')).toBeDefined());
    await user.click(screen.getByRole('button', { name: 'Save focal set' }));
    await waitFor(() => expect(screen.getByText('Saved')).toBeDefined());

    // `-` removes from the working copy, through the Python matcher.
    await user.click(screen.getByRole('radio', { name: /remove mode/i }));
    await user.type(screen.getByLabelText(/remove every matching focal entry/i), 'focal_2');
    await user.click(screen.getByRole('button', { name: 'Apply remove' }));

    await waitFor(() => expect(backend.calls).toContain('matchFocalHeaders'));
    await waitFor(() => expect(screen.getByText('Unsaved changes')).toBeDefined());
    await settlePresence();

    expect(
      (screen.getByRole('button', { name: /run molecular diagnosis/i }) as HTMLButtonElement)
        .disabled,
    ).toBe(true);
  });

  it('asks the backend about presence rather than matching in the renderer', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);

    await user.type(screen.getByLabelText(/add every matching header/i), 'Target');
    await user.click(screen.getByRole('button', { name: 'Apply add' }));
    await waitFor(() => expect(screen.getByText('focal_1|AU|Target')).toBeDefined());
    await settlePresence();

    await waitFor(() => expect(backend.calls).toContain('headerPresence'));
    // The legacy substring validator is not part of this path any more.
    expect(backend.api.analysis.validateFocalStrings).not.toHaveBeenCalled();
    expect(backend.api.analysis.runMolecularDiagnosis).not.toHaveBeenCalled();
  });

  it('offers the FASTA pool, defaulting to the first linked file', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);

    const pool = screen.getByLabelText('FASTA POOL') as HTMLSelectElement;
    expect(pool.value).toBe('f1');

    await user.selectOptions(pool, 'all');
    expect((screen.getByLabelText('FASTA POOL') as HTMLSelectElement).value).toBe('all');
  });
});

/*
 * The undo granularity a user actually experiences.
 *
 * The unit rules live in state/typingBurst.test.tsx; what this proves is the
 * WIRING — that the workspace really restarts the idle timer on each document
 * change, and that the pause closes the group without rewriting the text.
 */
describe('typing bursts in the focal editor', () => {
  async function openWorkspace(user: ReturnType<typeof userEvent.setup>) {
    await openProjectPage(user);
    await user.click(screen.getByRole('tab', { name: 'Molecular Diagnosis' }));
  }

  const editorText = () =>
    (document.querySelector('.cm-content') as HTMLElement).textContent ?? '';

  /** A pause long enough for the burst to expire, with no document change. */
  async function pause() {
    await act(async () => {
      vi.advanceTimersByTime(TYPING_BURST_IDLE_MS + 100);
      await Promise.resolve();
    });
  }

  it('starts a new undo step after a pause, and one Undo takes back only that', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);

    const content = document.querySelector('.cm-content') as HTMLElement;
    await user.click(content);
    await user.type(content, 'aa');
    await pause();

    // The pause closed the group; it must NOT have rewritten the text.
    expect(editorText()).toContain('aa');

    await user.type(content, 'bb');
    await pause();
    // Both bursts are in the document. (Where the caret lands between them is
    // the editor's business; what matters is that they are two steps.)
    expect(editorText()).toContain('aa');
    expect(editorText()).toContain('bb');

    await user.click(screen.getByRole('button', { name: 'Undo focal set edit' }));

    // Back to the end of the FIRST burst, not to an empty field.
    await waitFor(() => expect(editorText()).not.toContain('bb'));
    expect(editorText()).toContain('aa');
  });

  it('leaves an unfinished entry exactly as typed while the timer fires', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);

    const content = document.querySelector('.cm-content') as HTMLElement;
    await user.click(content);
    // A half-written second entry: the separator is there, the header is not.
    await user.type(content, 'one ; tw');
    await pause();

    // Canonicalising here would drop the trailing fragment or re-space it.
    expect(editorText()).toContain('one ; tw');
  });
});
