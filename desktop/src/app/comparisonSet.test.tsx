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
  SaveFocalSetRequest,
  SourceStatusPayload,
} from '../backendContract';

/**
 * The Comparison Set in the real workspace, against a stubbed backend.
 *
 * What these prove is the wiring: the comparison row's `+`/`−` edit the
 * comparison list and nothing else, Save carries both lists, a header in both
 * sets is painted in BOTH editors with a persistent explanation, and Run waits
 * for overlap and missing comparison entries to be resolved.
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
  sequenceCount: 4,
  alignmentLength: 8,
  duplicateHeaderCount: 0,
  locked: false,
  message: null,
};

const OPEN_RESULT: OpenProjectResult = {
  projectDir: '/projects/p',
  outputsDir: '/projects/p/outputs',
  metadata: { projectUuid: 'uuid-1', title: 'Comparison project' },
  capabilities: { sqliteVersion: '3.43.1', fts5: true, trigram: true, acceleratedSearch: true },
  sources: [SOURCE],
};

const ok = <T,>(result: T): BackendResult<T> => ({ ok: true, result });

/** What the fixture FASTA contains. Anything else is missing. */
const ALL_HEADERS = [
  'focal_1|AU|Target',
  'focal_2|GB|Target',
  'other_1|FR|Contrast',
  'other_2|DE|Contrast',
];

function makeBackend(initialSets: FocalSetPayload[] = []) {
  const calls: string[] = [];
  const saves: SaveFocalSetRequest[] = [];
  let focalSets = initialSets;

  const record =
    <A extends unknown[], R>(name: string, impl: (...args: A) => R) =>
    (...args: A): R => {
      calls.push(name);
      return impl(...args);
    };

  const presenceFor = (headers: readonly string[]) =>
    headers.map<HeaderPresencePayload>((header) => {
      const present = ALL_HEADERS.includes(header);
      return {
        header,
        state: present ? 'present_current' : 'missing',
        occurrences: (present ? { f1: 1 } : {}) as Record<string, number>,
      };
    });

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
      create: record('create', () => Promise.resolve(ok(OPEN_RESULT))),
      open: record('open', () => Promise.resolve(ok(OPEN_RESULT))),
      close: () => Promise.resolve(ok({ closed: true })),
      setTitle: vi.fn(),
      refreshSources: record('refreshSources', () => Promise.resolve(ok({ sources: [SOURCE] }))),
      validateFastaCandidate: vi.fn(),
      linkFasta: vi.fn(),
      setFastaFileLocked: vi.fn(),
      unlinkFasta: vi.fn(),
      relinkFasta: vi.fn(),
      reindexFasta: vi.fn(),
      searchHeaders: vi.fn(),
      resolveFocalAddQuery: record('resolveFocalAddQuery', (query: string) =>
        Promise.resolve(
          ok({
            query,
            headers: ALL_HEADERS.filter((header) =>
              header.toLowerCase().includes(query.toLowerCase()),
            ),
          }),
        ),
      ),
      listFocalSets: record('listFocalSets', () => Promise.resolve(ok({ focalSets }))),
      getFocalSet: vi.fn(),
      createFocalSet: vi.fn(),
      renameFocalSet: vi.fn(),
      setFocalSetLocked: vi.fn(),
      deleteFocalSet: vi.fn(),
      replaceFocalEntries: vi.fn(),
      replaceFocalEntryHeaders: vi.fn(),
      saveFocalSet: record('saveFocalSet', (request: SaveFocalSetRequest) => {
        saves.push(request);
        const saved: FocalSetPayload = {
          id: request.focalSetId ?? 'set-1',
          title: request.title,
          locked: false,
          entries: request.headers.map((header, index) => ({ id: `e${index}`, header })),
          comparisonHeaders: request.comparisonHeaders ?? [],
        };
        focalSets = [saved];
        return Promise.resolve(ok({ focalSet: saved }));
      }),
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
      onDiagnosisProgress: () => () => undefined,
    },
    projectDialog: {
      selectDirectory: () => Promise.resolve('/projects/p'),
    },
  };

  return { api, calls, saves };
}

let backend: ReturnType<typeof makeBackend>;

function install(initialSets: FocalSetPayload[] = []) {
  backend = makeBackend(initialSets);
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

async function settlePresence() {
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
  await waitFor(() => expect(screen.getByText('Comparison project')).toBeDefined());
  await user.click(screen.getByRole('tab', { name: 'Molecular Diagnosis' }));
}

/** [focal editor text, comparison editor text], without the empty-field placeholder. */
function editorTexts(): [string, string] {
  const text = (element: Element | undefined) => {
    if (!element) return '';
    const copy = element.cloneNode(true) as Element;
    copy.querySelectorAll('.cm-placeholder').forEach((node) => node.remove());
    return copy.textContent ?? '';
  };
  const [focal, comparison] = Array.from(document.querySelectorAll('.cm-content'));
  return [text(focal), text(comparison)];
}

const runButton = () =>
  screen.getByRole('button', { name: /run molecular diagnosis/i }) as HTMLButtonElement;

async function addFocal(user: ReturnType<typeof userEvent.setup>, query: string) {
  await user.type(screen.getByLabelText(/add every matching header/i), query);
  await user.click(screen.getByRole('button', { name: 'Apply add' }));
}

async function addComparison(user: ReturnType<typeof userEvent.setup>, query: string) {
  await user.type(
    screen.getByLabelText('Comparison string: add matching headers to the comparison set'),
    query,
  );
  await user.click(screen.getByRole('button', { name: 'Apply comparison add' }));
}

describe('the Comparison Set section', () => {
  it('is present, explained, and empty for a new focal set', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);

    expect(screen.getByText('COMPARISON SET')).toBeDefined();
    expect(
      screen.getByText('Leave blank to compare against all non-focal specimens in the FASTA pool.'),
    ).toBeDefined();
    expect(screen.getByRole('group', { name: 'Comparison set entries, separated by semicolons' }))
      .toBeDefined();
    expect(editorTexts()[1]).toBe('');
  });

  it('`+` in the comparison row fills only the comparison editor, and Save sends both lists', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);

    await user.type(screen.getByLabelText('Focal set title'), 'Targets');
    await addFocal(user, 'Target');
    await addComparison(user, 'Contrast');
    await waitFor(() => expect(editorTexts()[1]).toContain('other_1|FR|Contrast'));

    const [focalText, comparisonText] = editorTexts();
    expect(focalText).not.toContain('Contrast');
    expect(comparisonText).toContain('other_2|DE|Contrast');
    expect(comparisonText).not.toContain('Target');

    await settlePresence();
    await user.click(screen.getByRole('button', { name: 'Save focal set' }));
    await waitFor(() => expect(backend.saves).toHaveLength(1));
    expect(backend.saves[0]).toMatchObject({
      headers: ['focal_1|AU|Target', 'focal_2|GB|Target'],
      comparisonHeaders: ['other_1|FR|Contrast', 'other_2|DE|Contrast'],
    });

    await settlePresence();
    await waitFor(() => expect(runButton().disabled).toBe(false));
  });

  it('`−` in the comparison row removes only comparison entries', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);

    await addFocal(user, 'Target');
    await addComparison(user, 'Contrast');
    await waitFor(() => expect(editorTexts()[1]).toContain('other_2|DE|Contrast'));

    await user.click(screen.getByRole('radio', { name: 'Remove from comparison set' }));
    await user.type(
      screen.getByLabelText('Comparison string: remove matching comparison entries'),
      'other_2',
    );
    await user.click(screen.getByRole('button', { name: 'Apply comparison remove' }));

    await waitFor(() => expect(editorTexts()[1]).not.toContain('other_2'));
    expect(editorTexts()[1]).toContain('other_1|FR|Contrast');
    expect(editorTexts()[0]).toContain('focal_2|GB|Target');
  });
});

describe('focal/comparison overlap', () => {
  it('is painted in both editors, explained persistently, and blocks Run until resolved', async () => {
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);

    await user.type(screen.getByLabelText('Focal set title'), 'Targets');
    await addFocal(user, 'Target');
    await addComparison(user, 'focal_1');
    await waitFor(() => expect(editorTexts()[1]).toContain('focal_1|AU|Target'));
    await settlePresence();

    // Kept in BOTH lists: nothing is silently removed.
    expect(editorTexts()[0]).toContain('focal_1|AU|Target');

    // Painted as a conflict in each editor, over the green presence tone.
    await waitFor(() => expect(document.querySelectorAll('.cm-focal--conflict')).toHaveLength(2));
    for (const span of Array.from(document.querySelectorAll('.cm-focal--conflict'))) {
      expect(span.textContent).toBe('focal_1|AU|Target');
    }
    // The non-overlapping focal entry keeps its normal tone.
    expect(document.querySelectorAll('.cm-focal--match').length).toBeGreaterThan(0);

    expect(screen.getByText(/cannot belong to both/i)).toBeDefined();
    expect(runButton().disabled).toBe(true);
    expect(screen.getByText(/focal_1\|AU\|Target is in both the focal and the comparison set/))
      .toBeDefined();

    // Resolve it from the comparison side.
    await user.click(screen.getByRole('radio', { name: 'Remove from comparison set' }));
    await user.type(
      screen.getByLabelText('Comparison string: remove matching comparison entries'),
      'focal_1',
    );
    await user.click(screen.getByRole('button', { name: 'Apply comparison remove' }));

    await waitFor(() => expect(document.querySelectorAll('.cm-focal--conflict')).toHaveLength(0));
    expect(screen.queryByText(/cannot belong to both/i)).toBeNull();
  });
});

describe('saved focal sets with a comparison list', () => {
  it('restores both lists when a saved set is opened', async () => {
    install([
      {
        id: 'set-1',
        title: 'Targets',
        locked: false,
        entries: [{ id: 'e1', header: 'focal_1|AU|Target' }],
        comparisonHeaders: ['other_2|DE|Contrast'],
      },
    ]);
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);
    await settlePresence();

    await waitFor(() => expect(editorTexts()[1]).toBe('other_2|DE|Contrast'));
    expect(editorTexts()[0]).toBe('focal_1|AU|Target');
    expect(screen.getByText('Saved')).toBeDefined();
    await waitFor(() => expect(runButton().disabled).toBe(false));
  });

  it('a set saved before comparison sets runs with the default comparison', async () => {
    install([
      {
        id: 'set-1',
        title: 'Old',
        locked: false,
        entries: [{ id: 'e1', header: 'focal_1|AU|Target' }],
      },
    ]);
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);
    await settlePresence();

    expect(editorTexts()[1]).toBe('');
    await waitFor(() => expect(runButton().disabled).toBe(false));
  });

  it('blocks Run when a saved comparison entry is not in the FASTA pool', async () => {
    install([
      {
        id: 'set-1',
        title: 'Targets',
        locked: false,
        entries: [{ id: 'e1', header: 'focal_1|AU|Target' }],
        comparisonHeaders: ['ghost|XX|Nowhere'],
      },
    ]);
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);
    await settlePresence();

    await waitFor(() =>
      expect(
        screen.getByText(/Comparison entry ghost\|XX\|Nowhere is not in the selected FASTA/),
      ).toBeDefined(),
    );
    expect(runButton().disabled).toBe(true);
    await waitFor(() =>
      expect(document.querySelector('.cm-focal--missing')?.textContent).toBe('ghost|XX|Nowhere'),
    );
  });

  it('a locked set locks its comparison list too', async () => {
    install([
      {
        id: 'set-1',
        title: 'Targets',
        locked: true,
        entries: [{ id: 'e1', header: 'focal_1|AU|Target' }],
        comparisonHeaders: ['other_2|DE|Contrast'],
      },
    ]);
    const user = userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
    await openWorkspace(user);

    expect(
      (screen.getByRole('button', { name: 'Apply comparison add' }) as HTMLButtonElement).disabled,
    ).toBe(true);
    const [, comparison] = Array.from(document.querySelectorAll('.cm-content'));
    expect(comparison?.getAttribute('contenteditable')).toBe('false');
    expect(
      (screen.getByRole('button', { name: 'Undo comparison set edit' }) as HTMLButtonElement)
        .disabled,
    ).toBe(true);
  });
});
