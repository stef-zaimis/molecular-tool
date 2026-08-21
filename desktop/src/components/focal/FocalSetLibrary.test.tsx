import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act, cleanup, render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { App } from '../../app/App';
import { ProjectProvider } from '../../app/state/ProjectContext';
import type {
  BackendResult,
  FocalSetPayload,
  OpenProjectResult,
  SourceStatusPayload,
} from '../../backendContract';

/**
 * The FOCAL SET LIBRARY, driven through the real workspace.
 *
 * The rules it has to keep are all about NOT writing behind the user's back:
 * selecting a row must not save the one being left, locking must refuse an
 * unsaved set rather than saving it first, and deleting must be confirmed and
 * must only reach the backend for something that is actually stored.
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
  sequenceCount: 2,
  alignmentLength: 8,
  duplicateHeaderCount: 0,
  locked: false,
  message: null,
};

const OPEN_RESULT: OpenProjectResult = {
  projectDir: '/p',
  outputsDir: '/p/outputs',
  metadata: { projectUuid: 'u', title: 'Library project' },
  capabilities: { sqliteVersion: '3', fts5: true, trigram: true, acceleratedSearch: true },
  sources: [SOURCE],
};

const SET_A: FocalSetPayload = {
  id: 'set-a',
  title: 'Leptacis tipulae',
  locked: false,
  entries: [
    { id: 'a1', header: 'h1' },
    { id: 'a2', header: 'h2' },
  ],
};

const SET_B: FocalSetPayload = {
  id: 'set-b',
  title: 'Leptacis lignicola',
  locked: false,
  entries: [{ id: 'b1', header: 'h3' }],
};

const ok = <T,>(result: T): BackendResult<T> => ({ ok: true, result });

function makeBackend(sets: FocalSetPayload[]) {
  const calls: string[] = [];
  const lockCalls: { id: string; locked: boolean }[] = [];
  const deleteCalls: string[] = [];
  let focalSets = sets;

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
      exportFocalSet: () => Promise.resolve({ ok: true as const, path: '/tmp/x' }),
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
      linkFasta: vi.fn(),
      unlinkFasta: vi.fn(),
      relinkFasta: vi.fn(),
      reindexFasta: vi.fn(),
      searchHeaders: vi.fn(),
      resolveFocalAddQuery: (query: string) => Promise.resolve(ok({ query, headers: [] })),
      listFocalSets: () => Promise.resolve(ok({ focalSets })),
      getFocalSet: vi.fn(),
      createFocalSet: vi.fn(),
      renameFocalSet: vi.fn(),
      setFocalSetLocked: (focalSetId: string, locked: boolean) => {
        calls.push('setFocalSetLocked');
        lockCalls.push({ id: focalSetId, locked });
        const current = focalSets.find((set) => set.id === focalSetId) ?? SET_A;
        return Promise.resolve(ok({ focalSet: { ...current, locked } }));
      },
      deleteFocalSet: (focalSetId: string) => {
        calls.push('deleteFocalSet');
        deleteCalls.push(focalSetId);
        focalSets = focalSets.filter((set) => set.id !== focalSetId);
        return Promise.resolve(ok({ deleted: focalSetId }));
      },
      replaceFocalEntries: vi.fn(),
      replaceFocalEntryHeaders: vi.fn(),
      saveFocalSet: (request: {
        focalSetId?: string | null;
        title: string;
        headers: readonly string[];
      }) => {
        calls.push('saveFocalSet');
        const set: FocalSetPayload = {
          id: request.focalSetId ?? 'set-new',
          title: request.title,
          locked: false,
          entries: request.headers.map((header, index) => ({ id: `n${index}`, header })),
        };
        return Promise.resolve(ok({ focalSet: set }));
      },
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
      matchFocalHeaders: (query: string, headers: readonly string[]) =>
        Promise.resolve(ok({ query, matched: headers })),
      runMolecularDiagnosis: vi.fn(),
    },
    projectDialog: { selectDirectory: () => Promise.resolve('/p') },
  };

  return { api, calls, lockCalls, deleteCalls };
}

let backend: ReturnType<typeof makeBackend>;

function install(sets: FocalSetPayload[]) {
  backend = makeBackend(sets);
  (window as unknown as { desktop: unknown }).desktop = backend.api;
}

beforeEach(() => {
  vi.useFakeTimers({ shouldAdvanceTime: true });
  install([SET_A, SET_B]);
});

afterEach(() => {
  cleanup();
  vi.useRealTimers();
  delete (window as unknown as { desktop?: unknown }).desktop;
});

function setup() {
  return userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
}

async function openLibrary(user: ReturnType<typeof userEvent.setup>) {
  render(
    <ProjectProvider>
      <App />
    </ProjectProvider>,
  );
  await user.click(screen.getByRole('button', { name: /open existing/i }));
  await waitFor(() => expect(screen.getByText('Library project')).toBeDefined());
  await user.click(screen.getByRole('tab', { name: 'Molecular Diagnosis' }));
  await user.click(screen.getByRole('button', { name: 'FOCAL SET LIBRARY' }));
  await waitFor(() => expect(screen.getByRole('region', { name: /focal set library/i })));
}

async function settle() {
  await act(async () => {
    vi.advanceTimersByTime(400);
    await Promise.resolve();
  });
}

describe('the library lists working drafts', () => {
  it('shows every draft with its WORKING membership count', async () => {
    const user = setup();
    await openLibrary(user);

    expect(screen.getByRole('button', { name: 'Select Leptacis tipulae' })).toBeDefined();
    expect(screen.getByRole('button', { name: 'Select Leptacis lignicola' })).toBeDefined();
    expect(screen.getByText('[n=2]')).toBeDefined();
    expect(screen.getByText('[n=1]')).toBeDefined();
  });

  it('marks an INACTIVE dirty draft with the unsaved asterisk', async () => {
    const user = setup();
    await openLibrary(user);

    // Edit the active draft: while it is active the Save control says so, and
    // the row does not need to.
    const content = document.querySelector('.cm-content') as HTMLElement;
    await user.click(content);
    await user.type(content, '; extra');
    await settle();
    expect(screen.queryByTitle('Unsaved changes')).toBeNull();

    // Move to the other one: now the edited row carries the asterisk.
    await user.click(screen.getByRole('button', { name: 'Select Leptacis lignicola' }));
    await waitFor(() => expect(screen.getByTitle('Unsaved changes')).toBeDefined());
  });

  it('selecting a row does not save the one being left', async () => {
    const user = setup();
    await openLibrary(user);

    const content = document.querySelector('.cm-content') as HTMLElement;
    await user.click(content);
    await user.type(content, '; extra');
    await settle();

    await user.click(screen.getByRole('button', { name: 'Select Leptacis lignicola' }));
    await settle();

    expect(backend.calls).not.toContain('saveFocalSet');
    // And the edit survives on the draft that was left.
    await user.click(screen.getByRole('button', { name: 'Select Leptacis tipulae' }));
    await waitFor(() =>
      expect(document.querySelector('.cm-content')?.textContent).toContain('extra'),
    );
  });
});

describe('locking', () => {
  it('refuses a dirty draft rather than auto-saving it', async () => {
    const user = setup();
    await openLibrary(user);

    const content = document.querySelector('.cm-content') as HTMLElement;
    await user.click(content);
    await user.type(content, '; extra');
    await settle();

    await user.click(screen.getByRole('button', { name: /^Lock Leptacis tipulae$/ }));

    await waitFor(() => expect(screen.getByText(/save the changes.*before locking/i)));
    expect(backend.calls).not.toContain('setFocalSetLocked');
    expect(backend.calls).not.toContain('saveFocalSet');
  });

  it('locks a saved, clean draft and adopts the canonical response', async () => {
    const user = setup();
    await openLibrary(user);

    await user.click(screen.getByRole('button', { name: /^Lock Leptacis tipulae$/ }));

    await waitFor(() => expect(backend.lockCalls).toEqual([{ id: 'set-a', locked: true }]));
    // Rename and delete disappear; the lock stays visible so it can be undone.
    await waitFor(() =>
      expect(screen.queryByRole('button', { name: /^Rename Leptacis tipulae$/ })).toBeNull(),
    );
    expect(screen.queryByRole('button', { name: /^Delete Leptacis tipulae$/ })).toBeNull();
    expect(screen.getByRole('button', { name: /^Unlock Leptacis tipulae$/ })).toBeDefined();
  });

  it('a locked set stays selectable', async () => {
    install([{ ...SET_A, locked: true }, SET_B]);
    const user = setup();
    await openLibrary(user);

    await user.click(screen.getByRole('button', { name: 'Select Leptacis tipulae' }));
    await waitFor(() => expect(screen.getByText('Locked')).toBeDefined());
  });
});

describe('deleting', () => {
  it('asks first, and only then calls the backend', async () => {
    const user = setup();
    await openLibrary(user);

    await user.click(screen.getByRole('button', { name: /^Delete Leptacis lignicola$/ }));
    expect(backend.calls).not.toContain('deleteFocalSet');
    expect(screen.getByRole('alertdialog', { name: /delete "leptacis lignicola"/i })).toBeDefined();

    await user.click(screen.getByRole('button', { name: 'Delete' }));
    await waitFor(() => expect(backend.deleteCalls).toEqual(['set-b']));
    await waitFor(() =>
      expect(screen.queryByRole('button', { name: 'Select Leptacis lignicola' })).toBeNull(),
    );
  });

  it('cancelling changes nothing', async () => {
    const user = setup();
    await openLibrary(user);

    await user.click(screen.getByRole('button', { name: /^Delete Leptacis lignicola$/ }));
    await user.click(screen.getByRole('button', { name: 'Cancel' }));

    expect(backend.deleteCalls).toEqual([]);
    expect(screen.getByRole('button', { name: 'Select Leptacis lignicola' })).toBeDefined();
  });

  it('warns that unsaved edits go too when the draft is dirty', async () => {
    const user = setup();
    await openLibrary(user);

    const content = document.querySelector('.cm-content') as HTMLElement;
    await user.click(content);
    await user.type(content, '; extra');
    await settle();

    await user.click(screen.getByRole('button', { name: /^Delete Leptacis tipulae$/ }));
    expect(
      screen.getByRole('alertdialog', { name: /unsaved changes will be discarded/i }),
    ).toBeDefined();
  });

  it('removes a never-saved draft locally, without a backend call', async () => {
    install([]);
    const user = setup();
    await openLibrary(user);

    await user.click(screen.getByRole('button', { name: /^Delete this focal set$/ }));
    await user.click(screen.getByRole('button', { name: 'Delete' }));

    await settle();
    expect(backend.deleteCalls).toEqual([]);
    // A project always has somewhere to type: a fresh blank draft replaces it.
    expect(screen.getByText('[n=0]')).toBeDefined();
  });
});

describe('renaming', () => {
  it('selects the draft and opens its title for editing, locally', async () => {
    const user = setup();
    await openLibrary(user);

    await user.click(screen.getByRole('button', { name: /^Rename Leptacis lignicola$/ }));

    const field = (await screen.findByLabelText('Focal set title')) as HTMLInputElement;
    expect(field.value).toBe('Leptacis lignicola');

    await user.clear(field);
    await user.type(field, 'Renamed');
    // Local until Save: no rename RPC, no save.
    expect(backend.calls).not.toContain('renameFocalSet');
    expect(backend.calls).not.toContain('saveFocalSet');
  });
});

/*
 * Row controls, as STRUCTURE.
 *
 * jsdom has no CSS, so the visible states (hidden at rest, grey on row hover,
 * lit under the pointer) are verified in a browser. What is pinned here is what
 * the CSS acts on and what a screen reader walks: the pencil belongs to the
 * name group, the active row is not a special case with controls of its own,
 * and no control is nested inside another button.
 */
describe('focal row controls', () => {
  it('puts the pencil in the name group, right after the select button', async () => {
    const user = setup();
    await openLibrary(user);

    const group = document.querySelector('.focal-row__name-group');
    expect(group).toBeTruthy();
    const children = Array.from(group!.children).map((element) => element.className);
    expect(children[0]).toContain('focal-row__select');
    expect(children[1]).toContain('focal-row__icon--rename');
  });

  it('never nests a button inside the row selection button', async () => {
    const user = setup();
    await openLibrary(user);

    for (const button of Array.from(document.querySelectorAll('.focal-row button'))) {
      expect(button.querySelector('button')).toBeNull();
    }
  });

  it('gives the ACTIVE row the same controls as any other row', async () => {
    const user = setup();
    await openLibrary(user);

    const active = document.querySelector('.focal-row.is-active');
    const inactive = document.querySelector('.focal-row:not(.is-active)');
    expect(active).toBeTruthy();
    expect(inactive).toBeTruthy();

    const controls = (row: Element) =>
      Array.from(row.querySelectorAll('.focal-row__icon')).map((element) =>
        element.className.replace('focal-row__icon ', ''),
      );
    expect(controls(active!)).toEqual(controls(inactive!));
  });

  it('leaves a locked row with its lock and nothing else', async () => {
    install([{ ...SET_A, locked: true }, SET_B]);
    const user = setup();
    await openLibrary(user);

    const locked = document.querySelector('.focal-row.is-locked');
    expect(locked).toBeTruthy();
    expect(locked!.querySelector('.focal-row__icon--rename')).toBeNull();
    expect(locked!.querySelector('.focal-row__icon--danger')).toBeNull();
    expect(locked!.querySelector('.focal-row__icon--lock')).toBeTruthy();
  });
});

/*
 * The right-hand column.
 *
 * Its top-to-bottom ORDER is the whole mechanism behind the viewer's dynamic
 * top: the viewer is the flex item after the pane, so it starts below whatever
 * the library actually renders instead of at a constant that assumed the
 * library's maximum size.
 */
describe('the workspace right-hand column', () => {
  it('renders the viewer as the item after the selectors and library', async () => {
    const user = setup();
    await openLibrary(user);

    const column = document.querySelector('.workspace-right');
    expect(column).toBeTruthy();
    const order = Array.from(column!.children).map((element) => element.className);
    expect(order[0]).toContain('right-pane');
    expect(order[1]).toContain('viewer');

    // The library is inside the pane, above the viewer — not inside the viewer,
    // which would tie its width to the drag.
    expect(column!.querySelector('.right-pane .focal-library')).toBeTruthy();
  });

  it('has exactly ONE left rail, and it is the drag handle', async () => {
    const user = setup();
    await openLibrary(user);

    const rails = document.querySelectorAll('.viewer__rail');
    expect(rails).toHaveLength(1);
    // The old detached edge and the placeholder's own rail are both gone.
    expect(document.querySelector('.viewer__edge')).toBeNull();
    expect(document.querySelector('.alignment__rail')).toBeNull();

    const handle = screen.getByRole('separator', { name: /resize the sequence visualizer/i });
    expect(handle.className).toContain('viewer__rail');
  });

  it('keeps the viewer mounted with no library, and drops it when hidden', async () => {
    const user = setup();
    await openLibrary(user);

    await user.click(screen.getByRole('button', { name: 'FOCAL SET LIBRARY' }));
    expect(document.querySelector('.focal-library')).toBeNull();
    expect(document.querySelector('.viewer')).toBeTruthy();

    await user.click(screen.getByRole('button', { name: 'SEQUENCE VISUALIZER' }));
    expect(document.querySelector('.viewer')).toBeNull();
  });
});
