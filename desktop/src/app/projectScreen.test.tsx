import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act, cleanup, render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { App } from './App';
import { ProjectProvider } from './state/ProjectContext';
import type {
  BackendResult,
  FocalSetPayload,
  OpenProjectResult,
  SourceStatusPayload,
} from '../backendContract';

/**
 * The ONE project screen, FASTA vetting, and the per-source lock.
 *
 * The behaviours under test are the ones a user would notice going wrong:
 * a malformed file quietly becoming part of the project, Create throwing them
 * onto a different-looking page, or a locked file still being removable.
 */

const GOOD = '/data/good.fasta';
const RAGGED = '/data/ragged.fasta';

function source(id: string, displayName: string, locked = false): SourceStatusPayload {
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
    sequenceCount: 375,
    alignmentLength: 712,
    duplicateHeaderCount: 0,
    locked,
    message: null,
  };
}

const ok = <T,>(result: T): BackendResult<T> => ({ ok: true, result });

interface Options {
  /** Paths the picker returns. */
  readonly picks?: readonly string[];
  /** Sources an already-open project has. */
  readonly sources?: readonly SourceStatusPayload[];
}

function makeBackend({ picks = [GOOD], sources = [] }: Options = {}) {
  const calls: string[] = [];
  const linked: string[] = [];
  const lockCalls: { id: string; locked: boolean }[] = [];
  let current = [...sources];

  const openResult = (): OpenProjectResult => ({
    projectDir: '/p',
    outputsDir: '/p/outputs',
    metadata: { projectUuid: 'u', title: 'European Leptacis' },
    capabilities: { sqliteVersion: '3', fts5: true, trigram: true, acceleratedSearch: true },
    sources: current,
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
      selectFastaFile: () => Promise.resolve(null),
      selectFastaFiles: () => Promise.resolve(picks),
      exportFocalSet: () => Promise.resolve({ ok: true as const, path: '/x' }),
    },
    analysis: {
      ping: vi.fn(),
      loadFasta: vi.fn(),
      validateFocalStrings: vi.fn(),
      runMolecularDiagnosis: vi.fn(),
    },
    project: {
      create: (_dir: string, _title?: string) => {
        calls.push('create');
        return Promise.resolve(ok(openResult()));
      },
      open: () => {
        calls.push('open');
        return Promise.resolve(ok(openResult()));
      },
      close: () => Promise.resolve(ok({ closed: true })),
      setTitle: vi.fn(),
      refreshSources: () => Promise.resolve(ok({ sources: current })),
      /** The real rule: judged by CONTENT, not by extension. */
      validateFastaCandidate: (path: string) => {
        calls.push('validateFastaCandidate');
        if (path === RAGGED) {
          return Promise.resolve({
            ok: false as const,
            error: {
              code: 'FASTA_NOT_ALIGNED',
              message: 'The FASTA file needs to be aligned.',
            },
          });
        }
        return Promise.resolve(
          ok({
            candidate: {
              path,
              displayName: path.split('/').pop() ?? path,
              sequenceCount: 375,
              alignmentLength: 712,
              duplicateHeaderCount: 0,
            },
          }),
        );
      },
      linkFasta: (path: string) => {
        calls.push('linkFasta');
        linked.push(path);
        const added = source(`f${linked.length}`, path.split('/').pop() ?? path);
        current = [...current, added];
        return Promise.resolve(ok({ fastaFileId: added.fastaFileId, source: added }));
      },
      setFastaFileLocked: (fastaFileId: string, locked: boolean) => {
        calls.push('setFastaFileLocked');
        lockCalls.push({ id: fastaFileId, locked });
        const updated = current.map((item) =>
          item.fastaFileId === fastaFileId ? { ...item, locked } : item,
        );
        current = updated;
        return Promise.resolve(
          ok({ source: updated.find((item) => item.fastaFileId === fastaFileId)! }),
        );
      },
      unlinkFasta: (fastaFileId: string) => {
        calls.push('unlinkFasta');
        current = current.filter((item) => item.fastaFileId !== fastaFileId);
        return Promise.resolve(ok({ removed: fastaFileId }));
      },
      relinkFasta: vi.fn(),
      reindexFasta: vi.fn(),
      searchHeaders: vi.fn(),
      resolveFocalAddQuery: (query: string) => Promise.resolve(ok({ query, headers: [] })),
      listFocalSets: () => Promise.resolve(ok({ focalSets: [] as FocalSetPayload[] })),
      getFocalSet: vi.fn(),
      createFocalSet: vi.fn(),
      renameFocalSet: vi.fn(),
      setFocalSetLocked: vi.fn(),
      deleteFocalSet: vi.fn(),
      replaceFocalEntries: vi.fn(),
      replaceFocalEntryHeaders: vi.fn(),
      saveFocalSet: vi.fn(),
      headerPresence: () => Promise.resolve(ok({ entries: [] })),
      matchFocalHeaders: vi.fn(),
      runMolecularDiagnosis: vi.fn(),
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
    projectDialog: { selectDirectory: () => Promise.resolve('/p') },
  };

  return { api, calls, linked, lockCalls, emitProgress };
}

let backend: ReturnType<typeof makeBackend>;

function install(options: Options = {}) {
  backend = makeBackend(options);
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

function setup() {
  return userEvent.setup({ advanceTimers: vi.advanceTimersByTime });
}

async function openNewProjectPage(user: ReturnType<typeof userEvent.setup>) {
  render(
    <ProjectProvider>
      <App />
    </ProjectProvider>,
  );
  await user.click(screen.getByRole('button', { name: /create new/i }));
}

async function settle() {
  await act(async () => {
    vi.advanceTimersByTime(400);
    await Promise.resolve();
  });
}

/** The FASTA table rows, whichever project state produced them. */
function rowNames(): string[] {
  return Array.from(document.querySelectorAll('.source-row__name')).map(
    (element) => element.textContent ?? '',
  );
}

describe('one screen, before and after Create', () => {
  it('shows the title field with Create beside it, and the same FASTA table', async () => {
    const user = setup();
    await openNewProjectPage(user);

    expect(screen.getByLabelText('Project title')).toBeDefined();
    expect(screen.getByRole('button', { name: 'Create project' })).toBeDefined();
    expect(screen.getByRole('button', { name: 'Browse' })).toBeDefined();
    // The SAME table component, empty for now.
    expect(screen.getByRole('region', { name: 'FASTA files' })).toBeDefined();
    // Analyses belong to a project, so they are not offered yet.
    expect(screen.queryByText('SELECT ANALYSES')).toBeNull();
  });

  it('keeps the chosen files on screen through Create, in the same table', async () => {
    const user = setup();
    await openNewProjectPage(user);

    await user.click(screen.getByRole('button', { name: 'Browse' }));
    await waitFor(() => expect(rowNames()).toEqual(['good.fasta']));

    await user.type(screen.getByLabelText('Project title'), 'European Leptacis');
    await user.click(screen.getByRole('button', { name: 'Create project' }));

    // The page transitions in place: same table, now the project's own.
    await waitFor(() => expect(screen.getByText('European Leptacis')).toBeDefined());
    expect(rowNames()).toEqual(['good.fasta']);
    expect(screen.getByText('SELECT ANALYSES')).toBeDefined();
    expect(backend.linked).toEqual([GOOD]);
    // Order matters: there is nothing to link into until the project exists.
    expect(backend.calls.indexOf('create')).toBeLessThan(backend.calls.indexOf('linkFasta'));
  });

  it('Browse links straight into an open project', async () => {
    install({ sources: [source('f0', 'existing.fasta')] });
    const user = setup();
    render(
      <ProjectProvider>
        <App />
      </ProjectProvider>,
    );
    await user.click(screen.getByRole('button', { name: /open existing/i }));
    await waitFor(() => expect(screen.getByText('European Leptacis')).toBeDefined());

    await user.click(screen.getByRole('button', { name: 'Browse' }));

    await waitFor(() => expect(backend.linked).toEqual([GOOD]));
    await waitFor(() => expect(rowNames()).toEqual(['existing.fasta', 'good.fasta']));
  });
});

describe('FASTA candidates are vetted before they are accepted', () => {
  it('refuses an unaligned file inline and adds nothing', async () => {
    install({ picks: [RAGGED] });
    const user = setup();
    await openNewProjectPage(user);

    await user.click(screen.getByRole('button', { name: 'Browse' }));

    await waitFor(() =>
      expect(screen.getByRole('alert').textContent).toMatch(/needs to be aligned/i),
    );
    expect(rowNames()).toEqual([]);
    expect(backend.calls).not.toContain('linkFasta');
  });

  it('keeps the valid files from a mixed selection and reports only the rejects', async () => {
    install({ picks: [GOOD, RAGGED] });
    const user = setup();
    await openNewProjectPage(user);

    await user.click(screen.getByRole('button', { name: 'Browse' }));

    await waitFor(() => expect(rowNames()).toEqual(['good.fasta']));
    const alert = screen.getByRole('alert');
    expect(alert.textContent).toMatch(/ragged\.fasta/);
    expect(alert.textContent).toMatch(/needs to be aligned/i);
    expect(alert.textContent).not.toMatch(/good\.fasta/);
  });

  it('vets before linking even when the project is already open', async () => {
    install({ picks: [RAGGED], sources: [source('f0', 'existing.fasta')] });
    const user = setup();
    render(
      <ProjectProvider>
        <App />
      </ProjectProvider>,
    );
    await user.click(screen.getByRole('button', { name: /open existing/i }));
    await waitFor(() => expect(screen.getByText('European Leptacis')).toBeDefined());

    await user.click(screen.getByRole('button', { name: 'Browse' }));

    await waitFor(() => expect(backend.calls).toContain('validateFastaCandidate'));
    expect(backend.calls).not.toContain('linkFasta');
    expect(rowNames()).toEqual(['existing.fasta']);
  });
});

describe('the FASTA table', () => {
  async function openWithSources(user: ReturnType<typeof userEvent.setup>) {
    install({
      sources: [
        source('f1', 'one.fasta'),
        source('f2', 'two.fasta'),
        source('f3', 'three.fasta'),
      ],
    });
    render(
      <ProjectProvider>
        <App />
      </ProjectProvider>,
    );
    await user.click(screen.getByRole('button', { name: /open existing/i }));
    await waitFor(() => expect(screen.getByText('European Leptacis')).toBeDefined());
  }

  it('labels the units on the FIRST row only', async () => {
    const user = setup();
    await openWithSources(user);

    const units = Array.from(document.querySelectorAll('.source-row__unit')).map(
      (element) => element.textContent,
    );
    // seq / bp / PIS, once.
    expect(units).toEqual(['seq', 'bp', 'PIS']);

    const rows = document.querySelectorAll('.source-row');
    expect(rows[0].querySelectorAll('.source-row__unit')).toHaveLength(3);
    expect(rows[1].querySelectorAll('.source-row__unit')).toHaveLength(0);
  });

  it('moves the labels to the new first row when row 1 is removed', async () => {
    const user = setup();
    await openWithSources(user);

    await user.click(screen.getByRole('button', { name: /remove one\.fasta/i }));
    await user.click(screen.getByRole('button', { name: 'Remove' }));

    await waitFor(() => expect(rowNames()).toEqual(['two.fasta', 'three.fasta']));
    const rows = document.querySelectorAll('.source-row');
    expect(rows[0].querySelectorAll('.source-row__unit')).toHaveLength(3);
    expect(rows[1].querySelectorAll('.source-row__unit')).toHaveLength(0);
  });

  it('gives every row an edge element that is part of the row grid', async () => {
    const user = setup();
    await openWithSources(user);

    /*
     * Structure only. jsdom has no layout engine and does not load the
     * stylesheet, so the edge's rendered height is checked in the browser, not
     * here; what this pins down is that it is a direct child of the row's grid
     * (which is what `align-self: stretch` acts on) rather than a box nested
     * inside some inner padding.
     */
    const rows = Array.from(document.querySelectorAll('.source-row__body'));
    expect(rows).toHaveLength(3);
    for (const row of rows) {
      const edge = row.firstElementChild;
      expect(edge?.className).toBe('source-row__edge');
    }
  });
});

describe('the per-source lock', () => {
  async function openLocked(user: ReturnType<typeof userEvent.setup>, locked: boolean) {
    install({ sources: [source('f1', 'one.fasta', locked)] });
    render(
      <ProjectProvider>
        <App />
      </ProjectProvider>,
    );
    await user.click(screen.getByRole('button', { name: /open existing/i }));
    await waitFor(() => expect(screen.getByText('European Leptacis')).toBeDefined());
  }

  it('locks through the backend and adopts the canonical answer', async () => {
    const user = setup();
    await openLocked(user, false);

    await user.click(screen.getByRole('button', { name: /^Lock one\.fasta$/ }));

    await waitFor(() => expect(backend.lockCalls).toEqual([{ id: 'f1', locked: true }]));
    await waitFor(() =>
      expect(screen.getByRole('button', { name: /^Unlock one\.fasta$/ })).toBeDefined(),
    );
  });

  it('a locked row loses delete and keeps its lock', async () => {
    const user = setup();
    await openLocked(user, true);

    expect(screen.queryByRole('button', { name: /^Remove one\.fasta$/ })).toBeNull();
    expect(screen.getByRole('button', { name: /^Unlock one\.fasta$/ })).toBeDefined();
    // The lock stays visible without hovering: it is the row's state.
    const lock = document.querySelector('.source-row.is-locked .source-row__icon--lock');
    expect(lock).toBeTruthy();
  });

  it('removal is confirmed in-row and only then reaches the backend', async () => {
    const user = setup();
    await openLocked(user, false);

    await user.click(screen.getByRole('button', { name: /^Remove one\.fasta$/ }));
    expect(backend.calls).not.toContain('unlinkFasta');
    expect(screen.getByRole('alertdialog', { name: /confirm removal/i })).toBeDefined();

    await user.click(screen.getByRole('button', { name: 'Cancel' }));
    expect(backend.calls).not.toContain('unlinkFasta');

    await user.click(screen.getByRole('button', { name: /^Remove one\.fasta$/ }));
    await user.click(screen.getByRole('button', { name: 'Remove' }));
    await waitFor(() => expect(backend.calls).toContain('unlinkFasta'));
    await settle();
  });
});

/*
 * A CANDIDATE has no `fastaFileId`, so its lock cannot travel to a database
 * that does not exist yet. It is still the same control with the same meaning,
 * and the row must behave identically on both sides of Create — including
 * arriving locked once the file is linked.
 */
describe('locking a FASTA before the project exists', () => {
  async function chooseOne(user: ReturnType<typeof userEvent.setup>) {
    await openNewProjectPage(user);
    await user.click(screen.getByRole('button', { name: 'Browse' }));
    await waitFor(() => expect(rowNames()).toEqual(['good.fasta']));
  }

  it('locks a pending row without calling the backend', async () => {
    const user = setup();
    await chooseOne(user);

    await user.click(screen.getByRole('button', { name: /^Lock good\.fasta$/ }));

    // The lock is held locally: there is nothing to lock yet.
    expect(backend.lockCalls).toEqual([]);
    expect(backend.calls).not.toContain('setFastaFileLocked');
    // ...and it reads exactly like a linked locked row.
    expect(screen.getByRole('button', { name: /^Unlock good\.fasta$/ })).toBeDefined();
    expect(document.querySelector('.source-row.is-locked')).toBeTruthy();
  });

  it('hides rename and delete while locked, and restores them on unlock', async () => {
    const user = setup();
    await chooseOne(user);

    await user.click(screen.getByRole('button', { name: /^Lock good\.fasta$/ }));
    expect(screen.queryByRole('button', { name: /^Remove good\.fasta$/ })).toBeNull();
    expect(screen.queryByRole('button', { name: /^Rename good\.fasta/ })).toBeNull();

    await user.click(screen.getByRole('button', { name: /^Unlock good\.fasta$/ }));
    expect(screen.getByRole('button', { name: /^Remove good\.fasta$/ })).toBeDefined();
    expect(screen.getByRole('button', { name: /^Rename good\.fasta/ })).toBeDefined();
    expect(document.querySelector('.source-row.is-locked')).toBeNull();
  });

  it('applies the pending lock to the real source once Create links it', async () => {
    const user = setup();
    await chooseOne(user);
    await user.click(screen.getByRole('button', { name: /^Lock good\.fasta$/ }));

    await user.type(screen.getByLabelText('Project title'), 'European Leptacis');
    await user.click(screen.getByRole('button', { name: 'Create project' }));

    // The lock reaches the backend for the id the link produced, and only then.
    await waitFor(() => expect(backend.lockCalls).toEqual([{ id: 'f1', locked: true }]));
    expect(backend.calls.indexOf('linkFasta')).toBeLessThan(
      backend.calls.indexOf('setFastaFileLocked'),
    );
    // The row transitions pending -> linked without changing what it offers.
    await waitFor(() =>
      expect(screen.getByRole('button', { name: /^Unlock good\.fasta$/ })).toBeDefined(),
    );
    expect(screen.queryByRole('button', { name: /^Remove good\.fasta$/ })).toBeNull();
  });

  it('leaves an unlocked candidate unlocked after Create', async () => {
    const user = setup();
    await chooseOne(user);

    await user.type(screen.getByLabelText('Project title'), 'European Leptacis');
    await user.click(screen.getByRole('button', { name: 'Create project' }));

    await waitFor(() => expect(backend.linked).toEqual([GOOD]));
    expect(backend.lockCalls).toEqual([]);
    await waitFor(() =>
      expect(screen.getByRole('button', { name: /^Lock good\.fasta$/ })).toBeDefined(),
    );
  });
});

/*
 * The pencil belongs to the NAME, not to a column of its own.
 */
describe('the rename pencil sits with the file name', () => {
  it('is inside the name group, immediately after the name', async () => {
    const user = setup();
    await openNewProjectPage(user);
    await user.click(screen.getByRole('button', { name: 'Browse' }));
    await waitFor(() => expect(rowNames()).toEqual(['good.fasta']));

    const group = document.querySelector('.source-row__name-group');
    expect(group).toBeTruthy();
    const children = Array.from(group!.children).map((element) => element.className);
    expect(children[0]).toBe('source-row__name');
    expect(children[1]).toContain('source-row__icon--rename');
  });
});
