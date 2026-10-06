import type { DesktopApi } from '../preload';

/**
 * Safe accessor for the preload bridge.
 *
 * In the packaged/dev Electron app `window.desktop` is always present. This
 * fallback exists so the same renderer bundle can also be loaded in a plain
 * browser at exact mockup dimensions for visual comparison against docs/,
 * and so component tests do not need to stub globals.
 *
 * It is a no-op, never a simulation: nothing here pretends to resize a window
 * or open a file dialog.
 */
/** Uniform "there is no backend here" answer for the browser fallback. */
const unavailable = {
  ok: false as const,
  error: {
    code: 'BACKEND_UNAVAILABLE',
    message: 'The analysis backend is only available in the desktop application.',
  },
};

const noopApi: DesktopApi = {
  window: {
    minimize: () => undefined,
    toggleMaximize: () => undefined,
    close: () => undefined,
    isMaximized: () => Promise.resolve(false),
    onMaximizedChanged: () => () => undefined,
  },
  shell: {
    enterWorkspaceLayout: () => undefined,
    enterLauncherLayout: () => undefined,
    showItemInFolder: () => Promise.resolve(false),
  },
  dialog: {
    selectFastaFile: () => Promise.resolve(null),
    selectFastaFiles: () => Promise.resolve([]),
    exportFocalSet: () => Promise.resolve({ ok: false as const, code: 'NO_DESKTOP_RUNTIME' }),
  },
  analysis: {
    ping: () => Promise.resolve(unavailable),
    loadFasta: () => Promise.resolve(unavailable),
    validateFocalStrings: () => Promise.resolve(unavailable),
    runMolecularDiagnosis: () => Promise.resolve(unavailable),
  },
  project: {
    create: () => Promise.resolve(unavailable),
    open: () => Promise.resolve(unavailable),
    close: () => Promise.resolve(unavailable),
    setTitle: () => Promise.resolve(unavailable),
    refreshSources: () => Promise.resolve(unavailable),
    validateFastaCandidate: () => Promise.resolve(unavailable),
    linkFasta: () => Promise.resolve(unavailable),
    setFastaFileLocked: () => Promise.resolve(unavailable),
    unlinkFasta: () => Promise.resolve(unavailable),
    relinkFasta: () => Promise.resolve(unavailable),
    reindexFasta: () => Promise.resolve(unavailable),
    searchHeaders: () => Promise.resolve(unavailable),
    listFocalSets: () => Promise.resolve(unavailable),
    getFocalSet: () => Promise.resolve(unavailable),
    createFocalSet: () => Promise.resolve(unavailable),
    renameFocalSet: () => Promise.resolve(unavailable),
    setFocalSetLocked: () => Promise.resolve(unavailable),
    deleteFocalSet: () => Promise.resolve(unavailable),
    replaceFocalEntries: () => Promise.resolve(unavailable),
    replaceFocalEntryHeaders: () => Promise.resolve(unavailable),
    saveFocalSet: () => Promise.resolve(unavailable),
    headerPresence: () => Promise.resolve(unavailable),
    matchFocalHeaders: () => Promise.resolve(unavailable),
    resolveFocalAddQuery: () => Promise.resolve(unavailable),
    addFocalEntries: () => Promise.resolve(unavailable),
    removeFocalEntries: () => Promise.resolve(unavailable),
    focalPresence: () => Promise.resolve(unavailable),
    runMolecularDiagnosis: () => Promise.resolve(unavailable),
    // No backend, so no progress will ever arrive; unsubscribing is a no-op.
    onDiagnosisProgress: () => () => undefined,
  },
  projectDialog: {
    selectDirectory: () => Promise.resolve(null),
  },
};

export function desktop(): DesktopApi {
  return typeof window !== 'undefined' && window.desktop ? window.desktop : noopApi;
}

/** True when running inside Electron with the preload bridge available. */
export function isDesktopRuntime(): boolean {
  return typeof window !== 'undefined' && Boolean(window.desktop);
}
