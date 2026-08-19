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
    exportFocalSet: () => Promise.resolve({ ok: false as const, code: 'NO_DESKTOP_RUNTIME' }),
  },
  analysis: {
    ping: () => Promise.resolve(unavailable),
    loadFasta: () => Promise.resolve(unavailable),
    validateFocalStrings: () => Promise.resolve(unavailable),
    runMolecularDiagnosis: () => Promise.resolve(unavailable),
  },
};

export function desktop(): DesktopApi {
  return typeof window !== 'undefined' && window.desktop ? window.desktop : noopApi;
}

/** True when running inside Electron with the preload bridge available. */
export function isDesktopRuntime(): boolean {
  return typeof window !== 'undefined' && Boolean(window.desktop);
}
