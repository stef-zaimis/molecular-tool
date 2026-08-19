import { contextBridge, ipcRenderer } from 'electron';

/**
 * The entire renderer -> main surface.
 *
 * Deliberately narrow: the renderer gets no Node, no `fs`, no `child_process`
 * and no generic `invoke`. Every capability is an explicit named method so the
 * boundary stays reviewable as backend work lands behind it.
 */
const desktopApi = {
  window: {
    minimize: (): void => ipcRenderer.send('window:minimize'),
    toggleMaximize: (): void => ipcRenderer.send('window:toggle-maximize'),
    close: (): void => ipcRenderer.send('window:close'),
    isMaximized: (): Promise<boolean> => ipcRenderer.invoke('window:is-maximized'),
    onMaximizedChanged: (listener: (maximized: boolean) => void): (() => void) => {
      const handler = (_event: unknown, maximized: boolean) => listener(maximized);
      ipcRenderer.on('window:maximized-changed', handler);
      return () => ipcRenderer.removeListener('window:maximized-changed', handler);
    },
  },
  shell: {
    enterWorkspaceLayout: (): void => ipcRenderer.send('shell:enter-workspace-layout'),
    enterLauncherLayout: (): void => ipcRenderer.send('shell:enter-launcher-layout'),
  },
  dialog: {
    /** Returns the chosen absolute path, or null if the user cancelled. */
    selectFastaFile: (): Promise<string | null> => ipcRenderer.invoke('dialog:select-fasta-file'),
    /** Writes the focal set to a user-chosen .txt file, one entry per line. */
    exportFocalSet: (payload: {
      suggestedName: string;
      lines: readonly string[];
    }): Promise<ExportResult> => ipcRenderer.invoke('dialog:export-focal-set', payload),
  },
  fasta: {
    /**
     * Header lines only. Reads no residues and performs no analysis — see the
     * handler in main.ts. Used to validate focal strings against the real file.
     */
    readHeaders: (filePath: string): Promise<ReadHeadersResult> =>
      ipcRenderer.invoke('fasta:read-headers', filePath),
  },
} as const;

export type ReadHeadersResult =
  | { readonly ok: true; readonly path: string; readonly headers: string[]; readonly duplicateCount: number }
  | { readonly ok: false; readonly code: string; readonly message?: string };

export type ExportResult =
  | { readonly ok: true; readonly path: string }
  | { readonly ok: false; readonly code: string; readonly message?: string };

export type DesktopApi = typeof desktopApi;

contextBridge.exposeInMainWorld('desktop', desktopApi);
