import { contextBridge, ipcRenderer } from 'electron';
import type {
  BackendResult,
  FastaLoadResult,
  FocalValidationResult,
  MolecularDiagnosisRequest,
  MolecularDiagnosisResult,
} from './backendContract';

/**
 * The entire renderer -> main surface.
 *
 * Deliberately narrow: the renderer gets no Node, no `fs`, no `child_process`
 * and no generic `invoke`. Every capability is an explicit named method so the
 * boundary stays reviewable.
 *
 * Everything under `analysis` is forwarded to the Python service by the main
 * process. The renderer never spawns a process and never touches the
 * scientific code; it sends plain data and receives plain data.
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
    /** Reveal a produced output file in the OS file manager. */
    showItemInFolder: (filePath: string): Promise<boolean> =>
      ipcRenderer.invoke('shell:show-item-in-folder', filePath),
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
  analysis: {
    /** Liveness probe for the Python backend. */
    ping: (): Promise<BackendResult<{ ok: boolean; python: string }>> =>
      ipcRenderer.invoke('backend:ping', {}),

    /** Parse + validate a FASTA through the Python reader. */
    loadFasta: (path: string): Promise<BackendResult<FastaLoadResult>> =>
      ipcRenderer.invoke('backend:load-fasta', { path }),

    /** Authoritative green/red state for each focal string. */
    validateFocalStrings: (
      path: string,
      focalStrings: readonly string[],
    ): Promise<BackendResult<FocalValidationResult>> =>
      ipcRenderer.invoke('backend:validate-focal-strings', { path, focalStrings }),

    /** Run (or continue) the Molecular Diagnosis pipeline. */
    runMolecularDiagnosis: (
      request: MolecularDiagnosisRequest,
    ): Promise<BackendResult<MolecularDiagnosisResult>> =>
      ipcRenderer.invoke('backend:run-molecular-diagnosis', request),
  },
} as const;

export type DesktopApi = typeof desktopApi;

export type ExportResult =
  | { readonly ok: true; readonly path: string }
  | { readonly ok: false; readonly code: string; readonly message?: string };

contextBridge.exposeInMainWorld('desktop', desktopApi);
