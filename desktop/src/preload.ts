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
  },
} as const;

export type DesktopApi = typeof desktopApi;

contextBridge.exposeInMainWorld('desktop', desktopApi);
