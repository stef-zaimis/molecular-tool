import { app, BrowserWindow, ipcMain, dialog, screen, shell } from 'electron';
import path from 'node:path';
import fs from 'node:fs';
import { PythonBridge } from './backend/pythonBridge';

// Injected by @electron-forge/plugin-vite.
declare const MAIN_WINDOW_VITE_DEV_SERVER_URL: string | undefined;
declare const MAIN_WINDOW_VITE_NAME: string;

/**
 * Window sizing.
 *
 * The launcher is deliberately a small "setup wizard" window matching
 * docs/01-launcher.png (540x289 content box). Entering project creation
 * grows the window to a workspace size, clamped to the display work area
 * so we never open larger than the user's screen.
 */
const LAUNCHER_SIZE = { width: 540, height: 289 };
const WORKSPACE_PREFERRED = { width: 1920, height: 1080 };
const WORKSPACE_MIN = { width: 1180, height: 760 };

let mainWindow: BrowserWindow | null = null;

/**
 * The Python analysis backend, started lazily on first use.
 *
 * Repo root is the parent of the Electron app directory in development, so the
 * service is importable as `molecular_diagnosis.service` and the project
 * virtualenv is found next to it. `MOLECULAR_TOOL_ROOT` overrides it.
 */
let backend: PythonBridge | null = null;

function getBackend(): PythonBridge {
  if (!backend) {
    const repoRoot = process.env.MOLECULAR_TOOL_ROOT ?? path.resolve(app.getAppPath(), '..');
    backend = new PythonBridge(repoRoot, (message) => console.log(`[backend] ${message}`));
  }
  return backend;
}

function createWindow(): void {
  mainWindow = new BrowserWindow({
    width: LAUNCHER_SIZE.width,
    height: LAUNCHER_SIZE.height,
    useContentSize: true,
    resizable: false,
    maximizable: false,
    fullscreenable: false,
    frame: false,
    show: false,
    backgroundColor: '#525F72',
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      nodeIntegration: false,
      sandbox: true,
      webSecurity: true,
    },
  });

  mainWindow.once('ready-to-show', () => mainWindow?.show());

  const emitMaximized = () => {
    if (!mainWindow || mainWindow.isDestroyed()) return;
    mainWindow.webContents.send('window:maximized-changed', mainWindow.isMaximized());
  };
  mainWindow.on('maximize', emitMaximized);
  mainWindow.on('unmaximize', emitMaximized);

  if (MAIN_WINDOW_VITE_DEV_SERVER_URL) {
    void mainWindow.loadURL(MAIN_WINDOW_VITE_DEV_SERVER_URL);
  } else {
    void mainWindow.loadFile(path.join(__dirname, `../renderer/${MAIN_WINDOW_VITE_NAME}/index.html`));
  }

  mainWindow.on('closed', () => {
    mainWindow = null;
  });
}

/** Grow the launcher into the main workspace window, centred and clamped to the display. */
function applyWorkspaceLayout(): void {
  if (!mainWindow || mainWindow.isDestroyed()) return;

  const area = screen.getDisplayNearestPoint(mainWindow.getBounds()).workAreaSize;
  const width = Math.max(WORKSPACE_MIN.width, Math.min(WORKSPACE_PREFERRED.width, area.width));
  const height = Math.max(WORKSPACE_MIN.height, Math.min(WORKSPACE_PREFERRED.height, area.height));

  mainWindow.setResizable(true);
  mainWindow.setMaximizable(true);
  mainWindow.setFullScreenable(true);
  mainWindow.setMinimumSize(WORKSPACE_MIN.width, WORKSPACE_MIN.height);
  mainWindow.setContentSize(width, height, false);
  mainWindow.center();
}

/** Shrink back to the compact launcher window. */
function applyLauncherLayout(): void {
  if (!mainWindow || mainWindow.isDestroyed()) return;

  if (mainWindow.isMaximized()) mainWindow.unmaximize();
  mainWindow.setMinimumSize(LAUNCHER_SIZE.width, LAUNCHER_SIZE.height);
  mainWindow.setContentSize(LAUNCHER_SIZE.width, LAUNCHER_SIZE.height, false);
  mainWindow.setResizable(false);
  mainWindow.setMaximizable(false);
  mainWindow.setFullScreenable(false);
  mainWindow.center();
}

function registerIpc(): void {
  ipcMain.on('window:minimize', () => mainWindow?.minimize());

  ipcMain.on('window:toggle-maximize', () => {
    if (!mainWindow) return;
    if (mainWindow.isMaximized()) mainWindow.unmaximize();
    else if (mainWindow.isMaximizable()) mainWindow.maximize();
  });

  ipcMain.on('window:close', () => mainWindow?.close());

  ipcMain.handle('window:is-maximized', () => mainWindow?.isMaximized() ?? false);

  ipcMain.on('shell:enter-workspace-layout', applyWorkspaceLayout);
  ipcMain.on('shell:enter-launcher-layout', applyLauncherLayout);

  /**
   * The only filesystem touchpoint in this phase: a native "choose a FASTA file"
   * dialog. It returns a path string to the renderer and nothing else — no
   * reading, no parsing. FASTA parsing stays in the Python package.
   */
  ipcMain.handle('dialog:select-fasta-file', async () => {
    if (!mainWindow) return null;

    const result = await dialog.showOpenDialog(mainWindow, {
      title: 'Select aligned FASTA file',
      properties: ['openFile'],
      filters: [
        { name: 'FASTA files', extensions: ['fasta', 'fa', 'fna', 'aln'] },
        { name: 'All files', extensions: ['*'] },
      ],
    });

    if (result.canceled || result.filePaths.length === 0) return null;
    return result.filePaths[0];
  });

  /**
   * Everything scientific goes through the Python service.
   *
   * These handlers are deliberately thin: they validate that a payload is the
   * right shape and forward it. No parsing, no matching and no analysis happens
   * in Electron — that all lives in `molecular_diagnosis`.
   */
  const forward = (method: string) => async (_event: unknown, params: unknown) => {
    if (params !== undefined && (typeof params !== 'object' || params === null)) {
      return {
        ok: false as const,
        error: { code: 'INVALID_PARAMETER', message: 'The request payload was not an object.' },
      };
    }
    return getBackend().call(method, (params ?? {}) as Record<string, unknown>);
  };

  ipcMain.handle('backend:load-fasta', forward('loadFasta'));
  ipcMain.handle('backend:validate-focal-strings', forward('validateFocalStrings'));
  ipcMain.handle('backend:run-molecular-diagnosis', forward('runMolecularDiagnosis'));
  ipcMain.handle('backend:ping', forward('ping'));

  /** Reveal a produced output file in the OS file manager. */
  ipcMain.handle('shell:show-item-in-folder', async (_event, filePath: unknown) => {
    if (typeof filePath !== 'string' || !filePath) return false;
    if (!fs.existsSync(filePath)) return false;
    shell.showItemInFolder(filePath);
    return true;
  });

  /** Writes the focal set to a user-chosen .txt file, one entry per line. */
  ipcMain.handle('dialog:export-focal-set', async (_event, payload: unknown) => {
    if (!mainWindow) return { ok: false as const, code: 'NO_WINDOW' };

    const { suggestedName, lines } = (payload ?? {}) as {
      suggestedName?: unknown;
      lines?: unknown;
    };
    if (!Array.isArray(lines) || !lines.every((line) => typeof line === 'string')) {
      return { ok: false as const, code: 'INVALID_PAYLOAD' };
    }

    const result = await dialog.showSaveDialog(mainWindow, {
      title: 'Export focal set',
      defaultPath: typeof suggestedName === 'string' && suggestedName ? suggestedName : 'focal-set.txt',
      filters: [
        { name: 'Text files', extensions: ['txt'] },
        { name: 'All files', extensions: ['*'] },
      ],
    });

    if (result.canceled || !result.filePath) return { ok: false as const, code: 'CANCELLED' };

    try {
      // Trailing newline so the file is well-formed for line-based tools.
      await fs.promises.writeFile(result.filePath, `${lines.join('\n')}\n`, 'utf8');
      return { ok: true as const, path: result.filePath };
    } catch (error) {
      return {
        ok: false as const,
        code: 'WRITE_FAILED',
        message: error instanceof Error ? error.message : 'Could not write the file.',
      };
    }
  });
}

app.whenReady().then(() => {
  registerIpc();
  createWindow();

  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) createWindow();
  });
});

app.on('before-quit', () => {
  backend?.dispose();
  backend = null;
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') app.quit();
});

// Block any attempt by renderer content to open new windows or navigate away.
app.on('web-contents-created', (_event, contents) => {
  contents.setWindowOpenHandler(() => ({ action: 'deny' }));
  contents.on('will-navigate', (event) => event.preventDefault());
});
