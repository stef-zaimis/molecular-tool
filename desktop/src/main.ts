import { app, BrowserWindow, ipcMain, dialog, screen } from 'electron';
import path from 'node:path';

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
}

app.whenReady().then(() => {
  registerIpc();
  createWindow();

  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) createWindow();
  });
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') app.quit();
});

// Block any attempt by renderer content to open new windows or navigate away.
app.on('web-contents-created', (_event, contents) => {
  contents.setWindowOpenHandler(() => ({ action: 'deny' }));
  contents.on('will-navigate', (event) => event.preventDefault());
});
