import { contextBridge, ipcRenderer } from 'electron';
import type {
  BackendResult,
  FastaLoadResult,
  FocalPresencePayload,
  FocalQueryResult,
  FocalReplacementResult,
  FocalSetPayload,
  FocalValidationResult,
  MolecularDiagnosisRequest,
  MolecularDiagnosisResult,
  OpenProjectResult,
  ProjectDiagnosisRequest,
  ProjectDiagnosisResult,
  ProjectMetadata,
  RelinkResult,
  SearchHeadersResult,
  SourceStatusPayload,
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

  /**
   * Persistent projects.
   *
   * The renderer never sees SQLite. It cannot open a database, cannot run a
   * query, and has no SQL string anywhere in it; every one of these calls is a
   * named operation that the Python service performs against the project it
   * exclusively owns.
   */
  project: {
    /**
     * Initialise a NEW project in a folder. Refused if one is already there.
     *
     * Distinct from `open` on purpose: these are the two launcher buttons, and
     * conflating them turns "I picked the wrong folder" into "my project is
     * empty".
     */
    create: (projectDir: string, title?: string): Promise<BackendResult<OpenProjectResult>> =>
      ipcRenderer.invoke('project:create', { projectDir, title }),

    /**
     * Open an EXISTING project. Never creates one: a folder with no
     * project.sqlite comes back as `PROJECT_NOT_FOUND`.
     */
    open: (projectDir: string): Promise<BackendResult<OpenProjectResult>> =>
      ipcRenderer.invoke('project:open', { projectDir }),

    close: (): Promise<BackendResult<{ closed: boolean }>> =>
      ipcRenderer.invoke('project:close', {}),

    /** Rename the open project. The only mutable project metadata. */
    setTitle: (title: string): Promise<BackendResult<{ metadata: ProjectMetadata }>> =>
      ipcRenderer.invoke('project:set-title', { title }),

    /**
     * Re-read the state of every linked file.
     *
     * `strong: false` (the default) is stat-only and cheap enough for
     * lifecycle points and window focus. `strong: true` hashes files not yet
     * proven this session, and must never be triggered by typing.
     */
    refreshSources: (
      strong = false,
    ): Promise<BackendResult<{ sources: readonly SourceStatusPayload[] }>> =>
      ipcRenderer.invoke('project:refresh-sources', { strong }),

    linkFasta: (
      path: string,
    ): Promise<BackendResult<{ fastaFileId: string; source: SourceStatusPayload }>> =>
      ipcRenderer.invoke('project:link-fasta', { path }),

    unlinkFasta: (fastaFileId: string): Promise<BackendResult<{ removed: string }>> =>
      ipcRenderer.invoke('project:unlink-fasta', { fastaFileId }),

    /** Point a linked file at a new path, e.g. after the user moved it. */
    relinkFasta: (fastaFileId: string, path: string): Promise<BackendResult<RelinkResult>> =>
      ipcRenderer.invoke('project:relink-fasta', { fastaFileId, path }),

    reindexFasta: (
      fastaFileId: string,
    ): Promise<BackendResult<{ source: SourceStatusPayload }>> =>
      ipcRenderer.invoke('project:reindex-fasta', { fastaFileId }),

    /**
     * Plain substring search. Non-mutating, so files it could not read come
     * back in `unavailable` instead of failing the call.
     *
     * `fastaFileIds` omitted searches the whole project; an EMPTY array is an
     * empty scope and matches nothing. The two are deliberately different.
     */
    searchHeaders: (
      query: string,
      fastaFileIds?: readonly string[],
      limit?: number,
    ): Promise<BackendResult<SearchHeadersResult>> =>
      ipcRenderer.invoke('project:search-headers', { query, fastaFileIds, limit }),

    listFocalSets: (): Promise<BackendResult<{ focalSets: readonly FocalSetPayload[] }>> =>
      ipcRenderer.invoke('project:list-focal-sets', {}),

    getFocalSet: (
      focalSetId: string,
    ): Promise<BackendResult<{ focalSet: FocalSetPayload }>> =>
      ipcRenderer.invoke('project:get-focal-set', { focalSetId }),

    createFocalSet: (
      title: string,
    ): Promise<
      BackendResult<{ focalSetId: string; title: string; focalSet: FocalSetPayload }>
    > => ipcRenderer.invoke('project:create-focal-set', { title }),

    /** Refused with `FOCAL_SET_LOCKED` while the set is locked. */
    renameFocalSet: (
      focalSetId: string,
      title: string,
    ): Promise<BackendResult<{ focalSet: FocalSetPayload }>> =>
      ipcRenderer.invoke('project:rename-focal-set', { focalSetId, title }),

    /** Locking is a data rule, not a UI state: the backend enforces it. */
    setFocalSetLocked: (
      focalSetId: string,
      locked: boolean,
    ): Promise<BackendResult<{ focalSet: FocalSetPayload }>> =>
      ipcRenderer.invoke('project:set-focal-set-locked', { focalSetId, locked }),

    deleteFocalSet: (focalSetId: string): Promise<BackendResult<{ deleted: string }>> =>
      ipcRenderer.invoke('project:delete-focal-set', { focalSetId }),

    /**
     * The large textbox: replace the set's explicit membership wholesale.
     *
     * `text` is the raw `a ; b;c` box. These are COMPLETE headers, not a
     * query — a header that matches no FASTA is kept so it can show red.
     * Applied as a diff, so it is safe to call on a debounce; the caller is
     * responsible for ignoring stale responses.
     */
    replaceFocalEntries: (
      focalSetId: string,
      text: string,
    ): Promise<BackendResult<FocalReplacementResult>> =>
      ipcRenderer.invoke('project:replace-focal-entries', { focalSetId, text }),

    /** Same operation, when the caller has already split the box itself. */
    replaceFocalEntryHeaders: (
      focalSetId: string,
      headers: readonly string[],
    ): Promise<BackendResult<FocalReplacementResult>> =>
      ipcRenderer.invoke('project:replace-focal-entries', { focalSetId, headers }),

    /**
     * `+` mode: a transient query expands to exact headers, which become
     * members. Refuses (rather than adding a subset) if any file in the scope
     * cannot be searched. An empty query is `FOCAL_QUERY_EMPTY`.
     */
    addFocalEntries: (
      focalSetId: string,
      query: string,
      fastaFileIds?: readonly string[],
    ): Promise<BackendResult<FocalQueryResult>> =>
      ipcRenderer.invoke('project:add-focal-entries', { focalSetId, query, fastaFileIds }),

    /**
     * `-` mode: drop members whose stored exact header contains the query.
     * Never re-searches the FASTA, and an empty query is refused rather than
     * matching — and deleting — every member.
     */
    removeFocalEntries: (
      focalSetId: string,
      query: string,
    ): Promise<BackendResult<{ query: string; removed: readonly string[] }>> =>
      ipcRenderer.invoke('project:remove-focal-entries', { focalSetId, query }),

    /**
     * Where each member currently is. Verified by byte offset, not by rehashing.
     *
     * `fastaFileIds` bounds which files may confer presence — pass the files a
     * run would actually read so the colours match what Run will allow.
     */
    focalPresence: (
      focalSetId: string,
      selectedFastaFileId?: string | null,
      fastaFileIds?: readonly string[],
    ): Promise<BackendResult<{ entries: readonly FocalPresencePayload[] }>> =>
      ipcRenderer.invoke('project:focal-presence', {
        focalSetId,
        selectedFastaFileId,
        fastaFileIds,
      }),

    runMolecularDiagnosis: (
      request: ProjectDiagnosisRequest,
    ): Promise<BackendResult<ProjectDiagnosisResult>> =>
      ipcRenderer.invoke('project:run-molecular-diagnosis', request),
  },

  /** Pick a project directory to create or open. Null when cancelled. */
  projectDialog: {
    selectDirectory: (): Promise<string | null> =>
      ipcRenderer.invoke('dialog:select-project-directory'),
  },
} as const;

export type DesktopApi = typeof desktopApi;

export type ExportResult =
  | { readonly ok: true; readonly path: string }
  | { readonly ok: false; readonly code: string; readonly message?: string };

contextBridge.exposeInMainWorld('desktop', desktopApi);
