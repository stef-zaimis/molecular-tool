import { createContext, useCallback, useContext, useEffect, useMemo, useReducer, useRef } from 'react';
import type { Dispatch, ReactNode } from 'react';
import { appReducer, createInitialState, focalHistoryFor } from './projectState';
import type { AppAction, AppState, FocalHistory } from './projectState';
import type { FocalSet } from '../../contract';
import type { DiagnosisResumeState, FocalStringValidation } from '../../backendContract';
import { desktop } from '../desktopApi';
import type { DesktopApi } from '../../preload';
import { describeSources, sourceBannerMessage } from './sourceStatus';
import type { SourceView } from './sourceStatus';

interface ProjectContextValue {
  readonly state: AppState;
  readonly dispatch: Dispatch<AppAction>;
  /** Convenience: the focal set the workspace is currently editing. */
  readonly activeFocalSet: FocalSet;
  /** Undo/redo stacks for the active focal set. */
  readonly activeFocalHistory: FocalHistory;
  /**
   * Backend verdict per focal string, aligned by index with the active set,
   * or null when nothing has been validated yet (renders neutral, not red).
   */
  readonly focalValidations: readonly FocalStringValidation[] | null;
  /** Pick a FASTA and load it through the Python parser/validator. */
  readonly chooseFastaFile: () => Promise<void>;
  /** Run, or continue, the Molecular Diagnosis pipeline. */
  readonly runDiagnosis: (resume?: DiagnosisResumeState | null) => Promise<void>;

  /* ---- persistent project ---------------------------------------- */

  /** Every linked FASTA, already turned into display semantics. */
  readonly sourceViews: readonly SourceView[];
  /** One line summarising anything the user must deal with, or null. */
  readonly sourceBanner: string | null;
  /**
   * Open an EXISTING project directory. Resolves to whether it opened.
   *
   * Never creates one: a folder with no `project.sqlite` comes back as
   * `PROJECT_NOT_FOUND` and is surfaced as a notice.
   */
  readonly openProject: (projectDir: string) => Promise<boolean>;
  /** Initialise a NEW project. Refused if the folder already holds one. */
  readonly createProject: (projectDir: string, title?: string) => Promise<boolean>;
  /** Prompt for a directory, then open the existing project in it. */
  readonly chooseProjectDirectory: () => Promise<boolean>;
  /** Prompt for a directory, then create a project in it. */
  readonly chooseNewProjectDirectory: (title?: string) => Promise<boolean>;
  /** Rename the open project. */
  readonly setProjectTitle: (title: string) => Promise<boolean>;
  /**
   * Re-read linked-file state.
   *
   * `strong` hashes files that have not been proven this session. It is for
   * deliberate user actions and for immediately before a run — never for
   * incidental refreshes, and never for typing.
   */
  readonly refreshSources: (strong?: boolean) => Promise<void>;
  readonly linkFasta: (path: string) => Promise<void>;
  readonly reindexSource: (fastaFileId: string) => Promise<void>;
  /** Prompt for the file's new location and repoint the link at it. */
  readonly relinkSource: (fastaFileId: string) => Promise<void>;
}

const ProjectContext = createContext<ProjectContextValue | null>(null);

/** `project.open` and `project.create` resolve to the same shape. */
type DesktopProjectOpen = DesktopApi['project']['open'];

/** Directory the outputs are written to: the folder holding the FASTA. */
export function defaultOutputDirectory(fastaPath: string): string {
  const separator = fastaPath.includes('\\') ? '\\' : '/';
  const cut = fastaPath.lastIndexOf(separator);
  return cut > 0 ? fastaPath.slice(0, cut) : fastaPath;
}

export function ProjectProvider({ children }: { children: ReactNode }): JSX.Element {
  const [state, dispatch] = useReducer(appReducer, undefined, createInitialState);

  const chooseFastaFile = useCallback(async () => {
    const path = await desktop().dialog.selectFastaFile();
    // A cancelled dialog returns null; leave any previous choice intact.
    if (!path) return;

    dispatch({ type: 'setFastaPath', path });
    dispatch({ type: 'alignmentLoading', path });

    // Parsing and alignment validation both happen in Python.
    const response = await desktop().analysis.loadFasta(path);
    if (response.ok) {
      dispatch({ type: 'alignmentLoaded', data: response.result });
    } else {
      dispatch({ type: 'alignmentFailed', path, error: response.error });
    }
  }, []);

  const activeFocalSet =
    state.focalSets.find((set) => set.id === state.activeFocalSetId) ?? state.focalSets[0];

  const fastaPath = state.alignment.status === 'loaded' ? state.alignment.data.path : null;
  const focalStrings = activeFocalSet.strings;

  /*
   * Green/red comes from the backend matcher, not from a second implementation
   * in TypeScript, so the colours cannot disagree with which sequences the
   * analysis will actually select. Re-validated whenever the set or the file
   * changes; a stale in-flight response is discarded by sequence number.
   */
  const validationSequence = useRef(0);

  useEffect(() => {
    if (!fastaPath || focalStrings.length === 0) return;

    const ticket = ++validationSequence.current;
    dispatch({ type: 'focalValidating' });

    void desktop()
      .analysis.validateFocalStrings(fastaPath, focalStrings)
      .then((response) => {
        if (ticket !== validationSequence.current) return;
        if (response.ok) {
          dispatch({
            type: 'focalValidated',
            results: response.result.results,
            unionMatchCount: response.result.unionMatchCount,
          });
        } else {
          dispatch({ type: 'focalValidationFailed', error: response.error });
        }
      });
  }, [fastaPath, focalStrings]);

  /* ------------------------------------------------------------------ *
   * Persistent project
   * ------------------------------------------------------------------ */

  /*
   * Opening and creating are separate calls to separate backend methods.
   *
   * They used to be one: `project.open` created a database if the folder did
   * not have one, which meant "Open Existing Project" on the wrong folder
   * produced an empty project that looked exactly like lost work.
   */
  const applyOpenResult = useCallback(
    (projectDir: string, response: Awaited<ReturnType<DesktopProjectOpen>>) => {
      if (response.ok) {
        // The response already includes a status sweep, so a file that went
        // missing while the app was closed is visible on the first paint
        // rather than appearing a moment later.
        dispatch({ type: 'projectOpened', result: response.result });
        return true;
      }
      dispatch({ type: 'projectOpenFailed', projectDir, error: response.error });
      dispatch({ type: 'showNotice', message: response.error.message });
      return false;
    },
    [],
  );

  const openProject = useCallback(
    async (projectDir: string) => {
      dispatch({ type: 'projectOpening', projectDir });
      return applyOpenResult(projectDir, await desktop().project.open(projectDir));
    },
    [applyOpenResult],
  );

  const createProject = useCallback(
    async (projectDir: string, title?: string) => {
      dispatch({ type: 'projectOpening', projectDir });
      return applyOpenResult(projectDir, await desktop().project.create(projectDir, title));
    },
    [applyOpenResult],
  );

  const chooseProjectDirectory = useCallback(async () => {
    const directory = await desktop().projectDialog.selectDirectory();
    // A cancelled dialog is not a failure; leave any open project alone.
    if (!directory) return false;
    return openProject(directory);
  }, [openProject]);

  const chooseNewProjectDirectory = useCallback(
    async (title?: string) => {
      const directory = await desktop().projectDialog.selectDirectory();
      if (!directory) return false;
      return createProject(directory, title);
    },
    [createProject],
  );

  const setProjectTitle = useCallback(async (title: string) => {
    const response = await desktop().project.setTitle(title);
    if (response.ok) return true;
    dispatch({ type: 'showNotice', message: response.error.message });
    return false;
  }, []);

  /*
   * Cheap by default. `refreshSources()` is a stat comparison per file, which
   * is what makes it safe to run on window focus and on every screen change;
   * `refreshSources(true)` escalates to SHA-256 for files not yet proven this
   * session, and is only ever called from an explicit user action or from the
   * run path.
   */
  const refreshSources = useCallback(async (strong = false) => {
    dispatch({ type: 'sourcesRefreshing' });
    const response = await desktop().project.refreshSources(strong);
    if (response.ok) {
      dispatch({ type: 'sourcesRefreshed', sources: response.result.sources });
    } else {
      dispatch({ type: 'sourcesRefreshFailed' });
    }
  }, []);

  const linkFasta = useCallback(async (path: string) => {
    const response = await desktop().project.linkFasta(path);
    if (response.ok) {
      dispatch({ type: 'sourceUpdated', source: response.result.source });
    } else {
      dispatch({ type: 'showNotice', message: response.error.message });
    }
  }, []);

  const reindexSource = useCallback(async (fastaFileId: string) => {
    // The status shows "Re-indexing" for the whole read, so a large file does
    // not sit there still labelled "Changed on disk" with nothing happening.
    dispatch({ type: 'sourceActivity', fastaFileId, activity: 'reindexing' });
    const response = await desktop().project.reindexFasta(fastaFileId);
    if (response.ok) {
      dispatch({ type: 'sourceUpdated', source: response.result.source });
    } else {
      dispatch({ type: 'sourceActivity', fastaFileId, activity: 'idle' });
      dispatch({ type: 'showNotice', message: response.error.message });
    }
  }, []);

  const relinkSource = useCallback(async (fastaFileId: string) => {
    const path = await desktop().dialog.selectFastaFile();
    if (!path) return;

    dispatch({ type: 'sourceActivity', fastaFileId, activity: 'relinking' });
    const response = await desktop().project.relinkFasta(fastaFileId, path);
    if (response.ok) {
      dispatch({ type: 'sourceUpdated', source: response.result.source });
      dispatch({
        type: 'showNotice',
        message: response.result.identical
          ? 'Relinked to an identical copy; nothing needed rebuilding.'
          : 'Relinked and re-indexed from the new file.',
      });
    } else {
      dispatch({ type: 'sourceActivity', fastaFileId, activity: 'idle' });
      dispatch({ type: 'showNotice', message: response.error.message });
    }
  }, []);

  /*
   * Lifecycle refreshes.
   *
   * A linked file can change while the app is in the background, so the state
   * on screen is only trustworthy if it is re-read at the moments the user
   * comes back to it: regaining window focus, and entering a screen that shows
   * or acts on sources. Both take the cheap path.
   *
   * Note what is NOT here: nothing in the focal editor triggers this. Focal
   * text changes are matched against already-loaded data, so typing never
   * causes a stat sweep, let alone a rehash.
   */
  const projectOpen = state.project.status === 'open';

  useEffect(() => {
    if (!projectOpen) return;
    const onFocus = () => void refreshSources(false);
    window.addEventListener('focus', onFocus);
    return () => window.removeEventListener('focus', onFocus);
  }, [projectOpen, refreshSources]);

  useEffect(() => {
    if (!projectOpen) return;
    void refreshSources(false);
  }, [projectOpen, state.screen, state.activeAnalysis, refreshSources]);

  const runDiagnosis = useCallback(
    async (resume: DiagnosisResumeState | null = null) => {
      const path = state.draft.fastaPath;
      if (!path) return;

      const config = state.draft.molecularDiagnosis;
      dispatch({ type: 'diagnosisStarted', continuing: resume !== null });

      const response = await desktop().analysis.runMolecularDiagnosis({
        fastaPath: path,
        focalStrings: activeFocalSet.strings,
        outputDirectory: defaultOutputDirectory(path),
        ignoreGaps: config.ignoreGaps,
        giveBenefitOfDoubtToAmbiguousBases: config.giveBenefitOfDoubtToAmbiguousBases,
        minCandidateSize: config.minCandidateSize,
        maxCandidateSize: config.maxCandidateSize,
        resume,
      });

      if (response.ok) {
        dispatch({ type: 'diagnosisSucceeded', result: response.result });
      } else {
        dispatch({ type: 'diagnosisFailed', error: response.error });
      }
    },
    [state.draft.fastaPath, state.draft.molecularDiagnosis, activeFocalSet.strings],
  );

  const sourceViews = useMemo(
    () => describeSources(state.sources.items, state.sources.activity),
    [state.sources.items, state.sources.activity],
  );

  const value = useMemo<ProjectContextValue>(() => {
    const validations =
      state.focalValidation.status === 'loaded' ? state.focalValidation.results : null;

    return {
      state,
      dispatch,
      activeFocalSet,
      activeFocalHistory: focalHistoryFor(state, activeFocalSet.id),
      // Only trust the verdicts if they describe the set currently on screen.
      focalValidations:
        validations && validations.length === activeFocalSet.strings.length ? validations : null,
      chooseFastaFile,
      runDiagnosis,
      sourceViews,
      sourceBanner: sourceBannerMessage(sourceViews),
      openProject,
      createProject,
      chooseProjectDirectory,
      chooseNewProjectDirectory,
      setProjectTitle,
      refreshSources,
      linkFasta,
      reindexSource,
      relinkSource,
    };
  }, [
    state,
    activeFocalSet,
    chooseFastaFile,
    runDiagnosis,
    sourceViews,
    openProject,
    createProject,
    chooseProjectDirectory,
    chooseNewProjectDirectory,
    setProjectTitle,
    refreshSources,
    linkFasta,
    reindexSource,
    relinkSource,
  ]);

  return <ProjectContext.Provider value={value}>{children}</ProjectContext.Provider>;
}

export function useProject(): ProjectContextValue {
  const context = useContext(ProjectContext);
  if (!context) throw new Error('useProject must be used inside <ProjectProvider>.');
  return context;
}
