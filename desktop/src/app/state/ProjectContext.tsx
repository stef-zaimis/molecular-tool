import {
  createContext,
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useReducer,
  useRef,
} from 'react';
import type { Dispatch, ReactNode } from 'react';
import {
  activeDraft as selectActiveDraft,
  appReducer,
  createInitialState,
  scopeFileIds,
  selectedFileId,
} from './projectState';
import type { AppAction, AppState } from './projectState';
import type { FocalDraft } from './focalDrafts';
import { lockBlockedReason, normaliseHeaders, saveRequestFor } from './focalDrafts';
import { evaluateRunGate } from './runGate';
import type { RunGate } from './runGate';
import type { DiagnosisResumeState, FastaCandidate } from '../../backendContract';
import { desktop } from '../desktopApi';
import type { DesktopApi } from '../../preload';
import { describeSources, sourceBannerMessage } from './sourceStatus';
import type { SourceView } from './sourceStatus';

/** How long typing must settle before the backend is asked about presence. */
const PRESENCE_DEBOUNCE_MS = 300;

interface ProjectContextValue {
  readonly state: AppState;
  readonly dispatch: Dispatch<AppAction>;

  /* ---- focal drafts ------------------------------------------------ */

  /** The working copy the workspace is editing. Never a database row. */
  readonly activeDraft: FocalDraft;
  /** Commit the working copy. The only focal write ordinary editing performs. */
  readonly saveActiveDraft: () => Promise<boolean>;
  /** `+`: expand a query to exact headers within the current scope and append them. */
  readonly addByQuery: (query: string) => Promise<boolean>;
  /** `−`: drop working headers containing the query. Never searches the FASTA. */
  readonly removeByQuery: (query: string) => Promise<boolean>;
  /**
   * Lock or unlock a SAVED, CLEAN draft.
   *
   * Refuses a new or dirty one rather than auto-saving: locking marks a set as
   * final, and there has to be a final version to mark.
   */
  readonly setDraftLocked: (key: string, locked: boolean) => Promise<boolean>;
  /**
   * Delete a draft. A saved set is deleted from the project; a never-saved one
   * only ever existed here, so it is dropped locally.
   */
  readonly deleteDraft: (key: string) => Promise<boolean>;

  /* ---- scope and run ----------------------------------------------- */

  /** The linked files the current scope covers, in project order. */
  readonly scopeSourceIds: readonly string[];
  readonly runGate: RunGate;
  /** Run, or continue, Molecular Diagnosis over the project scope. */
  readonly runDiagnosis: (resume?: DiagnosisResumeState | null) => Promise<void>;

  /* ---- persistent project ------------------------------------------ */

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
  /**
   * Prompt for FASTA files, VET each one, and either link them (project open)
   * or hold them as pending candidates (project not created yet).
   *
   * Valid files are always kept; rejected ones are reported inline and never
   * discard the rest of the selection.
   */
  readonly chooseFastaFiles: () => Promise<void>;
  /**
   * Link vetted candidates gathered before Create.
   *
   * Takes them as an argument rather than reading state: opening the new
   * project clears `pending`, so by the time this runs the list it was meant to
   * link is already gone. The caller captures it first.
   */
  readonly linkCandidates: (
    candidates: readonly FastaCandidate[],
    /** Paths locked before Create, re-applied through the backend once linked. */
    lockedPaths?: readonly string[],
  ) => Promise<void>;
  readonly unlinkSource: (fastaFileId: string) => Promise<void>;
  /** Lock or unlock a linked source. A locked source stays analysable. */
  readonly setSourceLocked: (fastaFileId: string, locked: boolean) => Promise<void>;
  readonly reindexSource: (fastaFileId: string) => Promise<void>;
  /** Prompt for the file's new location and repoint the link at it. */
  readonly relinkSource: (fastaFileId: string) => Promise<void>;
}

const ProjectContext = createContext<ProjectContextValue | null>(null);

/** `project.open` and `project.create` resolve to the same shape. */
type DesktopProjectOpen = DesktopApi['project']['open'];

export function ProjectProvider({ children }: { children: ReactNode }): JSX.Element {
  const [state, dispatch] = useReducer(appReducer, undefined, createInitialState);

  const activeDraft = selectActiveDraft(state);
  /*
   * THREE scopes, deliberately distinct.
   *
   *  1. run scope     — the files an analysis reads
   *  2. `+` scope     — the files a query searches
   *  3. presence scope — the files a header is COMPARED against
   *
   * (1) and (2) are the same thing and follow the FASTA pool selector. (3) is
   * always every linked file, because "present, but in a different file" is
   * only sayable if the other files were consulted. Widening (2) to match (3)
   * would make `+` add headers the run will not read; narrowing (3) to match
   * (2) would paint a header red that the project plainly contains.
   */
  const scopeSourceIds = useMemo(() => scopeFileIds(state), [state]);
  const presenceComparisonIds = useMemo(
    () => state.sources.items.map((item) => item.fastaFileId),
    [state.sources.items],
  );
  const selectedFasta = selectedFileId(state);

  /*
   * Read through refs so the FASTA callbacks stay stable. They are passed to
   * the picker, which is rendered on a page that re-renders as the title is
   * typed; rebuilding them each keystroke would churn every consumer.
   */
  const projectIsOpen = useRef(state.project.status === 'open');
  projectIsOpen.current = state.project.status === 'open';

  /* ------------------------------------------------------------------ *
   * Persistent project
   * ------------------------------------------------------------------ */

  const loadFocalSets = useCallback(async () => {
    const response = await desktop().project.listFocalSets();
    if (response.ok) {
      dispatch({ type: 'focalSetsLoaded', focalSets: response.result.focalSets });
    } else {
      dispatch({ type: 'showNotice', message: response.error.message });
    }
  }, []);

  /*
   * Opening and creating are separate calls to separate backend methods.
   *
   * They used to be one: `project.open` created a database if the folder did
   * not have one, which meant "Open Existing Project" on the wrong folder
   * produced an empty project that looked exactly like lost work.
   */
  const applyOpenResult = useCallback(
    async (projectDir: string, response: Awaited<ReturnType<DesktopProjectOpen>>) => {
      if (response.ok) {
        // The response already includes a status sweep, so a file that went
        // missing while the app was closed is visible on the first paint
        // rather than appearing a moment later.
        dispatch({ type: 'projectOpened', result: response.result });
        await loadFocalSets();
        return true;
      }
      dispatch({ type: 'projectOpenFailed', projectDir, error: response.error });
      dispatch({ type: 'showNotice', message: response.error.message });
      return false;
    },
    [loadFocalSets],
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

  const setProjectTitle = useCallback(async (title: string) => {
    const response = await desktop().project.setTitle(title);
    if (response.ok) {
      dispatch({ type: 'projectTitleChanged', title: response.result.metadata.title });
      return true;
    }
    dispatch({ type: 'showNotice', message: response.error.message });
    return false;
  }, []);

  /*
   * Cheap by default. `refreshSources()` is a stat comparison per file, which
   * is what makes it safe to run on window focus and on every screen change;
   * `refreshSources(true)` escalates to SHA-256 for files not yet proven this
   * session, and is only ever called from an explicit user action.
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

  /**
   * Vet a selection, one file at a time.
   *
   * Per file rather than all-or-nothing: one malformed file among five must not
   * throw away the other four, and the user has to be told WHICH one was
   * refused and why. The check is the backend's own scanner, so nothing is
   * judged by its extension.
   */
  const vetCandidates = useCallback(async (paths: readonly string[]) => {
    const accepted: FastaCandidate[] = [];
    const rejected: string[] = [];

    for (const path of paths) {
      const response = await desktop().project.validateFastaCandidate(path);
      if (response.ok) {
        accepted.push(response.result.candidate);
      } else {
        rejected.push(`${fileName(path)}: ${response.error.message}`);
      }
    }

    return { accepted, rejected };
  }, []);

  /**
   * Link paths, optionally carrying a lock the user set before the project
   * existed.
   *
   * The lock is applied through `project.setFastaFileLocked` the moment the
   * file HAS an id — the same call the row makes when it is clicked on a linked
   * source — so there is one lock mechanism, not a pending one and a real one.
   * A file that links but fails to lock stays linked and says so.
   */
  const linkFastaFiles = useCallback(
    async (paths: readonly string[], lockedPaths: readonly string[] = []) => {
      if (paths.length === 0) return [] as string[];

      const shouldLock = new Set(lockedPaths);
      const failures: string[] = [];
      for (const path of paths) {
        const response = await desktop().project.linkFasta(path);
        if (!response.ok) {
          failures.push(`${fileName(path)}: ${response.error.message}`);
          continue;
        }

        dispatch({ type: 'sourceUpdated', source: response.result.source });

        if (shouldLock.has(path)) {
          const locked = await desktop().project.setFastaFileLocked(
            response.result.fastaFileId,
            true,
          );
          if (locked.ok) {
            dispatch({ type: 'sourceUpdated', source: locked.result.source });
          } else {
            failures.push(`${fileName(path)}: ${locked.error.message}`);
          }
        }
      }
      return failures;
    },
    [],
  );

  const chooseFastaFiles = useCallback(async () => {
    const paths = await desktop().dialog.selectFastaFiles();
    if (paths.length === 0) return;

    dispatch({ type: 'setSourceError', message: null });
    const { accepted, rejected } = await vetCandidates(paths);

    if (projectIsOpen.current) {
      // The project exists: vetted files go straight in.
      const failures = await linkFastaFiles(accepted.map((candidate) => candidate.path));
      const problems = [...rejected, ...failures];
      if (problems.length > 0) dispatch({ type: 'setSourceError', message: problems.join(' · ') });
      return;
    }

    // No project yet: hold them, showing the same table the project will.
    if (accepted.length > 0) dispatch({ type: 'addPendingCandidates', candidates: accepted });
    if (rejected.length > 0) dispatch({ type: 'setSourceError', message: rejected.join(' · ') });
  }, [vetCandidates, linkFastaFiles]);

  const linkCandidates = useCallback(
    async (candidates: readonly FastaCandidate[], lockedPaths: readonly string[] = []) => {
      const failures = await linkFastaFiles(
        candidates.map((candidate) => candidate.path),
        lockedPaths,
      );
      if (failures.length > 0) dispatch({ type: 'setSourceError', message: failures.join(' · ') });
    },
    [linkFastaFiles],
  );

  const setSourceLocked = useCallback(async (fastaFileId: string, locked: boolean) => {
    const response = await desktop().project.setFastaFileLocked(fastaFileId, locked);
    if (response.ok) {
      dispatch({ type: 'sourceUpdated', source: response.result.source });
    } else {
      dispatch({ type: 'showNotice', message: response.error.message });
    }
  }, []);

  const unlinkSource = useCallback(async (fastaFileId: string) => {
    const response = await desktop().project.unlinkFasta(fastaFileId);
    if (response.ok) {
      dispatch({ type: 'sourceRemoved', fastaFileId });
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

  /* ------------------------------------------------------------------ *
   * Focal drafts
   * ------------------------------------------------------------------ */

  const saveActiveDraft = useCallback(async () => {
    const request = saveRequestFor(activeDraft);
    if (!request.title) {
      dispatch({ type: 'showNotice', message: 'Give the focal set a title before saving it.' });
      return false;
    }

    const response = await desktop().project.saveFocalSet(request);
    if (!response.ok) {
      dispatch({ type: 'showNotice', message: response.error.message });
      return false;
    }

    // The backend response is canonical — it carries the ids, ordering and
    // deduplication the database actually applied — so it replaces the working
    // copy rather than merely marking it clean.
    dispatch({ type: 'focalDraftSaved', key: activeDraft.key, focalSet: response.result.focalSet });
    return true;
  }, [activeDraft]);

  const setDraftLocked = useCallback(
    async (key: string, locked: boolean) => {
      const draft = state.focalDrafts.find((candidate) => candidate.key === key);
      if (!draft) return false;

      // Locking is only meaningful for something that exists and is settled.
      // Auto-saving first would commit a version the user never approved.
      const blocked = locked ? lockBlockedReason(draft) : null;
      if (blocked) {
        dispatch({ type: 'showNotice', message: blocked });
        return false;
      }
      if (draft.persistedId === null) {
        dispatch({ type: 'showNotice', message: 'Save this focal set before locking it.' });
        return false;
      }

      const response = await desktop().project.setFocalSetLocked(draft.persistedId, locked);
      if (!response.ok) {
        dispatch({ type: 'showNotice', message: response.error.message });
        return false;
      }

      // The canonical response decides, not the click.
      dispatch({ type: 'focalDraftLockChanged', key, focalSet: response.result.focalSet });
      return true;
    },
    [state.focalDrafts],
  );

  const deleteDraft = useCallback(
    async (key: string) => {
      const draft = state.focalDrafts.find((candidate) => candidate.key === key);
      if (!draft) return false;

      if (draft.locked) {
        dispatch({
          type: 'showNotice',
          message: 'This focal set is locked. Unlock it before deleting it.',
        });
        return false;
      }

      // A never-saved draft has no database row to remove; dropping it locally
      // IS the whole deletion.
      if (draft.persistedId !== null) {
        const response = await desktop().project.deleteFocalSet(draft.persistedId);
        if (!response.ok) {
          dispatch({ type: 'showNotice', message: response.error.message });
          return false;
        }
      }

      dispatch({ type: 'removeFocalDraft', key });
      return true;
    },
    [state.focalDrafts],
  );

  const addByQuery = useCallback(
    async (query: string) => {
      if (!query.trim()) {
        dispatch({ type: 'showNotice', message: 'Type something to search for.' });
        return false;
      }

      /*
       * Resolved through `resolveFocalAddQuery`, NOT `searchHeaders`.
       *
       * `searchHeaders` is a capped preview: expanding `+` through it silently
       * truncated a query matching more headers than the cap, so the focal set
       * ended up smaller than the user asked for. This operation is uncapped,
       * and refuses outright (SEARCH_SCOPE_UNAVAILABLE) if a file in the scope
       * could not be searched, rather than returning a partial answer.
       *
       * The scope here is the ANALYSIS/search scope — one file when one is
       * selected — not the wider presence comparison scope.
       */
      const response = await desktop().project.resolveFocalAddQuery(query, scopeSourceIds);
      if (!response.ok) {
        dispatch({ type: 'showNotice', message: response.error.message });
        return false;
      }

      const { headers } = response.result;
      if (headers.length === 0) {
        dispatch({ type: 'showNotice', message: `No header contains "${query.trim()}".` });
        return false;
      }

      dispatch({
        type: 'setDraftHeaders',
        key: activeDraft.key,
        // Deterministic hit order, deduplicated against what is already there
        // by `normaliseHeaders`. The query itself is never stored.
        // No `coalesce`: `+` is one deliberate action and one undo step.
        headers: [...activeDraft.headers, ...headers],
      });
      return true;
    },
    [activeDraft, scopeSourceIds],
  );

  const removeByQuery = useCallback(
    async (query: string) => {
      if (!query.trim()) {
        dispatch({ type: 'showNotice', message: 'Type something to remove.' });
        return false;
      }

      // Matched in Python so `−` uses the same casefold as everything else,
      // and against the WORKING copy, not the FASTA.
      const response = await desktop().project.matchFocalHeaders(query, activeDraft.headers);
      if (!response.ok) {
        dispatch({ type: 'showNotice', message: response.error.message });
        return false;
      }

      const doomed = new Set(response.result.matched);
      if (doomed.size === 0) {
        dispatch({
          type: 'showNotice',
          message: `No focal entry contains "${query.trim()}".`,
        });
        return false;
      }

      dispatch({
        type: 'setDraftHeaders',
        key: activeDraft.key,
        headers: activeDraft.headers.filter((header) => !doomed.has(header)),
      });
      return true;
    },
    [activeDraft],
  );

  /* ------------------------------------------------------------------ *
   * Presence
   * ------------------------------------------------------------------ */

  /*
   * Debounced, generation-tagged, and index-only.
   *
   * `project.headerPresence` opens no file and hashes nothing, so there is no
   * filesystem work to cancel — the debounce exists to avoid pointless round
   * trips, and the generation token exists so an answer about an older draft
   * can never repaint the current one.
   */
  const presenceGeneration = useRef(0);
  const draftHeaders = activeDraft.headers;
  const projectOpen = state.project.status === 'open';
  const comparisonKey = presenceComparisonIds.join(' ');
  const scopeKey = scopeSourceIds.join(' ');

  useEffect(() => {
    if (!projectOpen) return;

    const headers = normaliseHeaders(draftHeaders);
    if (headers.length === 0) {
      dispatch({ type: 'presenceCleared' });
      return;
    }

    /*
     * TWO questions, so up to two answers.
     *
     *  display — compared against every linked file, so "absent here, present
     *            over there" is sayable and can read orange;
     *  run     — compared against the analysis scope alone, because that is
     *            what a run will actually read.
     *
     * Gating on the display answer was a correctness bug: an unrelated
     * unavailable FASTA turned a header's display state UNKNOWN, and the gate
     * then blamed "a selected FASTA is unavailable" while the selected file was
     * healthy. Under All files the two scopes coincide, so one call serves both.
     */
    const sameScope = selectedFasta === null && comparisonKey === scopeKey;

    const ticket = ++presenceGeneration.current;
    const timer = setTimeout(() => {
      dispatch({ type: 'presenceRequested', generation: ticket });

      const display = desktop().project.headerPresence(
        headers,
        selectedFasta,
        presenceComparisonIds,
      );
      const run = sameScope
        ? display
        : desktop().project.headerPresence(headers, selectedFasta, scopeSourceIds);

      void Promise.all([display, run]).then(([displayResponse, runResponse]) => {
        if (!displayResponse.ok) {
          dispatch({ type: 'presenceFailed', generation: ticket, error: displayResponse.error });
          return;
        }
        if (!runResponse.ok) {
          dispatch({ type: 'presenceFailed', generation: ticket, error: runResponse.error });
          return;
        }
        dispatch({
          type: 'presenceLoaded',
          generation: ticket,
          entries: displayResponse.result.entries,
          runEntries: runResponse.result.entries,
        });
      });
    }, PRESENCE_DEBOUNCE_MS);

    return () => clearTimeout(timer);
    // The joined keys stand in for the id arrays so a re-render with equal
    // arrays does not re-ask.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [projectOpen, draftHeaders, selectedFasta, comparisonKey, scopeKey]);

  /*
   * Lifecycle refreshes.
   *
   * A linked file can change while the app is in the background, so the state
   * on screen is only trustworthy if it is re-read at the moments the user
   * comes back to it: regaining window focus, and entering a screen that shows
   * or acts on sources. Both take the cheap path.
   *
   * Note what is NOT here: nothing in the focal editor triggers this. Focal
   * edits are answered from the header index, so typing never causes a stat
   * sweep, let alone a rehash.
   */
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

  /* ------------------------------------------------------------------ *
   * Runs
   * ------------------------------------------------------------------ */

  const scopeSources = useMemo(
    () => state.sources.items.filter((item) => scopeSourceIds.includes(item.fastaFileId)),
    [state.sources.items, scopeSourceIds],
  );

  const runGate = useMemo(
    () =>
      evaluateRunGate({
        draft: activeDraft,
        // The RUN-scope answer, never the display one.
        presence: state.runPresence,
        singleFileScope: state.fastaScope.kind === 'file',
        scopeSources,
        running: state.diagnosisRun.status === 'running',
      }),
    [
      activeDraft,
      state.runPresence,
      state.fastaScope.kind,
      scopeSources,
      state.diagnosisRun.status,
    ],
  );

  const runDiagnosis = useCallback(
    async (resume: DiagnosisResumeState | null = null) => {
      /*
       * The WHOLE gate is enforced here, not just the saved-set half.
       *
       * A disabled button is a hint, not a guarantee: a keyboard path, a stale
       * render, or a future caller can all reach this function while the gate
       * is closed. Re-checking it here means every refusal — dirty set, unknown
       * presence, an entry missing from the scope, no linked file — is caught
       * with the same sentence the button was showing, instead of becoming a
       * backend error after the fact.
       */
      if (!runGate.canRun) {
        if (runGate.reason) dispatch({ type: 'showNotice', message: runGate.reason });
        return;
      }

      const focalSetId = activeDraft.persistedId;
      // Belt and braces: the gate already refuses a new draft, and there would
      // be no saved focal set for a report to name.
      if (focalSetId === null) {
        dispatch({
          type: 'showNotice',
          message: 'Save this focal set before running an analysis with it.',
        });
        return;
      }

      const config = state.molecularDiagnosis;
      const singleFile = state.fastaScope.kind === 'file';

      dispatch({ type: 'diagnosisStarted', continuing: resume !== null });

      const response = await desktop().project.runMolecularDiagnosis({
        focalSetId,
        fastaFileIds: scopeSourceIds,
        singleFile,
        options: {
          ignoreGaps: config.ignoreGaps,
          giveBenefitOfDoubtToAmbiguousBases: config.giveBenefitOfDoubtToAmbiguousBases,
          minCandidateSize: config.minCandidateSize,
          maxCandidateSize: config.maxCandidateSize,
        },
        resume,
      });

      if (response.ok) {
        dispatch({ type: 'diagnosisSucceeded', result: response.result });
      } else {
        dispatch({ type: 'diagnosisFailed', error: response.error });
      }
    },
    [
      runGate,
      activeDraft.persistedId,
      state.molecularDiagnosis,
      state.fastaScope.kind,
      scopeSourceIds,
    ],
  );

  const sourceViews = useMemo(
    () => describeSources(state.sources.items, state.sources.activity),
    [state.sources.items, state.sources.activity],
  );

  const value = useMemo<ProjectContextValue>(
    () => ({
      state,
      dispatch,
      activeDraft,
      saveActiveDraft,
      addByQuery,
      removeByQuery,
      setDraftLocked,
      deleteDraft,
      scopeSourceIds,
      runGate,
      runDiagnosis,
      sourceViews,
      sourceBanner: sourceBannerMessage(sourceViews),
      openProject,
      createProject,
      chooseProjectDirectory,
      setProjectTitle,
      refreshSources,
      chooseFastaFiles,
      linkCandidates,
      unlinkSource,
      setSourceLocked,
      reindexSource,
      relinkSource,
    }),
    [
      state,
      activeDraft,
      saveActiveDraft,
      addByQuery,
      removeByQuery,
      setDraftLocked,
      deleteDraft,
      scopeSourceIds,
      runGate,
      runDiagnosis,
      sourceViews,
      openProject,
      createProject,
      chooseProjectDirectory,
      setProjectTitle,
      refreshSources,
      chooseFastaFiles,
      linkCandidates,
      unlinkSource,
      setSourceLocked,
      reindexSource,
      relinkSource,
    ],
  );

  return <ProjectContext.Provider value={value}>{children}</ProjectContext.Provider>;
}

/** Show just the file name; the full path goes in a title attribute. */
export function fileName(filePath: string): string {
  const parts = filePath.split(/[\\/]/);
  return parts[parts.length - 1] || filePath;
}

export function useProject(): ProjectContextValue {
  const context = useContext(ProjectContext);
  if (!context) throw new Error('useProject must be used inside <ProjectProvider>.');
  return context;
}
