import type { AnalysisKind, AnalysisSelection, MolecularDiagnosisConfig } from '../../contract';
import type {
  BackendError,
  DiagnosisResumeState,
  FastaCandidate,
  FocalSetPayload,
  HeaderPresencePayload,
  OpenProjectResult,
  ProjectCapabilities,
  ProjectDiagnosisResult,
  SourceStatusPayload,
} from '../../backendContract';
import type { FocalDraft } from './focalDrafts';
import {
  blankDraft,
  draftFromPayload,
  draftSaved,
  endBurst,
  redoDraft,
  undoDraft,
  withHeaders,
} from './focalDrafts';
import type { SourceActivity } from './sourceStatus';

/**
 * One coherent state model for the whole application.
 *
 * Screens read from here and dispatch into here; there are no ad-hoc click
 * handlers holding their own copies of project data.
 *
 * The shape follows one architecture, end to end:
 *
 *   open/create project -> linked FASTA sources -> selected FASTA scope
 *     -> working focal-set draft -> saved persistent focal set -> run
 *
 * Everything durable lives in the project database and arrives here as backend
 * payloads. What this file owns is the SESSION: which screen is showing, which
 * scope and draft are selected, and what the user has typed but not yet saved.
 */

/**
 * ONE project screen, not two.
 *
 * Describing a project and working on one are the same page in the design: the
 * title field gains a Create button before the project exists and becomes a
 * heading afterwards, and the SAME FASTA table is used throughout. Routing them
 * as separate screens produced a visible jump that the design does not have.
 */
export type ScreenId = 'launcher' | 'project' | 'workspace';

/**
 * Focal-set editing is modal, matching the design: `+` and `−` are mutually
 * exclusive, and Enter applies whichever mode is selected.
 */
export type FocalEditMode = 'add' | 'remove';

/**
 * Which linked FASTAs an analysis and its presence colours apply to.
 *
 * A single file is the default because that is the common case and the only
 * one where "present in another file" (orange) means anything.
 */
export type FastaScope = { readonly kind: 'all' } | { readonly kind: 'file'; readonly fastaFileId: string };

/**
 * Where the sequence visualizer's left edge sits, as a fraction of the
 * workspace width — indexed by `viewerSnap`.
 *
 * Snap positions rather than a free-floating pixel value: the viewer overlays
 * the analysis column, so an arbitrary edge would let it be parked half across
 * a control it also covers. A short list keeps every arrangement one the
 * layout was designed for.
 *
 * `0` is full expansion — the edge reaches the left boundary of the workspace
 * content. The previous list stopped at 34%, which left the viewer unable to
 * cover the width it is eventually meant to fill. It never travels outside the
 * window: the last snap keeps the edge, and therefore the drag target, well
 * inside the right-hand side.
 */
export const VIEWER_SNAPS: readonly number[] = [0, 0.34, 0.546, 0.7, 0.85];

/** The design's resting position: left edge at x=1049 of 1920. */
export const DEFAULT_VIEWER_SNAP = 2;

/**
 * A project being described before it exists.
 *
 * The new-project screen runs before there is any database to write to, so the
 * title and the chosen FASTA paths are held here and committed by Create.
 */
export interface PendingProject {
  readonly title: string;
  /**
   * FASTA files already VETTED but not yet linked.
   *
   * Candidates, not paths: each was scanned by the backend before it got here,
   * so the pending rows show the same sequence/alignment figures the linked
   * rows will, and a file that cannot be used never reaches this list.
   */
  readonly candidates: readonly FastaCandidate[];
  /**
   * Paths the user has locked BEFORE the project exists, held by path because
   * a candidate has no `fastaFileId` yet.
   *
   * The lock is the same idea either side of Create — "stop editing this
   * link" — so the row behaves identically before and after. It is a local
   * intent while it lives here; Create hands each one to
   * `project.setFastaFileLocked` as soon as its file has a row to lock.
   */
  readonly lockedPaths: readonly string[];
  readonly creating: boolean;
}

/**
 * Presence of the working draft's headers, keyed BY HEADER rather than by
 * position.
 *
 * Reordering entries must not invalidate what is already known about them, and
 * an index-keyed map would silently mis-colour every entry after a move.
 */
export interface FocalPresenceState {
  readonly status: 'idle' | 'checking' | 'loaded' | 'failed';
  readonly byHeader: Readonly<Record<string, HeaderPresencePayload>>;
  /**
   * Request token. Async responses that do not carry the current generation are
   * stale and dropped — the user has typed since they were asked for.
   */
  readonly generation: number;
  readonly error: BackendError | null;
}

/** A Molecular Diagnosis run over the project scope, including the continuation offer. */
export type DiagnosisRunState =
  | { readonly status: 'idle' }
  | { readonly status: 'running'; readonly continuing: boolean }
  | {
      readonly status: 'succeeded';
      readonly result: ProjectDiagnosisResult;
      /** Set when the search stopped only because it hit the maximum size. */
      readonly pendingContinuation: DiagnosisResumeState | null;
    }
  | { readonly status: 'failed'; readonly error: BackendError };

/**
 * The persistent project this session has open.
 *
 * `sources` is a live reading of the filesystem refreshed at lifecycle points,
 * NOT a stored project property: nothing there is ever written back.
 */
export interface OpenProject {
  readonly projectDir: string;
  readonly outputsDir: string;
  readonly title: string;
  readonly capabilities: ProjectCapabilities;
}

export type ProjectOpenState =
  | { readonly status: 'closed' }
  | { readonly status: 'opening'; readonly projectDir: string }
  | { readonly status: 'open'; readonly project: OpenProject }
  | { readonly status: 'failed'; readonly projectDir: string; readonly error: BackendError };

export interface SourcesState {
  readonly items: readonly SourceStatusPayload[];
  /** Per-file transient work, so the UI can say "Re-indexing" instead of "Changed". */
  readonly activity: Readonly<Record<string, SourceActivity>>;
  /** True while a project-wide refresh is in flight. */
  readonly refreshing: boolean;
  /** Timestamp of the last completed refresh, for the "checked just now" line. */
  readonly checkedAt: number | null;
}

export interface AppState {
  readonly screen: ScreenId;
  /** Which workspace analysis tab is showing. */
  readonly activeAnalysis: AnalysisKind;
  /**
   * Analyses selected for this project.
   *
   * Session-local: the schema has nowhere to put this, and inventing a
   * migration for a checkbox would be the wrong trade. It resets when the app
   * restarts, which is visible and harmless.
   */
  readonly analyses: AnalysisSelection;
  readonly molecularDiagnosis: MolecularDiagnosisConfig;
  /** The project-to-be, until Create brings it into existence. */
  readonly pending: PendingProject;
  /**
   * Why the last FASTA selection was refused, shown inline by the picker.
   *
   * Kept out of the general notice bar because it belongs beside the control
   * that produced it — "The FASTA file needs to be aligned." means nothing
   * floating at the bottom of the window.
   */
  readonly sourceError: string | null;
  /** Working copies of the project's focal sets. Editing never writes. */
  readonly focalDrafts: readonly FocalDraft[];
  readonly activeFocalKey: string;
  readonly focalMode: FocalEditMode;
  /**
   * Presence for DISPLAY: compared against every linked file, so a header
   * absent from the selected one but present elsewhere can read orange.
   */
  readonly focalPresence: FocalPresenceState;
  /**
   * Presence for the RUN GATE: compared against the analysis scope only.
   *
   * Separate because the two answer different questions. Gating on the display
   * result made an unrelated unavailable FASTA report "a selected FASTA is
   * unavailable" while the actually-selected file was perfectly healthy — the
   * unknown came from a file the run would never read.
   */
  readonly runPresence: FocalPresenceState;
  /** Which linked FASTAs the workspace is working against. */
  readonly fastaScope: FastaScope;
  /** Whether the focal-set library occupies the area above the visualizer. */
  readonly libraryVisible: boolean;
  /**
   * Whether the visualizer is shown.
   *
   * Independent of the library: neither, either or both may be visible. Hiding
   * the viewer keeps `viewerSnap`, so reopening it restores the width the user
   * last chose rather than resetting to the default.
   */
  readonly viewerVisible: boolean;
  /** Index into `VIEWER_SNAPS`. */
  readonly viewerSnap: number;
  /**
   * The draft whose title the pencil just put into edit state, or null.
   *
   * Session-only and transient: the rename itself is an ordinary working-copy
   * edit, committed by Save like any other.
   */
  readonly titleEditingKey: string | null;
  /** The most recent Molecular Diagnosis run. */
  readonly diagnosisRun: DiagnosisRunState;
  /** Transient banner text. */
  readonly notice: string | null;
  readonly project: ProjectOpenState;
  /** Live state of every linked FASTA. */
  readonly sources: SourcesState;
}

export const ANALYSIS_ORDER: readonly AnalysisKind[] = [
  'molecularDiagnosis',
  'sequencePunishmentTest',
  'consensusSequenceGeneration',
];

export const ANALYSIS_LABELS: Readonly<Record<AnalysisKind, string>> = {
  molecularDiagnosis: 'Molecular Diagnosis',
  sequencePunishmentTest: 'Sequence Punishment Test',
  consensusSequenceGeneration: 'Consensus Sequence Generation',
};

/** Analyses that actually have a workspace in this build. */
export const IMPLEMENTED_ANALYSES: readonly AnalysisKind[] = ['molecularDiagnosis'];

/**
 * Defaults mirror the Python signatures documented in REPO_MAP.md §7 so the
 * frontend does not silently introduce different starting values:
 *   include_gappy_consensus_dmc_sites = False
 *   include_ambiguous_dmc_bd          = False
 *   min_combination_length            = 1
 *   max_combination_length            = 2
 */
export const DEFAULT_MOLECULAR_DIAGNOSIS_CONFIG: MolecularDiagnosisConfig = {
  ignoreGaps: false,
  giveBenefitOfDoubtToAmbiguousBases: false,
  minCandidateSize: 1,
  maxCandidateSize: 2,
  focalSetId: null,
};

const DEFAULT_ANALYSES: AnalysisSelection = {
  molecularDiagnosis: true,
  sequencePunishmentTest: false,
  consensusSequenceGeneration: false,
};

const EMPTY_SOURCES: SourcesState = {
  items: [],
  activity: {},
  refreshing: false,
  checkedAt: null,
};

const EMPTY_PRESENCE: FocalPresenceState = {
  status: 'idle',
  byHeader: {},
  generation: 0,
  error: null,
};

const EMPTY_PENDING: PendingProject = {
  title: '',
  candidates: [],
  lockedPaths: [],
  creating: false,
};

let idCounter = 0;
/** Session-local id, for things that never reach the database. */
export function nextId(prefix: string): string {
  idCounter += 1;
  return `${prefix}-${idCounter}`;
}

export function createInitialState(): AppState {
  const draft = blankDraft();

  return {
    screen: 'launcher',
    activeAnalysis: 'molecularDiagnosis',
    analyses: { ...DEFAULT_ANALYSES },
    molecularDiagnosis: { ...DEFAULT_MOLECULAR_DIAGNOSIS_CONFIG },
    pending: EMPTY_PENDING,
    sourceError: null,
    focalDrafts: [draft],
    activeFocalKey: draft.key,
    focalMode: 'add',
    focalPresence: EMPTY_PRESENCE,
    runPresence: EMPTY_PRESENCE,
    fastaScope: { kind: 'all' },
    libraryVisible: false,
    viewerVisible: true,
    viewerSnap: DEFAULT_VIEWER_SNAP,
    titleEditingKey: null,
    diagnosisRun: { status: 'idle' },
    notice: null,
    project: { status: 'closed' },
    sources: EMPTY_SOURCES,
  };
}

export type AppAction =
  | { type: 'navigate'; screen: ScreenId }
  /* new-project screen -------------------------------------------------- */
  | { type: 'setPendingTitle'; title: string }
  | { type: 'addPendingCandidates'; candidates: readonly FastaCandidate[] }
  | { type: 'removePendingCandidate'; path: string }
  | { type: 'setPendingCandidateLocked'; path: string; locked: boolean }
  | { type: 'projectCreating'; creating: boolean }
  | { type: 'setSourceError'; message: string | null }
  /* analyses ------------------------------------------------------------ */
  | { type: 'toggleAnalysis'; analysis: AnalysisKind }
  | { type: 'setActiveAnalysis'; analysis: AnalysisKind }
  | { type: 'updateMolecularDiagnosis'; patch: Partial<MolecularDiagnosisConfig> }
  /* focal drafts -------------------------------------------------------- */
  | { type: 'focalSetsLoaded'; focalSets: readonly FocalSetPayload[] }
  | { type: 'addFocalDraft' }
  | { type: 'selectFocalDraft'; key: string }
  | { type: 'setDraftTitle'; key: string; title: string }
  | {
      type: 'setDraftHeaders';
      key: string;
      headers: readonly string[];
      /** True for a keystroke inside a continuous manual typing burst. */
      coalesce?: boolean;
    }
  | { type: 'endDraftBurst'; key: string }
  | { type: 'undoDraftEdit'; key: string }
  | { type: 'redoDraftEdit'; key: string }
  | { type: 'focalDraftSaved'; key: string; focalSet: FocalSetPayload }
  | { type: 'focalDraftLockChanged'; key: string; focalSet: FocalSetPayload }
  | { type: 'removeFocalDraft'; key: string }
  | { type: 'setFocalMode'; mode: FocalEditMode }
  | { type: 'setTitleEditing'; key: string | null }
  /* library and visualizer ---------------------------------------------- */
  | { type: 'setLibraryVisible'; visible: boolean }
  | { type: 'setViewerVisible'; visible: boolean }
  | { type: 'setViewerSnap'; snap: number }
  /* presence ------------------------------------------------------------ */
  | { type: 'presenceRequested'; generation: number }
  | {
      type: 'presenceLoaded';
      generation: number;
      entries: readonly HeaderPresencePayload[];
      /** The run-scope answer, when it differs from the display one. */
      runEntries: readonly HeaderPresencePayload[];
    }
  | { type: 'presenceFailed'; generation: number; error: BackendError }
  | { type: 'presenceCleared' }
  /* scope --------------------------------------------------------------- */
  | { type: 'setFastaScope'; scope: FastaScope }
  /* run ----------------------------------------------------------------- */
  | { type: 'diagnosisStarted'; continuing: boolean }
  | { type: 'diagnosisSucceeded'; result: ProjectDiagnosisResult }
  | { type: 'diagnosisFailed'; error: BackendError }
  | { type: 'dismissContinuation' }
  | { type: 'clearDiagnosisRun' }
  /* notices and project lifecycle --------------------------------------- */
  | { type: 'showNotice'; message: string }
  | { type: 'dismissNotice' }
  | { type: 'projectOpening'; projectDir: string }
  | { type: 'projectOpened'; result: OpenProjectResult }
  | { type: 'projectOpenFailed'; projectDir: string; error: BackendError }
  | { type: 'projectClosed' }
  | { type: 'projectTitleChanged'; title: string }
  | { type: 'sourcesRefreshing' }
  | { type: 'sourcesRefreshed'; sources: readonly SourceStatusPayload[] }
  | { type: 'sourcesRefreshFailed' }
  | { type: 'sourceActivity'; fastaFileId: string; activity: SourceActivity }
  | { type: 'sourceUpdated'; source: SourceStatusPayload }
  | { type: 'sourceRemoved'; fastaFileId: string };

/** The first enabled analysis, used when the active tab becomes unavailable. */
export function firstEnabledAnalysis(analyses: AnalysisSelection): AnalysisKind | null {
  return ANALYSIS_ORDER.find((kind) => analyses[kind]) ?? null;
}

/** Analyses selected for this project, in canonical order. */
export function selectedAnalyses(analyses: AnalysisSelection): readonly AnalysisKind[] {
  return ANALYSIS_ORDER.filter((kind) => analyses[kind]);
}

export function activeDraft(state: AppState): FocalDraft {
  return state.focalDrafts.find((draft) => draft.key === state.activeFocalKey) ?? state.focalDrafts[0];
}

/**
 * The linked files the current scope covers, in the project's own order.
 *
 * `all` deliberately resolves to every linked file rather than to "no filter":
 * the run and the presence check must both be told exactly which files they
 * cover, so nothing can be judged against a file the run will not read.
 */
export function scopeFileIds(state: AppState): readonly string[] {
  const all = state.sources.items.map((item) => item.fastaFileId);
  const scope = state.fastaScope;
  if (scope.kind === 'all') return all;
  return all.filter((id) => id === scope.fastaFileId);
}

/** The selected file id, or null for the All-files scope. */
export function selectedFileId(state: AppState): string | null {
  return state.fastaScope.kind === 'file' ? state.fastaScope.fastaFileId : null;
}

/**
 * Keep the scope pointing at a file that still exists.
 *
 * Defaults to the first linked file — the design's resting state — and falls
 * back to it when the selected file is unlinked. It never silently promotes a
 * vanished selection to "All files", which would widen an analysis without the
 * user asking.
 */
function reconcileScope(scope: FastaScope, sources: readonly SourceStatusPayload[]): FastaScope {
  if (sources.length === 0) return { kind: 'all' };
  if (scope.kind === 'all') return scope;
  const stillLinked = sources.some((item) => item.fastaFileId === scope.fastaFileId);
  return stillLinked ? scope : { kind: 'file', fastaFileId: sources[0].fastaFileId };
}

/** The scope a freshly opened project starts in: its first linked file. */
function initialScope(sources: readonly SourceStatusPayload[]): FastaScope {
  return sources.length > 0 ? { kind: 'file', fastaFileId: sources[0].fastaFileId } : { kind: 'all' };
}

function mapDraft(
  state: AppState,
  key: string,
  update: (draft: FocalDraft) => FocalDraft,
): AppState {
  const current = state.focalDrafts.find((draft) => draft.key === key);
  if (!current) return state;
  const next = update(current);
  if (next === current) return state;
  return {
    ...state,
    focalDrafts: state.focalDrafts.map((draft) => (draft.key === key ? next : draft)),
  };
}

/**
 * A locked set is read-only, enforced here as well as in Python.
 *
 * The backend refusal is the one that matters; this stops the working copy
 * from drifting into a state that could only ever be rejected.
 */
function editable(draft: FocalDraft): boolean {
  return !draft.locked;
}

export function appReducer(state: AppState, action: AppAction): AppState {
  switch (action.type) {
    case 'navigate':
      return { ...state, screen: action.screen, notice: null, sourceError: null };

    /* ---------------------------------------------------------------- */
    /* The project-to-be                                                */
    /* ---------------------------------------------------------------- */

    case 'setPendingTitle':
      return { ...state, pending: { ...state.pending, title: action.title } };

    case 'addPendingCandidates': {
      // The picker can return a file already in the list, and linking the same
      // path twice is refused by the backend anyway.
      const known = new Set(state.pending.candidates.map((candidate) => candidate.path));
      const added = action.candidates.filter((candidate) => !known.has(candidate.path));
      if (added.length === 0) return state;
      return {
        ...state,
        pending: { ...state.pending, candidates: [...state.pending.candidates, ...added] },
      };
    }

    case 'removePendingCandidate':
      return {
        ...state,
        pending: {
          ...state.pending,
          candidates: state.pending.candidates.filter(
            (candidate) => candidate.path !== action.path,
          ),
          // Nothing may hold a lock for a row that is gone: re-adding the same
          // file must come back unlocked, like any other fresh choice.
          lockedPaths: state.pending.lockedPaths.filter((path) => path !== action.path),
        },
      };

    case 'setPendingCandidateLocked': {
      const held = state.pending.lockedPaths.includes(action.path);
      if (held === action.locked) return state;
      return {
        ...state,
        pending: {
          ...state.pending,
          lockedPaths: action.locked
            ? [...state.pending.lockedPaths, action.path]
            : state.pending.lockedPaths.filter((path) => path !== action.path),
        },
      };
    }

    case 'projectCreating':
      return { ...state, pending: { ...state.pending, creating: action.creating } };

    case 'setSourceError':
      return { ...state, sourceError: action.message };

    /* ---------------------------------------------------------------- */
    /* Analyses                                                          */
    /* ---------------------------------------------------------------- */

    case 'toggleAnalysis': {
      const analyses: AnalysisSelection = {
        ...state.analyses,
        [action.analysis]: !state.analyses[action.analysis],
      };

      // Keep the workspace tab valid: if the active analysis was just turned
      // off, fall back to the first one that is still enabled.
      const activeAnalysis = analyses[state.activeAnalysis]
        ? state.activeAnalysis
        : (firstEnabledAnalysis(analyses) ?? state.activeAnalysis);

      return { ...state, activeAnalysis, analyses };
    }

    case 'setActiveAnalysis':
      // Guarded at the reducer, not just the button, so an unselected analysis
      // can never be activated by any code path.
      if (!state.analyses[action.analysis]) return state;
      return { ...state, activeAnalysis: action.analysis };

    case 'updateMolecularDiagnosis':
      return {
        ...state,
        molecularDiagnosis: { ...state.molecularDiagnosis, ...action.patch },
      };

    /* ---------------------------------------------------------------- */
    /* Focal drafts — working copies, never writes                       */
    /* ---------------------------------------------------------------- */

    case 'focalSetsLoaded': {
      // A project with no focal sets gets one blank local draft. No database
      // row is created for it: an empty, untitled set the user never asked for
      // would be indistinguishable from one they made and abandoned.
      const drafts =
        action.focalSets.length > 0 ? action.focalSets.map(draftFromPayload) : [blankDraft()];
      const active = drafts.some((draft) => draft.key === state.activeFocalKey)
        ? state.activeFocalKey
        : drafts[0].key;
      return {
        ...state,
        focalDrafts: drafts,
        activeFocalKey: active,
        focalPresence: EMPTY_PRESENCE,
        runPresence: EMPTY_PRESENCE,
      };
    }

    case 'addFocalDraft': {
      const draft = blankDraft();
      return {
        ...state,
        focalDrafts: [...state.focalDrafts, draft],
        activeFocalKey: draft.key,
        focalPresence: EMPTY_PRESENCE,
        runPresence: EMPTY_PRESENCE,
        titleEditingKey: draft.key,
      };
    }

    case 'selectFocalDraft': {
      // Switching away from a dirty draft keeps it dirty. Nothing is saved on
      // the way out: an implicit save would commit a set the user was still
      // deciding about.
      if (!state.focalDrafts.some((draft) => draft.key === action.key)) return state;
      /*
       * Re-selecting the draft that is ALREADY active is a strict no-op.
       *
       * Falling through would clear presence, drop the rename state and make
       * the token colours flicker through neutral on the way back to the same
       * answer. Guarded here as well as in the row handler, so no path can do
       * it.
       */
      if (action.key === state.activeFocalKey) return state;
      const closed = mapDraft(state, state.activeFocalKey, endBurst);
      return {
        ...closed,
        activeFocalKey: action.key,
        focalPresence: EMPTY_PRESENCE,
        runPresence: EMPTY_PRESENCE,
        titleEditingKey: null,
      };
    }

    case 'setDraftTitle':
      return mapDraft(state, action.key, (draft) =>
        editable(draft) ? { ...draft, title: action.title } : draft,
      );

    case 'setDraftHeaders':
      return mapDraft(state, action.key, (draft) =>
        editable(draft)
          ? withHeaders(draft, action.headers, { coalesce: action.coalesce ?? false })
          : draft,
      );

    case 'endDraftBurst':
      return mapDraft(state, action.key, endBurst);

    case 'undoDraftEdit':
      return mapDraft(state, action.key, (draft) =>
        editable(draft) ? undoDraft(draft) : draft,
      );

    case 'redoDraftEdit':
      return mapDraft(state, action.key, (draft) =>
        editable(draft) ? redoDraft(draft) : draft,
      );

    case 'focalDraftSaved':
      return {
        ...mapDraft(state, action.key, (draft) => draftSaved(draft, action.focalSet)),
        // A save closes the rename the pencil opened.
        titleEditingKey: state.titleEditingKey === action.key ? null : state.titleEditingKey,
      };

    case 'focalDraftLockChanged':
      // ONLY the lock is adopted from the response. Locking does not change the
      // stored membership, and overwriting the working copy here would silently
      // discard edits the user had not saved.
      return mapDraft(state, action.key, (draft) => ({
        ...draft,
        locked: action.focalSet.locked,
      }));

    case 'removeFocalDraft': {
      const index = state.focalDrafts.findIndex((draft) => draft.key === action.key);
      if (index === -1) return state;

      const remaining = state.focalDrafts.filter((draft) => draft.key !== action.key);
      // A project always has somewhere to type. If the last draft goes, a blank
      // LOCAL one takes its place — no database row is created for it.
      const drafts = remaining.length > 0 ? remaining : [blankDraft()];

      const activeFocalKey =
        state.activeFocalKey === action.key
          ? // The neighbour that slid into its place, or the new last one.
            (drafts[Math.min(index, drafts.length - 1)] ?? drafts[0]).key
          : state.activeFocalKey;

      return {
        ...state,
        focalDrafts: drafts,
        activeFocalKey,
        focalPresence:
          activeFocalKey === state.activeFocalKey ? state.focalPresence : EMPTY_PRESENCE,
        runPresence:
          activeFocalKey === state.activeFocalKey ? state.runPresence : EMPTY_PRESENCE,
        titleEditingKey: state.titleEditingKey === action.key ? null : state.titleEditingKey,
      };
    }

    case 'setFocalMode':
      return { ...state, focalMode: action.mode };

    case 'setTitleEditing':
      return { ...state, titleEditingKey: action.key };

    case 'setLibraryVisible':
      return { ...state, libraryVisible: action.visible };

    case 'setViewerVisible':
      // `viewerSnap` is deliberately untouched: the width survives hiding.
      return { ...state, viewerVisible: action.visible };

    case 'setViewerSnap':
      return {
        ...state,
        viewerSnap: Math.max(0, Math.min(VIEWER_SNAPS.length - 1, action.snap)),
      };

    /* ---------------------------------------------------------------- */
    /* Presence                                                          */
    /* ---------------------------------------------------------------- */

    case 'presenceRequested':
      return {
        ...state,
        focalPresence: { ...state.focalPresence, status: 'checking', generation: action.generation },
        runPresence: { ...state.runPresence, status: 'checking', generation: action.generation },
      };

    case 'presenceLoaded': {
      // Stale responses are dropped rather than applied: the draft has moved on.
      if (action.generation !== state.focalPresence.generation) return state;
      const index = (entries: readonly HeaderPresencePayload[]) => {
        const byHeader: Record<string, HeaderPresencePayload> = {};
        for (const entry of entries) byHeader[entry.header] = entry;
        return byHeader;
      };
      return {
        ...state,
        focalPresence: {
          status: 'loaded',
          byHeader: index(action.entries),
          generation: action.generation,
          error: null,
        },
        runPresence: {
          status: 'loaded',
          byHeader: index(action.runEntries),
          generation: action.generation,
          error: null,
        },
      };
    }

    case 'presenceFailed':
      if (action.generation !== state.focalPresence.generation) return state;
      return {
        ...state,
        focalPresence: { ...state.focalPresence, status: 'failed', error: action.error },
        runPresence: { ...state.runPresence, status: 'failed', error: action.error },
      };

    case 'presenceCleared':
      return {
        ...state,
        focalPresence: { ...EMPTY_PRESENCE, generation: state.focalPresence.generation },
        runPresence: { ...EMPTY_PRESENCE, generation: state.runPresence.generation },
      };

    /* ---------------------------------------------------------------- */
    /* Scope                                                             */
    /* ---------------------------------------------------------------- */

    case 'setFastaScope': {
      // Re-selecting what is already selected must change nothing. Clearing
      // presence here would drop the colours while leaving the presence
      // effect's inputs identical, so nothing would ask for them again and the
      // set would sit neutral until something unrelated moved.
      const same =
        state.fastaScope.kind === action.scope.kind &&
        (state.fastaScope.kind !== 'file' ||
          action.scope.kind !== 'file' ||
          state.fastaScope.fastaFileId === action.scope.fastaFileId);
      if (same) return state;

      // Presence answers are scope-specific: what is orange for one file is
      // green for another, so they are dropped rather than shown stale.
      return {
        ...state,
        fastaScope: action.scope,
        focalPresence: { ...EMPTY_PRESENCE, generation: state.focalPresence.generation },
        runPresence: { ...EMPTY_PRESENCE, generation: state.runPresence.generation },
      };
    }

    /* ---------------------------------------------------------------- */
    /* Runs                                                              */
    /* ---------------------------------------------------------------- */

    case 'diagnosisStarted':
      return { ...state, diagnosisRun: { status: 'running', continuing: action.continuing } };

    case 'diagnosisSucceeded':
      return {
        ...state,
        diagnosisRun: {
          status: 'succeeded',
          result: action.result,
          pendingContinuation: action.result.canContinue ? action.result.resume : null,
        },
      };

    case 'diagnosisFailed':
      return { ...state, diagnosisRun: { status: 'failed', error: action.error } };

    case 'dismissContinuation':
      return state.diagnosisRun.status === 'succeeded'
        ? { ...state, diagnosisRun: { ...state.diagnosisRun, pendingContinuation: null } }
        : state;

    case 'clearDiagnosisRun':
      return { ...state, diagnosisRun: { status: 'idle' } };

    /* ---------------------------------------------------------------- */
    /* Notices                                                           */
    /* ---------------------------------------------------------------- */

    case 'showNotice':
      return { ...state, notice: action.message };

    case 'dismissNotice':
      return { ...state, notice: null };

    /* ---------------------------------------------------------------- */
    /* Persistent project and linked-source status                       */
    /* ---------------------------------------------------------------- */

    case 'projectOpening':
      return {
        ...state,
        project: { status: 'opening', projectDir: action.projectDir },
        sources: EMPTY_SOURCES,
      };

    case 'projectOpened': {
      // Opening already carries a full status sweep, so the project page can
      // show a missing file immediately rather than after a first refresh tick.
      const draft = blankDraft();
      return {
        ...state,
        project: {
          status: 'open',
          project: {
            projectDir: action.result.projectDir,
            outputsDir: action.result.outputsDir,
            title: action.result.metadata.title,
            capabilities: action.result.capabilities,
          },
        },
        sources: {
          items: action.result.sources,
          activity: {},
          refreshing: false,
          checkedAt: Date.now(),
        },
        // Focal sets belong to the project that was just closed; the new
        // project's own are loaded straight after opening.
        focalDrafts: [draft],
        activeFocalKey: draft.key,
        focalPresence: EMPTY_PRESENCE,
        runPresence: EMPTY_PRESENCE,
        fastaScope: initialScope(action.result.sources),
        diagnosisRun: { status: 'idle' },
        pending: EMPTY_PENDING,
        titleEditingKey: null,
      };
    }

    case 'projectOpenFailed':
      return {
        ...state,
        project: { status: 'failed', projectDir: action.projectDir, error: action.error },
        sources: EMPTY_SOURCES,
      };

    case 'projectClosed':
      return { ...state, project: { status: 'closed' }, sources: EMPTY_SOURCES };

    case 'projectTitleChanged':
      return state.project.status === 'open'
        ? {
            ...state,
            project: {
              status: 'open',
              project: { ...state.project.project, title: action.title },
            },
          }
        : state;

    case 'sourcesRefreshing':
      return { ...state, sources: { ...state.sources, refreshing: true } };

    case 'sourcesRefreshed':
      return {
        ...state,
        sources: {
          items: action.sources,
          // A completed sweep supersedes any per-file activity it observed.
          activity: {},
          refreshing: false,
          checkedAt: Date.now(),
        },
        fastaScope: reconcileScope(state.fastaScope, action.sources),
      };

    case 'sourcesRefreshFailed':
      // Keep the last known statuses: a failed refresh is not evidence that
      // anything changed, and blanking the list would look like data loss.
      return { ...state, sources: { ...state.sources, refreshing: false } };

    case 'sourceActivity': {
      const activity = { ...state.sources.activity };
      if (action.activity === 'idle') {
        delete activity[action.fastaFileId];
      } else {
        activity[action.fastaFileId] = action.activity;
      }
      return { ...state, sources: { ...state.sources, activity } };
    }

    case 'sourceUpdated': {
      const activity = { ...state.sources.activity };
      delete activity[action.source.fastaFileId];
      const known = state.sources.items.some(
        (item) => item.fastaFileId === action.source.fastaFileId,
      );
      const items = known
        ? state.sources.items.map((item) =>
            item.fastaFileId === action.source.fastaFileId ? action.source : item,
          )
        : [...state.sources.items, action.source];

      return {
        ...state,
        sources: { ...state.sources, items, activity },
        // The first file linked into an empty project becomes the scope, so
        // the workspace opens on something rather than on "All files".
        fastaScope: state.sources.items.length === 0 ? initialScope(items) : state.fastaScope,
      };
    }

    case 'sourceRemoved': {
      const items = state.sources.items.filter(
        (item) => item.fastaFileId !== action.fastaFileId,
      );
      const activity = { ...state.sources.activity };
      delete activity[action.fastaFileId];
      return {
        ...state,
        sources: { ...state.sources, items, activity },
        fastaScope: reconcileScope(state.fastaScope, items),
      };
    }

    default:
      return state;
  }
}
