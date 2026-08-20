import type {
  AnalysisKind,
  AnalysisSelection,
  FocalSet,
  MolecularDiagnosisConfig,
  ProjectDraft,
} from '../../contract';
import type {
  BackendError,
  DiagnosisResumeState,
  FastaLoadResult,
  FocalStringValidation,
  MolecularDiagnosisResult,
  OpenProjectResult,
  ProjectCapabilities,
  SourceStatusPayload,
} from '../../backendContract';
import type { SourceActivity } from './sourceStatus';

/**
 * One coherent state model for the whole application.
 *
 * Screens read from here and dispatch into here; there are no ad-hoc click
 * handlers holding their own copies of project data. This is what makes the
 * project draft survive Launcher -> Creation -> Workspace navigation.
 *
 * Deliberately a plain useReducer + Context: the state is small and entirely
 * session-local, so a state-management dependency would not earn its weight.
 */

export type ScreenId = 'launcher' | 'projectCreation' | 'workspace';

/**
 * Focal-set editing is modal, matching the design: `+` and `−` are mutually
 * exclusive, and Enter applies whichever mode is selected.
 */
export type FocalEditMode = 'add' | 'remove';

/**
 * Undo/redo for focal-set edits. Each entry is a complete snapshot of the
 * set — the sets are small (tens of strings), so snapshots are simpler and
 * safer than replaying inverse operations.
 */
export interface FocalHistory {
  readonly past: readonly (readonly string[])[];
  readonly future: readonly (readonly string[])[];
}

/**
 * The selected FASTA, as parsed and validated by the Python backend.
 * The frontend never parses FASTA itself.
 */
export type AlignmentLoadState =
  | { readonly status: 'idle' }
  | { readonly status: 'loading'; readonly path: string }
  | { readonly status: 'loaded'; readonly data: FastaLoadResult }
  | { readonly status: 'failed'; readonly path: string; readonly error: BackendError };

/**
 * Per-focal-string green/red, as decided by the backend matcher.
 * Null means "not validated yet", which renders neutral rather than red.
 */
export type FocalValidationState =
  | { readonly status: 'idle' }
  | { readonly status: 'validating' }
  | { readonly status: 'loaded'; readonly results: readonly FocalStringValidation[]; readonly unionMatchCount: number }
  | { readonly status: 'failed'; readonly error: BackendError };

/** A Molecular Diagnosis run, including the continuation offer. */
export type DiagnosisRunState =
  | { readonly status: 'idle' }
  | { readonly status: 'running'; readonly continuing: boolean }
  | {
      readonly status: 'succeeded';
      readonly result: MolecularDiagnosisResult;
      /** Set when the search stopped only because it hit the maximum size. */
      readonly pendingContinuation: DiagnosisResumeState | null;
    }
  | { readonly status: 'failed'; readonly error: BackendError };

/**
 * The persistent project this session has open, if any.
 *
 * `sources` is a live reading of the filesystem refreshed at lifecycle points,
 * NOT a stored project property: nothing here is ever written back. It is
 * replaced wholesale by whatever the backend last reported.
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
  readonly draft: ProjectDraft;
  /**
   * Focal sets for the current draft. The UI currently edits exactly one;
   * modelled as a list because the design implies more than one may exist.
   */
  readonly focalSets: readonly FocalSet[];
  readonly activeFocalSetId: string;
  /** Per-focal-set undo/redo stacks, keyed by focal set id. */
  readonly focalHistory: Readonly<Record<string, FocalHistory>>;
  readonly focalMode: FocalEditMode;
  /** The parsed/validated FASTA, from the Python backend. */
  readonly alignment: AlignmentLoadState;
  /** Backend verdict on each focal string of the active set. */
  readonly focalValidation: FocalValidationState;
  /** The most recent Molecular Diagnosis run. */
  readonly diagnosisRun: DiagnosisRunState;
  /** Transient banner text, e.g. for actions that are not implemented yet. */
  readonly notice: string | null;
  /** The persistent project, when one is open. */
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

const EMPTY_HISTORY: FocalHistory = { past: [], future: [] };

const EMPTY_SOURCES: SourcesState = {
  items: [],
  activity: {},
  refreshing: false,
  checkedAt: null,
};

let idCounter = 0;
/** Session-local id. Not a persistence key — there is no persistence yet. */
export function nextId(prefix: string): string {
  idCounter += 1;
  return `${prefix}-${idCounter}`;
}

function createFocalSet(): FocalSet {
  return { id: nextId('focal'), title: '', strings: [] };
}

export function createInitialState(): AppState {
  const focalSet = createFocalSet();

  return {
    screen: 'launcher',
    activeAnalysis: 'molecularDiagnosis',
    draft: {
      id: nextId('project'),
      name: '',
      fastaPath: null,
      analyses: { ...DEFAULT_ANALYSES },
      molecularDiagnosis: { ...DEFAULT_MOLECULAR_DIAGNOSIS_CONFIG, focalSetId: focalSet.id },
      createdAt: new Date().toISOString(),
    },
    focalSets: [focalSet],
    activeFocalSetId: focalSet.id,
    focalHistory: { [focalSet.id]: EMPTY_HISTORY },
    focalMode: 'add',
    alignment: { status: 'idle' },
    focalValidation: { status: 'idle' },
    diagnosisRun: { status: 'idle' },
    notice: null,
    project: { status: 'closed' },
    sources: EMPTY_SOURCES,
  };
}

export type AppAction =
  | { type: 'navigate'; screen: ScreenId }
  | { type: 'setProjectName'; name: string }
  | { type: 'setFastaPath'; path: string | null }
  | { type: 'alignmentLoading'; path: string }
  | { type: 'alignmentLoaded'; data: FastaLoadResult }
  | { type: 'alignmentFailed'; path: string; error: BackendError }
  | { type: 'focalValidating' }
  | {
      type: 'focalValidated';
      results: readonly FocalStringValidation[];
      unionMatchCount: number;
    }
  | { type: 'focalValidationFailed'; error: BackendError }
  | { type: 'diagnosisStarted'; continuing: boolean }
  | { type: 'diagnosisSucceeded'; result: MolecularDiagnosisResult }
  | { type: 'diagnosisFailed'; error: BackendError }
  | { type: 'dismissContinuation' }
  | { type: 'clearDiagnosisRun' }
  | { type: 'toggleAnalysis'; analysis: AnalysisKind }
  | { type: 'setActiveAnalysis'; analysis: AnalysisKind }
  | { type: 'updateMolecularDiagnosis'; patch: Partial<MolecularDiagnosisConfig> }
  | { type: 'setFocalTitle'; focalSetId: string; title: string }
  | { type: 'setFocalMode'; mode: FocalEditMode }
  | { type: 'addFocalString'; focalSetId: string; value: string }
  | { type: 'removeFocalString'; focalSetId: string; value: string }
  | { type: 'undoFocalEdit'; focalSetId: string }
  | { type: 'redoFocalEdit'; focalSetId: string }
  | { type: 'addFocalSet' }
  | { type: 'removeActiveFocalSet' }
  | { type: 'showNotice'; message: string }
  | { type: 'dismissNotice' }
  | { type: 'resetDraft' }
  | { type: 'projectOpening'; projectDir: string }
  | { type: 'projectOpened'; result: OpenProjectResult }
  | { type: 'projectOpenFailed'; projectDir: string; error: BackendError }
  | { type: 'projectClosed' }
  | { type: 'sourcesRefreshing' }
  | { type: 'sourcesRefreshed'; sources: readonly SourceStatusPayload[] }
  | { type: 'sourcesRefreshFailed' }
  | { type: 'sourceActivity'; fastaFileId: string; activity: SourceActivity }
  | { type: 'sourceUpdated'; source: SourceStatusPayload };

/** The first enabled analysis, used when the active tab becomes unavailable. */
export function firstEnabledAnalysis(analyses: AnalysisSelection): AnalysisKind | null {
  return ANALYSIS_ORDER.find((kind) => analyses[kind]) ?? null;
}

/** Analyses selected for this project, in canonical order. */
export function selectedAnalyses(analyses: AnalysisSelection): readonly AnalysisKind[] {
  return ANALYSIS_ORDER.filter((kind) => analyses[kind]);
}

/**
 * Can the user leave project creation for an analysis workspace?
 *
 * Three gates: a FASTA must be chosen, that FASTA must have loaded and
 * validated in the Python backend (a ragged file is not analysable, so
 * entering the workspace with one would only fail later), and at least one
 * analysis that actually has a screen must be selected.
 */
export function canEnterWorkspace(state: AppState): boolean {
  return (
    state.draft.fastaPath !== null &&
    state.alignment.status === 'loaded' &&
    state.draft.analyses.molecularDiagnosis
  );
}

export function focalHistoryFor(state: AppState, focalSetId: string): FocalHistory {
  return state.focalHistory[focalSetId] ?? EMPTY_HISTORY;
}

/** Applies a focal-set change and records one undo step. */
function commitFocalStrings(
  state: AppState,
  focalSetId: string,
  next: readonly string[],
): AppState {
  const current = state.focalSets.find((set) => set.id === focalSetId);
  if (!current) return state;

  const history = focalHistoryFor(state, focalSetId);

  return {
    ...state,
    focalSets: state.focalSets.map((set) =>
      set.id === focalSetId ? { ...set, strings: next } : set,
    ),
    focalHistory: {
      ...state.focalHistory,
      // A fresh edit always invalidates the redo branch.
      [focalSetId]: { past: [...history.past, current.strings], future: [] },
    },
  };
}

export function appReducer(state: AppState, action: AppAction): AppState {
  switch (action.type) {
    case 'navigate':
      return { ...state, screen: action.screen, notice: null };

    case 'setProjectName':
      return { ...state, draft: { ...state.draft, name: action.name } };

    case 'setFastaPath':
      return {
        ...state,
        draft: { ...state.draft, fastaPath: action.path },
        // Results belong to the previous file; drop them until the new one loads.
        alignment: action.path === null ? { status: 'idle' } : state.alignment,
        focalValidation: { status: 'idle' },
        diagnosisRun: { status: 'idle' },
      };

    case 'alignmentLoading':
      return { ...state, alignment: { status: 'loading', path: action.path } };

    case 'alignmentLoaded':
      return {
        ...state,
        alignment: { status: 'loaded', data: action.data },
        // A new file invalidates the previous verdicts and the previous run.
        focalValidation: { status: 'idle' },
        diagnosisRun: { status: 'idle' },
      };

    case 'alignmentFailed':
      return {
        ...state,
        alignment: { status: 'failed', path: action.path, error: action.error },
        focalValidation: { status: 'idle' },
        diagnosisRun: { status: 'idle' },
      };

    case 'focalValidating':
      return { ...state, focalValidation: { status: 'validating' } };

    case 'focalValidated':
      return {
        ...state,
        focalValidation: {
          status: 'loaded',
          results: action.results,
          unionMatchCount: action.unionMatchCount,
        },
      };

    case 'focalValidationFailed':
      return { ...state, focalValidation: { status: 'failed', error: action.error } };

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

    case 'toggleAnalysis': {
      const analyses: AnalysisSelection = {
        ...state.draft.analyses,
        [action.analysis]: !state.draft.analyses[action.analysis],
      };

      // Keep the workspace tab valid: if the active analysis was just turned
      // off, fall back to the first one that is still enabled.
      const activeAnalysis = analyses[state.activeAnalysis]
        ? state.activeAnalysis
        : (firstEnabledAnalysis(analyses) ?? state.activeAnalysis);

      return { ...state, activeAnalysis, draft: { ...state.draft, analyses } };
    }

    case 'setActiveAnalysis':
      // Guard at the reducer, not just the button, so an unselected analysis
      // can never be activated by any code path.
      if (!state.draft.analyses[action.analysis]) return state;
      return { ...state, activeAnalysis: action.analysis };

    case 'updateMolecularDiagnosis':
      return {
        ...state,
        draft: {
          ...state.draft,
          molecularDiagnosis: { ...state.draft.molecularDiagnosis, ...action.patch },
        },
      };

    case 'setFocalTitle':
      return {
        ...state,
        focalSets: state.focalSets.map((set) =>
          set.id === action.focalSetId ? { ...set, title: action.title } : set,
        ),
      };

    case 'setFocalMode':
      return { ...state, focalMode: action.mode };

    case 'addFocalString': {
      const value = action.value.trim();
      const set = state.focalSets.find((s) => s.id === action.focalSetId);
      if (!set) return state;
      // Empty and duplicate adds change nothing, so they record no history.
      if (value.length === 0 || set.strings.includes(value)) return state;
      return commitFocalStrings(state, action.focalSetId, [...set.strings, value]);
    }

    case 'removeFocalString': {
      const value = action.value.trim();
      const set = state.focalSets.find((s) => s.id === action.focalSetId);
      if (!set) return state;
      // A removal that matches nothing is a no-op and must not push history.
      if (value.length === 0 || !set.strings.includes(value)) return state;
      return commitFocalStrings(
        state,
        action.focalSetId,
        set.strings.filter((entry) => entry !== value),
      );
    }

    case 'undoFocalEdit': {
      const history = focalHistoryFor(state, action.focalSetId);
      const set = state.focalSets.find((s) => s.id === action.focalSetId);
      if (!set || history.past.length === 0) return state;

      const previous = history.past[history.past.length - 1];
      return {
        ...state,
        focalSets: state.focalSets.map((s) =>
          s.id === action.focalSetId ? { ...s, strings: previous } : s,
        ),
        focalHistory: {
          ...state.focalHistory,
          [action.focalSetId]: {
            past: history.past.slice(0, -1),
            future: [set.strings, ...history.future],
          },
        },
      };
    }

    case 'redoFocalEdit': {
      const history = focalHistoryFor(state, action.focalSetId);
      const set = state.focalSets.find((s) => s.id === action.focalSetId);
      if (!set || history.future.length === 0) return state;

      const [next, ...rest] = history.future;
      return {
        ...state,
        focalSets: state.focalSets.map((s) =>
          s.id === action.focalSetId ? { ...s, strings: next } : s,
        ),
        focalHistory: {
          ...state.focalHistory,
          [action.focalSetId]: { past: [...history.past, set.strings], future: rest },
        },
      };
    }

    case 'addFocalSet': {
      const focalSet = createFocalSet();
      return {
        ...state,
        focalSets: [...state.focalSets, focalSet],
        activeFocalSetId: focalSet.id,
        focalHistory: { ...state.focalHistory, [focalSet.id]: EMPTY_HISTORY },
      };
    }

    case 'removeActiveFocalSet': {
      if (state.focalSets.length <= 1) return state;
      const remaining = state.focalSets.filter((set) => set.id !== state.activeFocalSetId);
      const history = { ...state.focalHistory };
      delete history[state.activeFocalSetId];
      return {
        ...state,
        focalSets: remaining,
        activeFocalSetId: remaining[0].id,
        focalHistory: history,
      };
    }

    case 'showNotice':
      return { ...state, notice: action.message };

    case 'dismissNotice':
      return { ...state, notice: null };

    case 'resetDraft':
      return { ...createInitialState(), screen: state.screen };

    /* ---------------------------------------------------------------- */
    /* Persistent project and linked-source status                      */
    /* ---------------------------------------------------------------- */

    case 'projectOpening':
      return {
        ...state,
        project: { status: 'opening', projectDir: action.projectDir },
        sources: EMPTY_SOURCES,
      };

    case 'projectOpened':
      // Opening already carries a full status sweep, so the workspace can show
      // a missing file immediately rather than after a first refresh tick.
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
      };

    case 'projectOpenFailed':
      return {
        ...state,
        project: { status: 'failed', projectDir: action.projectDir, error: action.error },
        sources: EMPTY_SOURCES,
      };

    case 'projectClosed':
      return { ...state, project: { status: 'closed' }, sources: EMPTY_SOURCES };

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
      return {
        ...state,
        sources: {
          ...state.sources,
          items: known
            ? state.sources.items.map((item) =>
                item.fastaFileId === action.source.fastaFileId ? action.source : item,
              )
            : [...state.sources.items, action.source],
          activity,
        },
      };
    }

    default:
      return state;
  }
}
