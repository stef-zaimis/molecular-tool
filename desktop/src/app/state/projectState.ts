import type {
  AnalysisKind,
  AnalysisSelection,
  FocalSet,
  MolecularDiagnosisConfig,
  ProjectDraft,
} from '../../contract';

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

/** Headers of the selected FASTA. Header lines only — nothing else is read. */
export interface AlignmentHeaders {
  readonly path: string;
  readonly headers: readonly string[];
  readonly duplicateCount: number;
}

export type AlignmentLoadState =
  | { readonly status: 'idle' }
  | { readonly status: 'loading'; readonly path: string }
  | { readonly status: 'loaded'; readonly data: AlignmentHeaders }
  | { readonly status: 'failed'; readonly path: string; readonly message: string };

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
  /** Headers of the selected FASTA, used to validate focal strings. */
  readonly alignment: AlignmentLoadState;
  /** Transient banner text, e.g. for actions that are not implemented yet. */
  readonly notice: string | null;
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
    notice: null,
  };
}

export type AppAction =
  | { type: 'navigate'; screen: ScreenId }
  | { type: 'setProjectName'; name: string }
  | { type: 'setFastaPath'; path: string | null }
  | { type: 'alignmentLoading'; path: string }
  | { type: 'alignmentLoaded'; data: AlignmentHeaders }
  | { type: 'alignmentFailed'; path: string; message: string }
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
  | { type: 'resetDraft' };

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
 * Two gates: a FASTA must be chosen (no analysis is meaningful without one),
 * and at least one analysis that actually has a screen must be selected.
 */
export function canEnterWorkspace(state: AppState): boolean {
  return state.draft.fastaPath !== null && state.draft.analyses.molecularDiagnosis;
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
        // Headers belong to the previous file; drop them until the new one loads.
        alignment: action.path === null ? { status: 'idle' } : state.alignment,
      };

    case 'alignmentLoading':
      return { ...state, alignment: { status: 'loading', path: action.path } };

    case 'alignmentLoaded':
      return { ...state, alignment: { status: 'loaded', data: action.data } };

    case 'alignmentFailed':
      return {
        ...state,
        alignment: { status: 'failed', path: action.path, message: action.message },
      };

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

    default:
      return state;
  }
}
