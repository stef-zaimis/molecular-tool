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
 *
 * NOTE: the mockup (docs/03) shows "Ignore gaps" ticked. That is a mockup
 * state, not a default; we keep the Python default of unticked. See UI_NOTES.
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
    notice: null,
  };
}

export type AppAction =
  | { type: 'navigate'; screen: ScreenId }
  | { type: 'setProjectName'; name: string }
  | { type: 'setFastaPath'; path: string | null }
  | { type: 'toggleAnalysis'; analysis: AnalysisKind }
  | { type: 'setActiveAnalysis'; analysis: AnalysisKind }
  | { type: 'updateMolecularDiagnosis'; patch: Partial<MolecularDiagnosisConfig> }
  | { type: 'setFocalTitle'; focalSetId: string; title: string }
  | { type: 'setFocalStrings'; focalSetId: string; strings: readonly string[] }
  | { type: 'addFocalSet' }
  | { type: 'removeActiveFocalSet' }
  | { type: 'showNotice'; message: string }
  | { type: 'dismissNotice' }
  | { type: 'resetDraft' };

/** The first enabled analysis, used when the active tab becomes unavailable. */
export function firstEnabledAnalysis(analyses: AnalysisSelection): AnalysisKind | null {
  return ANALYSIS_ORDER.find((kind) => analyses[kind]) ?? null;
}

export function appReducer(state: AppState, action: AppAction): AppState {
  switch (action.type) {
    case 'navigate':
      return { ...state, screen: action.screen, notice: null };

    case 'setProjectName':
      return { ...state, draft: { ...state.draft, name: action.name } };

    case 'setFastaPath':
      return { ...state, draft: { ...state.draft, fastaPath: action.path } };

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
      // Guard at the reducer, not just the button, so disabled tabs can never
      // be activated by any code path.
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

    case 'setFocalStrings':
      return {
        ...state,
        focalSets: state.focalSets.map((set) =>
          set.id === action.focalSetId ? { ...set, strings: action.strings } : set,
        ),
      };

    case 'addFocalSet': {
      const focalSet = createFocalSet();
      return {
        ...state,
        focalSets: [...state.focalSets, focalSet],
        activeFocalSetId: focalSet.id,
      };
    }

    case 'removeActiveFocalSet': {
      if (state.focalSets.length <= 1) return state;
      const remaining = state.focalSets.filter((set) => set.id !== state.activeFocalSetId);
      return { ...state, focalSets: remaining, activeFocalSetId: remaining[0].id };
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
