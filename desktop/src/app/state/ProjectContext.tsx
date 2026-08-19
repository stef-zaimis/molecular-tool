import { createContext, useCallback, useContext, useEffect, useMemo, useReducer, useRef } from 'react';
import type { Dispatch, ReactNode } from 'react';
import { appReducer, createInitialState, focalHistoryFor } from './projectState';
import type { AppAction, AppState, FocalHistory } from './projectState';
import type { FocalSet } from '../../contract';
import type { DiagnosisResumeState, FocalStringValidation } from '../../backendContract';
import { desktop } from '../desktopApi';

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
}

const ProjectContext = createContext<ProjectContextValue | null>(null);

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
    };
  }, [state, activeFocalSet, chooseFastaFile, runDiagnosis]);

  return <ProjectContext.Provider value={value}>{children}</ProjectContext.Provider>;
}

export function useProject(): ProjectContextValue {
  const context = useContext(ProjectContext);
  if (!context) throw new Error('useProject must be used inside <ProjectProvider>.');
  return context;
}
