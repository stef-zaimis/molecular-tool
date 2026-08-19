import { createContext, useCallback, useContext, useMemo, useReducer } from 'react';
import type { Dispatch, ReactNode } from 'react';
import { appReducer, createInitialState, focalHistoryFor } from './projectState';
import type { AppAction, AppState, FocalHistory } from './projectState';
import type { FocalSet } from '../../contract';
import { desktop } from '../desktopApi';

interface ProjectContextValue {
  readonly state: AppState;
  readonly dispatch: Dispatch<AppAction>;
  /** Convenience: the focal set the workspace is currently editing. */
  readonly activeFocalSet: FocalSet;
  /** Undo/redo stacks for the active focal set. */
  readonly activeFocalHistory: FocalHistory;
  /**
   * Headers of the loaded FASTA, or null when none are available yet.
   * Null means "cannot validate", which the editor renders as neutral rather
   * than as a failed match.
   */
  readonly fastaHeaders: readonly string[] | null;
  /** Pick a FASTA and load its headers. Resolves once loading has settled. */
  readonly chooseFastaFile: () => Promise<void>;
}

const ProjectContext = createContext<ProjectContextValue | null>(null);

export function ProjectProvider({ children }: { children: ReactNode }): JSX.Element {
  const [state, dispatch] = useReducer(appReducer, undefined, createInitialState);

  const chooseFastaFile = useCallback(async () => {
    const path = await desktop().dialog.selectFastaFile();
    // A cancelled dialog returns null; leave any previous choice intact.
    if (!path) return;

    dispatch({ type: 'setFastaPath', path });
    dispatch({ type: 'alignmentLoading', path });

    const result = await desktop().fasta.readHeaders(path);
    if (result.ok) {
      dispatch({
        type: 'alignmentLoaded',
        data: { path: result.path, headers: result.headers, duplicateCount: result.duplicateCount },
      });
    } else {
      dispatch({
        type: 'alignmentFailed',
        path,
        message: result.message ?? 'The file could not be read.',
      });
    }
  }, []);

  const value = useMemo<ProjectContextValue>(() => {
    const activeFocalSet =
      state.focalSets.find((set) => set.id === state.activeFocalSetId) ?? state.focalSets[0];

    return {
      state,
      dispatch,
      activeFocalSet,
      activeFocalHistory: focalHistoryFor(state, activeFocalSet.id),
      fastaHeaders: state.alignment.status === 'loaded' ? state.alignment.data.headers : null,
      chooseFastaFile,
    };
  }, [state, chooseFastaFile]);

  return <ProjectContext.Provider value={value}>{children}</ProjectContext.Provider>;
}

export function useProject(): ProjectContextValue {
  const context = useContext(ProjectContext);
  if (!context) throw new Error('useProject must be used inside <ProjectProvider>.');
  return context;
}
