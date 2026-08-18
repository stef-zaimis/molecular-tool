import { createContext, useContext, useMemo, useReducer } from 'react';
import type { Dispatch, ReactNode } from 'react';
import { appReducer, createInitialState } from './projectState';
import type { AppAction, AppState } from './projectState';
import type { FocalSet } from '../../contract';

interface ProjectContextValue {
  readonly state: AppState;
  readonly dispatch: Dispatch<AppAction>;
  /** Convenience: the focal set the workspace is currently editing. */
  readonly activeFocalSet: FocalSet;
}

const ProjectContext = createContext<ProjectContextValue | null>(null);

export function ProjectProvider({ children }: { children: ReactNode }): JSX.Element {
  const [state, dispatch] = useReducer(appReducer, undefined, createInitialState);

  const value = useMemo<ProjectContextValue>(() => {
    const activeFocalSet =
      state.focalSets.find((set) => set.id === state.activeFocalSetId) ?? state.focalSets[0];

    return { state, dispatch, activeFocalSet };
  }, [state]);

  return <ProjectContext.Provider value={value}>{children}</ProjectContext.Provider>;
}

export function useProject(): ProjectContextValue {
  const context = useContext(ProjectContext);
  if (!context) throw new Error('useProject must be used inside <ProjectProvider>.');
  return context;
}
