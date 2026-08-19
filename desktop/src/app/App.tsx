import { TitleBar } from '../components/chrome/TitleBar';
import { TabStrip } from '../components/chrome/TabStrip';
import { LauncherScreen } from '../screens/LauncherScreen';
import { ProjectCreationScreen } from '../screens/ProjectCreationScreen';
import { MolecularDiagnosisScreen } from '../screens/MolecularDiagnosisScreen';
import { useProject } from './state/ProjectContext';
import { ANALYSIS_LABELS, selectedAnalyses } from './state/projectState';
import type { AnalysisKind } from '../contract';

/**
 * Placeholder for an analysis that is selected but has no workspace yet.
 *
 * Molecular Diagnosis is the only analysis screen in this build. Since tabs are
 * now generated from the project's selection, the others are reachable, and a
 * blank page would read as a bug rather than as unfinished work.
 */
function NotBuiltYet({ label }: { label: string }): JSX.Element {
  return (
    <div className="not-built">
      <p className="not-built__title">{label}</p>
      <p className="not-built__body">
        This analysis workspace has not been built yet. It is selected for the project, so it
        appears here, but only Molecular Diagnosis can be configured in this build.
      </p>
    </div>
  );
}

/**
 * The workspace shell: window chrome + analysis tabs + the active analysis page.
 *
 * The tab list IS the project's analysis selection — analyses that were not
 * selected have no tab, rather than a disabled one. Home is the project page,
 * so it navigates to project creation rather than acting as a history "back".
 */
function WorkspaceScreen(): JSX.Element {
  const { state, dispatch } = useProject();

  const tabs = selectedAnalyses(state.draft.analyses).map((analysis) => ({
    id: analysis,
    label: ANALYSIS_LABELS[analysis],
  }));

  return (
    <div className="app-shell">
      <TitleBar title="Molecular Diagnosis Tool" variant="workspace" />

      <TabStrip
        tabs={tabs}
        activeId={state.activeAnalysis}
        onHome={() => dispatch({ type: 'navigate', screen: 'projectCreation' })}
        onSelect={(id) => dispatch({ type: 'setActiveAnalysis', analysis: id as AnalysisKind })}
        homeLabel="Project page"
      />

      <main className="app-body">
        {state.activeAnalysis === 'molecularDiagnosis' ? (
          <MolecularDiagnosisScreen />
        ) : (
          <NotBuiltYet label={ANALYSIS_LABELS[state.activeAnalysis]} />
        )}
      </main>
    </div>
  );
}

export function App(): JSX.Element {
  const { state } = useProject();

  switch (state.screen) {
    case 'launcher':
      return <LauncherScreen />;
    case 'projectCreation':
      return <ProjectCreationScreen />;
    case 'workspace':
      return <WorkspaceScreen />;
    default:
      return <LauncherScreen />;
  }
}
