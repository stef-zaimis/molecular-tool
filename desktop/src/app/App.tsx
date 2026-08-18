import { TitleBar } from '../components/chrome/TitleBar';
import { TabStrip } from '../components/chrome/TabStrip';
import { LauncherScreen } from '../screens/LauncherScreen';
import { ProjectCreationScreen } from '../screens/ProjectCreationScreen';
import { MolecularDiagnosisScreen } from '../screens/MolecularDiagnosisScreen';
import { useProject } from './state/ProjectContext';
import { ANALYSIS_LABELS, ANALYSIS_ORDER } from './state/projectState';
import type { AnalysisKind } from '../contract';

/**
 * The workspace shell: window chrome + analysis tabs + the active analysis page.
 *
 * Tab availability is driven entirely by the analyses chosen during project
 * creation, so the two screens cannot disagree.
 */
function WorkspaceScreen(): JSX.Element {
  const { state, dispatch } = useProject();

  const tabs = ANALYSIS_ORDER.map((analysis) => ({
    id: analysis,
    label: ANALYSIS_LABELS[analysis],
    enabled: state.draft.analyses[analysis],
  }));

  return (
    <div className="app-shell">
      <TitleBar title="Molecular Diagnosis Tool" variant="workspace" />

      <TabStrip
        tabs={tabs}
        activeId={state.activeAnalysis}
        onHome={() => dispatch({ type: 'navigate', screen: 'projectCreation' })}
        onSelect={(id) => dispatch({ type: 'setActiveAnalysis', analysis: id as AnalysisKind })}
        homeLabel="Back to project settings"
      />

      <main className="app-body">
        {state.activeAnalysis === 'molecularDiagnosis' ? (
          <MolecularDiagnosisScreen />
        ) : (
          // Unreachable while the other tabs stay disabled; kept as an explicit
          // branch so adding an analysis screen is a one-line change.
          <div className="app-body" />
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
