import { TitleBar } from '../components/chrome/TitleBar';
import { TabStrip } from '../components/chrome/TabStrip';
import { LauncherScreen } from '../screens/LauncherScreen';
import { ProjectScreen } from '../screens/ProjectScreen';
import { MolecularDiagnosisScreen } from '../screens/MolecularDiagnosisScreen';
import { useProject } from './state/ProjectContext';
import { useWorkspaceZoom } from './uiScale';
import { ANALYSIS_LABELS, IMPLEMENTED_ANALYSES, selectedAnalyses } from './state/projectState';
import type { AnalysisKind } from '../contract';

/**
 * Placeholder for an analysis that is selected but has no workspace yet.
 *
 * Molecular Diagnosis is the only analysis screen in this build. Since tabs are
 * generated from the project's selection, the others are reachable, and a blank
 * page would read as a bug rather than as unfinished work.
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
 * The project shell: window chrome, the analysis tab strip, and a page.
 *
 * Home IS the project page, and the analysis tabs are the project's own
 * selection — so moving between the project page and a workspace is one tab
 * strip, not a wizard with a Continue button.
 */
function ProjectShell({ screen }: { screen: 'project' | 'workspace' }): JSX.Element {
  const { state, dispatch } = useProject();

  // Analyses belong to a project, so there are no tabs until one exists.
  const tabs =
    state.project.status === 'open'
      ? selectedAnalyses(state.analyses).map((analysis) => ({
          id: analysis,
          label: ANALYSIS_LABELS[analysis],
        }))
      : [];

  const onProjectPage = screen === 'project';

  return (
    <div className="app-shell">
      <TitleBar title="Molecular Diagnosis Tool" variant="workspace" />

      <TabStrip
        tabs={tabs}
        activeId={onProjectPage ? null : state.activeAnalysis}
        homeActive={onProjectPage}
        onHome={() => dispatch({ type: 'navigate', screen: 'project' })}
        onSelect={(id) => {
          dispatch({ type: 'setActiveAnalysis', analysis: id as AnalysisKind });
          dispatch({ type: 'navigate', screen: 'workspace' });
        }}
        homeLabel="Project page"
      />

      {onProjectPage ? (
        <ProjectScreen />
      ) : (
        <main className="app-body">
          {IMPLEMENTED_ANALYSES.includes(state.activeAnalysis) ? (
            <MolecularDiagnosisScreen />
          ) : (
            <NotBuiltYet label={ANALYSIS_LABELS[state.activeAnalysis]} />
          )}
        </main>
      )}
    </div>
  );
}

export function App(): JSX.Element {
  const { state } = useProject();

  // Keeps the 1920x1080 design frame usable at 125% Windows scaling. See
  // uiScale.ts: at 1920 it is exactly 1 and changes nothing, and the launcher —
  // which has its own small window and its own frame — is never scaled.
  useWorkspaceZoom(state.screen !== 'launcher');

  switch (state.screen) {
    case 'launcher':
      return <LauncherScreen />;
    case 'project':
      return <ProjectShell screen="project" />;
    case 'workspace':
      return <ProjectShell screen="workspace" />;
    default:
      return <LauncherScreen />;
  }
}
