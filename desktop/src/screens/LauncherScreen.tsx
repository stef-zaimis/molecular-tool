import { TitleBar } from '../components/chrome/TitleBar';
import { useProject } from '../app/state/ProjectContext';
import { desktop } from '../app/desktopApi';
import { NoticeBar } from '../components/controls/NoticeBar';
import './LauncherScreen.css';

/**
 * Screen 1 — docs/01-launcher.png.
 *
 * A compact opening window in the spirit of a setup wizard. Two stacked
 * buttons, each a bold uppercase action over a lighter noun.
 */
export function LauncherScreen(): JSX.Element {
  const { state, dispatch, chooseProjectDirectory } = useProject();

  const createNewProject = () => {
    // Grow the OS window first so the larger screen does not paint into the
    // compact launcher frame for a beat.
    desktop().shell.enterWorkspaceLayout();
    dispatch({ type: 'navigate', screen: 'projectCreation' });
  };

  /*
   * A project is a directory holding `project.sqlite`. Opening one restores
   * its linked FASTA files and focal sets, and immediately reports the current
   * state of those files on disk — including any that have gone missing since
   * the project was last used.
   */
  const openExistingProject = async () => {
    const opened = await chooseProjectDirectory();
    if (!opened) return;
    desktop().shell.enterWorkspaceLayout();
    dispatch({ type: 'navigate', screen: 'projectCreation' });
  };

  return (
    <div className="app-shell launcher">
      <TitleBar title="Molecular Diagnosis Tool" variant="launcher" canMaximize={false} />

      <main className="launcher__body">
        <button type="button" className="launcher__button" onClick={createNewProject}>
          <span className="launcher__action">Create New</span>
          <span className="launcher__object">Project</span>
        </button>

        <button
          type="button"
          className="launcher__button"
          onClick={() => void openExistingProject()}
          aria-describedby={state.notice ? 'launcher-notice' : undefined}
        >
          <span className="launcher__action">Open Existing</span>
          <span className="launcher__object">Project</span>
        </button>

        {state.notice && (
          <div className="notice--launcher-wrap" id="launcher-notice">
            <NoticeBar
              message={state.notice}
              onDismiss={() => dispatch({ type: 'dismissNotice' })}
            />
          </div>
        )}
      </main>
    </div>
  );
}
