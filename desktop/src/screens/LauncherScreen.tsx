import { TitleBar } from '../components/chrome/TitleBar';
import { useProject } from '../app/state/ProjectContext';
import { desktop } from '../app/desktopApi';
import { NoticeBar } from '../components/controls/NoticeBar';
import { UNFINISHED_TEXT } from '../copy/helpText';
import './LauncherScreen.css';

/**
 * Screen 1 — docs/01-launcher.png.
 *
 * A compact opening window in the spirit of a setup wizard. Two stacked
 * buttons, each a bold uppercase action over a lighter noun.
 */
export function LauncherScreen(): JSX.Element {
  const { state, dispatch } = useProject();

  const createNewProject = () => {
    // Grow the OS window first so the larger screen does not paint into the
    // compact launcher frame for a beat.
    desktop().shell.enterWorkspaceLayout();
    dispatch({ type: 'navigate', screen: 'projectCreation' });
  };

  const openExistingProject = () => {
    // TODO(backend): there is no project file format in the repository yet, so
    // there is nothing to open. The control is reproduced from the design and
    // explains itself only when actually activated — see UI_NOTES Q2.
    dispatch({ type: 'showNotice', message: UNFINISHED_TEXT.openExistingProject });
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
          onClick={openExistingProject}
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
