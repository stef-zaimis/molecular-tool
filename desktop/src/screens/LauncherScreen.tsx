import { TitleBar } from '../components/chrome/TitleBar';
import { useProject } from '../app/state/ProjectContext';
import { desktop } from '../app/desktopApi';
import { NoticeBar } from '../components/controls/NoticeBar';
import './LauncherScreen.css';

/**
 * Screen 1 — the `page1` layer of docs/new/all_pages-20Aug2026.svg.
 *
 * A compact opening window in the spirit of a setup wizard: two stacked
 * actions, each a bold uppercase verb over a lighter noun, with Recents and
 * Browse as the two ways of reaching an existing project.
 *
 * The two primary actions are genuinely different backend operations, not two
 * routes to one. Create initialises a project; Open refuses a folder that does
 * not already contain one.
 */
export function LauncherScreen(): JSX.Element {
  const { state, dispatch, chooseProjectDirectory } = useProject();

  const createNewProject = () => {
    // Grow the OS window first so the larger screen does not paint into the
    // compact launcher frame for a beat.
    desktop().shell.enterWorkspaceLayout();
    dispatch({ type: 'navigate', screen: 'project' });
  };

  /*
   * A project is a directory holding `project.sqlite`. Opening one restores its
   * linked FASTA files and focal sets, and immediately reports the state of
   * those files on disk — including any that went missing since it was last
   * used. A folder without one is reported, never silently initialised.
   */
  const openExistingProject = async () => {
    const opened = await chooseProjectDirectory();
    if (!opened) return;
    desktop().shell.enterWorkspaceLayout();
    dispatch({ type: 'navigate', screen: 'project' });
  };

  return (
    <div className="app-shell launcher">
      <TitleBar title="Molecular Diagnosis Tool" variant="launcher" canMaximize={false} />

      <main className="launcher__body">
        <button
          type="button"
          className="launcher__button launcher__button--primary"
          onClick={createNewProject}
        >
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

        <div className="launcher__sub">
          {/*
            Recents needs a store of previously opened project directories,
            which nothing writes yet. It is drawn because the design has it,
            and it says so rather than doing nothing when clicked.
          */}
          <button
            type="button"
            className="launcher__sub-button launcher__sub-button--accented"
            onClick={() =>
              dispatch({
                type: 'showNotice',
                message: 'Recently opened projects are not recorded yet. Use Browse.',
              })
            }
          >
            Recents
          </button>

          <button
            type="button"
            className="launcher__sub-button"
            onClick={() => void openExistingProject()}
          >
            Browse
          </button>
        </div>

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
