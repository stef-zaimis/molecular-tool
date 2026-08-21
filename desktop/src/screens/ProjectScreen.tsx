import { useEffect, useRef, useState } from 'react';
import { Checkbox, TextField } from '../components/controls/Controls';
import { HelpButton } from '../components/controls/HelpButton';
import { NoticeBar } from '../components/controls/NoticeBar';
import { PencilIcon } from '../components/icons/Icons';
import { SourceTable } from '../components/sources/SourceTable';
import { desktop } from '../app/desktopApi';
import { useProject } from '../app/state/ProjectContext';
import {
  ANALYSIS_LABELS,
  ANALYSIS_ORDER,
  IMPLEMENTED_ANALYSES,
} from '../app/state/projectState';
import { HELP_TEXT } from '../copy/helpText';
import type { AnalysisKind } from '../contract';
import './ProjectScreen.css';

const HELP_BY_ANALYSIS: Record<AnalysisKind, string> = {
  molecularDiagnosis: HELP_TEXT.molecularDiagnosis,
  sequencePunishmentTest: HELP_TEXT.sequencePunishmentTest,
  consensusSequenceGeneration: HELP_TEXT.consensusSequenceGeneration,
};

/**
 * The project title once the project exists: a heading with the design's pencil.
 *
 * Editing is a real backend operation (`project.setTitle`), not local state:
 * the title names a project that outlives this session, so a rename has to
 * reach the database or not happen at all.
 */
function SavedTitle({ title }: { title: string }): JSX.Element {
  const { setProjectTitle } = useProject();
  const [editing, setEditing] = useState(false);
  const [value, setValue] = useState(title);
  const inputRef = useRef<HTMLInputElement>(null);

  useEffect(() => {
    if (!editing) setValue(title);
  }, [title, editing]);

  useEffect(() => {
    if (editing) inputRef.current?.focus();
  }, [editing]);

  const commit = async () => {
    const next = value.trim();
    // An unchanged or emptied title is a cancel, not a rename: a project with
    // no name is not a state the database accepts, and saying so would be
    // noise for what is obviously an abandoned edit.
    if (next && next !== title) await setProjectTitle(next);
    setEditing(false);
  };

  if (!editing) {
    return (
      <div className="project__title-row">
        <h2 className="project__title" title={title}>
          {title}
        </h2>
        <button
          type="button"
          className="project__title-edit"
          onClick={() => setEditing(true)}
          aria-label="Rename this project"
          title="Rename this project"
        >
          <PencilIcon size={22} />
        </button>
      </div>
    );
  }

  return (
    <div className="project__title-row">
      <TextField
        inputRef={inputRef}
        ariaLabel="Project title"
        value={value}
        onChange={setValue}
        onSubmit={() => void commit()}
        width="var(--w-field-wide)"
        action={{ label: 'Save', ariaLabel: 'Save project title', onClick: () => void commit() }}
      />
      <button
        type="button"
        className="project__title-cancel"
        onClick={() => {
          setValue(title);
          setEditing(false);
        }}
      >
        Cancel
      </button>
    </div>
  );
}

/**
 * Page 3 — ONE screen for describing a project and for working on one.
 *
 * Before Create it shows a PROJECT TITLE field with Create attached to its
 * right, the same way Browse attaches to the FASTA picker. After Create the
 * page does not navigate anywhere: the title becomes a heading, SELECT ANALYSES
 * appears, and the FASTA table — the SAME component, already populated with the
 * files that were chosen — starts linking directly into the project.
 *
 * That continuity is the point. The previous two-screen flow made creating a
 * project look like leaving the page you had just filled in.
 */
export function ProjectScreen(): JSX.Element {
  const {
    state,
    dispatch,
    createProject,
    linkCandidates,
    chooseFastaFiles,
  } = useProject();

  const { pending } = state;
  const projectOpen = state.project.status === 'open';
  const canCreate = pending.title.trim().length > 0 && !pending.creating;

  const create = async () => {
    if (!canCreate) return;

    // Captured before creating: opening the new project resets `pending`.
    const chosen = pending.candidates;
    // Locks the user set on those rows while they were still candidates. They
    // are re-applied to the real sources as each file is linked, so a row that
    // was locked before Create stays locked after it.
    const lockedPaths = pending.lockedPaths;

    const directory = await desktop().projectDialog.selectDirectory();
    // A cancelled folder picker is not a failure; the form is still standing.
    if (!directory) return;

    dispatch({ type: 'projectCreating', creating: true });
    const created = await createProject(directory, pending.title.trim());
    if (!created) {
      dispatch({ type: 'projectCreating', creating: false });
      return;
    }

    // The project exists from here on. Linking reports its own per-file
    // failures inline and leaves the project intact either way. No navigation:
    // this page simply becomes the saved-project page.
    await linkCandidates(chosen, lockedPaths);
    dispatch({ type: 'projectCreating', creating: false });
  };

  const createBlockedReason = pending.creating
    ? 'Creating the project...'
    : pending.title.trim().length === 0
      ? 'Give the project a title to continue.'
      : undefined;

  return (
    <main className="app-body project__body themed-scroll">
      {projectOpen ? (
        <SavedTitle title={state.project.project.title} />
      ) : (
        <div className="project__row">
          <label className="project__label" htmlFor="project-title">
            PROJECT TITLE
          </label>
          <TextField
            id="project-title"
            ariaLabel="Project title"
            value={pending.title}
            onChange={(title) => dispatch({ type: 'setPendingTitle', title })}
            onSubmit={() => void create()}
            width="var(--w-field-wide)"
            action={{
              label: pending.creating ? 'Creating...' : 'Create',
              ariaLabel: 'Create project',
              onClick: () => void create(),
              disabled: !canCreate,
              title: createBlockedReason ?? 'Choose a folder and create this project',
            }}
          />
          {createBlockedReason && !pending.creating && (
            <span className="project__hint">{createBlockedReason}</span>
          )}
        </div>
      )}

      <div className="project__row">
        <label className="project__label" htmlFor="select-fasta">
          SELECT FASTA FILES
        </label>
        <TextField
          id="select-fasta"
          ariaLabel="Add FASTA files"
          value=""
          onChange={() => undefined}
          readOnly
          compact
          width="var(--w-field-wide)"
          title="Choose one or more aligned FASTA files"
          action={{ label: 'Browse', onClick: () => void chooseFastaFiles() }}
        />
      </div>

      {/*
        Inline, beside the control that produced it. A refusal like "The FASTA
        file needs to be aligned." is about the file just chosen, and belongs
        where the choosing happened rather than in a banner at the bottom.
      */}
      {state.sourceError && (
        <p className="project__source-error" role="alert">
          {state.sourceError}
        </p>
      )}

      <SourceTable />

      {/* Analyses are a property of a project, so they appear once there is one. */}
      {projectOpen && (
        <div className="project__analyses">
          <span className="project__analyses-label">SELECT ANALYSES</span>

          <ul className="project__analyses-list">
            {ANALYSIS_ORDER.map((analysis) => {
              const implemented = IMPLEMENTED_ANALYSES.includes(analysis);
              return (
                <li key={analysis} className="project__analysis-row">
                  <Checkbox
                    checked={state.analyses[analysis]}
                    onChange={() => dispatch({ type: 'toggleAnalysis', analysis })}
                    label={ANALYSIS_LABELS[analysis]}
                  />
                  <HelpButton
                    label={`About ${ANALYSIS_LABELS[analysis]}`}
                    text={HELP_BY_ANALYSIS[analysis]}
                  />
                  {/* Selecting an unbuilt analysis gives it a tab; the tab says
                      what it is rather than opening an empty screen. */}
                  {!implemented && state.analyses[analysis] && (
                    <span className="project__analysis-note">no workspace in this build</span>
                  )}
                </li>
              );
            })}
          </ul>
        </div>
      )}

      {state.notice && (
        <div className="notice--workspace">
          <NoticeBar
            message={state.notice}
            onDismiss={() => dispatch({ type: 'dismissNotice' })}
          />
        </div>
      )}
    </main>
  );
}
