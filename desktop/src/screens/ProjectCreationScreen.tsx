import { TitleBar } from '../components/chrome/TitleBar';
import { TabStrip } from '../components/chrome/TabStrip';
import { Checkbox, FieldRow, TextField } from '../components/controls/Controls';
import { HelpButton } from '../components/controls/HelpButton';
import { useProject } from '../app/state/ProjectContext';
import { desktop } from '../app/desktopApi';
import { ANALYSIS_LABELS, ANALYSIS_ORDER } from '../app/state/projectState';
import { HELP_TEXT } from '../copy/helpText';
import type { AnalysisKind } from '../contract';
import './ProjectCreationScreen.css';

const HELP_BY_ANALYSIS: Record<AnalysisKind, string> = {
  molecularDiagnosis: HELP_TEXT.molecularDiagnosis,
  sequencePunishmentTest: HELP_TEXT.sequencePunishmentTest,
  consensusSequenceGeneration: HELP_TEXT.consensusSequenceGeneration,
};

/** Show just the file name; the full path goes in the title attribute. */
function baseName(path: string): string {
  const parts = path.split(/[\\/]/);
  return parts[parts.length - 1] || path;
}

/**
 * Screen 2 — docs/02-project-creation.png.
 *
 * A project creation form, not a dashboard: two fields and an analysis
 * selection, on an otherwise empty canvas.
 */
export function ProjectCreationScreen(): JSX.Element {
  const { state, dispatch } = useProject();
  const { draft } = state;

  const chooseFasta = async () => {
    const path = await desktop().dialog.selectFastaFile();
    // A cancelled dialog returns null; leave any previous choice intact.
    if (path) dispatch({ type: 'setFastaPath', path });
  };

  const canContinue = draft.analyses.molecularDiagnosis;

  const goToWorkspace = () => {
    if (!canContinue) return;
    dispatch({ type: 'setActiveAnalysis', analysis: 'molecularDiagnosis' });
    dispatch({ type: 'navigate', screen: 'workspace' });
  };

  return (
    <div className="app-shell creation">
      <TitleBar title="Molecular Diagnosis Tool" variant="workspace" />

      <TabStrip
        homeActive
        onHome={() => {
          desktop().shell.enterLauncherLayout();
          dispatch({ type: 'navigate', screen: 'launcher' });
        }}
        homeLabel="Back to launcher"
      />

      <main className="app-body creation__body">
        <div className="creation__form">
          <FieldRow label="PROJECT TITLE" htmlFor="project-title">
            <TextField
              id="project-title"
              ariaLabel="Project title"
              value={draft.name}
              onChange={(name) => dispatch({ type: 'setProjectName', name })}
              width="var(--w-field-wide)"
            />
          </FieldRow>

          <FieldRow label="SELECT FASTA FILE" htmlFor="fasta-path">
            <TextField
              id="fasta-path"
              ariaLabel="Selected FASTA file"
              value={draft.fastaPath ? baseName(draft.fastaPath) : ''}
              onChange={() => undefined}
              readOnly
              compact
              width="var(--w-field-wide)"
              title={draft.fastaPath ?? undefined}
              action={{ label: 'Browse', onClick: () => void chooseFasta() }}
            />
          </FieldRow>

          <div className="creation__analyses">
            <span className="creation__analyses-label">SELECT ANALYSES</span>

            <ul className="creation__analyses-list">
              {ANALYSIS_ORDER.map((analysis) => (
                <li key={analysis} className="creation__analysis-row">
                  <Checkbox
                    checked={draft.analyses[analysis]}
                    onChange={() => dispatch({ type: 'toggleAnalysis', analysis })}
                    label={ANALYSIS_LABELS[analysis]}
                  />
                  <HelpButton
                    label={`About ${ANALYSIS_LABELS[analysis]}`}
                    text={HELP_BY_ANALYSIS[analysis]}
                  />
                </li>
              ))}
            </ul>
          </div>
        </div>

        {/*
          The reference mockup shows no Continue control, but the flow requires
          one to reach the workspace. Styled in the same language as the other
          action buttons and parked bottom-right. See UI_NOTES §4.
        */}
        <div className="creation__footer">
          <button
            type="button"
            className="creation__continue"
            onClick={goToWorkspace}
            disabled={!canContinue}
            title={
              canContinue
                ? undefined
                : 'Select Molecular Diagnosis to continue — it is the only analysis workspace in this build.'
            }
          >
            Continue
          </button>
        </div>
      </main>
    </div>
  );
}
