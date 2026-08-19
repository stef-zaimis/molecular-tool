import { useState } from 'react';
import { AlignmentPlaceholder } from '../components/alignment/AlignmentPlaceholder';
import { FocalStringEditor } from '../components/focal/FocalStringEditor';
import { HelpButton } from '../components/controls/HelpButton';
import { Checkbox, NumberSpinner, TextField } from '../components/controls/Controls';
import { MinusCircleIcon, PlusCircleIcon } from '../components/icons/Icons';
import { NoticeBar } from '../components/controls/NoticeBar';
import { useProject } from '../app/state/ProjectContext';
import { desktop } from '../app/desktopApi';
import { HELP_TEXT, UNFINISHED_TEXT } from '../copy/helpText';
import {
  describeFocalStringProblem,
  validateFocalString,
} from '../components/focal/focalMatching';
import type { FocalEditMode } from '../app/state/projectState';
import './MolecularDiagnosisScreen.css';

/** Filename-safe stem for the exported focal set. */
function exportFileName(title: string): string {
  const stem = title.trim().replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '');
  return `${stem || 'focal-set'}.txt`;
}

/**
 * Screen 3 — the Molecular Diagnosis workspace.
 *
 * Focal-set editing is modal: `+` and `−` are mutually exclusive, and Enter
 * applies the selected mode to whatever is in the STRING field. Every applied
 * edit is one undo step.
 */
export function MolecularDiagnosisScreen(): JSX.Element {
  const { state, dispatch, activeFocalSet, activeFocalHistory, fastaHeaders } = useProject();
  const config = state.draft.molecularDiagnosis;
  const mode = state.focalMode;

  const [pendingString, setPendingString] = useState('');
  const [inlineError, setInlineError] = useState<string | null>(null);

  const setMode = (next: FocalEditMode) => {
    dispatch({ type: 'setFocalMode', mode: next });
    setInlineError(null);
  };

  /** Enter applies the current mode. */
  const applyPendingString = () => {
    const value = pendingString.trim();

    // ';' is rejected in both modes: it is the display separator, so allowing
    // it inside an entry would make the rendered set ambiguous.
    const problem = validateFocalString(
      value,
      mode === 'add' ? activeFocalSet.strings : [],
    );
    if (problem) {
      setInlineError(describeFocalStringProblem(problem));
      return;
    }

    if (mode === 'add') {
      dispatch({ type: 'addFocalString', focalSetId: activeFocalSet.id, value });
      setPendingString('');
      setInlineError(null);
      return;
    }

    // Remove mode: only an exact entry is removed, and a miss is not destructive.
    if (!activeFocalSet.strings.includes(value)) {
      setInlineError('No focal string exactly matches that text.');
      return;
    }

    dispatch({ type: 'removeFocalString', focalSetId: activeFocalSet.id, value });
    setPendingString('');
    setInlineError(null);
  };

  const exportFocalSet = async () => {
    const result = await desktop().dialog.exportFocalSet({
      suggestedName: exportFileName(activeFocalSet.title),
      lines: activeFocalSet.strings,
    });

    if (!result.ok && result.code !== 'CANCELLED') {
      dispatch({
        type: 'showNotice',
        message: result.message ?? 'The focal set could not be exported.',
      });
    }
  };

  return (
    <div className="diagnosis">
      <div className="diagnosis__left">
        <div className="diagnosis__row">
          <label className="diagnosis__label" htmlFor="focal-set-title">
            FOCAL SET TITLE
          </label>
          <TextField
            id="focal-set-title"
            ariaLabel="Focal set title"
            value={activeFocalSet.title}
            onChange={(title) =>
              dispatch({ type: 'setFocalTitle', focalSetId: activeFocalSet.id, title })
            }
            width="var(--w-field)"
            action={{
              label: 'Save',
              // TODO(backend): no persistence exists; the title lives in session
              // state only. See UI_NOTES.
              onClick: () =>
                dispatch({ type: 'showNotice', message: UNFINISHED_TEXT.saveFocalSet }),
            }}
          />
        </div>

        <div className="diagnosis__row">
          <div className="diagnosis__set-controls" role="radiogroup" aria-label="Focal string edit mode">
            <button
              type="button"
              role="radio"
              aria-checked={mode === 'add'}
              className={`mode-button${mode === 'add' ? ' is-selected' : ''}`}
              title="Add mode: Enter adds the string to the focal set"
              onClick={() => setMode('add')}
            >
              <PlusCircleIcon size={30} />
              <span className="sr-only">Add mode</span>
            </button>
            <button
              type="button"
              role="radio"
              aria-checked={mode === 'remove'}
              className={`mode-button${mode === 'remove' ? ' is-selected' : ''}`}
              title="Remove mode: Enter removes the matching string from the focal set"
              onClick={() => setMode('remove')}
            >
              <MinusCircleIcon size={30} />
              <span className="sr-only">Remove mode</span>
            </button>
          </div>

          <label className="diagnosis__label" htmlFor="focal-string">
            STRING
          </label>

          <TextField
            id="focal-string"
            ariaLabel={
              mode === 'add'
                ? 'Focal search string to add'
                : 'Focal search string to remove'
            }
            value={pendingString}
            onChange={(value) => {
              setPendingString(value);
              if (inlineError) setInlineError(null);
            }}
            onSubmit={applyPendingString}
            compact
            width="var(--w-field)"
            action={{ label: 'Enter', onClick: applyPendingString }}
          />

          {inlineError && (
            <p className="diagnosis__inline-error" role="alert">
              {inlineError}
            </p>
          )}
        </div>

        <div className="diagnosis__editor">
          <FocalStringEditor
            strings={activeFocalSet.strings}
            headers={fastaHeaders}
            onUndo={() => dispatch({ type: 'undoFocalEdit', focalSetId: activeFocalSet.id })}
            onRedo={() => dispatch({ type: 'redoFocalEdit', focalSetId: activeFocalSet.id })}
            canUndo={activeFocalHistory.past.length > 0}
            canRedo={activeFocalHistory.future.length > 0}
            onExport={() => void exportFocalSet()}
            canExport={activeFocalSet.strings.length > 0}
            emptyHint={
              mode === 'add'
                ? 'No focal strings yet. Type one above and press Enter.'
                : 'No focal strings to remove.'
            }
          />
        </div>

        <div className="diagnosis__parameters">
          <span className="diagnosis__label diagnosis__label--section">SELECT PARAMETERS</span>

          <div className="diagnosis__parameter-list">
            <div className="diagnosis__parameter-row">
              <Checkbox
                checked={config.ignoreGaps}
                onChange={(ignoreGaps) =>
                  dispatch({ type: 'updateMolecularDiagnosis', patch: { ignoreGaps } })
                }
                label="Ignore gaps"
              />
              <HelpButton label="About Ignore gaps" text={HELP_TEXT.ignoreGaps} />
            </div>

            <div className="diagnosis__parameter-row">
              <Checkbox
                checked={config.giveBenefitOfDoubtToAmbiguousBases}
                onChange={(value) =>
                  dispatch({
                    type: 'updateMolecularDiagnosis',
                    patch: { giveBenefitOfDoubtToAmbiguousBases: value },
                  })
                }
                label="Give BoTD to ambiguous bases"
              />
              <HelpButton
                label="About benefit of the doubt for ambiguous bases"
                text={HELP_TEXT.benefitOfDoubt}
              />
            </div>

            <div className="diagnosis__parameter-row diagnosis__parameter-row--numeric">
              <span className="diagnosis__parameter-name">Min. candidate DNC size</span>
              <NumberSpinner
                ariaLabel="Minimum candidate DNC size"
                value={config.minCandidateSize}
                min={1}
                max={20}
                onChange={(minCandidateSize) =>
                  dispatch({
                    type: 'updateMolecularDiagnosis',
                    patch: {
                      minCandidateSize,
                      maxCandidateSize: Math.max(minCandidateSize, config.maxCandidateSize),
                    },
                  })
                }
              />
              <HelpButton label="About minimum candidate size" text={HELP_TEXT.minCandidateSize} />
            </div>

            <div className="diagnosis__parameter-row diagnosis__parameter-row--numeric">
              <span className="diagnosis__parameter-name">Max. candidate DNC size</span>
              <NumberSpinner
                ariaLabel="Maximum candidate DNC size"
                value={config.maxCandidateSize}
                min={1}
                max={20}
                onChange={(maxCandidateSize) =>
                  dispatch({
                    type: 'updateMolecularDiagnosis',
                    patch: {
                      maxCandidateSize,
                      minCandidateSize: Math.min(maxCandidateSize, config.minCandidateSize),
                    },
                  })
                }
              />
              <HelpButton label="About maximum candidate size" text={HELP_TEXT.maxCandidateSize} />
            </div>
          </div>
        </div>
      </div>

      {state.notice && (
        <div className="notice--workspace">
          <NoticeBar
            message={state.notice}
            onDismiss={() => dispatch({ type: 'dismissNotice' })}
          />
        </div>
      )}

      <AlignmentPlaceholder
        caption={
          state.draft.fastaPath
            ? state.draft.fastaPath.split(/[\\/]/).pop()
            : 'No alignment loaded'
        }
      />
    </div>
  );
}
