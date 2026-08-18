import { useState } from 'react';
import { AlignmentPlaceholder } from '../components/alignment/AlignmentPlaceholder';
import { FocalStringEditor } from '../components/focal/FocalStringEditor';
import { HelpButton } from '../components/controls/HelpButton';
import { Checkbox, IconButton, NumberSpinner, TextField } from '../components/controls/Controls';
import { MinusCircleIcon, PlusCircleIcon } from '../components/icons/Icons';
import { useProject } from '../app/state/ProjectContext';
import { MOCK_FASTA_HEADERS } from '../fixtures/mockHeaders';
import { NoticeBar } from '../components/controls/NoticeBar';
import { HELP_TEXT, UNFINISHED_TEXT } from '../copy/helpText';
import { serializeFocalTokens, tokenizeFocalInput } from '../components/focal/focalMatching';
import './MolecularDiagnosisScreen.css';

/**
 * Screen 3 — docs/03-molecular-diagnosis.png.
 *
 * Left: focal-set configuration and analysis parameters.
 * Right: the reserved alignment workspace.
 *
 * MOCK: focal strings are validated against MOCK_FASTA_HEADERS, not against
 * the file chosen on the previous screen. Nothing here parses FASTA.
 */
export function MolecularDiagnosisScreen(): JSX.Element {
  const { state, dispatch, activeFocalSet } = useProject();
  const config = state.draft.molecularDiagnosis;

  // Draft text for the STRING field; committed into the focal set on Enter.
  const [pendingString, setPendingString] = useState('');

  const editorValue = serializeFocalTokens(activeFocalSet.strings);

  const commitEditorValue = (raw: string) => {
    dispatch({
      type: 'setFocalStrings',
      focalSetId: activeFocalSet.id,
      strings: tokenizeFocalInput(raw),
    });
  };

  const appendPendingString = () => {
    const additions = tokenizeFocalInput(pendingString);
    if (additions.length === 0) return;

    // Preserve order, drop duplicates already present.
    const existing = new Set(activeFocalSet.strings);
    const merged = [...activeFocalSet.strings, ...additions.filter((t) => !existing.has(t))];

    dispatch({ type: 'setFocalStrings', focalSetId: activeFocalSet.id, strings: merged });
    setPendingString('');
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
              // TODO(backend): no persistence exists. The title already lives in
              // session state, so this only explains itself. See UI_NOTES Q3.
              onClick: () =>
                dispatch({ type: 'showNotice', message: UNFINISHED_TEXT.saveFocalSet }),
            }}
          />
        </div>

        <div className="diagnosis__row">
          <div className="diagnosis__set-controls">
            <IconButton label="Add focal set" onClick={() => dispatch({ type: 'addFocalSet' })}>
              <PlusCircleIcon size={30} />
            </IconButton>
            <IconButton
              label="Remove focal set"
              onClick={() => dispatch({ type: 'removeActiveFocalSet' })}
              disabled={state.focalSets.length <= 1}
            >
              <MinusCircleIcon size={30} />
            </IconButton>
          </div>

          <label className="diagnosis__label diagnosis__label--inline" htmlFor="focal-string">
            STRING
          </label>

          <TextField
            id="focal-string"
            ariaLabel="Add a focal search string"
            value={pendingString}
            onChange={setPendingString}
            onSubmit={appendPendingString}
            compact
            width="var(--w-field)"
            action={{ label: 'Enter', onClick: appendPendingString }}
          />
        </div>

        <div className="diagnosis__editor">
          <FocalStringEditor
            value={editorValue}
            onChange={commitEditorValue}
            headers={MOCK_FASTA_HEADERS}
            onImport={() =>
              dispatch({ type: 'showNotice', message: UNFINISHED_TEXT.loadFocalStrings })
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
