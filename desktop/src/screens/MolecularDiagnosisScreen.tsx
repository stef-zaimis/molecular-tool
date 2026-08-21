import { useCallback, useEffect, useRef, useState } from 'react';
import { FocalSetEditor } from '../components/focal/FocalSetEditor';
import { WorkspaceRightPane } from '../components/alignment/WorkspaceRightPane';
import { HelpButton } from '../components/controls/HelpButton';
import { Checkbox, NumberSpinner, TextField } from '../components/controls/Controls';
import { MinusCircleIcon, PlusCircleIcon } from '../components/icons/Icons';
import { NoticeBar } from '../components/controls/NoticeBar';
import { DiagnosisRunPanel } from '../components/analysis/DiagnosisRunPanel';
import { useProject } from '../app/state/ProjectContext';
import { desktop } from '../app/desktopApi';
import { HELP_TEXT } from '../copy/helpText';
import {
  canRedo,
  canUndo,
  isDirty,
  isNewDraft,
  saveBlockedReason,
} from '../app/state/focalDrafts';
import { useTypingBurst } from '../app/state/typingBurst';
import type { FocalEditMode } from '../app/state/projectState';
import './MolecularDiagnosisScreen.css';

/** Filename-safe stem for the exported focal set. */
function exportFileName(title: string): string {
  const stem = title.trim().replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '');
  return `${stem || 'focal-set'}.txt`;
}

/**
 * The Molecular Diagnosis workspace — screens 4 and 6 of the new design.
 *
 * Screen 4 is a never-saved set: the title is an input with Enter, and the
 * right pane is the visualizer alone. Screen 6 is a saved one: the title
 * becomes a heading with SAVE beside it, and the focal-set library can occupy
 * the area above the visualizer. Both are this component; which one you see is
 * a property of the draft, not a separate screen.
 *
 * The chain underneath is unchanged from the pass that landed it:
 *
 *   linked sources -> selected scope -> working draft -> saved set -> run
 *
 * Editing writes nothing. Save is the only focal write, and Run refuses a
 * draft that is new or dirty rather than quietly saving it.
 */
export function MolecularDiagnosisScreen(): JSX.Element {
  const {
    state,
    dispatch,
    activeDraft,
    saveActiveDraft,
    addByQuery,
    removeByQuery,
    runGate,
    runDiagnosis,
  } = useProject();

  const config = state.molecularDiagnosis;
  const mode = state.focalMode;
  const scope = state.fastaScope;

  const [pendingString, setPendingString] = useState('');
  const [busy, setBusy] = useState(false);
  const [viewerDragging, setViewerDragging] = useState(false);
  const titleInputRef = useRef<HTMLInputElement>(null);

  /*
   * Bumped at each CANONICALISATION BOUNDARY, which is what tells the editor to
   * square its text up with the stored membership: a save, a `+`/`−`, or moving
   * to a different focal set. Between boundaries the text is left exactly as
   * typed, so a half-finished entry survives the next keystroke.
   */
  const [canonicalToken, setCanonicalToken] = useState(0);
  const canonicalise = useCallback(() => setCanonicalToken((token) => token + 1), []);

  useEffect(canonicalise, [activeDraft.key, canonicalise]);

  /*
   * The idle timer that ends a typing burst.
   *
   * It closes the undo GROUP and nothing else — no canonicalisation — so a
   * pause in the middle of typing a header never rewrites what is on screen.
   * The boundaries below (blur, `+`/`−`, Save, undo/redo, switching sets) close
   * it explicitly as well; the timer only covers the case where the user keeps
   * the field focused and simply stops.
   */
  const burst = useTypingBurst<string>(
    useCallback((key: string) => dispatch({ type: 'endDraftBurst', key }), [dispatch]),
  );

  // A different focal set is a different history. Nothing pending applies to
  // it, and the reducer has already closed the burst on the one being left.
  useEffect(() => burst.cancel, [activeDraft.key, burst.cancel]);

  const save = async () => {
    burst.endNow(activeDraft.key);
    const saved = await saveActiveDraft();
    // After a save the field must show exactly what was stored — no struck-out
    // duplicate left behind for an entry the database does not hold.
    if (saved) canonicalise();
    return saved;
  };

  const saveBlocked = saveBlockedReason(activeDraft);
  const dirty = isDirty(activeDraft);
  const isNew = isNewDraft(activeDraft);
  const locked = activeDraft.locked;

  /* The library's pencil selects a draft and asks for its title; honour that. */
  const renaming = state.titleEditingKey === activeDraft.key;
  useEffect(() => {
    if (renaming) titleInputRef.current?.focus();
  }, [renaming]);

  const setMode = (next: FocalEditMode) => dispatch({ type: 'setFocalMode', mode: next });

  /** Enter applies the current mode, through the backend. */
  const applyPendingString = async () => {
    const value = pendingString.trim();
    if (!value || busy) return;

    // `+`/`−` are one deliberate action and one undo step each, so whatever was
    // being typed is a finished step before either runs.
    burst.endNow(activeDraft.key);

    setBusy(true);
    const applied = mode === 'add' ? await addByQuery(value) : await removeByQuery(value);
    setBusy(false);
    // Cleared only when something happened, so a query that matched nothing
    // can be corrected rather than retyped.
    if (applied) {
      setPendingString('');
      // A discrete `+`/`−` is a safe moment to square the text up.
      canonicalise();
    }
  };

  /*
   * Manual edits arrive here as complete membership.
   *
   * `typing` marks a continuous burst, which the reducer folds into ONE undo
   * step — undoing a typed header one character at a time would be useless.
   * Nothing here expands a substring: this field holds exact headers, and `+`
   * is the thing that searches.
   */
  const onEditorChange = useCallback(
    (headers: readonly string[], typing: boolean) => {
      dispatch({ type: 'setDraftHeaders', key: activeDraft.key, headers, coalesce: typing });
      // Every document change restarts the countdown; a change that is not
      // typing closes the burst in the reducer, so no timer should survive it.
      if (typing) burst.touch(activeDraft.key);
      else burst.cancel();
    },
    [dispatch, activeDraft.key, burst],
  );

  const onEditorBlur = useCallback(() => {
    burst.endNow(activeDraft.key);
  }, [burst, activeDraft.key]);

  const exportFocalSet = async () => {
    const result = await desktop().dialog.exportFocalSet({
      suggestedName: exportFileName(activeDraft.title),
      lines: activeDraft.headers,
    });

    if (!result.ok && result.code !== 'CANCELLED') {
      dispatch({
        type: 'showNotice',
        message: result.message ?? 'The focal set could not be exported.',
      });
    }
  };

  const scopeCaption =
    scope.kind === 'all'
      ? `All files (${state.sources.items.length})`
      : (state.sources.items.find((item) => item.fastaFileId === scope.fastaFileId)
          ?.displayName ?? 'No FASTA selected');

  return (
    <div className={`diagnosis${viewerDragging ? ' is-viewer-dragging' : ''}`}>
      <div className="diagnosis__left themed-scroll">
        {/*
          Screen 4 vs screen 6: an unsaved set gets the labelled input with
          Enter; once it has been saved its name is a heading with SAVE beside
          it. Renaming from the library re-opens the input form.
        */}
        {isNew || renaming ? (
          <div className="diagnosis__row">
            <label className="diagnosis__label" htmlFor="focal-set-title">
              FOCAL SET TITLE
            </label>
            <TextField
              id="focal-set-title"
              inputRef={titleInputRef}
              ariaLabel="Focal set title"
              value={activeDraft.title}
              readOnly={locked}
              onChange={(title) =>
                dispatch({ type: 'setDraftTitle', key: activeDraft.key, title })
              }
              onSubmit={() => void save()}
              width="var(--w-field)"
              action={{
                label: 'Enter',
                ariaLabel: 'Save focal set',
                onClick: () => void save(),
                disabled: saveBlocked !== null,
                title: saveBlocked ?? 'Save this focal set to the project',
              }}
            />
            {/*
              No help dot here. "FOCAL SET TITLE" is a name for a focal set;
              a tooltip explaining that would be noise beside the one control
              on this page that genuinely needs explaining.
            */}
          </div>
        ) : (
          <div className="diagnosis__title-row">
            <h2 className="diagnosis__title" title={activeDraft.title}>
              {activeDraft.title}
            </h2>
            <button
              type="button"
              className="diagnosis__save"
              onClick={() => void save()}
              disabled={saveBlocked !== null}
              title={saveBlocked ?? 'Save this focal set to the project'}
            >
              SAVE
            </button>
            {/*
              Stated, not implied by an enabled button: "is what I am looking at
              what a run would use?" must never be a guess.
            */}
            <span className={`diagnosis__save-state${dirty ? ' is-dirty' : ''}`} role="status">
              {locked ? 'Locked' : dirty ? 'Unsaved changes' : 'Saved'}
            </span>
          </div>
        )}

        <div className="diagnosis__row diagnosis__row--string">
          <div
            className="diagnosis__set-controls"
            role="radiogroup"
            aria-label="Focal entry edit mode"
          >
            <button
              type="button"
              role="radio"
              aria-checked={mode === 'add'}
              className={`mode-button${mode === 'add' ? ' is-selected' : ''}`}
              title="Add mode: Enter adds every matching FASTA header to the focal set"
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
              title="Remove mode: Enter removes every focal entry containing the string"
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
                ? 'Search string: add every matching header'
                : 'Search string: remove every matching focal entry'
            }
            value={pendingString}
            readOnly={locked}
            onChange={setPendingString}
            onSubmit={() => void applyPendingString()}
            compact
            width="var(--w-field)"
            action={{
              label: 'Enter',
              ariaLabel: mode === 'add' ? 'Apply add' : 'Apply remove',
              onClick: () => void applyPendingString(),
              disabled: busy || locked,
            }}
          />
          <HelpButton label="About the focal search string" text={HELP_TEXT.focalString} />
        </div>

        {/*
          FASTA POOL — the analysis and search scope. It decides which files a
          run reads AND which files `+` searches, so the two can never describe
          different things. Presence colours compare against every linked file
          regardless, which is what makes orange possible.
        */}
        <div className="diagnosis__row">
          <label className="diagnosis__label" htmlFor="fasta-pool">
            FASTA POOL
          </label>
          <div className="diagnosis__pool-wrap">
            <select
              id="fasta-pool"
              className="diagnosis__pool"
              value={scope.kind === 'all' ? 'all' : scope.fastaFileId}
              onChange={(event) =>
                dispatch({
                  type: 'setFastaScope',
                  scope:
                    event.target.value === 'all'
                      ? { kind: 'all' }
                      : { kind: 'file', fastaFileId: event.target.value },
                })
              }
            >
              <option value="all">All files</option>
              {state.sources.items.map((source) => (
                <option key={source.fastaFileId} value={source.fastaFileId}>
                  {source.displayName}
                </option>
              ))}
            </select>
          </div>
          <HelpButton label="About the FASTA pool" text={HELP_TEXT.fastaPool} />
        </div>

        <div className="diagnosis__editor">
          <FocalSetEditor
            headers={activeDraft.headers}
            presence={state.focalPresence.byHeader}
            readOnly={locked}
            onChange={onEditorChange}
            onEditingEnd={onEditorBlur}
            canonicalToken={canonicalToken}
            /* Undo/redo are boundaries: the burst they would walk into ends first. */
            onUndo={() => {
              burst.endNow(activeDraft.key);
              dispatch({ type: 'undoDraftEdit', key: activeDraft.key });
            }}
            onRedo={() => {
              burst.endNow(activeDraft.key);
              dispatch({ type: 'redoDraftEdit', key: activeDraft.key });
            }}
            canUndo={canUndo(activeDraft) && !locked}
            canRedo={canRedo(activeDraft) && !locked}
            onExport={() => void exportFocalSet()}
            canExport={activeDraft.headers.length > 0}
          />

          {/*
            The Search row beneath the box, drawn in the design but not yet
            wired to anything. Disabled rather than fake: a field that looks
            live and does nothing is worse than one that says it is not ready.
          */}
          <div className="diagnosis__search" aria-hidden="true">
            <input
              className="diagnosis__search-input"
              type="text"
              placeholder="Search"
              disabled
              tabIndex={-1}
            />
            <span className="diagnosis__search-enter">ENTER</span>
          </div>
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

            {/*
              In the design but not computed in this build. An em dash says
              "not available"; a number would be a claim.
            */}
            <div className="diagnosis__parameter-row diagnosis__parameter-row--numeric">
              <span className="diagnosis__parameter-name">Estimated run time</span>
              <span className="diagnosis__estimate">—</span>
            </div>
          </div>
        </div>

        <DiagnosisRunPanel
          run={state.diagnosisRun}
          gate={runGate}
          onRun={() => void runDiagnosis(null)}
          onContinue={(resume) => void runDiagnosis(resume)}
          onDismissContinuation={() => dispatch({ type: 'dismissContinuation' })}
          onReveal={(filePath) => void desktop().shell.showItemInFolder(filePath)}
        />
      </div>

      {state.notice && (
        <div className="notice--workspace">
          <NoticeBar
            message={state.notice}
            onDismiss={() => dispatch({ type: 'dismissNotice' })}
          />
        </div>
      )}

      <WorkspaceRightPane caption={scopeCaption} onDraggingChange={setViewerDragging} />
    </div>
  );
}
