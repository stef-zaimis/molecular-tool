import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
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
  overlappingHeaders,
  saveBlockedReason,
} from '../app/state/focalDrafts';
import type { HeaderListKind } from '../app/state/focalDrafts';
import { useTypingBurst } from '../app/state/typingBurst';
import type { FocalEditMode } from '../app/state/projectState';
import './MolecularDiagnosisScreen.css';

/** Filename-safe stem for the exported focal set. */
function exportFileName(title: string, suffix = ''): string {
  const stem = title.trim().replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '');
  return `${stem || 'focal-set'}${suffix}.txt`;
}

interface HeaderQueryRowProps {
  readonly id: string;
  /** "focal" or "comparison": what an entry is called in labels. */
  readonly entryNoun: string;
  readonly mode: FocalEditMode;
  readonly onModeChange: (mode: FocalEditMode) => void;
  readonly value: string;
  readonly onChange: (value: string) => void;
  readonly onSubmit: () => void;
  readonly busy: boolean;
  readonly locked: boolean;
  readonly help: { readonly label: string; readonly text: string };
}

/**
 * The STRING row: the +/− mode pair, the search string and its Enter.
 *
 * One component for both the focal and the comparison list, so `+` and `−`
 * cannot come to mean different things in the two places. The focal row's
 * labels are exactly what they were before the comparison row existed.
 */
function HeaderQueryRow({
  id,
  entryNoun,
  mode,
  onModeChange,
  value,
  onChange,
  onSubmit,
  busy,
  locked,
  help,
}: HeaderQueryRowProps): JSX.Element {
  const Noun = entryNoun.charAt(0).toUpperCase() + entryNoun.slice(1);
  const isFocal = entryNoun === 'focal';
  return (
    <div className="diagnosis__row diagnosis__row--string">
      <div
        className="diagnosis__set-controls"
        role="radiogroup"
        aria-label={`${Noun} entry edit mode`}
      >
        <button
          type="button"
          role="radio"
          aria-checked={mode === 'add'}
          className={`mode-button${mode === 'add' ? ' is-selected' : ''}`}
          title={`Add mode: Enter adds every matching FASTA header to the ${entryNoun} set`}
          onClick={() => onModeChange('add')}
        >
          <PlusCircleIcon size={30} />
          <span className="sr-only">{isFocal ? 'Add mode' : 'Add to comparison set'}</span>
        </button>
        <button
          type="button"
          role="radio"
          aria-checked={mode === 'remove'}
          className={`mode-button${mode === 'remove' ? ' is-selected' : ''}`}
          title={`Remove mode: Enter removes every ${entryNoun} entry containing the string`}
          onClick={() => onModeChange('remove')}
        >
          <MinusCircleIcon size={30} />
          <span className="sr-only">{isFocal ? 'Remove mode' : 'Remove from comparison set'}</span>
        </button>
      </div>

      <label className="diagnosis__label" htmlFor={id}>
        STRING
      </label>

      <TextField
        id={id}
        ariaLabel={
          isFocal
            ? mode === 'add'
              ? 'Search string: add every matching header'
              : 'Search string: remove every matching focal entry'
            : mode === 'add'
              ? 'Comparison string: add matching headers to the comparison set'
              : 'Comparison string: remove matching comparison entries'
        }
        value={value}
        readOnly={locked}
        onChange={onChange}
        onSubmit={onSubmit}
        compact
        width="var(--w-field)"
        action={{
          label: 'Enter',
          ariaLabel: isFocal
            ? mode === 'add'
              ? 'Apply add'
              : 'Apply remove'
            : mode === 'add'
              ? 'Apply comparison add'
              : 'Apply comparison remove',
          onClick: onSubmit,
          disabled: busy || locked,
        }}
      />
      <HelpButton label={help.label} text={help.text} />
    </div>
  );
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
    stopDiagnosis,
    continuationCurrent,
  } = useProject();

  const config = state.molecularDiagnosis;
  const mode = state.focalMode;
  const scope = state.fastaScope;

  const [pendingString, setPendingString] = useState('');
  /* The comparison row keeps its own string and mode: two rows, two queries. */
  const [comparisonString, setComparisonString] = useState('');
  const [comparisonMode, setComparisonMode] = useState<FocalEditMode>('add');
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

  /* The comparison list has its own undo stack, so its own burst timer. */
  const comparisonBurst = useTypingBurst<string>(
    useCallback(
      (key: string) => dispatch({ type: 'endDraftBurst', key, list: 'comparison' }),
      [dispatch],
    ),
  );
  useEffect(() => comparisonBurst.cancel, [activeDraft.key, comparisonBurst.cancel]);

  /*
   * Headers in BOTH lists. Both editors paint them, and a persistent line says
   * why Run is waiting. Nothing is removed from either side automatically.
   */
  const overlapKey = overlappingHeaders(activeDraft).join('\n');
  const conflicts = useMemo(
    () => new Set(overlapKey ? overlapKey.split('\n') : []),
    [overlapKey],
  );

  const save = async () => {
    burst.endNow(activeDraft.key);
    comparisonBurst.endNow(activeDraft.key);
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

  /** Enter applies the row's current mode to its list, through the backend. */
  const applyQuery = async (list: HeaderListKind) => {
    const isFocal = list === 'focal';
    const value = (isFocal ? pendingString : comparisonString).trim();
    const rowMode = isFocal ? mode : comparisonMode;
    if (!value || busy) return;

    // `+`/`−` are one deliberate action and one undo step each, so whatever was
    // being typed is a finished step before either runs.
    (isFocal ? burst : comparisonBurst).endNow(activeDraft.key);

    setBusy(true);
    const applied =
      rowMode === 'add' ? await addByQuery(value, list) : await removeByQuery(value, list);
    setBusy(false);
    // Cleared only when something happened, so a query that matched nothing
    // can be corrected rather than retyped.
    if (applied) {
      (isFocal ? setPendingString : setComparisonString)('');
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

  const onComparisonChange = useCallback(
    (headers: readonly string[], typing: boolean) => {
      dispatch({
        type: 'setDraftHeaders',
        key: activeDraft.key,
        headers,
        coalesce: typing,
        list: 'comparison',
      });
      if (typing) comparisonBurst.touch(activeDraft.key);
      else comparisonBurst.cancel();
    },
    [dispatch, activeDraft.key, comparisonBurst],
  );

  const onComparisonBlur = useCallback(() => {
    comparisonBurst.endNow(activeDraft.key);
  }, [comparisonBurst, activeDraft.key]);

  const exportFocalSet = async (list: HeaderListKind = 'focal') => {
    const result = await desktop().dialog.exportFocalSet({
      suggestedName: exportFileName(activeDraft.title, list === 'focal' ? '' : '-comparison'),
      lines: list === 'focal' ? activeDraft.headers : activeDraft.comparison,
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

        <HeaderQueryRow
          id="focal-string"
          entryNoun="focal"
          mode={mode}
          onModeChange={setMode}
          value={pendingString}
          onChange={setPendingString}
          onSubmit={() => void applyQuery('focal')}
          busy={busy}
          locked={locked}
          help={{ label: 'About the focal search string', text: HELP_TEXT.focalString }}
        />

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
            onExport={() => void exportFocalSet('focal')}
            canExport={activeDraft.headers.length > 0}
            conflicts={conflicts}
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

        {/*
          COMPARISON SET — optional, saved WITH the focal set. Blank means the
          default every run has always had: compare against every non-focal
          specimen in the FASTA pool. The same editor, the same +/−, the same
          presence colours; headers also in the focal set are painted as
          conflicts in both boxes.
        */}
        <section className="diagnosis__comparison" aria-labelledby="comparison-set-heading">
          <div className="diagnosis__comparison-heading">
            <span
              id="comparison-set-heading"
              className="diagnosis__label diagnosis__label--section"
            >
              COMPARISON SET
            </span>
            <span className="diagnosis__note">
              Leave blank to compare against all non-focal specimens in the FASTA pool.
            </span>
            <HelpButton label="About the comparison set" text={HELP_TEXT.comparisonSet} />
          </div>

          <HeaderQueryRow
            id="comparison-string"
            entryNoun="comparison"
            mode={comparisonMode}
            onModeChange={setComparisonMode}
            value={comparisonString}
            onChange={setComparisonString}
            onSubmit={() => void applyQuery('comparison')}
            busy={busy}
            locked={locked}
            help={{ label: 'About the comparison search string', text: HELP_TEXT.focalString }}
          />

          <div className="diagnosis__editor diagnosis__editor--comparison">
            <FocalSetEditor
              headers={activeDraft.comparison}
              presence={state.focalPresence.byHeader}
              readOnly={locked}
              onChange={onComparisonChange}
              onEditingEnd={onComparisonBlur}
              canonicalToken={canonicalToken}
              onUndo={() => {
                comparisonBurst.endNow(activeDraft.key);
                dispatch({ type: 'undoDraftEdit', key: activeDraft.key, list: 'comparison' });
              }}
              onRedo={() => {
                comparisonBurst.endNow(activeDraft.key);
                dispatch({ type: 'redoDraftEdit', key: activeDraft.key, list: 'comparison' });
              }}
              canUndo={canUndo(activeDraft, 'comparison') && !locked}
              canRedo={canRedo(activeDraft, 'comparison') && !locked}
              onExport={() => void exportFocalSet('comparison')}
              canExport={activeDraft.comparison.length > 0}
              conflicts={conflicts}
              noun="comparison set"
              entryNoun="comparison"
              placeholderText="Optional: type headers separated by ; or use + above"
            />
          </div>

          {/*
            Persistent, not a toast: it stays exactly as long as the problem
            does, and says what to do about it.
          */}
          {conflicts.size > 0 && (
            <p className="diagnosis__conflict" role="status">
              {conflicts.size === 1
                ? '1 specimen is in both the focal and the comparison set (highlighted). '
                : `${conflicts.size} specimens are in both the focal and the comparison set (highlighted). `}
              A specimen cannot belong to both. Remove it from one set to run.
            </p>
          )}
        </section>

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
          onStop={() => void stopDiagnosis()}
          onContinue={(resume) => void runDiagnosis(resume)}
          continuationCurrent={continuationCurrent}
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
