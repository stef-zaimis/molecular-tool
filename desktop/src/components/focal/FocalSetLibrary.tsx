import { useEffect, useRef, useState } from 'react';
import { useProject } from '../../app/state/ProjectContext';
import { isDirty, isNewDraft } from '../../app/state/focalDrafts';
import type { FocalDraft } from '../../app/state/focalDrafts';
import { LockClosedIcon, LockOpenIcon, PencilIcon, TrashIcon } from '../icons/Icons';
import './FocalSetLibrary.css';

/**
 * The FOCAL SET LIBRARY panel — screen 6 of the design.
 *
 * It lists WORKING DRAFTS, not database rows. A set the user started and has
 * not saved is as real to them as one that is stored, so it appears here with
 * an asterisk rather than being invisible until it happens to be committed.
 *
 * Selecting a row never saves the row being left. That is the same rule the
 * rest of the workspace follows: leaving a draft is not a decision to commit
 * it, and an implicit save here would silently write a set the user was still
 * deciding about.
 */

interface RowProps {
  readonly draft: FocalDraft;
  readonly active: boolean;
}

function LibraryRow({ draft, active }: RowProps): JSX.Element {
  const { dispatch, setDraftLocked, deleteDraft } = useProject();
  const [confirmingDelete, setConfirmingDelete] = useState(false);
  const confirmRef = useRef<HTMLButtonElement>(null);

  useEffect(() => {
    if (confirmingDelete) confirmRef.current?.focus();
  }, [confirmingDelete]);

  const dirty = isDirty(draft);
  /*
   * The active draft's unsaved state is already shown by the Save control next
   * to its title, so the asterisk here would say the same thing twice. It
   * appears once the user moves on and that control is describing a different
   * set.
   */
  const showAsterisk = dirty && !active;

  /*
   * Clicking the row that is ALREADY active is a strict no-op.
   *
   * It used to re-dispatch a selection, which cleared presence, reset the
   * rename state and made the token colours flicker through neutral on the way
   * back to the same answer. Guarded here as well as in the reducer, so
   * neither path can do it.
   */
  const select = () => {
    if (active) return;
    dispatch({ type: 'selectFocalDraft', key: draft.key });
  };

  const rename = () => {
    // Renaming is an edit to a working copy, so the row is selected first and
    // the title field takes focus. Save commits it like any other change.
    select();
    dispatch({ type: 'setTitleEditing', key: draft.key });
  };

  const deleteMessage = isNewDraft(draft)
    ? `Delete "${draft.title.trim() || 'Untitled'}"? It has never been saved, so it exists only here.`
    : dirty
      ? `Delete "${draft.savedTitle}" from the project? Its unsaved changes will be discarded as well.`
      : `Delete "${draft.savedTitle}" from the project?`;

  return (
    <li
      className={`focal-row${active ? ' is-active' : ''}${draft.locked ? ' is-locked' : ''}`}
    >
      <span className="focal-row__edge" aria-hidden="true" />

      {/*
        Selection button and pencil are SIBLINGS inside one inline group, never
        nested: a button inside a button is invalid and unpredictable for
        assistive tech. The group is what sizes to the title, so the pencil
        begins where the rendered title ends rather than out at the [n=x]
        column, and a long title truncates with the pencil still against it.
      */}
      <span className="focal-row__name-group">
        <button
          type="button"
          className="focal-row__select"
          onClick={select}
          aria-current={active ? 'true' : undefined}
          // The row's own controls repeat the title in their names; this one says
          // what activating the row itself does.
          aria-label={`Select ${draft.title.trim() || 'untitled focal set'}`}
        >
          <span className="focal-row__marker" aria-hidden="true">
            {active ? '▶' : ''}
          </span>
          <span className="focal-row__title">
            {showAsterisk && (
              <span className="focal-row__unsaved" title="Unsaved changes">
                *
              </span>
            )}
            {draft.title.trim() || <span className="focal-row__untitled">Untitled</span>}
          </span>
        </button>

        {!draft.locked && (
          <button
            type="button"
            className="focal-row__icon focal-row__icon--rename"
            onClick={rename}
            title="Rename this focal set"
            aria-label={`Rename ${draft.title.trim() || 'this focal set'}`}
          >
            <PencilIcon size={18} />
          </button>
        )}
      </span>

      <span className="focal-row__count">[n={draft.headers.length}]</span>

      <span className="focal-row__actions">
        {confirmingDelete ? (
          /*
           * In-row, replacing the action area. The previous panel expanded
           * below the row, which pushed every later row down and could put the
           * buttons out of view in a scrolled list.
           */
          <span
            className="focal-row__confirm"
            role="alertdialog"
            // The row has no space for the sentence, so it IS the accessible
            // name and the tooltip rather than being dropped.
            aria-label={deleteMessage}
            title={deleteMessage}
            onKeyDown={(event) => {
              if (event.key === 'Escape') {
                event.stopPropagation();
                setConfirmingDelete(false);
              }
            }}
          >
            <button
              ref={confirmRef}
              type="button"
              className="focal-row__confirm-action focal-row__confirm-action--danger"
              onClick={() => {
                setConfirmingDelete(false);
                void deleteDraft(draft.key);
              }}
            >
              Delete
            </button>
            <button
              type="button"
              className="focal-row__confirm-action"
              onClick={() => setConfirmingDelete(false)}
            >
              Cancel
            </button>
          </span>
        ) : (
          <>
            <button
              type="button"
              className="focal-row__icon focal-row__icon--lock"
              onClick={() => void setDraftLocked(draft.key, !draft.locked)}
              title={draft.locked ? 'Unlock this focal set' : 'Lock this focal set'}
              aria-label={`${draft.locked ? 'Unlock' : 'Lock'} ${
                draft.title.trim() || 'this focal set'
              }`}
              aria-pressed={draft.locked}
            >
              {draft.locked ? <LockClosedIcon size={18} /> : <LockOpenIcon size={18} />}
            </button>

            {/* A locked set keeps its lock and loses delete, per the design. */}
            {!draft.locked && (
              <button
                type="button"
                className="focal-row__icon focal-row__icon--danger"
                onClick={() => setConfirmingDelete(true)}
                title="Delete this focal set"
                aria-label={`Delete ${draft.title.trim() || 'this focal set'}`}
              >
                <TrashIcon size={18} />
              </button>
            )}
          </>
        )}
      </span>
    </li>
  );
}

export function FocalSetLibrary(): JSX.Element {
  const { state, dispatch } = useProject();

  return (
    <section className="focal-library" aria-label="Focal set library">
      <ul className="focal-library__list themed-scroll">
        {state.focalDrafts.map((draft) => (
          <LibraryRow
            key={draft.key}
            draft={draft}
            active={draft.key === state.activeFocalKey}
          />
        ))}
      </ul>

      {/*
        The design shows no dedicated "new focal set" control, so this stays a
        modest row action rather than an invented panel. See UI_NOTES §18.
      */}
      <div className="focal-library__footer">
        <button
          type="button"
          className="focal-library__add"
          onClick={() => dispatch({ type: 'addFocalDraft' })}
        >
          + New focal set
        </button>
      </div>
    </section>
  );
}
