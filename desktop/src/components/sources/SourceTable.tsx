import { useState } from 'react';
import { useProject } from '../../app/state/ProjectContext';
import type { SourceView } from '../../app/state/sourceStatus';
import type { FastaCandidate } from '../../backendContract';
import { HelpButton } from '../controls/HelpButton';
import { LockClosedIcon, LockOpenIcon, PencilIcon, TrashIcon } from '../icons/Icons';
import { HELP_TEXT } from '../../copy/helpText';
import './SourceTable.css';

/**
 * The FASTA list, from page 3 of the design:
 *
 *   ▌[1] Leptacis_tipulae_only.fasta ✎   375 seq ?   712 bp ?   591 PISs ?   🔓 🗑
 *   ▌[2] Leptacis_allSequences-BOLD…       2,123         756         596
 *
 * ONE component for both project states. Before the project exists these rows
 * are vetted CANDIDATES; afterwards they are linked sources with live status.
 * The user is looking at the same list either way, so it is the same list —
 * a second, differently-styled "pending files" box was the thing to remove.
 *
 * Unit labels and their help dots appear on the FIRST row only, exactly as the
 * reference draws them: they explain the columns once, and the numbers below
 * line up underneath because the columns are fixed rather than content-sized.
 */

/** What a row needs, whichever kind it is. */
interface RowModel {
  readonly key: string;
  readonly displayName: string;
  readonly path: string;
  readonly sequenceCount: number | null;
  readonly alignmentLength: number | null;
  readonly locked: boolean;
  /** Live filesystem status; absent for a candidate that is not linked yet. */
  readonly view: SourceView | null;
  /** The linked file's id, or null while it is only a candidate. */
  readonly fastaFileId: string | null;
}

function fromCandidate(candidate: FastaCandidate, locked: boolean): RowModel {
  return {
    key: candidate.path,
    displayName: candidate.displayName,
    path: candidate.path,
    sequenceCount: candidate.sequenceCount,
    alignmentLength: candidate.alignmentLength,
    locked,
    view: null,
    fastaFileId: null,
  };
}

function fromView(view: SourceView): RowModel {
  const { payload } = view;
  return {
    key: payload.fastaFileId,
    displayName: payload.displayName,
    path: payload.sourcePath,
    sequenceCount: payload.sequenceCount,
    alignmentLength: payload.alignmentLength,
    locked: payload.locked,
    view,
    fastaFileId: payload.fastaFileId,
  };
}

function formatCount(value: number | null): string {
  return value === null ? '—' : value.toLocaleString();
}

interface SourceRowProps {
  readonly row: RowModel;
  readonly index: number;
  /** Only the first row carries the unit labels and their help dots. */
  readonly showUnits: boolean;
  readonly onRemove: () => void;
  /** Locks or unlocks this row, wherever its lock happens to live. */
  readonly onSetLocked: (locked: boolean) => void;
}

function SourceRow({ row, index, showUnits, onRemove, onSetLocked }: SourceRowProps): JSX.Element {
  const { reindexSource, relinkSource } = useProject();
  const [confirming, setConfirming] = useState(false);

  const view = row.view;
  const tone = view?.tone ?? 'ok';

  return (
    <li className={`source-row source-row--${tone}${row.locked ? ' is-locked' : ''}`}>
      <span className="source-row__body">
        {/*
          The full-height light edge from the reference. It is the row's first
          GRID CELL, not a box floating inside it, which is what lets
          `align-self: stretch` scale it with the row height.
        */}
        <span className="source-row__edge" aria-hidden="true" />

        <span className="source-row__index">[{index + 1}]</span>

        {/*
          Name and pencil are ONE inline group, so the pencil begins where the
          rendered name ends rather than at a column edge far to its right. The
          group is what truncates: a long name ellipsises and the pencil stays
          against it, and neither can reach the seq/bp/PIS columns because the
          group's track ends before them.
        */}
        <span className="source-row__name-group">
          <span className="source-row__name" title={row.path}>
            {row.displayName}
          </span>

          {/*
            A locked row loses the pencil entirely: hiding it while leaving it
            clickable would be a lock that only looks like one.
          */}
          {!row.locked && (
            <button
              type="button"
              className="source-row__icon source-row__icon--rename"
              disabled
              title="Renaming a linked file is not available yet."
              aria-label={`Rename ${row.displayName} (not available yet)`}
            >
              <PencilIcon size={18} />
            </button>
          )}
        </span>

        <span className="source-row__metric">
          <span className="source-row__value">{formatCount(row.sequenceCount)}</span>
          {showUnits && (
            <>
              <span className="source-row__unit">seq</span>
              <HelpButton
                variant="on-row"
                label="About the sequence count"
                text={HELP_TEXT.sequenceCount}
              />
            </>
          )}
        </span>

        <span className="source-row__metric">
          <span className="source-row__value">{formatCount(row.alignmentLength)}</span>
          {showUnits && (
            <>
              <span className="source-row__unit">bp</span>
              <HelpButton
                variant="on-row"
                label="About the alignment length"
                text={HELP_TEXT.alignmentLength}
              />
            </>
          )}
        </span>

        {/* Nothing computes parsimony-informative sites yet; a dash says so,
            a 0 would be a claim. */}
        <span className="source-row__metric source-row__metric--empty">
          <span className="source-row__value">—</span>
          {showUnits && (
            <>
              <span className="source-row__unit">PIS</span>
              <HelpButton
                variant="on-row"
                label="About parsimony-informative sites"
                text={HELP_TEXT.pis}
              />
            </>
          )}
        </span>

        <span className="source-row__actions">
          {confirming ? (
            /* In-row, over the actions: the list must not grow or scroll to
               ask a yes/no question about one of its rows. */
            <span className="source-row__confirm" role="alertdialog" aria-label="Confirm removal">
              <button
                type="button"
                className="source-row__confirm-action source-row__confirm-action--danger"
                onClick={() => {
                  setConfirming(false);
                  onRemove();
                }}
              >
                Remove
              </button>
              <button
                type="button"
                className="source-row__confirm-action"
                onClick={() => setConfirming(false)}
                onKeyDown={(event) => {
                  if (event.key === 'Escape') setConfirming(false);
                }}
              >
                Cancel
              </button>
            </span>
          ) : (
            <>
              {/*
                A locked row keeps only its lock, and keeps it visible.

                It works before the project exists too: a candidate has no
                `fastaFileId` to lock through the backend yet, so the lock is
                held locally by path and applied for real the moment Create
                links the file. Same button, same states, same meaning.
              */}
              <button
                type="button"
                className="source-row__icon source-row__icon--lock"
                onClick={() => onSetLocked(!row.locked)}
                title={
                  row.locked
                    ? 'Unlock this file'
                    : 'Lock this file so it cannot be removed'
                }
                aria-label={`${row.locked ? 'Unlock' : 'Lock'} ${row.displayName}`}
                aria-pressed={row.locked}
              >
                {row.locked ? <LockClosedIcon size={18} /> : <LockOpenIcon size={18} />}
              </button>

              {!row.locked && (
                <button
                  type="button"
                  className="source-row__icon source-row__icon--danger"
                  onClick={() => setConfirming(true)}
                  title="Remove this file from the project"
                  aria-label={`Remove ${row.displayName}`}
                >
                  <TrashIcon size={18} />
                </button>
              )}
            </>
          )}
        </span>
      </span>

      {/* Only a row that needs attention grows a second line. */}
      {view?.needsAttention && (
        <div className="source-row__status" role="status">
          <span className="source-row__state">{view.label}</span>
          <span className="source-row__detail">{view.detail}</span>
          {(view.payload.state === 'missing' || view.payload.state === 'unreadable') && (
            <button
              type="button"
              className="source-row__action"
              onClick={() => void relinkSource(view.payload.fastaFileId)}
              disabled={!view.canRelink}
            >
              Relink
            </button>
          )}
          {(view.payload.state === 'stale' || view.payload.state === 'never_indexed') && (
            <button
              type="button"
              className="source-row__action"
              onClick={() => void reindexSource(view.payload.fastaFileId)}
              disabled={!view.canReindex}
            >
              Re-index
            </button>
          )}
        </div>
      )}
    </li>
  );
}

export function SourceTable(): JSX.Element {
  const {
    state,
    sourceViews,
    sourceBanner,
    refreshSources,
    unlinkSource,
    setSourceLocked,
    dispatch,
  } = useProject();

  const projectOpen = state.project.status === 'open';
  const lockedPaths = new Set(state.pending.lockedPaths);
  const rows: readonly RowModel[] = projectOpen
    ? sourceViews.map(fromView)
    : state.pending.candidates.map((candidate) =>
        fromCandidate(candidate, lockedPaths.has(candidate.path)),
      );

  return (
    <section className="source-table" aria-label="FASTA files">
      {sourceBanner && projectOpen && (
        <p className="source-table__banner" role="alert">
          {sourceBanner}
        </p>
      )}

      {rows.length === 0 ? (
        <p className="source-table__empty">
          {projectOpen
            ? 'No FASTA files are linked to this project yet. Use Browse to add one.'
            : 'No FASTA files chosen yet. They can also be added after the project is created.'}
        </p>
      ) : (
        <ul className="source-table__list">
          {rows.map((row, index) => (
            <SourceRow
              key={row.key}
              row={row}
              index={index}
              showUnits={index === 0}
              onRemove={() =>
                row.fastaFileId
                  ? void unlinkSource(row.fastaFileId)
                  : dispatch({ type: 'removePendingCandidate', path: row.path })
              }
              /*
               * One lock, two places to keep it. A linked source locks through
               * the backend, which owns the answer; a candidate locks locally
               * by path until Create gives it a row to lock.
               */
              onSetLocked={(locked) =>
                row.fastaFileId
                  ? void setSourceLocked(row.fastaFileId, locked)
                  : dispatch({ type: 'setPendingCandidateLocked', path: row.path, locked })
              }
            />
          ))}
        </ul>
      )}

      {projectOpen && rows.length > 0 && (
        <div className="source-table__footer">
          <button
            type="button"
            className="source-row__action"
            onClick={() => void refreshSources(true)}
            disabled={state.sources.refreshing}
          >
            {/* An explicit user action, so this one is allowed to hash. */}
            {state.sources.refreshing ? 'Checking...' : 'Verify now'}
          </button>
        </div>
      )}
    </section>
  );
}
