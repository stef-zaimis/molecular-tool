import { useProject } from '../../app/state/ProjectContext';
import type { SourceView } from '../../app/state/sourceStatus';
import './SourceStatusPanel.css';

/**
 * Live state of the project's linked FASTA files.
 *
 * The panel exists because a project stores a PATH, not a copy: between one
 * session and the next a file can be moved, edited, or deleted, and the user
 * needs to find that out here rather than from a failed analysis. Everything
 * shown is re-read from disk; nothing is a remembered status.
 *
 * Two recovery actions are offered inline, matched to the two failure modes:
 * a missing file needs a new path (Relink), a changed file needs its index
 * rebuilt (Re-index).
 */
function SourceRow({ view }: { view: SourceView }): JSX.Element {
  const { reindexSource, relinkSource } = useProject();
  const { payload } = view;

  return (
    <li className={`source-row source-row--${view.tone}`}>
      <div className="source-row__main">
        <span className="source-row__name" title={payload.sourcePath}>
          {payload.displayName}
        </span>
        <span className="source-row__state" role="status">
          {view.label}
        </span>
      </div>

      <p className="source-row__detail">{view.detail}</p>

      {payload.state === 'current' && payload.sequenceCount !== null && (
        <p className="source-row__facts">
          {payload.sequenceCount} sequences, alignment length {payload.alignmentLength}
          {payload.duplicateHeaderCount > 0 &&
            ` (${payload.duplicateHeaderCount} duplicate headers)`}
        </p>
      )}

      <div className="source-row__actions">
        {(payload.state === 'missing' || payload.state === 'unreadable') && (
          <button
            type="button"
            className="source-row__action"
            onClick={() => void relinkSource(payload.fastaFileId)}
            disabled={!view.canRelink}
          >
            Relink
          </button>
        )}
        {(payload.state === 'stale' || payload.state === 'never_indexed') && (
          <button
            type="button"
            className="source-row__action"
            onClick={() => void reindexSource(payload.fastaFileId)}
            disabled={!view.canReindex}
          >
            Re-index
          </button>
        )}
      </div>
    </li>
  );
}

export function SourceStatusPanel(): JSX.Element | null {
  const { state, sourceViews, sourceBanner, refreshSources } = useProject();

  if (state.project.status !== 'open') return null;

  const { capabilities } = state.project.project;

  return (
    <section className="sources" aria-label="Linked FASTA files">
      <header className="sources__header">
        <span className="sources__title">LINKED FASTA FILES</span>
        <button
          type="button"
          className="sources__refresh"
          onClick={() => void refreshSources(true)}
          disabled={state.sources.refreshing}
        >
          {/* Explicit user action, so this one is allowed to hash. */}
          {state.sources.refreshing ? 'Checking...' : 'Verify now'}
        </button>
      </header>

      {sourceBanner && (
        <p className="sources__banner" role="alert">
          {sourceBanner}
        </p>
      )}

      {sourceViews.length === 0 ? (
        <p className="sources__empty">No FASTA files are linked to this project yet.</p>
      ) : (
        <ul className="sources__list">
          {sourceViews.map((view) => (
            <SourceRow key={view.payload.fastaFileId} view={view} />
          ))}
        </ul>
      )}

      {!capabilities.acceleratedSearch && (
        <p className="sources__note">
          This SQLite build ({capabilities.sqliteVersion}) has no trigram full-text index, so header
          search runs unaccelerated. Results are identical, only slower on large files.
        </p>
      )}
    </section>
  );
}
