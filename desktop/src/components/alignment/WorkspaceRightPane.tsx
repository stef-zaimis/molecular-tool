import { useCallback, useEffect, useState } from 'react';
import { AlignmentPlaceholder } from './AlignmentPlaceholder';
import { FocalSetLibrary } from '../focal/FocalSetLibrary';
import { useProject } from '../../app/state/ProjectContext';
import { VIEWER_SNAPS, DEFAULT_VIEWER_SNAP } from '../../app/state/projectState';
import './WorkspaceRightPane.css';

/**
 * The workspace's right-hand side: display selectors, an optional focal-set
 * library, and an independently resizable sequence visualizer beneath it.
 *
 * The three are SEPARATE elements. The library used to be rendered inside the
 * visualizer, which tied it to the viewer's width — dragging the viewer moved
 * and resized the list, which is not what the design shows and not what the
 * two controls mean. Now:
 *
 *  - the selectors and the library are right-anchored at a fixed width and do
 *    not move when the viewer is dragged;
 *  - the viewer sits BELOW them in the same column, so its top follows the
 *    library's real height instead of a hardcoded position;
 *  - the viewer's visible LEFT RAIL is the drag handle — one rail, which is
 *    also the light edge the design draws;
 *  - either may be shown or hidden on its own, and hiding the viewer keeps its
 *    width for when it comes back.
 *
 * SHELL ONLY: there is no sequence transport or renderer here, and the viewer's
 * interior is the existing deliberate placeholder.
 */

/** Which snap a fraction of the workspace width is nearest. */
function nearestSnap(fraction: number): number {
  let best = 0;
  let bestDistance = Infinity;
  VIEWER_SNAPS.forEach((snap, index) => {
    const distance = Math.abs(snap - fraction);
    if (distance < bestDistance) {
      bestDistance = distance;
      best = index;
    }
  });
  return best;
}

interface WorkspaceRightPaneProps {
  /** Caption for the placeholder interior: the current FASTA scope. */
  readonly caption: string;
  /** Raised while the edge is being dragged, so tooltips stay out of the way. */
  readonly onDraggingChange: (dragging: boolean) => void;
}

export function WorkspaceRightPane({
  caption,
  onDraggingChange,
}: WorkspaceRightPaneProps): JSX.Element {
  const { state, dispatch } = useProject();
  const [dragging, setDragging] = useState(false);
  /** Live fraction while dragging, so the edge tracks the pointer smoothly. */
  const [previewFraction, setPreviewFraction] = useState<number | null>(null);

  const snapFraction = VIEWER_SNAPS[state.viewerSnap] ?? VIEWER_SNAPS[DEFAULT_VIEWER_SNAP];
  const left = previewFraction ?? snapFraction;

  useEffect(() => onDraggingChange(dragging), [dragging, onDraggingChange]);

  const beginDrag = useCallback(
    (event: React.PointerEvent<HTMLDivElement>) => {
      event.preventDefault();
      const host = event.currentTarget.closest('.diagnosis');
      if (!host) return;

      const bounds = host.getBoundingClientRect();
      setDragging(true);
      (event.target as HTMLElement).setPointerCapture?.(event.pointerId);

      // Clamped to the allowed range, which now reaches 0 — the viewer can be
      // pulled out to the full width of the workspace content.
      const clamp = (value: number) =>
        Math.min(VIEWER_SNAPS[VIEWER_SNAPS.length - 1], Math.max(VIEWER_SNAPS[0], value));

      const fractionAt = (clientX: number) => clamp((clientX - bounds.left) / bounds.width);

      const move = (moveEvent: PointerEvent) => {
        // Preview follows the pointer; the release decides the snap.
        setPreviewFraction(fractionAt(moveEvent.clientX));
      };

      const end = (endEvent: PointerEvent) => {
        window.removeEventListener('pointermove', move);
        window.removeEventListener('pointerup', end);
        window.removeEventListener('pointercancel', end);
        setDragging(false);
        setPreviewFraction(null);
        dispatch({ type: 'setViewerSnap', snap: nearestSnap(fractionAt(endEvent.clientX)) });
      };

      window.addEventListener('pointermove', move);
      window.addEventListener('pointerup', end);
      window.addEventListener('pointercancel', end);
    },
    [dispatch],
  );

  const onEdgeKeyDown = (event: React.KeyboardEvent<HTMLDivElement>) => {
    if (event.key === 'ArrowLeft') {
      event.preventDefault();
      dispatch({ type: 'setViewerSnap', snap: state.viewerSnap - 1 });
    } else if (event.key === 'ArrowRight') {
      event.preventDefault();
      dispatch({ type: 'setViewerSnap', snap: state.viewerSnap + 1 });
    } else if (event.key === 'Home') {
      event.preventDefault();
      dispatch({ type: 'setViewerSnap', snap: DEFAULT_VIEWER_SNAP });
    }
  };

  return (
    /*
     * ONE column: the selectors and library, then the viewer beneath them.
     *
     * The viewer's top is no longer a number. It used to be one of two
     * constants — 73 without the library, 302 with it — so a library of two
     * rows pushed the viewer as far down as a library of eight would, and the
     * band between them stayed permanently empty. As the next flex item after
     * the pane, the viewer begins a fixed 16px below whatever the pane actually
     * renders, and climbs back up when rows are removed. When the list reaches
     * its own max-height it scrolls internally, so the viewer stops descending
     * there.
     *
     * Horizontal independence is unchanged: the pane is right-anchored at its
     * own width, and the viewer's left edge is a margin percentage of this
     * column, so dragging the viewer moves neither the pane nor the analysis
     * controls underneath it.
     */
    <div className="workspace-right">
      {/*
        Selectors and library: right-anchored, fixed width, ABOVE the viewer and
        entirely independent of it.
      */}
      <div className="right-pane">
        <div className="right-pane__selectors" role="group" aria-label="Right pane">
          <button
            type="button"
            className={`right-pane__selector${state.libraryVisible ? ' is-active' : ''}`}
            aria-pressed={state.libraryVisible}
            onClick={() =>
              dispatch({ type: 'setLibraryVisible', visible: !state.libraryVisible })
            }
          >
            FOCAL SET LIBRARY
          </button>
          <button
            type="button"
            className={`right-pane__selector${state.viewerVisible ? ' is-active' : ''}`}
            aria-pressed={state.viewerVisible}
            onClick={() => dispatch({ type: 'setViewerVisible', visible: !state.viewerVisible })}
          >
            SEQUENCE VISUALIZER
          </button>
        </div>

        {state.libraryVisible && (
          <div className="right-pane__library">
            <FocalSetLibrary />
          </div>
        )}
      </div>

      {state.viewerVisible && (
        <section
          className={`viewer${dragging ? ' is-dragging' : ''}`}
          style={{
            marginLeft: `${left * 100}%`,
          }}
          aria-label="Sequence visualizer"
        >
          {/*
            THE rail — one left edge, not two.
            
            The design draws a single full-height light rail down the viewer's
            left boundary with the panel starting immediately after it. That
            rail used to live inside the placeholder while a second, separate
            drag edge was drawn beside it, so the workspace showed two thin
            edges where the reference has one. The rail now belongs to the
            viewer, which is what owns the drag, and the placeholder draws only
            what is inside the frame.
          */}
          <div
            className="viewer__rail"
            role="separator"
            aria-orientation="vertical"
            aria-label="Resize the sequence visualizer"
            aria-valuenow={state.viewerSnap + 1}
            aria-valuemin={1}
            aria-valuemax={VIEWER_SNAPS.length}
            tabIndex={0}
            onPointerDown={beginDrag}
            onDoubleClick={() => dispatch({ type: 'setViewerSnap', snap: DEFAULT_VIEWER_SNAP })}
            onKeyDown={onEdgeKeyDown}
            title="Drag to resize · double-click to reset"
          >
            <span className="viewer__rail-thumb" aria-hidden="true" />
          </div>

          <div className="viewer__stage">
            <AlignmentPlaceholder caption={caption} />
          </div>
        </section>
      )}
    </div>
  );
}
