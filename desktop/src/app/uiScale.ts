import { useEffect } from 'react';

/**
 * Windows DPI robustness, by scaling the workspace rather than reflowing it.
 *
 * The design is authored at 1920x1080 and every measurement in `tokens.css`
 * comes from that frame. At 125% Windows scaling the SAME physical screen
 * reports a 1536x864 CSS viewport, which is 0.8 of the design in both axes —
 * the window did not get narrower, the pixels got bigger. Reflowing the layout
 * for that case would produce a second design nobody drew, and the 1920 frame
 * would then have to be maintained against it.
 *
 * So the frame is kept and the whole shell is zoomed by the ratio instead. At
 * 1920 the factor is exactly 1 and nothing is touched; at 1536 it is 0.8, so
 * the renderer lays out at 1920x1080 and paints it into 1536x864 — which, on a
 * 125% display, is the same physical size the user would see at 100%. Text
 * included: 27px at 0.8 zoom on a 1.25 device ratio is 27 device-independent
 * pixels again.
 *
 * `zoom` rather than `transform: scale()` because zoom participates in layout:
 * viewport units, fixed positioning, scrolling and hit testing all follow it,
 * where a transform would leave them describing the unscaled box.
 */

/** The frame every token in `tokens.css` was measured against. */
export const DESIGN_WIDTH = 1920;

/**
 * Never shrink past this. Below it the window is so small that scaling further
 * would produce unreadable text; internal scrolling is the better answer, and
 * a window that small is outside what this desktop UI targets.
 */
export const MIN_ZOOM = 0.55;

/**
 * The zoom for a given viewport width.
 *
 * Width only. A short window (1920x900) is a different problem — the analysis
 * column already scrolls internally — and scaling for height would shrink a
 * layout that has no horizontal problem, leaving a wide empty margin down the
 * right-hand side.
 *
 * Never above 1: a wider window keeps the reference size and absorbs the slack
 * the way the design does, rather than magnifying everything.
 */
export function workspaceZoom(viewportWidth: number): number {
  if (!Number.isFinite(viewportWidth) || viewportWidth <= 0) return 1;
  const fitted = viewportWidth / DESIGN_WIDTH;
  return Math.min(1, Math.max(MIN_ZOOM, Number(fitted.toFixed(4))));
}

/** The zoom currently applied, for code that converts pointer deltas. */
export function currentZoom(): number {
  if (typeof document === 'undefined') return 1;
  const value = Number.parseFloat(
    getComputedStyle(document.documentElement).getPropertyValue('--ui-zoom'),
  );
  return Number.isFinite(value) && value > 0 ? value : 1;
}

/**
 * Keep `--ui-zoom` on the document element in step with the window.
 *
 * One listener for the whole app: the value is read by CSS (`#root`) and by the
 * two places that translate pointer pixels into layout pixels.
 *
 * `active` is false on the LAUNCHER, and that is not a detail. The launcher has
 * its own 540x340 window sized to its own 540x289 design; measuring it against
 * the 1920 workspace frame would read a small window as a scaled-down large one
 * and shrink it to a third of its size. Only the workspace frame is scaled.
 */
export function useWorkspaceZoom(active: boolean): void {
  useEffect(() => {
    const apply = () => {
      document.documentElement.style.setProperty(
        '--ui-zoom',
        String(active ? workspaceZoom(window.innerWidth) : 1),
      );
    };

    apply();
    if (!active) return;

    window.addEventListener('resize', apply);
    return () => window.removeEventListener('resize', apply);
  }, [active]);
}
