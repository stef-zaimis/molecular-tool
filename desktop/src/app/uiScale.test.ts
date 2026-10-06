import { describe, expect, it } from 'vitest';
import { DESIGN_WIDTH, MIN_ZOOM, workspaceZoom } from './uiScale';

/**
 * The 125% Windows case, as arithmetic.
 *
 * The rule this file protects: the 1920x1080 reference must be pixel-for-pixel
 * untouched, and the same physical screen at 125% scaling (1536x864 CSS) must
 * get exactly 0.8 — which reproduces the reference layout at the same physical
 * size rather than reflowing it into a second design.
 */
describe('workspaceZoom', () => {
  it('is exactly 1 at the design width', () => {
    expect(workspaceZoom(DESIGN_WIDTH)).toBe(1);
  });

  it('is 0.8 at the 125%-scaled viewport of the same screen', () => {
    expect(workspaceZoom(1536)).toBe(0.8);
  });

  it('never magnifies a wider window', () => {
    expect(workspaceZoom(2560)).toBe(1);
    expect(workspaceZoom(1921)).toBe(1);
  });

  it('scales continuously between, so a window a few pixels short is not a cliff', () => {
    // An 1904px window (1920 minus a scrollbar or a border) must not drop to a
    // breakpoint's worth of shrinkage.
    expect(workspaceZoom(1904)).toBeGreaterThan(0.99);
    expect(workspaceZoom(1600)).toBeCloseTo(0.8333, 3);
  });

  it('stops shrinking at the floor rather than becoming unreadable', () => {
    expect(workspaceZoom(400)).toBe(MIN_ZOOM);
  });

  it('answers 1 for a viewport it cannot believe', () => {
    expect(workspaceZoom(0)).toBe(1);
    expect(workspaceZoom(Number.NaN)).toBe(1);
  });

  /*
   * The launcher is NOT measured against this frame — it has its own 540x340
   * window and its own 540x289 design, and treating that as a scaled-down
   * workspace would shrink it to a third of its size. `useWorkspaceZoom(false)`
   * is what keeps it at 1; this is the value it would otherwise get.
   */
  it('would badly mis-scale the launcher window, which is why it is excluded', () => {
    expect(workspaceZoom(540)).toBeLessThan(0.6);
  });
});
