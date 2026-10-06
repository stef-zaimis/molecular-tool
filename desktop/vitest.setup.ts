/**
 * jsdom gaps that the editor needs.
 *
 * CodeMirror measures the document to decide what to draw, which means calling
 * `Range.getClientRects()` and `Element.getBoundingClientRect()`. jsdom has no
 * layout engine, so it implements neither meaningfully: `getClientRects` is
 * missing from `Range` entirely, and CodeMirror's measure pass throws inside a
 * `requestAnimationFrame` callback where the failure cannot be caught.
 *
 * These stubs return empty/zero geometry. That is honest for a headless DOM —
 * nothing HAS a size here — and it is enough for CodeMirror to skip its
 * viewport optimisation and render the whole document, which is exactly what a
 * test wants to assert against.
 *
 * TEST-ONLY. Nothing in the app depends on these; a real browser supplies the
 * genuine implementations.
 */

if (typeof Range !== 'undefined' && !Range.prototype.getClientRects) {
  Range.prototype.getClientRects = function getClientRects() {
    const list: DOMRect[] = [];
    return Object.assign(list, {
      item: (index: number) => list[index] ?? null,
    }) as unknown as DOMRectList;
  };
}

if (typeof Range !== 'undefined') {
  const zeroRect = () =>
    ({
      x: 0,
      y: 0,
      top: 0,
      left: 0,
      right: 0,
      bottom: 0,
      width: 0,
      height: 0,
      toJSON: () => ({}),
    }) as DOMRect;

  // jsdom's own implementation returns a rect of zeros already; this only
  // guarantees one exists on the Range prototype for CodeMirror to call.
  if (!Range.prototype.getBoundingClientRect) {
    Range.prototype.getBoundingClientRect = zeroRect;
  }
}
