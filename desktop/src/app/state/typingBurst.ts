import { useCallback, useEffect, useMemo, useRef } from 'react';

/**
 * When a manual typing burst ends by itself.
 *
 * A burst is the unit Undo walks back, so its length decides how much work one
 * Undo can destroy. Closing it only at a focus or action boundary meant a burst
 * could stay open for an entire editing session: click into the field, type for
 * two minutes without leaving it, and a single Undo threw all of it away.
 *
 * 700ms is the usual "the user stopped typing" pause — long enough not to split
 * a header mid-word while someone is looking at the keyboard, short enough that
 * a deliberate pause between two edits records them as two steps.
 */
export const TYPING_BURST_IDLE_MS = 700;

export interface TypingBurst<T> {
  /**
   * A document change happened. Restarts the idle countdown, so a burst ends
   * only after a genuine pause rather than a fixed time after it began.
   */
  readonly touch: (payload: T) => void;
  /** Close the burst now — blur, `+`, `−`, Save, undo/redo, a canonical rewrite. */
  readonly endNow: (payload: T) => void;
  /** Drop a pending countdown WITHOUT closing anything. */
  readonly cancel: () => void;
}

/**
 * An idle timer that closes a typing burst.
 *
 * It ONLY closes the history group. It must never canonicalise, reformat or
 * otherwise touch the document: a user who pauses halfway through typing a
 * header has not asked for their unfinished text to be rewritten, and the
 * pause is not a decision about anything except undo granularity.
 *
 * The payload is carried through the timer so the burst that ends is the one
 * that was being typed in, not whichever draft happens to be active when the
 * timer fires.
 */
export function useTypingBurst<T>(
  end: (payload: T) => void,
  idleMs: number = TYPING_BURST_IDLE_MS,
): TypingBurst<T> {
  const timer = useRef<ReturnType<typeof setTimeout> | null>(null);
  const endRef = useRef(end);
  endRef.current = end;

  const cancel = useCallback(() => {
    if (timer.current === null) return;
    clearTimeout(timer.current);
    timer.current = null;
  }, []);

  const touch = useCallback(
    (payload: T) => {
      cancel();
      timer.current = setTimeout(() => {
        timer.current = null;
        endRef.current(payload);
      }, idleMs);
    },
    [cancel, idleMs],
  );

  const endNow = useCallback(
    (payload: T) => {
      cancel();
      endRef.current(payload);
    },
    [cancel],
  );

  // Unmounting must not leave a timer that dispatches into a gone component.
  useEffect(() => cancel, [cancel]);

  // Stable identity, so callers can list it as a hook dependency.
  return useMemo(() => ({ touch, endNow, cancel }), [touch, endNow, cancel]);
}
