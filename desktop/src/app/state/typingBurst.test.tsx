import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act, cleanup, renderHook } from '@testing-library/react';
import { useReducer } from 'react';
import { TYPING_BURST_IDLE_MS, useTypingBurst } from './typingBurst';
import { appReducer, createInitialState } from './projectState';
import type { AppAction, AppState } from './projectState';
import { activeDraft } from './projectState';

/**
 * How much work ONE Undo may destroy.
 *
 * A burst that only closed at a focus boundary could stay open for an entire
 * editing session, so Undo threw away everything typed since the field was
 * clicked into. These tests pin the opposite failure down as well: one step per
 * keystroke, which makes Undo useless in the other direction.
 *
 * The timer and the reducer are driven together, because the guarantee is about
 * the two of them: the timer must close the group, and closing the group must
 * make the next keystroke record a new step.
 */

/** The editor's wiring, minus CodeMirror: a change, then the idle countdown. */
function useEditingSession() {
  const [state, dispatch] = useReducer(appReducer, undefined, createInitialState);
  const burst = useTypingBurst<string>((key) => dispatch({ type: 'endDraftBurst', key }));

  const key = activeDraft(state).key;

  return {
    state,
    key,
    /** One document change from manual typing. */
    type: (headers: readonly string[]) => {
      dispatch({ type: 'setDraftHeaders', key, headers, coalesce: true });
      burst.touch(key);
    },
    /** Anything that is a boundary: blur, `+`, `−`, Save, undo/redo. */
    boundary: () => burst.endNow(key),
    dispatch: (action: AppAction) => dispatch(action),
  };
}

function draftOf(state: AppState) {
  return activeDraft(state);
}

beforeEach(() => {
  vi.useFakeTimers();
});

afterEach(() => {
  cleanup();
  vi.useRealTimers();
});

describe('a typing burst is one undo step', () => {
  it('coalesces rapid keystrokes into a single step', () => {
    const { result } = renderHook(() => useEditingSession());

    act(() => {
      result.current.type(['h']);
      result.current.type(['he']);
      result.current.type(['hea']);
      result.current.type(['head']);
    });

    expect(draftOf(result.current.state).headers).toEqual(['head']);
    // One step, not four: the history holds only the state the burst began at.
    expect(draftOf(result.current.state).history.past).toEqual([[]]);
  });

  it('does NOT record a step per keystroke', () => {
    const { result } = renderHook(() => useEditingSession());

    act(() => {
      for (const text of ['a', 'ab', 'abc', 'abcd', 'abcde']) result.current.type([text]);
    });

    expect(draftOf(result.current.state).history.past).toHaveLength(1);
  });

  it('keeps the burst open while the user is still typing', () => {
    const { result } = renderHook(() => useEditingSession());

    act(() => {
      result.current.type(['a']);
    });
    // Each change restarts the countdown, so a burst never expires mid-flow.
    act(() => {
      vi.advanceTimersByTime(TYPING_BURST_IDLE_MS - 100);
      result.current.type(['ab']);
    });
    act(() => {
      vi.advanceTimersByTime(TYPING_BURST_IDLE_MS - 100);
      result.current.type(['abc']);
    });

    expect(draftOf(result.current.state).burstOpen).toBe(true);
    expect(draftOf(result.current.state).history.past).toEqual([[]]);
  });
});

describe('the idle timeout ends the burst', () => {
  it('closes it after the pause, without touching the document', () => {
    const { result } = renderHook(() => useEditingSession());

    act(() => {
      result.current.type(['abc']);
    });
    expect(draftOf(result.current.state).burstOpen).toBe(true);

    act(() => {
      vi.advanceTimersByTime(TYPING_BURST_IDLE_MS + 1);
    });

    expect(draftOf(result.current.state).burstOpen).toBe(false);
    // The pause must never rewrite what is on screen: a half-typed header is
    // still the user's text, and canonicalising it here would edit their work.
    expect(draftOf(result.current.state).headers).toEqual(['abc']);
  });

  it('makes the next keystroke a NEW step', () => {
    const { result } = renderHook(() => useEditingSession());

    act(() => {
      result.current.type(['one']);
    });
    act(() => {
      vi.advanceTimersByTime(TYPING_BURST_IDLE_MS + 1);
    });
    act(() => {
      result.current.type(['one', 'two']);
      result.current.type(['one', 'twoo']);
    });

    // Two bursts, two steps.
    expect(draftOf(result.current.state).history.past).toEqual([[], ['one']]);
  });

  it('gives back only the most recent burst on one Undo', () => {
    const { result } = renderHook(() => useEditingSession());

    act(() => {
      result.current.type(['first']);
    });
    act(() => {
      vi.advanceTimersByTime(TYPING_BURST_IDLE_MS + 1);
    });
    act(() => {
      result.current.type(['first', 'second']);
    });
    act(() => {
      vi.advanceTimersByTime(TYPING_BURST_IDLE_MS + 1);
    });

    act(() => {
      result.current.dispatch({ type: 'undoDraftEdit', key: result.current.key });
    });

    // The first burst survives; only the second is undone.
    expect(draftOf(result.current.state).headers).toEqual(['first']);

    act(() => {
      result.current.dispatch({ type: 'undoDraftEdit', key: result.current.key });
    });
    expect(draftOf(result.current.state).headers).toEqual([]);
  });
});

describe('explicit boundaries still close the burst', () => {
  it('ends it immediately and cancels the pending timer', () => {
    const { result } = renderHook(() => useEditingSession());

    act(() => {
      result.current.type(['abc']);
      result.current.boundary();
    });

    expect(draftOf(result.current.state).burstOpen).toBe(false);

    // The cancelled timer must not fire later and close a burst that a new
    // edit has since opened.
    act(() => {
      result.current.type(['abc', 'd']);
    });
    act(() => {
      vi.advanceTimersByTime(TYPING_BURST_IDLE_MS - 50);
    });
    expect(draftOf(result.current.state).burstOpen).toBe(true);
  });
});

describe('the timer does not outlive what it belongs to', () => {
  it('is cleared on unmount', () => {
    const ended: string[] = [];
    const { result, unmount } = renderHook(() =>
      useTypingBurst<string>((key) => ended.push(key)),
    );

    act(() => result.current.touch('draft-1'));
    unmount();
    act(() => {
      vi.advanceTimersByTime(TYPING_BURST_IDLE_MS * 2);
    });

    expect(ended).toEqual([]);
  });

  it('ends the burst that was typed in, not whichever is active later', () => {
    const ended: string[] = [];
    const { result } = renderHook(() => useTypingBurst<string>((key) => ended.push(key)));

    act(() => result.current.touch('draft-1'));
    act(() => {
      vi.advanceTimersByTime(TYPING_BURST_IDLE_MS + 1);
    });

    expect(ended).toEqual(['draft-1']);
  });
});
