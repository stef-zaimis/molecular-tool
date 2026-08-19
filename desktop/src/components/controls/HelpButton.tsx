import { useCallback, useEffect, useId, useLayoutEffect, useRef, useState } from 'react';
import { createPortal } from 'react-dom';
import { HelpIcon } from '../icons/Icons';
import './HelpButton.css';

interface HelpButtonProps {
  /** Tooltip body. Source the text from src/copy/helpText.ts, never inline. */
  readonly text: string;
  /** Accessible name, e.g. "About Molecular Diagnosis". */
  readonly label: string;
}

const TOOLTIP_WIDTH = 317;
const GAP = 14;
/** Distance from the dot's vertical centre to the tooltip's top edge. */
const TOP_OFFSET = 14;
const VIEWPORT_MARGIN = 12;

interface Position {
  readonly left: number;
  readonly top: number;
  /** Which side of the dot the panel ended up on. */
  readonly side: 'right' | 'left';
  /** Notch offset from the panel's top edge. */
  readonly notchTop: number;
}

/**
 * The single tooltip implementation for the whole app.
 *
 * Rendered through a PORTAL into <body> with fixed positioning. It used to be
 * absolutely positioned inside the control, which meant any ancestor with
 * `overflow` clipped it and the alignment pane could paint over it. A portal
 * removes it from those stacking and clipping contexts entirely, so the panel
 * is always on top of the workspace regardless of where the dot lives.
 *
 * Opens on hover and on keyboard focus, toggles on click, closes on Escape or
 * on pointer/focus leaving.
 */
export function HelpButton({ text, label }: HelpButtonProps): JSX.Element {
  const [open, setOpen] = useState(false);
  const [pinned, setPinned] = useState(false);
  const [position, setPosition] = useState<Position | null>(null);
  const tooltipId = useId();
  const wrapperRef = useRef<HTMLSpanElement>(null);
  const buttonRef = useRef<HTMLButtonElement>(null);
  const panelRef = useRef<HTMLDivElement>(null);

  const place = useCallback(() => {
    const button = buttonRef.current;
    if (!button) return;

    const rect = button.getBoundingClientRect();
    const panelHeight = panelRef.current?.offsetHeight ?? 0;
    const centreY = rect.top + rect.height / 2;

    // Prefer the design's placement: to the right of the dot. Flip only when
    // the panel would run off the window.
    const fitsRight = rect.right + GAP + TOOLTIP_WIDTH <= window.innerWidth - VIEWPORT_MARGIN;
    const side: 'right' | 'left' = fitsRight ? 'right' : 'left';
    const left = fitsRight ? rect.right + GAP : rect.left - GAP - TOOLTIP_WIDTH;

    let top = centreY - TOP_OFFSET;
    const maxTop = window.innerHeight - VIEWPORT_MARGIN - panelHeight;
    if (panelHeight > 0 && top > maxTop) top = Math.max(VIEWPORT_MARGIN, maxTop);
    if (top < VIEWPORT_MARGIN) top = VIEWPORT_MARGIN;

    // Keep the notch aimed at the dot even after the panel has been nudged.
    const notchTop = Math.max(6, Math.min(centreY - top - 10, Math.max(6, panelHeight - 26)));

    setPosition({ left, top, side, notchTop });
  }, []);

  useLayoutEffect(() => {
    if (!open) return;
    place();
  }, [open, place, text]);

  useEffect(() => {
    if (!open) return;

    const close = () => {
      setOpen(false);
      setPinned(false);
    };
    const onKeyDown = (event: KeyboardEvent) => {
      if (event.key === 'Escape') close();
    };
    const onPointerDown = (event: PointerEvent) => {
      if (!wrapperRef.current?.contains(event.target as Node)) close();
    };

    document.addEventListener('keydown', onKeyDown);
    document.addEventListener('pointerdown', onPointerDown);
    window.addEventListener('resize', place);
    window.addEventListener('scroll', place, true);
    return () => {
      document.removeEventListener('keydown', onKeyDown);
      document.removeEventListener('pointerdown', onPointerDown);
      window.removeEventListener('resize', place);
      window.removeEventListener('scroll', place, true);
    };
  }, [open, place]);

  return (
    <span
      ref={wrapperRef}
      className="help"
      onMouseEnter={() => setOpen(true)}
      onMouseLeave={() => !pinned && setOpen(false)}
    >
      <button
        ref={buttonRef}
        type="button"
        className={`help__dot${open ? ' is-open' : ''}`}
        aria-label={label}
        aria-expanded={open}
        aria-describedby={open ? tooltipId : undefined}
        onClick={() => {
          setPinned((wasPinned) => !wasPinned);
          setOpen((wasOpen) => (pinned ? !wasOpen : true));
        }}
        onFocus={() => setOpen(true)}
        onBlur={() => !pinned && setOpen(false)}
      >
        <HelpIcon size={30} />
      </button>

      {open &&
        createPortal(
          <div
            ref={panelRef}
            className={`help__tooltip help__tooltip--${position?.side ?? 'right'}`}
            id={tooltipId}
            role="tooltip"
            style={{
              left: position?.left ?? -9999,
              top: position?.top ?? -9999,
              // Keep the panel invisible for the first frame, before it is measured.
              visibility: position ? 'visible' : 'hidden',
            }}
          >
            <span
              className="help__notch"
              style={{ top: position?.notchTop ?? 6 }}
              aria-hidden="true"
            />
            {text}
          </div>,
          document.body,
        )}
    </span>
  );
}
