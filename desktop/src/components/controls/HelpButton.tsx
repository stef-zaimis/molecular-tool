import { useEffect, useId, useRef, useState } from 'react';
import './HelpButton.css';

interface HelpButtonProps {
  /** Tooltip body. Source the text from src/copy/helpText.ts, never inline. */
  readonly text: string;
  /** Accessible name, e.g. "About Molecular Diagnosis". */
  readonly label: string;
}

/**
 * The single tooltip implementation for the whole app.
 *
 * Behaviour: opens on hover and on keyboard focus, toggles on click, closes on
 * Escape or on pointer/focus leaving. The panel is anchored to the right of the
 * dot with a pointer notch, matching docs/02-project-creation.png, so it never
 * covers the control it explains.
 */
export function HelpButton({ text, label }: HelpButtonProps): JSX.Element {
  const [open, setOpen] = useState(false);
  const [pinned, setPinned] = useState(false);
  const tooltipId = useId();
  const wrapperRef = useRef<HTMLSpanElement>(null);

  useEffect(() => {
    if (!open) return;

    const onKeyDown = (event: KeyboardEvent) => {
      if (event.key === 'Escape') {
        setOpen(false);
        setPinned(false);
      }
    };
    const onPointerDown = (event: PointerEvent) => {
      if (!wrapperRef.current?.contains(event.target as Node)) {
        setOpen(false);
        setPinned(false);
      }
    };

    document.addEventListener('keydown', onKeyDown);
    document.addEventListener('pointerdown', onPointerDown);
    return () => {
      document.removeEventListener('keydown', onKeyDown);
      document.removeEventListener('pointerdown', onPointerDown);
    };
  }, [open]);

  return (
    <span
      ref={wrapperRef}
      className="help"
      onMouseEnter={() => setOpen(true)}
      onMouseLeave={() => !pinned && setOpen(false)}
    >
      <button
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
        ?
      </button>

      {open && (
        <span className="help__tooltip" id={tooltipId} role="tooltip">
          {text}
        </span>
      )}
    </span>
  );
}
