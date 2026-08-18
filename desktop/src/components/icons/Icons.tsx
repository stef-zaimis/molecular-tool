/**
 * Icon set.
 *
 * docs/ contained only the three reference screenshots — no exported SVG
 * assets — so every glyph here is hand-drawn to match the reference shapes.
 * House rules, so the set stays coherent:
 *   - single-colour, `currentColor` only
 *   - stroked geometry with round caps/joins, uniform 2px stroke on a 24 grid
 *   - square 24x24 viewBox, visually centred
 *   - no emoji, no third-party icon font
 * Divergences from the mockup are listed in UI_NOTES.md §4.
 */

export interface IconProps {
  readonly size?: number;
  readonly className?: string;
  /** Stroke width on the 24-unit grid. */
  readonly strokeWidth?: number;
}

function svgProps(size: number, className?: string) {
  return {
    width: size,
    height: size,
    viewBox: '0 0 24 24',
    fill: 'none' as const,
    xmlns: 'http://www.w3.org/2000/svg',
    'aria-hidden': true,
    focusable: false,
    className,
  };
}

export function MenuIcon({ size = 26, className }: IconProps): JSX.Element {
  // Non-square by design: the reference glyph measures 35x26 with three 6px bars.
  const width = Math.round((size * 35) / 26);
  return (
    <svg
      width={width}
      height={size}
      viewBox="0 0 35 26"
      fill="none"
      xmlns="http://www.w3.org/2000/svg"
      aria-hidden
      focusable={false}
      className={className}
    >
      <g fill="currentColor">
        <rect x="0" y="0" width="35" height="6" />
        <rect x="0" y="10" width="35" height="6" />
        <rect x="0" y="20" width="35" height="6" />
      </g>
    </svg>
  );
}

/**
 * Solid house with an overhanging roof and a door notch.
 * Traced from the reference glyph, which measures 40x36 with a full-width
 * roof, so this one keeps its own aspect rather than the 24 grid.
 * `size` is the height.
 */
export function HomeIcon({ size = 36, className }: IconProps): JSX.Element {
  const width = Math.round((size * 40) / 36);
  return (
    <svg
      width={width}
      height={size}
      viewBox="0 0 40 36"
      fill="none"
      xmlns="http://www.w3.org/2000/svg"
      aria-hidden
      focusable={false}
      className={className}
    >
      <g fill="currentColor">
        <path d="M20 0 40 19H0Z" />
        <path d="M5 19h30v17H23V24h-7v12H5Z" />
      </g>
    </svg>
  );
}

/** Window minimise. The mockups draw this as a wide chevron/"V", not a dash. */
export function MinimizeIcon({ size = 24, className, strokeWidth = 2 }: IconProps): JSX.Element {
  return (
    <svg {...svgProps(size, className)}>
      <polyline
        points="3,4 12,20 21,4"
        stroke="currentColor"
        strokeWidth={strokeWidth}
        strokeLinecap="round"
        strokeLinejoin="round"
      />
    </svg>
  );
}

export function MaximizeIcon({ size = 24, className, strokeWidth = 2 }: IconProps): JSX.Element {
  return (
    <svg {...svgProps(size, className)}>
      <rect
        x="4.5"
        y="4.5"
        width="15"
        height="15"
        stroke="currentColor"
        strokeWidth={strokeWidth}
        strokeLinejoin="round"
      />
    </svg>
  );
}

/** Restore-down: two offset frames, shown when the window is maximised. */
export function RestoreIcon({ size = 24, className, strokeWidth = 2 }: IconProps): JSX.Element {
  return (
    <svg {...svgProps(size, className)}>
      <g stroke="currentColor" strokeWidth={strokeWidth} strokeLinejoin="round">
        <rect x="4" y="7.5" width="12" height="12" />
        <polyline points="8,7.5 8,4 20,4 20,16 16.5,16" fill="none" />
      </g>
    </svg>
  );
}

export function CloseIcon({ size = 24, className, strokeWidth = 2 }: IconProps): JSX.Element {
  return (
    <svg {...svgProps(size, className)}>
      <g stroke="currentColor" strokeWidth={strokeWidth} strokeLinecap="round">
        <line x1="4" y1="4" x2="20" y2="20" />
        <line x1="20" y1="4" x2="4" y2="20" />
      </g>
    </svg>
  );
}

/** Filled disc with a plus knocked out — matches the STRING row controls. */
export function PlusCircleIcon({ size = 24, className }: IconProps): JSX.Element {
  return (
    <svg {...svgProps(size, className)}>
      <circle cx="12" cy="12" r="10" fill="currentColor" />
      <g stroke="var(--icon-knockout, var(--c-bg))" strokeWidth="2.6" strokeLinecap="round">
        <line x1="12" y1="6.6" x2="12" y2="17.4" />
        <line x1="6.6" y1="12" x2="17.4" y2="12" />
      </g>
    </svg>
  );
}

export function MinusCircleIcon({ size = 24, className }: IconProps): JSX.Element {
  return (
    <svg {...svgProps(size, className)}>
      <circle cx="12" cy="12" r="10" fill="currentColor" />
      <line
        x1="6.6"
        y1="12"
        x2="17.4"
        y2="12"
        stroke="var(--icon-knockout, var(--c-bg))"
        strokeWidth="2.6"
        strokeLinecap="round"
      />
    </svg>
  );
}

/** Import/load: arrow descending into a tray. */
export function ImportIcon({ size = 24, className }: IconProps): JSX.Element {
  return (
    <svg {...svgProps(size, className)}>
      <path d="M10.6 2.6h2.8v8.1h3.5L12 16.7 7.1 10.7h3.5Z" fill="currentColor" />
      <path
        d="M3.4 15.4v4.1a1.6 1.6 0 0 0 1.6 1.6h14a1.6 1.6 0 0 0 1.6-1.6v-4.1"
        stroke="currentColor"
        strokeWidth="2.4"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
    </svg>
  );
}

export function UndoIcon({ size = 24, className }: IconProps): JSX.Element {
  return (
    <svg {...svgProps(size, className)}>
      <path
        d="M8.6 7.4H14a5.6 5.6 0 0 1 0 11.2h-2.4"
        stroke="currentColor"
        strokeWidth="2.6"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
      <path d="M9.6 3.2 4.2 7.4l5.4 4.2Z" fill="currentColor" />
    </svg>
  );
}

export function RedoIcon({ size = 24, className }: IconProps): JSX.Element {
  return (
    <svg {...svgProps(size, className)}>
      <path
        d="M15.4 7.4H10a5.6 5.6 0 0 0 0 11.2h2.4"
        stroke="currentColor"
        strokeWidth="2.6"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
      <path d="M14.4 3.2l5.4 4.2-5.4 4.2Z" fill="currentColor" />
    </svg>
  );
}

/**
 * Checkbox tick. In the mockups the stroke is heavy and deliberately breaks
 * out past the top-right corner of the box, so it is drawn on an oversized
 * canvas and positioned by the Checkbox component rather than clipped.
 */
export function CheckIcon({ size = 24, className }: IconProps): JSX.Element {
  return (
    <svg
      width={size}
      height={size}
      viewBox="0 0 24 24"
      fill="none"
      xmlns="http://www.w3.org/2000/svg"
      aria-hidden
      focusable={false}
      className={className}
    >
      <polyline
        points="3.4,12.4 9.2,18.6 21.4,2.6"
        stroke="currentColor"
        strokeWidth="2.8"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
    </svg>
  );
}

/** Stacked triangles for the number spinner. */
export function SpinnerArrowsIcon({ size = 24, className }: IconProps): JSX.Element {
  return (
    <svg {...svgProps(size, className)}>
      <path d="M12 3.6 17 9.6H7Z" fill="currentColor" />
      <path d="M12 20.4 7 14.4h10Z" fill="currentColor" />
    </svg>
  );
}
