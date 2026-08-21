import { useId } from 'react';
import type { ChangeEvent, KeyboardEvent, ReactNode, RefObject } from 'react';
import { CheckIcon, SpinnerArrowsIcon } from '../icons/Icons';
import './Controls.css';

/* ------------------------------------------------------------------ */
/* FieldRow — the shared "right-aligned label + control" layout         */
/* ------------------------------------------------------------------ */

interface FieldRowProps {
  readonly label: ReactNode;
  readonly htmlFor?: string;
  readonly children: ReactNode;
  /** Aligns the label with the first line when the control is tall. */
  readonly alignTop?: boolean;
}

export function FieldRow({ label, htmlFor, children, alignTop = false }: FieldRowProps): JSX.Element {
  return (
    <div className={`field-row${alignTop ? ' field-row--top' : ''}`}>
      <label className="field-row__label" htmlFor={htmlFor}>
        {label}
      </label>
      <div className="field-row__control">{children}</div>
    </div>
  );
}

/* ------------------------------------------------------------------ */
/* TextField — input with an optional attached action button           */
/* ------------------------------------------------------------------ */

interface TextFieldProps {
  readonly value: string;
  readonly onChange: (value: string) => void;
  readonly placeholder?: string;
  readonly ariaLabel: string;
  readonly id?: string;
  readonly readOnly?: boolean;
  /** For callers that need to focus the field, e.g. a rename started elsewhere. */
  readonly inputRef?: RefObject<HTMLInputElement>;
  readonly width?: number | string;
  readonly compact?: boolean;
  /** Native tooltip, used to surface a full path behind a shortened value. */
  readonly title?: string;
  readonly onSubmit?: () => void;
  readonly action?: {
    readonly label: string;
    readonly onClick: () => void;
    readonly disabled?: boolean;
    readonly title?: string;
    /**
     * Accessible name, when the visible label is not distinguishing.
     *
     * The design puts "Enter" on several unrelated buttons; the label stays as
     * drawn, and this says which one it is.
     */
    readonly ariaLabel?: string;
  };
}

/**
 * The recurring input pattern from the mockups: a flat dark field with a
 * darker action button butted directly against its right edge, no gap and no
 * border between them.
 */
export function TextField({
  value,
  onChange,
  placeholder,
  ariaLabel,
  id,
  readOnly = false,
  inputRef,
  width,
  compact = false,
  title,
  onSubmit,
  action,
}: TextFieldProps): JSX.Element {
  const handleKeyDown = (event: KeyboardEvent<HTMLInputElement>) => {
    if (event.key === 'Enter' && onSubmit) {
      event.preventDefault();
      onSubmit();
    }
  };

  return (
    <div className="textfield">
      <input
        id={id}
        ref={inputRef}
        style={width === undefined ? undefined : { width, flex: '0 0 auto' }}
        className={`textfield__input${compact ? ' is-compact' : ''}`}
        type="text"
        value={value}
        placeholder={placeholder}
        readOnly={readOnly}
        aria-label={ariaLabel}
        title={title}
        spellCheck={false}
        autoComplete="off"
        onChange={(event: ChangeEvent<HTMLInputElement>) => onChange(event.target.value)}
        onKeyDown={handleKeyDown}
      />
      {action && (
        <button
          type="button"
          className="textfield__action"
          onClick={action.onClick}
          disabled={action.disabled}
          title={action.title}
          aria-label={action.ariaLabel}
        >
          {action.label}
        </button>
      )}
    </div>
  );
}

/* ------------------------------------------------------------------ */
/* Checkbox                                                            */
/* ------------------------------------------------------------------ */

interface CheckboxProps {
  readonly checked: boolean;
  readonly onChange: (checked: boolean) => void;
  readonly label: string;
  readonly disabled?: boolean;
}

/**
 * Square checkbox with an oversized tick that breaks out past the top-right
 * corner, as drawn in the mockups. A real <input> carries the semantics; the
 * visual box is a sibling span.
 */
export function Checkbox({ checked, onChange, label, disabled = false }: CheckboxProps): JSX.Element {
  const id = useId();

  return (
    <div className={`checkbox${disabled ? ' is-disabled' : ''}`}>
      <input
        id={id}
        className="checkbox__input"
        type="checkbox"
        checked={checked}
        disabled={disabled}
        onChange={(event) => onChange(event.target.checked)}
      />
      <label className="checkbox__label" htmlFor={id}>
        <span className="checkbox__box" aria-hidden="true">
          {checked && <CheckIcon size={33} className="checkbox__check" />}
        </span>
        <span className="checkbox__text">{label}</span>
      </label>
    </div>
  );
}

/* ------------------------------------------------------------------ */
/* NumberSpinner                                                       */
/* ------------------------------------------------------------------ */

interface NumberSpinnerProps {
  readonly value: number;
  readonly onChange: (value: number) => void;
  readonly min?: number;
  readonly max?: number;
  readonly ariaLabel: string;
}

/**
 * Numeric field with the stepper from the design.
 *
 * ONE arrow-pair glyph, exactly as drawn: the up and down triangles are a
 * single piece of artwork, with two transparent hit areas laid over its halves.
 * A previous version rendered the glyph twice and cropped each copy, which read
 * as two stacked spinners.
 */
export function NumberSpinner({
  value,
  onChange,
  min = 1,
  max = 20,
  ariaLabel,
}: NumberSpinnerProps): JSX.Element {
  const clamp = (next: number) => Math.min(max, Math.max(min, next));

  return (
    <div className="spinner">
      <input
        className="spinner__input"
        type="number"
        inputMode="numeric"
        value={Number.isFinite(value) ? value : ''}
        min={min}
        max={max}
        aria-label={ariaLabel}
        onChange={(event) => {
          const parsed = Number.parseInt(event.target.value, 10);
          if (!Number.isNaN(parsed)) onChange(clamp(parsed));
        }}
      />
      <span className="spinner__stepper">
        <SpinnerArrowsIcon size={25} className="spinner__glyph" />
        <button
          type="button"
          className="spinner__hit spinner__hit--up"
          aria-label={`Increase ${ariaLabel}`}
          tabIndex={-1}
          disabled={value >= max}
          onClick={() => onChange(clamp(value + 1))}
        />
        <button
          type="button"
          className="spinner__hit spinner__hit--down"
          aria-label={`Decrease ${ariaLabel}`}
          tabIndex={-1}
          disabled={value <= min}
          onClick={() => onChange(clamp(value - 1))}
        />
      </span>
    </div>
  );
}

/* ------------------------------------------------------------------ */
/* IconButton                                                          */
/* ------------------------------------------------------------------ */

interface IconButtonProps {
  readonly children: ReactNode;
  readonly label: string;
  readonly onClick?: () => void;
  readonly disabled?: boolean;
  readonly variant?: 'bare' | 'panel';
}

export function IconButton({
  children,
  label,
  onClick,
  disabled = false,
  variant = 'bare',
}: IconButtonProps): JSX.Element {
  return (
    <button
      type="button"
      className={`icon-button icon-button--${variant}`}
      aria-label={label}
      title={label}
      disabled={disabled}
      onClick={onClick}
    >
      {children}
    </button>
  );
}
