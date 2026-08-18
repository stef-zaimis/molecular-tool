import { useEffect, useState } from 'react';
import { CloseIcon, MaximizeIcon, MenuIcon, MinimizeIcon, RestoreIcon } from '../icons/Icons';
import { desktop } from '../../app/desktopApi';
import './TitleBar.css';

export type TitleBarVariant = 'launcher' | 'workspace';

interface TitleBarProps {
  readonly title: string;
  readonly variant: TitleBarVariant;
  /** The launcher window cannot be maximised, so it shows no maximise control. */
  readonly canMaximize?: boolean;
  readonly onMenu?: () => void;
}

/**
 * Custom window chrome.
 *
 * The whole bar is an OS drag region except the interactive controls, which
 * opt out via `.no-drag`. Without that, the buttons would move the window
 * instead of firing.
 */
export function TitleBar({
  title,
  variant,
  canMaximize = true,
  onMenu,
}: TitleBarProps): JSX.Element {
  const [maximized, setMaximized] = useState(false);

  useEffect(() => {
    let cancelled = false;
    desktop()
      .window.isMaximized()
      .then((value) => {
        if (!cancelled) setMaximized(value);
      })
      // Never let window chrome take the UI down: fall back to "not maximised".
      .catch(() => undefined);
    const unsubscribe = desktop().window.onMaximizedChanged(setMaximized);
    return () => {
      cancelled = true;
      unsubscribe();
    };
  }, []);

  return (
    <header className={`titlebar titlebar--${variant} drag-region`}>
      {variant === 'workspace' && (
        <button
          type="button"
          className="titlebar__menu no-drag"
          onClick={onMenu}
          aria-label="Application menu"
        >
          <MenuIcon size={26} />
        </button>
      )}

      <h1 className="titlebar__title">{title}</h1>

      <div className="titlebar__controls no-drag">
        <button
          type="button"
          className="titlebar__control"
          onClick={() => desktop().window.minimize()}
          aria-label="Minimize"
        >
          <MinimizeIcon size={35} />
        </button>

        {canMaximize && (
          <button
            type="button"
            className="titlebar__control"
            onClick={() => desktop().window.toggleMaximize()}
            aria-label={maximized ? 'Restore' : 'Maximize'}
          >
            {maximized ? <RestoreIcon size={35} /> : <MaximizeIcon size={35} />}
          </button>
        )}

        <button
          type="button"
          className="titlebar__control titlebar__control--close"
          onClick={() => desktop().window.close()}
          aria-label="Close"
        >
          <CloseIcon size={35} />
        </button>
      </div>
    </header>
  );
}
