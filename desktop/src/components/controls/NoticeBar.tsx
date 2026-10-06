import './NoticeBar.css';

interface NoticeBarProps {
  readonly message: string;
  readonly onDismiss: () => void;
}

/**
 * Explains an action that exists in the design but has no implementation yet.
 *
 * Only rendered once the user actually activates such a control, so the
 * resting UI never shows placeholder text. Shared by every screen so
 * unfinished actions behave the same way wherever they appear.
 */
export function NoticeBar({ message, onDismiss }: NoticeBarProps): JSX.Element {
  return (
    <div className="notice" role="status">
      <p className="notice__text">{message}</p>
      <button type="button" className="notice__dismiss" onClick={onDismiss}>
        Dismiss
      </button>
    </div>
  );
}
