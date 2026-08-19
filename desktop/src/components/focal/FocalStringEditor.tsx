import { DownloadIcon, RedoIcon, UndoIcon } from '../icons/Icons';
import { IconButton } from '../controls/Controls';
import { TOKEN_SEPARATOR, evaluateFocalTokens } from './focalMatching';
import './FocalStringEditor.css';

interface FocalStringEditorProps {
  /** The focal set. An ARRAY is the source of truth — never a joined string. */
  readonly strings: readonly string[];
  /** Known FASTA headers, or null when none are loaded (=> neutral tokens). */
  readonly headers: readonly string[] | null;
  readonly onUndo: () => void;
  readonly onRedo: () => void;
  readonly canUndo: boolean;
  readonly canRedo: boolean;
  readonly onExport: () => void;
  readonly canExport: boolean;
  /** Shown in place of tokens when the set is empty. */
  readonly emptyHint?: string;
}

/**
 * The focal-set display from docs/02-project-creation.svg.
 *
 * READ-ONLY BY DESIGN. Entries are added and removed through the `+`/`−` modes
 * and the STRING field, and reverted through undo/redo. It previously wrapped a
 * textarea whose text was tokenised on every keystroke, which made a joined
 * string the de-facto source of truth and left half-typed entries in the set.
 * The array is now authoritative and this component only renders it, so there
 * is no string-splitting path that can corrupt the set.
 *
 * The ';' between entries is a separator glyph, not data.
 */
export function FocalStringEditor({
  strings,
  headers,
  onUndo,
  onRedo,
  canUndo,
  canRedo,
  onExport,
  canExport,
  emptyHint = 'No focal strings yet.',
}: FocalStringEditorProps): JSX.Element {
  const tokens = evaluateFocalTokens(strings, headers);
  const matched = tokens.filter((token) => token.state === 'match').length;
  const unmatched = tokens.filter((token) => token.state === 'noMatch').length;

  return (
    <div className="focal-editor">
      <div className="focal-editor__gutter">
        <IconButton label="Export focal set as a text file" onClick={onExport} disabled={!canExport}>
          <DownloadIcon size={30} />
        </IconButton>
      </div>

      <div className="focal-editor__box">
        <div className="focal-editor__content themed-scroll">
          {tokens.length === 0 ? (
            <span className="focal-editor__empty">{emptyHint}</span>
          ) : (
            tokens.map((token, index) => (
              <span key={`${token.text}-${index}`} className="focal-editor__entry">
                <span
                  className={`focal-token focal-token--${token.state}`}
                  title={
                    token.state === 'neutral'
                      ? 'Select a FASTA file to validate this string'
                      : `${token.matchCount} matching header${token.matchCount === 1 ? '' : 's'}`
                  }
                >
                  {token.text}
                </span>
                {index < tokens.length - 1 && (
                  <span className="focal-editor__separator">{TOKEN_SEPARATOR} </span>
                )}
              </span>
            ))
          )}
        </div>
      </div>

      <div className="focal-editor__history">
        <IconButton label="Undo focal set edit" onClick={onUndo} disabled={!canUndo} variant="panel">
          <UndoIcon size={26} />
        </IconButton>
        <IconButton label="Redo focal set edit" onClick={onRedo} disabled={!canRedo} variant="panel">
          <RedoIcon size={26} />
        </IconButton>
      </div>

      <p className="sr-only" aria-live="polite">
        {headers === null
          ? `${strings.length} focal strings. No FASTA loaded, so strings cannot be validated.`
          : `${strings.length} focal strings. ${matched} match at least one header, ${unmatched} match none.`}
      </p>
    </div>
  );
}
