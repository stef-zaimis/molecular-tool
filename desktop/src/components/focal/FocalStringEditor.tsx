import { DownloadIcon, RedoIcon, UndoIcon } from '../icons/Icons';
import { IconButton } from '../controls/Controls';
import { TOKEN_SEPARATOR } from './focalMatching';
import type { FocalStringValidation } from '../../backendContract';
import './FocalStringEditor.css';

interface FocalStringEditorProps {
  /** The focal set. An ARRAY is the source of truth — never a joined string. */
  readonly strings: readonly string[];
  /**
   * Backend verdict per entry, aligned by index, or null when nothing has been
   * validated yet. Null renders neutral: "cannot say", not "no match".
   */
  readonly validations: readonly FocalStringValidation[] | null;
  readonly onUndo: () => void;
  readonly onRedo: () => void;
  readonly canUndo: boolean;
  readonly canRedo: boolean;
  readonly onExport: () => void;
  readonly canExport: boolean;
  /** Shown in place of entries when the set is empty. */
  readonly emptyHint?: string;
}

type EntryState = 'match' | 'noMatch' | 'neutral';

function entryState(validation: FocalStringValidation | undefined): EntryState {
  if (!validation) return 'neutral';
  return validation.matched ? 'match' : 'noMatch';
}

function entryTitle(entry: string, validation: FocalStringValidation | undefined): string {
  if (!validation) return `${entry} — select a FASTA file to validate this string`;
  if (validation.matched) {
    const count = validation.matchCount;
    return `${entry} — matches ${count} header${count === 1 ? '' : 's'}`;
  }
  return `${entry} — matches no header in this FASTA`;
}

/**
 * The focal-set display.
 *
 * READ-ONLY BY DESIGN. Entries are added and removed through the `+`/`−` modes
 * and the STRING field, and reverted through undo/redo.
 *
 * The `;` between entries is a RENDERED SEPARATOR, not stored data and not an
 * editable character. That is what lets an entry legitimately contain a `;`
 * without the display becoming ambiguous — the separator is a distinct element
 * in its own colour, and each entry is one span carrying its own state.
 *
 * Colours come from the backend (`validateFocalStrings`), which runs the same
 * matcher the analysis uses, so they cannot disagree with the sequences the
 * run will actually select.
 */
export function FocalStringEditor({
  strings,
  validations,
  onUndo,
  onRedo,
  canUndo,
  canRedo,
  onExport,
  canExport,
  emptyHint = 'No focal strings yet.',
}: FocalStringEditorProps): JSX.Element {
  const matched = validations?.filter((entry) => entry.matched).length ?? 0;
  const unmatched = validations ? validations.length - matched : 0;

  return (
    <div className="focal-editor">
      <div className="focal-editor__gutter">
        <IconButton label="Export focal set as a text file" onClick={onExport} disabled={!canExport}>
          <DownloadIcon size={30} />
        </IconButton>
      </div>

      <div className="focal-editor__box">
        <div className="focal-editor__content themed-scroll">
          {strings.length === 0 ? (
            <span className="focal-editor__empty">{emptyHint}</span>
          ) : (
            strings.map((entry, index) => {
              const validation = validations?.[index];
              return (
                <span key={`${entry}-${index}`} className="focal-editor__entry">
                  <span
                    className={`focal-token focal-token--${entryState(validation)}`}
                    title={entryTitle(entry, validation)}
                  >
                    {entry}
                  </span>
                  {index < strings.length - 1 && (
                    <span className="focal-editor__separator" aria-hidden="true">
                      {TOKEN_SEPARATOR}{' '}
                    </span>
                  )}
                </span>
              );
            })
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
        {validations === null
          ? `${strings.length} focal strings. Not validated against a FASTA yet.`
          : `${strings.length} focal strings. ${matched} match at least one header, ${unmatched} match none.`}
      </p>
    </div>
  );
}
