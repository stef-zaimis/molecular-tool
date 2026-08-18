import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import type { KeyboardEvent } from 'react';
import { ImportIcon, RedoIcon, UndoIcon } from '../icons/Icons';
import { IconButton } from '../controls/Controls';
import {
  TOKEN_SEPARATOR,
  evaluateFocalTokens,
  tokenizeFocalInput,
} from './focalMatching';
import type { FocalTokenState } from './focalMatching';
import './FocalStringEditor.css';

interface FocalStringEditorProps {
  readonly value: string;
  readonly onChange: (value: string) => void;
  /** Known FASTA headers, or null when nothing is loaded (=> neutral tokens). */
  readonly headers: readonly string[] | null;
  readonly onImport?: () => void;
}

interface Segment {
  readonly text: string;
  readonly kind: 'token' | 'separator' | 'space';
  readonly state: FocalTokenState;
}

/**
 * Split raw text into render segments WITHOUT losing a single character.
 *
 * The highlight layer must be glyph-for-glyph identical to the textarea it
 * sits behind, otherwise the colours drift away from the caret. So separators
 * and whitespace are emitted as their own segments rather than trimmed away,
 * and only the token body carries a match colour.
 */
function segmentValue(value: string, headers: readonly string[] | null): Segment[] {
  const tokens = tokenizeFocalInput(value);
  const evaluated = evaluateFocalTokens(tokens, headers);
  const stateByToken = new Map<string, FocalTokenState>();
  evaluated.forEach((token) => stateByToken.set(token.text, token.state));

  const segments: Segment[] = [];

  value.split(TOKEN_SEPARATOR).forEach((chunk, index, all) => {
    if (chunk.length > 0) {
      const leading = chunk.length - chunk.trimStart().length;
      const trailing = chunk.length - chunk.trimEnd().length;
      const body = chunk.slice(leading, chunk.length - trailing);

      if (leading > 0) {
        segments.push({ text: chunk.slice(0, leading), kind: 'space', state: 'neutral' });
      }
      if (body.length > 0) {
        segments.push({
          text: body,
          kind: 'token',
          state: stateByToken.get(body) ?? 'neutral',
        });
      }
      if (trailing > 0) {
        segments.push({
          text: chunk.slice(chunk.length - trailing),
          kind: 'space',
          state: 'neutral',
        });
      }
    }

    if (index < all.length - 1) {
      segments.push({ text: TOKEN_SEPARATOR, kind: 'separator', state: 'neutral' });
    }
  });

  return segments;
}

const HISTORY_LIMIT = 100;

/**
 * The STRING area from docs/03-molecular-diagnosis.png.
 *
 * Implementation: a transparent <textarea> layered over a mirrored highlight
 * <div>. This keeps native caret placement, selection, IME and clipboard
 * behaviour — all of which a chips/contenteditable implementation would have to
 * re-create — while still colouring each token green/red/neutral. No editor
 * dependency (CodeMirror et al.) is warranted for a single-line-ish token list.
 *
 * Undo/redo is an explicit stack rather than the browser's, so that the on-screen
 * arrows and the keyboard shortcuts share one history and stay in sync.
 */
export function FocalStringEditor({
  value,
  onChange,
  headers,
  onImport,
}: FocalStringEditorProps): JSX.Element {
  const textareaRef = useRef<HTMLTextAreaElement>(null);
  const highlightRef = useRef<HTMLDivElement>(null);
  const [focused, setFocused] = useState(false);

  const [past, setPast] = useState<string[]>([]);
  const [future, setFuture] = useState<string[]>([]);
  // Suppresses history capture while we are ourselves applying an undo/redo.
  const applyingHistory = useRef(false);

  const segments = useMemo(() => segmentValue(value, headers), [value, headers]);

  const commit = useCallback(
    (next: string) => {
      if (next === value) return;
      if (!applyingHistory.current) {
        setPast((stack) => [...stack, value].slice(-HISTORY_LIMIT));
        setFuture([]);
      }
      onChange(next);
    },
    [onChange, value],
  );

  const undo = useCallback(() => {
    setPast((stack) => {
      if (stack.length === 0) return stack;
      const previous = stack[stack.length - 1];
      applyingHistory.current = true;
      setFuture((forward) => [value, ...forward].slice(0, HISTORY_LIMIT));
      onChange(previous);
      applyingHistory.current = false;
      return stack.slice(0, -1);
    });
  }, [onChange, value]);

  const redo = useCallback(() => {
    setFuture((stack) => {
      if (stack.length === 0) return stack;
      const [next, ...rest] = stack;
      applyingHistory.current = true;
      setPast((backward) => [...backward, value].slice(-HISTORY_LIMIT));
      onChange(next);
      applyingHistory.current = false;
      return rest;
    });
  }, [onChange, value]);

  // Keep the highlight layer scrolled in lockstep with the textarea.
  const syncScroll = useCallback(() => {
    if (!textareaRef.current || !highlightRef.current) return;
    highlightRef.current.scrollTop = textareaRef.current.scrollTop;
    highlightRef.current.scrollLeft = textareaRef.current.scrollLeft;
  }, []);

  useEffect(() => {
    syncScroll();
  }, [value, syncScroll]);

  const handleKeyDown = (event: KeyboardEvent<HTMLTextAreaElement>) => {
    const mod = event.ctrlKey || event.metaKey;
    if (!mod) return;

    const key = event.key.toLowerCase();
    if (key === 'z' && !event.shiftKey) {
      event.preventDefault();
      undo();
    } else if ((key === 'z' && event.shiftKey) || key === 'y') {
      event.preventDefault();
      redo();
    }
  };

  const tokenCount = tokenizeFocalInput(value).length;

  return (
    <div className="focal-editor">
      <div className="focal-editor__gutter">
        <IconButton
          label="Load focal strings from file"
          onClick={onImport}
          // TODO(backend): needs a project/focal-set file format. See UI_NOTES Q3.
          disabled={!onImport}
        >
          <ImportIcon size={24} />
        </IconButton>
      </div>

      <div className={`focal-editor__box${focused ? ' is-focused' : ''}`}>
        <div className="focal-editor__highlight themed-scroll" ref={highlightRef} aria-hidden="true">
          {segments.map((segment, index) => (
            <span
              key={index}
              className={
                segment.kind === 'token'
                  ? `focal-token focal-token--${segment.state}`
                  : `focal-editor__${segment.kind}`
              }
            >
              {segment.text}
            </span>
          ))}
          {/* Trailing newline keeps the mirror's height in step with the textarea. */}
          {'\n'}
        </div>

        <textarea
          ref={textareaRef}
          className="focal-editor__input themed-scroll"
          value={value}
          spellCheck={false}
          autoComplete="off"
          wrap="soft"
          aria-label="Focal search strings, separated by semicolons"
          aria-describedby="focal-editor-status"
          onChange={(event) => commit(event.target.value)}
          onScroll={syncScroll}
          onKeyDown={handleKeyDown}
          onFocus={() => setFocused(true)}
          onBlur={() => setFocused(false)}
        />
      </div>

      <div className="focal-editor__history">
        <IconButton label="Undo" onClick={undo} disabled={past.length === 0} variant="panel">
          <UndoIcon size={30} />
        </IconButton>
        <IconButton label="Redo" onClick={redo} disabled={future.length === 0} variant="panel">
          <RedoIcon size={30} />
        </IconButton>
      </div>

      <p id="focal-editor-status" className="sr-only" aria-live="polite">
        {headers === null
          ? `${tokenCount} focal strings entered. No alignment loaded, so strings cannot be validated.`
          : `${tokenCount} focal strings entered. ` +
            `${segments.filter((s) => s.state === 'match').length} match at least one header, ` +
            `${segments.filter((s) => s.state === 'noMatch').length} match none.`}
      </p>
    </div>
  );
}
