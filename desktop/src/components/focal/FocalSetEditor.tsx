import { useCallback, useEffect, useRef, useState } from 'react';
import { Compartment, EditorState, StateEffect, StateField } from '@codemirror/state';
import { Decoration, EditorView, ViewPlugin, keymap, placeholder } from '@codemirror/view';
import type { DecorationSet, ViewUpdate } from '@codemirror/view';
import { DownloadIcon, RedoIcon, UndoIcon } from '../icons/Icons';
import { IconButton } from '../controls/Controls';
import { scanFocalText, textMatchesHeaders, serialiseFocalText } from './focalText';
import { currentZoom } from '../../app/uiScale';
import type { FocalPresenceState, HeaderPresencePayload } from '../../backendContract';
import './FocalSetEditor.css';

/**
 * The focal-set field: a real editor, coloured per entry.
 *
 * Built on CodeMirror 6 rather than a `contenteditable` div. Per-token colour
 * on top of genuinely editable text needs correct caret, selection, IME,
 * clipboard and undo behaviour, and hand-rolling that on `contenteditable`
 * produces exactly the class of bug that is impossible to finish. CodeMirror's
 * decoration model gives the colours without any of it.
 *
 * What it is NOT: a query field. Text typed here is COMPLETE EXACT HEADERS
 * separated by `;`. Nothing typed into this box expands a substring — that is
 * what `+` is for, and conflating the two is how a focal set ends up with
 * members nobody chose.
 *
 * The document is the source of truth while the user is typing; the membership
 * array is derived from it. Incoming arrays (from `+`, `−`, undo/redo or a
 * save response) are pushed back in only when they express DIFFERENT
 * membership, so ordinary typing is never interrupted by a re-serialisation
 * that would move the caret.
 */

/** The design's box height, and the range the grip may drag it through. */
const DEFAULT_HEIGHT = 143;
const MIN_HEIGHT = 96;
const MAX_HEIGHT = 420;

type EntryTone = 'match' | 'elsewhere' | 'missing' | 'unknown' | 'pending' | 'duplicate';

const TONE_BY_PRESENCE: Record<FocalPresenceState, EntryTone> = {
  present_current: 'match',
  present_other: 'elsewhere',
  missing: 'missing',
  unknown: 'unknown',
};

const MARKS: Record<EntryTone, Decoration> = {
  match: Decoration.mark({ class: 'cm-focal cm-focal--match' }),
  elsewhere: Decoration.mark({ class: 'cm-focal cm-focal--elsewhere' }),
  missing: Decoration.mark({ class: 'cm-focal cm-focal--missing' }),
  unknown: Decoration.mark({ class: 'cm-focal cm-focal--unknown' }),
  pending: Decoration.mark({ class: 'cm-focal cm-focal--pending' }),
  duplicate: Decoration.mark({ class: 'cm-focal cm-focal--duplicate' }),
};

/** Read-only is reconfigured in place, so locking a set never remounts the editor. */
const editableCompartment = new Compartment();

/** Presence lives in editor state so a re-check repaints without a reconfigure. */
type PresenceMap = Readonly<Record<string, HeaderPresencePayload>>;

const setPresence = StateEffect.define<PresenceMap>();

const presenceField = StateField.define<PresenceMap>({
  create: () => ({}),
  update(value, transaction) {
    for (const effect of transaction.effects) {
      if (effect.is(setPresence)) return effect.value;
    }
    return value;
  },
});

function toneFor(text: string, isMember: boolean, presence: PresenceMap): EntryTone {
  if (!isMember) return 'duplicate';
  const verdict = presence[text];
  // No verdict yet is NOT "absent": it renders neutral, so a set never flashes
  // red on its way to being checked.
  return verdict ? TONE_BY_PRESENCE[verdict.state] : 'pending';
}

function buildDecorations(view: EditorView): DecorationSet {
  const presence = view.state.field(presenceField);
  const text = view.state.doc.toString();
  return Decoration.set(
    scanFocalText(text)
      .filter((span) => span.to > span.from)
      .map((span) => MARKS[toneFor(span.text, span.isMember, presence)].range(span.from, span.to)),
    true,
  );
}

const focalHighlighting = ViewPlugin.fromClass(
  class {
    decorations: DecorationSet;

    constructor(view: EditorView) {
      this.decorations = buildDecorations(view);
    }

    update(update: ViewUpdate) {
      const presenceChanged = update.transactions.some((transaction) =>
        transaction.effects.some((effect) => effect.is(setPresence)),
      );
      if (update.docChanged || update.viewportChanged || presenceChanged) {
        this.decorations = buildDecorations(update.view);
      }
    }
  },
  { decorations: (plugin) => plugin.decorations },
);

interface FocalSetEditorProps {
  /** Working membership. The editor derives its document from this. */
  readonly headers: readonly string[];
  readonly presence: PresenceMap;
  readonly readOnly: boolean;
  /**
   * Called with the membership a manual edit produced.
   *
   * `typing` marks a continuous manual burst, which the caller coalesces into
   * ONE history step. `+`/`−` arrive through their own actions and stay one
   * step each.
   */
  readonly onChange: (headers: readonly string[], typing: boolean) => void;
  /** Fired when the field loses focus, so a typing burst can be closed off. */
  readonly onEditingEnd: () => void;
  /**
   * Bumped by the caller at a CANONICALISATION BOUNDARY — a save, a `+`/`−`,
   * or switching focal sets.
   *
   * Between boundaries the document is left exactly as typed, so a half-written
   * entry or a trailing `; ` survives a keystroke. At a boundary the text is
   * rewritten from the canonical membership, which is what removes a struck-out
   * duplicate that is not actually stored. Rewriting on every keystroke instead
   * would fight the user mid-word.
   */
  readonly canonicalToken: number;
  readonly onUndo: () => void;
  readonly onRedo: () => void;
  readonly canUndo: boolean;
  readonly canRedo: boolean;
  readonly onExport: () => void;
  readonly canExport: boolean;
}

export function FocalSetEditor({
  headers,
  presence,
  readOnly,
  onChange,
  onEditingEnd,
  canonicalToken,
  onUndo,
  onRedo,
  canUndo,
  canRedo,
  onExport,
  canExport,
}: FocalSetEditorProps): JSX.Element {
  const hostRef = useRef<HTMLDivElement>(null);
  const boxRef = useRef<HTMLDivElement>(null);
  const viewRef = useRef<EditorView | null>(null);
  /** Height is component state; see the CSS note about the native grip. */
  const [height, setHeight] = useState(DEFAULT_HEIGHT);
  // Kept in refs so the editor is constructed once: rebuilding it on every
  // prop change would destroy the caret and the selection.
  const onChangeRef = useRef(onChange);
  const onEditingEndRef = useRef(onEditingEnd);
  onChangeRef.current = onChange;
  onEditingEndRef.current = onEditingEnd;

  /**
   * Rewrite the document from the canonical membership.
   *
   * Only ever called at a boundary. A no-op when the text already reads that
   * way, so a caret sitting in settled text is not disturbed.
   */
  const canonicalise = useCallback(() => {
    const view = viewRef.current;
    if (!view) return;
    const canonical = serialiseFocalText(headersRef.current);
    if (view.state.doc.toString() === canonical) return;
    view.dispatch({
      changes: { from: 0, to: view.state.doc.length, insert: canonical },
    });
  }, []);

  const headersRef = useRef(headers);
  headersRef.current = headers;
  const canonicaliseRef = useRef(canonicalise);
  canonicaliseRef.current = canonicalise;

  useEffect(() => {
    const host = hostRef.current;
    if (!host) return;

    const view = new EditorView({
      state: EditorState.create({
        doc: serialiseFocalText(headers),
        extensions: [
          presenceField,
          focalHighlighting,
          editableCompartment.of([
            EditorView.editable.of(!readOnly),
            EditorState.readOnly.of(readOnly),
          ]),
          EditorView.lineWrapping,
          placeholder('Type headers separated by ; or use + above'),
          // CodeMirror's own history is deliberately absent: undo/redo belong
          // to the draft, so that `+`, `−` and manual edits share one stack
          // and the toolbar buttons drive all three.
          keymap.of([
            { key: 'Mod-z', run: () => true, preventDefault: true },
            { key: 'Mod-y', run: () => true, preventDefault: true },
            { key: 'Mod-Shift-z', run: () => true, preventDefault: true },
          ]),
          EditorView.updateListener.of((update) => {
            if (!update.docChanged) return;
            onChangeRef.current(
              scanFocalText(update.state.doc.toString())
                .filter((span) => span.isMember)
                .map((span) => span.text),
              true,
            );
          }),
          EditorView.domEventHandlers({
            blur: () => {
              // Leaving the field is a safe boundary: whatever was being typed
              // is finished, so the text is squared up with what is stored.
              onEditingEndRef.current();
              canonicaliseRef.current();
              return false;
            },
          }),
        ],
      }),
      parent: host,
    });

    viewRef.current = view;
    return () => {
      view.destroy();
      viewRef.current = null;
    };
    // Constructed once. `headers` seeds the initial document only; later
    // changes arrive through the reconciliation effect below.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  /*
   * Reconcile incoming membership into the document.
   *
   * Only when the membership genuinely differs: if the document already
   * expresses these headers, whatever the user has typed — a trailing `; `, an
   * entry half-finished — is left exactly as it is.
   */
  useEffect(() => {
    const view = viewRef.current;
    if (!view) return;

    const current = view.state.doc.toString();
    if (textMatchesHeaders(current, headers)) return;

    view.dispatch({
      changes: { from: 0, to: view.state.doc.length, insert: serialiseFocalText(headers) },
    });
  }, [headers]);

  /*
   * A caller-signalled boundary: save, `+`/`−`, or a different focal set.
   *
   * This is what guarantees the promise that after Save the visible text is
   * exactly the saved membership — no struck-out duplicate left behind for an
   * entry the database does not hold.
   */
  useEffect(() => {
    canonicalise();
    // Deliberately keyed on the token alone: `canonicalise` reads current
    // headers through a ref, so a header change on its own must not fire it.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [canonicalToken]);

  /* Presence is pushed in as an effect so decorations repaint in place. */
  useEffect(() => {
    viewRef.current?.dispatch({ effects: setPresence.of(presence) });
  }, [presence]);

  /* Read-only is a live property: locking a set must not remount the editor. */
  useEffect(() => {
    viewRef.current?.dispatch({
      effects: editableCompartment.reconfigure([
        EditorView.editable.of(!readOnly),
        EditorState.readOnly.of(readOnly),
      ]),
    });
  }, [readOnly]);

  const spans = scanFocalText(serialiseFocalText(headers));
  const answered = spans.filter((span) => presence[span.text] !== undefined).length;
  const present = spans.filter(
    (span) => presence[span.text]?.state === 'present_current',
  ).length;

  return (
    <div className="focal-editor">
      <div className="focal-editor__gutter">
        <IconButton
          label="Export focal set as a text file"
          onClick={onExport}
          disabled={!canExport}
        >
          <DownloadIcon size={30} />
        </IconButton>
      </div>

      <div className="focal-editor__box" ref={boxRef} style={{ height }}>
        <div
          ref={hostRef}
          className="focal-editor__cm"
          role="group"
          aria-label="Focal set entries, separated by semicolons"
        />

        {/*
          Our own grip, drawn in the design's idiom. VERTICAL ONLY: the box
          width is part of the measured layout, and the native CSS `resize`
          both allowed horizontal dragging and sat on top of the scrollbar.
        */}
        <div
          className="focal-editor__grip"
          role="separator"
          aria-orientation="horizontal"
          aria-label="Resize the focal set field"
          tabIndex={0}
          onPointerDown={(event) => {
            event.preventDefault();
            const startY = event.clientY;
            /*
             * Pointer coordinates are in PAINTED pixels and `height` is in
             * layout pixels, and the shell may be zoomed for a high-DPI window
             * (see app/uiScale.ts). Dividing by the zoom keeps the box under
             * the pointer instead of drifting away from it at 1.25x the rate.
             */
            const zoom = currentZoom();
            const startHeight = boxRef.current
              ? boxRef.current.getBoundingClientRect().height / zoom
              : height;
            (event.target as HTMLElement).setPointerCapture?.(event.pointerId);

            const move = (moveEvent: PointerEvent) => {
              const next = startHeight + (moveEvent.clientY - startY) / zoom;
              setHeight(Math.min(MAX_HEIGHT, Math.max(MIN_HEIGHT, next)));
            };
            const end = () => {
              window.removeEventListener('pointermove', move);
              window.removeEventListener('pointerup', end);
              window.removeEventListener('pointercancel', end);
            };
            window.addEventListener('pointermove', move);
            window.addEventListener('pointerup', end);
            window.addEventListener('pointercancel', end);
          }}
          onKeyDown={(event) => {
            if (event.key === 'ArrowUp') {
              event.preventDefault();
              setHeight((current) => Math.max(MIN_HEIGHT, current - 16));
            } else if (event.key === 'ArrowDown') {
              event.preventDefault();
              setHeight((current) => Math.min(MAX_HEIGHT, current + 16));
            }
          }}
          onDoubleClick={() => setHeight(DEFAULT_HEIGHT)}
          title="Drag to resize · double-click to reset"
        />
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
        {headers.length === 0
          ? 'The focal set is empty.'
          : answered < headers.length
            ? `${headers.length} focal entries. Checking where they are.`
            : `${headers.length} focal entries. ${present} present in the selected scope.`}
      </p>
    </div>
  );
}
