import './AlignmentPlaceholder.css';

interface AlignmentPlaceholderProps {
  /** Shown beneath the title, e.g. the selected file name. */
  readonly caption?: string;
}

/**
 * Reserved space for the alignment viewer.
 *
 * DELIBERATELY NOT A RENDERER. No canvas, no WebGL, no FASTA parsing and no
 * fabricated nucleotide matrix — the viewer is its own task, and filling this
 * area with fake residues now would make it harder to tell real progress from
 * scaffolding later.
 *
 * What it does provide is the correct visual footprint: the panel occupies the
 * exact region the viewer will occupy, including the vertical scrollbar rail
 * that sits to its left in the reference, so the workspace composition is
 * already right when the real component drops in.
 */
export function AlignmentPlaceholder({ caption }: AlignmentPlaceholderProps): JSX.Element {
  return (
    <section className="alignment" aria-label="Alignment viewer">
      {/* Rail only — the real scrollbar belongs to the future viewer. */}
      <div className="alignment__rail" aria-hidden="true">
        <div className="alignment__rail-thumb" />
      </div>

      <div className="alignment__panel">
        <div className="alignment__empty">
          <p className="alignment__title">Alignment Viewer</p>
          {caption && <p className="alignment__caption">{caption}</p>}
        </div>
      </div>
    </section>
  );
}
