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
 * exact region the viewer will occupy, inside the frame the viewer draws around
 * it, so the workspace composition is already right when the real component
 * drops in.
 */
export function AlignmentPlaceholder({ caption }: AlignmentPlaceholderProps): JSX.Element {
  return (
    <section className="alignment" aria-label="Alignment viewer">
      {/*
        No rail here. The viewer's left rail is ONE element, owned by
        WorkspaceRightPane because it is also the drag handle; drawing a second
        one here produced the doubled edge the reference does not have.
      */}
      <div className="alignment__panel">
        <div className="alignment__empty">
          <p className="alignment__title">Alignment Viewer</p>
          {caption && <p className="alignment__caption">{caption}</p>}
        </div>
      </div>
    </section>
  );
}
