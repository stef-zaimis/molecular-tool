# UI_NOTES.md

Notes for the Electron/React desktop UI in `desktop/`.

This file is where the uncertainty lives. If something in the code looks like a
decision, the reasoning is here; if something is unresolved, it is in §10 rather
than hidden behind a plausible-looking implementation.

Scope of this phase: **UI only**. No Python was modified, no backend transport
exists, and no analysis can be executed. `python main.py` still launches the
existing Tkinter application, unchanged.

Reference material: `docs/02-project-creation.svg` (newest, authoritative for
the project-creation screen), plus `docs/01-launcher.png`,
`docs/02-project-creation.png`, `docs/03-molecular-diagnosis.png`. Behavioural
reference for the existing Python: `REPO_MAP.md` at the repo root.

> **Pass 2 (interaction semantics)** reworked navigation gating, focal-set
> editing, the DNC steppers, tooltip stacking, and replaced every hand-drawn
> icon with artwork extracted from the design SVG. Sections marked *(pass 2)*
> describe current behaviour; anything contradicting them is stale.

---

## 1. Visual assumptions

**Measured, not estimated.** `docs/` contained only the three PNGs — no exported
SVGs, no CSS, no font files. Every token in `src/styles/tokens.css` was obtained
by sampling pixels in the references (run-length scans along rows/columns and
ink bounding boxes), then verified by screenshotting the built app at the exact
reference dimensions and diffing the same measurements.

| Assumption | Evidence |
|---|---|
| The palette is exactly three slates plus white and black. | `#525F72` (body / active tab), `#2D3541` (title bar, dark buttons, tooltips, inactive tabs), `#3D4757` (inputs, panels, checkboxes, tab strip, alignment area) account for 97–99% of every reference image. |
| One type size everywhere. | Character advance measures 16.2px on every label, tab, input and button in screens 2 and 3. JetBrains Mono advance is 0.6em, so 16.2 / 0.6 = **27px**. Launcher buttons measure 14.9px advance → **25px**. |
| Reference frames are 1×, not 2×. | The launcher is 540×289; at 2× that would be a 270×144 window, which is not a plausible window size. |
| The hard shadow is offset **left and down**. | Launcher buttons: 4px of pure black at x157–160 (left of the face) and y160–163 (below it). Not the usual bottom-right drop shadow. |
| Window controls and the hamburger are `#525F72`, not white. | Sampled directly. They are the body colour on the dark title bar — the same low-contrast treatment used for inactive tab labels. |
| Inactive/disabled tabs use body-coloured text on `#2D3541`. | Sampled at x=600 on screen 3. This is why unselected analyses read as unavailable. |
| The workspace title is centred between the hamburger and the window controls, not on the window. | Title ink centre is x=915 on a 1920 frame; the midpoint of the gap between the hamburger (ends x=53) and the controls (start x=1778) is x=915.5. |
| Token colours in the STRING area. | Green `#82B987`, red `#D39DA7`, sampled from the two tokens drawn in screen 3. |
| The alignment panel bleeds off the right window edge. | It has a 67px top and 68px bottom inset but touches x=1920. Reproduced as-is; presumably the future viewer is wider than the window. |

**Verification method.** `desktop/` has no screenshot tooling committed. Captures
were taken with a throwaway Electron harness using offscreen rendering and
`force-device-scale-factor=1`, because this display runs at 125% scaling and an
on-screen 1920×1080 window neither fits nor composites. If you want this
repeatable, see §10 Q12.

---

## 2. Navigation and state assumptions

- One reducer + Context (`src/app/state/`), no state-management dependency. The
  state is small and session-local, so Redux/Zustand would not earn their weight.
- `Launcher → Project Creation → Workspace` are screens within a single window,
  not routes. There is no router.
- **The project draft survives all navigation in a session.** Going back to the
  launcher and forward again keeps the name, FASTA path, analyses, parameters and
  focal strings. Nothing is persisted to disk.
- Home is the Project page. From the workspace it navigates to Project Creation;
  on Project Creation it is the current location and does nothing. **Superseded
  pass 1 behaviour:** home used to fall back to the launcher from the project
  page, which made it read as a Back button. The launcher is now unreachable
  once a project exists — see §10 Q15.
- `setActiveAnalysis` is guarded **in the reducer**, not just on the button, so an
  unselected analysis cannot be activated by any code path.
- If you untick the analysis you are currently viewing, the active tab falls back
  to the first still-enabled analysis rather than showing an empty workspace.

---

## 2b. Interaction semantics (pass 2)

**Navigation**
- The house icon is the **Project page**, not a Back button. On project
  creation it is the current location and is inert (`aria-current="page"`, no
  click handler). From an analysis page it navigates to project creation.
- **Analysis tabs are generated from the project's analysis selection.** An
  unselected analysis has no tab at all — the previous build rendered all three
  and disabled two. `selectedAnalyses()` in `projectState.ts` is the single
  source, consumed by the workspace shell.
- **An analysis page cannot be entered without a FASTA.** `canEnterWorkspace()`
  requires both a chosen file and Molecular Diagnosis selected; Continue is
  disabled and its tooltip names whichever gate is unmet.

**Focal set**
- Stored as an **array of strings**. The `;` between entries is a rendered
  separator, never the storage format, and is rejected inside an entry.
- `+` and `-` are **mutually exclusive modes** (`role="radiogroup"`). The design
  already distinguishes them by tone, so the selected disc is chrome with a
  white symbol and the unselected one recedes into the surface colour.
- **Enter applies the current mode**: add appends; remove deletes an exact
  entry. A removal that matches nothing shows an inline message and changes
  nothing.
- **Undo/redo cover focal-set edits**, one add/remove per step, held in app
  state rather than in the editor component, because the controls that mutate
  the set live outside it. A new edit clears the redo branch. No-ops (empty add,
  duplicate add, missed removal) record no history entry.
- The **display is read-only**, which the brief explicitly permits. It was a
  textarea whose text was re-tokenised on every keystroke, which made a joined
  string the de-facto source of truth; entries now change only through `+`, `-`,
  undo and redo. Text can still be selected and copied.
- The **download icon exports** the set to a user-chosen `.txt`, one entry per
  line, through a save dialog in the main process.

**Focal validation against the real FASTA**
- Colours no longer come from a fixture. On choosing a file, the main process
  reads **only its header lines** (`fasta:read-headers`) and the renderer
  validates each entry against them: green = the literal string occurs in at
  least one header, red = it occurs in none, neutral = no file loaded yet.
- Matching is unchanged: **case-sensitive literal substring**
  (`header.includes(entry)`), matching Python's `target_string in header`. No
  fuzzy matching was introduced.
- This is a header read, not a FASTA parser: no residues are read, no alignment
  is validated, nothing is computed, and no Python was touched. See
  `FastaHeaderRequest` in `contract.ts`.

---

## 3. Figma controls whose behaviour was unclear

Implemented visually; behaviour is either minimal or deferred. None of these
invent semantics.

| Control | Where | What it does now |
|---|---|---|
| Hamburger (title bar) | Screens 2, 3 | Renders and has hover/focus states. **No menu** — there is no menu content in the design. See Q7. |
| `Save` next to FOCAL SET TITLE | Screen 3 | Shows a "not available yet" notice. The title is already held in session state; there is nothing to save to. |
| Download icon (left of the focal box) | Screen 3 | **Resolved in pass 2**: exports the focal set as a plain `.txt`, one entry per line. |
| `+` / `-` (STRING row) | Screen 3 | **Resolved in pass 2**: mutually exclusive add/remove *modes* for focal strings, not focal-set management. Adding or removing a whole focal set now has no UI — see Q16. |
| Undo / redo arrows | Screen 3 | **Resolved in pass 2**: undo/redo for focal-set edits, one add/remove at a time. Redo is dimmed until something has been undone, matching the design. |
| Resize grip (focal box) | Screen 3 | Real: the box is vertically resizable. |
| Scrollbar rail left of the alignment panel | Screen 3 | **Presentational only.** The real scrollbar belongs to the future viewer. |
| Multiple FASTA files (new SVG) | Screen 2 | **Not implemented, deliberately.** The design shows a list of FASTA files with its own undo/redo; the app stays single-FASTA because the scientific meaning of several files is undefined. See Q13. |
| Saved-title state with CREATE / RENAME (new SVG) | Screen 2 | Not implemented — it implies project persistence. See Q14. |
| RECENTS / BROWSE under Open Existing (new SVG) | Screen 1 | Not implemented; still blocked on there being no project format. |

---

## 4. Deviations from the mockups, and why

1. **A `Continue` button was added to Project Creation.** The mockup has no way
   to leave the screen, but the flow requires one. It is bottom-right, styled in
   the same language as the other dark action buttons, and disabled unless
   Molecular Diagnosis is selected (it is the only workspace that exists). If the
   design has a different intended affordance, this is the first thing to replace.

2. **"Ignore gaps" is unticked by default; the mockup shows it ticked.** The
   mockup is a presentation state, not a specification of defaults. The frontend
   default matches the Python default (`include_gappy_consensus_dmc_sites=False`,
   REPO_MAP §7) so the UI cannot silently change analysis behaviour. Say the word
   and it flips.

3. **Interaction states are invented.** The mockups only show resting states
   (plus one open tooltip). Hover, active, focus, and disabled colours are derived
   from the palette by small lightness steps and live in `tokens.css`. The
   keyboard focus ring (`--c-focus-ring`) has no counterpart in the design at all
   — it is required for a usable desktop app.

4. **Icons are the design's own artwork (pass 2).** `src/components/icons/Icons.tsx`
   is generated from `docs/02-project-creation.svg`: each glyph is pulled out by
   element id, its inherited group transforms are composed onto the shapes, and
   its viewBox is set to the artwork's own bounding box. Nothing is redrawn, and
   the file says not to hand-edit the path data. Covered: home, menu, tick,
   question mark, undo, redo, plus, minus, DNC stepper, download, minimise,
   maximise, close, resize grip.
   Two-tone glyphs (help, plus, minus, download) paint their disc with
   `--icon-bg` and their symbol with `--icon-fg`, which is how selected vs
   unselected mode is expressed without swapping artwork.
   The design has **no separate restore glyph**, so a maximised window shows the
   same maximise artwork. Window-control hover colours are taken from the
   design's own "negative" variants (close `#A35365`, maximise `#507454`,
   minimise `#886750`).

5. **The tooltip is 317px wide with a left-pointing notch**, matching the one
   open tooltip in screen 2. Only one tooltip style exists, so this is applied to
   every help dot.

6. **The alignment placeholder shows the text "Alignment Viewer"** plus the
   selected file name. The reference panel is empty. Per the brief the area
   should "look intentional, not like unfinished white space", so it carries a
   recessive (42% opacity) label. Remove it when the viewer lands.

7. **Numbers in the spinners are visible in our build.** The design's spinners
   appear empty. Ours show the actual values (1 and 2, the Python defaults),
   because an empty numeric control with no value would be misleading.
   Each field has exactly one stepper: the design's up/down pair is a single
   piece of artwork with two transparent hit areas over its halves. The previous
   build rendered that glyph twice and cropped each copy, which read as a
   duplicated control.

8. **A notice bar appears for unfinished actions.** Not in the design. It only
   renders after the user activates such a control, so the resting UI stays clean.

9. **Additions with no counterpart in the design (pass 2):** a one-line FASTA
   header-read status under the file field ("14 sequence headers read", or the
   read error); a floating inline validation message under the STRING field; and
   a centred "not built yet" panel for a selected analysis that has no
   workspace. All three exist because the new gating would otherwise fail
   silently.

10. **Tooltips are portalled into `<body>`** with fixed positioning, so no
    ancestor `overflow` can clip them and the alignment pane cannot paint over
    them. They flip to the left of the dot when they would run off the window.

**Measured fidelity after iteration** (reference vs. build, same measurement
method): title-bar and tab-strip bands exact; input positions and heights within
1px; checkbox rows within 1–2px; spinner rows exact; window-control ink within
1px; home glyph within 1px. Largest remaining discrepancy is text run width —
our strings render 5–6px wider over a 380px run (~1.5%), which is font rasteriser
difference, not layout.

---

## 5. Frontend naming for backend options

The brief called out one naming problem specifically. Resolution:

| UI label (from mockup) | Frontend state field | Python keyword | Polarity |
|---|---|---|---|
| `Ignore gaps` | `ignoreGaps` | `include_gappy_consensus_dmc_sites` | **Same** — no inversion |
| `Give BoTD to ambiguous bases` | `giveBenefitOfDoubtToAmbiguousBases` | `include_ambiguous_dmc_bd` | **Same** |
| `Min. candidate DNC size` | `minCandidateSize` | `min_combination_length` | Same |
| `Max. candidate DNC size` | `maxCandidateSize` | `max_combination_length` | Same |

The two names "ignore gaps" and "include gappy consensus DMC sites" sound like
opposites but describe the same switch from different ends: setting it true means
gap-containing focal columns are *not* disqualified, i.e. gaps are ignored when
deciding the consensus, i.e. gappy consensus sites are included. **`ignoreGaps ===
include_gappy_consensus_dmc_sites`.** There is no `!` anywhere in the mapping, and
no double negative is carried in frontend state. The adapter should be a
straight rename.

**`min_combination_length` is a trap and the frontend does not paper over it.**
Per REPO_MAP §4.4/§10.1 it is *not* a floor on which combination sizes are
searched — the search always starts at size 1, and this value only gates whether
the search may stop early. The tooltip says so explicitly, and `contract.ts`
carries a warning comment. Do not "fix" this in the frontend.

---

## 6. Backend gaps

Full detail is in `src/contract.ts`, marked `// BACKEND GAP:`. Summary:

1. **No project persistence.** No format, no reader, no writer anywhere in the
   repo. Blocks "Open Existing Project" and `Save`.
2. **One focal string, not many.** Python takes a single `target_string`. The
   frontend models a list because the design shows multiple tokens. Unresolved —
   see Q4. **Not** to be fixed by editing Python in this task.
3. **Results are files, not data.** `run_pipeline_core` fuses computation and
   output writing and returns paths; `FiveSiteResult` and the consensus result are
   computed, written, then dropped (REPO_MAP §6.6). Showing results in the UI
   requires splitting computation from serialisation in Python.
4. **No progress, no cancellation.** The Python core is one synchronous call with
   no callbacks and runs on the Tk main thread (REPO_MAP §12). The workspace has
   nowhere to show progress and no way to stop a run. Note the 5-site search is
   `C(n,5) × sequence count` with no cap (REPO_MAP §10.4), so this matters.
5. **Untyped errors.** Python raises bare `ValueError`s with prose messages; the
   existing Python tests match on message text. `AppErrorCode` in `contract.ts`
   is the target shape.
6. **No alignment transport.** Nothing exposes parsed FASTA. The real alignment is
   2354 × 736 ≈ 1.7M residues, so this needs a binary channel, not JSON.
7. **No consensus entry point.** Consensus runs only as a side effect of
   `run_pipeline_core`. Selecting it alone is not currently executable.
8. **Punishment takes no parameters.** Its thresholds are module constants.
9. **`output_dir` has no UI.** It is a required Python argument and appears
   nowhere in the mockups. See Q8.

---

## 7. Selected-analysis tabs — explicit assumption

**As instructed, recorded explicitly:** the analyses ticked during project
creation determine which workspace tabs are enabled. Unselected analyses render
as disabled tabs using the reference's low-contrast treatment, carry
`aria-selected`/`disabled`, and show a title explaining why.

Consequences and open edges:
- Only Molecular Diagnosis has a screen in this phase. `Continue` is disabled
  unless it is selected.
- No execution ordering is implied or implemented for multiple selected analyses.
- Whether analyses should be changeable *after* project creation is undecided —
  currently they cannot be, because the creation screen is reachable again via
  home, but nothing prevents editing there and returning.

---

## 8. Open Existing Project

Reproduced visually. Activating it shows a notice explaining that this build has
no project file format and offering "create a new project" instead. It does not
open a file picker, because there is no file type to pick and offering one would
imply a format that does not exist.

`OpenProjectRequest` / `OpenProjectResponse` / `PersistedProject` in
`contract.ts` describe what a real implementation needs, including two
unresolved questions: whether the FASTA is referenced or copied, and whether a
checksum is stored to detect the alignment changing under a saved project.

---

## 9. Focal-string semantics

Implemented for this phase:
- Tokens are separated by `;` in the editor, matching the reference.
- Empty tokens are dropped at tokenisation. This is deliberate: an empty token
  would substring-match *every* header.
- Validation is **case-sensitive substring containment** (`header.includes(token)`),
  which is exactly what the Python does (REPO_MAP §4.2). It is **not** exact-ID
  matching.
- Green = matches ≥1 known header, red = matches 0, neutral = cannot be validated.
- **Neutral is distinct from red on purpose.** "No alignment loaded" must not look
  like "your string is wrong".
- All matching lives in `focalMatching.ts` behind a `FocalMatchMode` parameter, so
  swapping to case-insensitive / exact-ID / field-aware / regex does not touch the
  editor. Those modes are implemented frontend-side but are **not** backed by
  Python.

Unresolved: see Q4 (multiple strings), Q5 (focal sets), Q6 (union vs intersection).

**Mock data:** validation runs against 18 representative headers in
`src/fixtures/mockHeaders.ts`, shaped like the real file
(`AACTA5253-20|AU|Leptacis`, with a `Leptacis_tipulae` subset). Nothing parses
FASTA anywhere in this app.

---

## 10. Questions for you, before backend integration

1. **`DNC` or `DMC`?** The mockup labels the spinners "Min./Max. candidate **DNC**
   size", but the codebase, REPO_MAP and the domain all use **DMC** (diagnostic
   molecular characters). I reproduced the mockup's spelling verbatim rather than
   guess, and it is a single string in `MolecularDiagnosisScreen.tsx`. Which is
   correct?

2. **Where does the home button go from the workspace?** I send it to Project
   Creation. Should it go to the launcher instead, or open a project overview
   that does not exist yet?

3. **What does the import icon next to the STRING box load?** I assumed "focal
   strings from a file". It could equally be "import a header list" or "load a
   saved focal set".

4. **Multiple focal strings — what do they mean scientifically?** The design shows
   several tokens, but Python accepts exactly one `target_string`. Do multiple
   strings mean (a) the union of matching headers as one focal set, (b) several
   independent focal sets analysed separately, or (c) something else? This changes
   the backend API and I did not want to guess at anything with scientific
   consequences.

5. **What are the `+` / `−` buttons adding and removing?** I assumed focal sets.
   If they add/remove *strings*, the STRING row and the token area overlap in
   purpose. Also: if multiple focal sets exist, how does the user switch between
   them? There is no visible switcher in the mockup.

6. **Union or intersection for multiple tokens?** If two strings both match a
   header, is it in the focal set once (union) or must a header match all tokens
   (intersection)? `unionMatchCount` currently assumes union.

7. **What is in the hamburger menu?** Nothing in the mockups indicates its
   contents.

8. **Where do output files go?** `output_dir` is required by Python and absent
   from the design. Options: a project directory chosen at creation, a fixed
   subfolder next to the FASTA, or a per-run picker.

9. **Should "Ignore gaps" default to ticked?** The mockup shows it ticked; I kept
   the Python default (unticked). This is a scientific default, so it is yours.

10. **Is `Save` on the focal set meant to persist to disk, or just to confirm the
    name within the session?**

11. **Should the analysis selection remain editable after project creation?**

12. **Do you want the screenshot harness committed?** It currently lives in a
    scratch directory. Committing it under `desktop/tools/` would make the
    visual-diff loop repeatable in CI, at the cost of a dev dependency on
    offscreen Electron. The same directory holds the icon-extraction script that
    generates `Icons.tsx` from the design SVG — that one is arguably more
    valuable to keep, so the icons can be regenerated when the design changes.

13. **What do multiple FASTA files mean?** The new SVG shows a list of them with
    its own undo/redo. Left unimplemented on purpose, per the brief. Are they
    merged into one alignment, analysed independently, or is one of them
    "active"? This is a scientific question, not a UI one.

14. **What does the project-title `Save` / CREATE / RENAME state do?** The new
    SVG's second project-creation frame shows a saved title with CREATE and
    RENAME buttons, which implies project persistence. Not implemented, since
    there is still no project format.

15. **Is there meant to be any route back to the launcher?** Home is now the
    project page, as instructed, so the launcher becomes unreachable once a
    project is created. No affordance for it exists in the design either.

16. **Should adding/removing whole focal sets still be possible?** `+`/`-` are
    now string modes, so the multi-focal-set actions they previously performed
    have no UI. The state model still supports several sets.

17. **Should a duplicate add be treated as a no-op?** It currently is, by
    symmetry with a missed removal: nothing changes, so no history entry is
    recorded and an inline message explains why.

---

## 11. Figma exports that would improve fidelity

Most valuable first:

1. **The three window-control glyphs, the hamburger, and the home icon as SVG.**
   The home glyph in particular has fine detail (the roof appears to be an
   outlined triangle with an inner fill and small eave gaps) that I approximated
   with a single filled triangle.
2. **The import, undo, redo and `+`/`−` glyphs as SVG.** These are the least
   pixel-faithful of the set.
3. **The checkbox tick as SVG**, including exactly how far it overhangs the box.
4. **Hover / active / focus / disabled states for buttons, tabs, checkboxes and
   the help dot.** Everything interactive is currently my interpretation.
5. **Exact letter-spacing and font weights per text style.** I measured advance
   width (which shows no extra tracking) but weights are inferred: bold for the
   title and the launcher's first line, regular elsewhere.
6. **The tooltip's full spec** — padding, max width, line height, arrow size, and
   what the *active* help dot looks like (in screen 2 it appears white and larger
   than the resting dots).
7. **A `Continue`/primary-button spec**, if one is intended (see §4.1).
8. **The intended empty state for the alignment panel.**
9. **Any spec for the alignment viewer's own scrollbar**, since the rail in
   screen 3 is currently presentational.
10. **A frame at a smaller width** (e.g. 1440 or 1280) showing what should
    compress. Right now the left panel is fixed at 1049px and the alignment area
    absorbs all slack; below ~1180px the viewer area gets very narrow.

---

## 12. Things deliberately NOT done

Per the brief: no alignment renderer, no Canvas, no WebGL, no FASTA parsing in
TypeScript, no mock nucleotide matrix, no FastAPI, no localhost HTTP, no Python
bridge, no project persistence, no analysis execution, no component library, and
no changes to any Python file. The known-stale Python test
(`tests/test_core.py::test_find_dmc_information`) was left failing exactly as
REPO_MAP §9 documents — Python baseline before and after this work is
**1 failed, 18 passed**.
