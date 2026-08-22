# UI_NOTES.md

Notes for the Electron/React desktop UI in `desktop/`.

This file is where the uncertainty lives. If something in the code looks like a
decision, the reasoning is here; if something is unresolved, it is in §10 rather
than hidden behind a plausible-looking implementation.

Scope: the UI **and** its integration with the Python backend. Molecular
Diagnosis now runs end to end through a child-process service; see §13.

Reference material: **`docs/new/all_pages-20Aug2026.svg` is authoritative** -
page 1 launcher, page 2 new project, page 3 project / FASTA list, page 4 a new
focal set, page 6 a saved focal set with the library. `docs/legacy/` (the three
PNGs, `02-project-creation.svg`, `all_pages.svg`) is HISTORY and must not drive
new work; sections written before pass 8 cite those paths and are left as
written. Behavioural reference for the Python: `REPO_MAP.md` at the repo root.

> **Pass 2 (interaction semantics)** reworked navigation gating, focal-set
> editing, the DNC steppers, tooltip stacking, and replaced every hand-drawn
> icon with artwork extracted from the design SVG.
>
> **Pass 3 (backend integration)** wired the frontend to Python over a
> child-process JSON service. Several pass-1 and pass-2 notes are now
> superseded — §13 is authoritative where they disagree.
>
> **Pass 10 (interaction corrections)** is §19. Where it disagrees with any
> earlier section, §19 is what the code does.
>
> **Pass 11 (run observability)** is §20: live progress during a run, and the
> diagnostics behind it.

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

**Implemented (pass 5).** A project is a directory holding `project.sqlite`
plus an `outputs/` folder, so the control opens a *directory* picker rather
than a file picker. Opening restores the project's linked FASTA files and
focal sets and immediately reports the live state of every linked file.

Both questions this section previously left open are now answered, by the
schema rather than by the UI:

* **Referenced, never copied.** `fasta_file` stores an absolute path. No FASTA
  content is ever written into SQLite, so a project stays small and the user's
  file remains the single source of truth.
* **Yes, a checksum is stored** — `indexed_sha256`, alongside
  `indexed_size_bytes` and `indexed_mtime_ns`. Size and mtime are the cheap
  detector; the hash is the proof. See §15.

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

---

## 13. Backend integration (pass 3)

**Architecture.** Electron main spawns `python -m molecular_diagnosis.service`
as a child process and talks to it in newline-delimited JSON over stdin/stdout.
The renderer never spawns anything and never sees the scientific code: it calls
`window.desktop.analysis.*`, which the preload forwards to main, which forwards
to the child. `desktop/src/backendContract.ts` holds the shapes that actually
cross that boundary (as opposed to `contract.ts`, which still describes the
wider future architecture).

stdout is protocol-only. The service redirects its own `sys.stdout` to stderr
so a stray `print` in any library cannot corrupt the stream; stderr is captured
by the bridge and attached to failures as diagnostics.

**Superseded notes**
- Focal green/red is no longer computed in TypeScript. It comes from
  `validateFocalStrings`, which runs the same matcher the analysis uses, so the
  colours cannot disagree with which sequences a run will select.
  `src/fixtures/mockHeaders.ts` is now a **test fixture only**.
- FASTA loading no longer uses a Node header reader. `loadFasta` parses and
  validates through `fasta_io.parse_fasta` / `validate_aligned_fasta`.
- `;` inside a focal entry is now **allowed** (pass 2 rejected it). The
  separator between entries is a rendered element, never a stored character, so
  an entry containing `;` is unambiguous. There is a test for it on both sides.
- Entering the workspace now also requires the FASTA to have **loaded and
  validated**, not merely to have been chosen. A ragged file blocks Continue.

**Focal-set semantics.** `molecular_diagnosis/focal.py` is the single matcher.
A header is focal when ANY selector occurs literally within it (OR), matching
is case-sensitive substring containment, selectors are trimmed and deduplicated
preserving order, and an empty selector is rejected rather than silently
matching everything. The four places that previously wrote
`target_string in header` by hand all route through it now.

**Additions with no counterpart in the design**
- A **Run Molecular Diagnosis** button and result panel. The design shows no way
  to start an analysis, and one is required to reach the pipeline at all. Styled
  in the existing language.
- A **continuation prompt**. When the search stops only because it hit the
  maximum candidate size, the panel offers to continue from the next size. The
  decision is the user's; the continuation state itself is an opaque blob
  produced by the backend and handed straight back, so a continuation resumes
  rather than restarts.

**Assumptions**
- **Output directory** is the folder containing the selected FASTA. Python
  requires one and the design has no control for it. Each run writes a fresh
  numbered set (`DMCs_output(2).txt`, ...) because `next_available_filename`
  never overwrites — so a three-step continuation leaves three sets of files.
  That is pre-existing pipeline behaviour, not new.
- **Interpreter resolution**: `MOLECULAR_TOOL_PYTHON`, else the project
  `.venv`, else `python`/`python3` on PATH. `MOLECULAR_TOOL_ROOT` overrides the
  repo root. Packaging a bundled interpreter is not addressed.
- **Multiple FASTA files** remain unimplemented, per the brief.

**Known gaps**
- No progress reporting and no cancellation: the service answers one request at
  a time and a long search cannot be interrupted. The UI stays responsive
  because the work is in another process, but it cannot show progress.
- `runSequencePunishment` exists in the service and is tested indirectly, but no
  UI reaches it.
- The legacy Tkinter `gui.py` had its two pipeline calls renamed to the new
  `focal_strings=` keyword so it keeps working while it is being sunset. No
  Tkinter architecture was ported into the new code.

---

## 14. Parity and workspace scrolling (pass 4)

**Parity verified.** `tests/test_parity.py` runs `tests/parity_driver.py` twice —
once against the current tree, once against a read-only `git worktree` of the
pre-integration commit — and diffs the results. With a single focal selector the
two agree on focal/non-focal selection, every `DMCResult` field, five-site
optimisation, consensus, the DMC and consensus text reports, the logical
contents of `comparison_output.xlsx`, and continuation/resume. Two further tests
drive `service.dispatch(...)` — the route Electron actually takes — against the
same baseline.

The workbook comparison is logical, not byte-wise: sheet names, cell values,
number formats, bold, fill colours, alignment, freeze panes, auto-filter and
column widths. The xlsx zip carries timestamps, so bytes are not stable.

**Search semantics were not touched.** `min_combination_length` still does not
raise the starting length, and the early stop still fires at the first
productive length at or above the minimum. There is a parity test asserting
exactly that (`test_parity_non_default_min_and_max`).

**Scrolling.** The Molecular Diagnosis page is now a two-pane workspace: the
analysis column (`.diagnosis__left`) scrolls on its own with 72px of bottom
padding, while the alignment pane beside it and the tab strip above it stay
put. The window itself never scrolls. Verified by measurement, not by eye: the
viewer's viewport position is unchanged at scroll top and scroll bottom, the
document never overflows, and a tall window produces no scrollbar at all.

**Scrollbar** is modelled on the `scroller` element in
`docs/02-project-creation.svg`: 8px wide (2.236 design units), fully rounded
pill, `#2D3541` track inset 6px top and bottom, `#7F92AE` thumb, with hover and
active tints on the same hue. `scrollbar-width` / `scrollbar-color` are
deliberately unset — in current Chromium they override the
`::-webkit-scrollbar` pseudo-elements and would restore the default scrollbar.

A 10px transparent right border insets the scrollbar from the pane edge so it
does not sit flush against the alignment viewer's own rail; the two are
inverses of each other and read as one bar when they touch. The viewer rail
itself was left alone — it belongs to the future viewer.

**Output directory remains different from Tkinter.** Tkinter asked for an output
directory (Browse, or "Use FASTA Location"); the new UI always uses the folder
containing the FASTA, because the current design has no output-directory
control. This is a deliberate deviation, not a parity failure, and was not
redesigned in this pass.

---

## 15. Persistent projects and source status (pass 6)

### The renderer holds no database

There is no SQL anywhere in `desktop/`. The renderer cannot open, query or
migrate the project database; it calls named operations
(`window.desktop.project.*`) that the Python service performs against the
database it exclusively owns. The main process is equally ignorant — it
forwards, and holds no project state of its own.

### Source state is a live reading, not a stored property

A project stores a path. Between sessions the file behind that path can be
moved, edited or deleted, so every status the UI shows is recomputed from the
filesystem and nothing about availability is ever persisted. `sourceStatus.ts`
keeps three ideas apart, because conflating them is what produces a UI that
lies:

| Idea | Question it answers |
| --- | --- |
| `available` | Can this file be read right now? |
| `indexUsable` | Does the stored header index provably describe those bytes? |
| `activity` | What is the app doing about it at this instant? |

A file can be available with an unusable index (it was edited), and it can have
a usable index while a check is in flight. Neither implies the other.

### The six states, and how each is phrased

| Backend state | Shown as | Tone | Analysable | Offered action |
| --- | --- | --- | --- | --- |
| `current` | Current | ok | yes | — |
| `unverified` | Unverified | pending | yes | — |
| `stale` | Changed on disk | warning | no | Re-index |
| `never_indexed` | Not indexed | warning | no | Re-index |
| `missing` | Missing | error | no | Relink |
| `unreadable` | Unreadable | error | no | Relink |

Two of these deserve their reasoning written down:

* **`unverified` is not a warning.** Size and mtime still match; the file
  simply has not been hashed during this session. Opening a project is
  deliberately stat-only so that a project linking a multi-gigabyte alignment
  opens instantly. Painting that as a problem would train users to ignore the
  colour that means something.
* **`unverified` is still analysable.** A run reads and verifies the file
  anyway, so refusing to start would only mean hashing it twice.

While a file is being checked, re-indexed or relinked, the row reports the
*activity* ("Re-indexing") rather than the underlying state, so a large file
does not sit there labelled "Changed on disk" with nothing apparently
happening. Activity never makes a file analysable: work in progress is not a
reason to let a run start.

### When status is refreshed, and when it is not

Refreshes take the cheap path — one `stat` per linked file — at:

* project open (the response carries a full sweep, so a file that vanished
  while the app was closed is visible on first paint);
* regaining window focus;
* entering a screen or analysis tab.

The strong path, which hashes anything not yet proven this session, runs only
for a deliberate user action ("Verify now") and on the run path.

**Nothing in the focal editor triggers either.** Focal text is matched against
data already in memory, so typing never causes a stat sweep, let alone a
rehash of the alignment.

### Errors

An expected refusal — "this file is missing", "that focal entry is not in this
file" — travels as a code plus a sentence written for the user, and carries no
traceback. A traceback would describe where the service chose to refuse, which
tells nobody anything. Unexpected exceptions still travel with theirs.


## 16. Backend hardening before the project-backed UI (pass 7)

This pass changed no visual design. It closed correctness holes and completed
the service contract the Molecular Diagnosis UI is about to be built against,
so the renderer has something stable and honest to compile and reason against.

### Create and Open are now two different operations

`project.open` used to create a `project.sqlite` in whatever folder it was
given. That made the launcher's two buttons the same button: "Open Existing
Project" on the wrong folder produced an empty project that looked exactly
like lost work.

* `project.create` initialises one, and refuses a folder that already holds a
  project (`PROJECT_ALREADY_EXISTS`).
* `project.open` opens an existing one and never creates (`PROJECT_NOT_FOUND`),
  and never renames it.
* `project.setTitle` is the one way the title changes.

A failed open leaves the currently open project open. The context exposes
`createProject` / `chooseNewProjectDirectory` alongside `openProject` /
`chooseProjectDirectory`; wiring the creation screen to them is the next pass.

### The rules live in Python, not in disabled buttons

A locked focal set may be selected, presence-checked and analysed, but every
mutation — rename, delete, `+`, `-`, textbox replace — is refused with
`FOCAL_SET_LOCKED` no matter who calls it. Grey the controls out for clarity;
do not rely on that greying for correctness. The same applies to the empty
`+`/`-` query (`FOCAL_QUERY_EMPTY`) and to run gating.

### The big textbox is membership, not a query

`project.replaceFocalEntries` takes the raw `a ; b;c` box and makes the set
exactly that: trimmed around the separators, empties ignored, duplicates
collapsed onto their first occurrence, order preserved. A header that matches
no FASTA is KEPT so it can show red — it is not silently dropped.

It is applied as a diff, so surviving entries keep their ids and their cached
locations. The editor may debounce and resend freely; it is responsible for
ignoring stale responses, and the backend is responsible for the call being
deterministic and transactional.

### Colours must match what Run will allow

`project.focalPresence` now takes `fastaFileIds`. Pass the files a run would
actually read, and the states returned are the states Run gates on: green when
the entry is in something the run will read, red when it is in none of them,
grey/unknown when a required file is unavailable. Omitting the field asks about
the whole project, which is a different question and a different answer.

### Refusals the UI needs sentences for

`DUPLICATE_HEADER_IN_FILE` (the file repeats a header, so analysing it would
silently discard sequences — the file stays linked and searchable),
`FOCAL_ENTRIES_NOT_IN_SCOPE`, `FOCAL_PRESENCE_UNKNOWN`,
`SEARCH_SCOPE_UNAVAILABLE` (a `+` refused rather than persisting a partial
expansion). Plain search does NOT refuse: it returns its hits plus an
`unavailable` list of the files it could not read, and the UI must show that
before the user treats the hit list as complete.

## 17. Screens 1–3 on the new design, and the project-backed workspace (pass 8)

The authoritative design is now `docs/new/all_pages-20Aug2026.svg`. `docs/legacy/`
is history and must not drive new work. Screens 1–3 (launcher, new project,
project page) are implemented from the new file; the Molecular Diagnosis screen
keeps its current layout but has been moved onto the persistent project, so the
next pass is predominantly visual.

### One architecture, end to end

    open/create project -> linked FASTA sources -> selected FASTA scope
      -> working focal-set draft -> saved persistent focal set
      -> project.runMolecularDiagnosis

Nothing in the workspace reads `draft.fastaPath` any more, and neither
`analysis.validateFocalStrings` nor `analysis.runMolecularDiagnosis` is reachable
from it. Those preload methods still exist for the non-project entry points.

### Saved data vs working copy

A focal set in the database is a scientific commitment: it names the sequences a
run will treat as focal, and a run's report names the set. So editing does not
write. `FocalDraft` (app/state/focalDrafts.ts) holds `persistedId`, the saved
title/headers, the working title/headers, `locked` and the undo history;
`dirty` is DERIVED by comparing normalised working values against saved ones,
never stored.

Consequences the UI must keep:

* `project.saveFocalSet` is the only write ordinary editing performs.
* Switching drafts does not save the one being left.
* Run is refused while a draft is new or dirty, and never auto-saves. The
  reason is stated on the button, not left to be discovered.
* A locked set is read-only in the reducer AND in Python, and stays selectable
  and runnable.

### `+`, `−` and presence all run in Python

`+` is `project.searchHeaders` over the current scope, and refuses outright if a
file in that scope could not be searched — a partial expansion looks complete
and is not reproducible. `−` is `project.matchFocalHeaders` against the WORKING
copy; it never searches the FASTA, and it lives in Python because
`str.casefold()` is not `String.prototype.toLowerCase()` and two
implementations of "the same" rule would eventually disagree.

Colours come from `project.headerPresence`, which answers from the header index
with one batched lookup. It opens no file, hashes nothing, reindexes nothing and
persists nothing, so typing costs no filesystem work. The renderer debounces it
(~300 ms) and tags each request with a generation, dropping stale answers.
Presence is keyed BY HEADER, so reordering entries cannot mis-colour them.

Four states, not two: green (in the selected scope), orange (in the project but
not in the SELECTED file), red (in none of them), and a dotted-underlined grey
for `unknown` — "could not check" and "is not there" send the user to different
problems.

### FASTA pool

`FastaScope` is either one file or All files, defaulting to the first linked
file. It drives BOTH the presence colours and the files a run reads, so the two
can never describe different scopes. A single file runs with
`singleFile: true`; All files passes every linked id with `singleFile: false`.

### Screens

* **Launcher** — Create New / Open Existing, plus Recents and Browse. Browse is
  the folder picker; Recents needs a store of previously opened directories that
  nothing writes yet, and says so.
* **New project** — title and chosen FASTA paths are LOCAL; there is no database
  to write to yet. Create requires a title, picks a directory, calls
  `project.create`, then links each FASTA. A file that fails to link is reported
  and skipped: the project that was just created is not torn down over one bad
  input.
* **Project page** — persistent title with a pencil (`project.setTitle`), Browse
  for multi-file linking, the linked-file rows from `page2-b`, and the analysis
  checkboxes. Selected analyses are the tabs, so the temporary Continue button is
  gone. Home returns here.

Deliberately deferred, and drawn as disabled rather than faked: source rename,
source lock, and drag reordering. PIS is shown as `—` because nothing computes
it — a 0 would be a claim. Analysis selection is session-local; the schema has
nowhere to put it and a migration for a checkbox would be the wrong trade.


## 18. Screens 4 and 6: the editing workspace (pass 9)

### Three scopes, deliberately distinct

There are three, and conflating any two produces a specific bug:

| scope | what it decides | follows |
|---|---|---|
| run | which files an analysis reads | the FASTA pool |
| `+` search | which files a query searches | the FASTA pool |
| presence comparison | which files a header is compared against | ALWAYS every linked file |

With one FASTA selected, `+` searches only that file — so a header it adds is
normally green, and `+` is NOT widened to the whole project just to make orange
possible. Orange comes from the comparison scope: a header typed by hand (or
inherited from another selection) that is absent from the selected file but
present elsewhere in the project reads orange, not red. Changing the selection
re-colours the existing entries accordingly.

Under **All files** the run, search and comparison scopes coincide and there is
no "outside the selected file", so `present_other` cannot occur and orange never
appears.

### `+` must never truncate

`project.searchHeaders` is a CAPPED preview. Expanding `+` through it silently
produced a smaller focal set than the user asked for. `+` now resolves through
`project.resolveFocalAddQuery`: same substring rule, same ordering, same scope
verification, no limit, no write. The capped preview stays for listing UI.

### The run gate lives in the handler

`runDiagnosis` re-checks `runGate.canRun` before calling the backend and
surfaces its reason. A disabled button is a hint; the guarantee has to survive a
keyboard path or a stale render.

### The focal field is a real editor

CodeMirror 6 (`@codemirror/state` + `@codemirror/view`, no language modes, no
history extension). Per-token colour on top of genuinely editable text needs
correct caret, selection, IME, clipboard and undo behaviour; hand-rolling that
on `contenteditable` is the class of bug that never finishes.

* entries are complete exact headers separated by `;`
* whitespace around a separator is dropped, whitespace INSIDE an entry is kept
* blanks vanish, duplicates collapse to the first occurrence (repeats are shown
  struck through — the text is real, the second membership is not)
* manual editing NEVER expands a substring; `+` is the thing that searches
* the document is the source of truth while typing, and an incoming array only
  replaces it when it expresses different membership, so the caret is never
  yanked mid-entry

Undo/redo belong to the DRAFT, not to the editor: CodeMirror's own history is
absent and `Mod-Z` is swallowed, so one stack covers manual typing, `+` and `−`.
A continuous typing burst coalesces into ONE step (`coalesce` on
`setDraftHeaders`, closed by blur, by `+`/`−`, by undo/redo and by switching
drafts). `+` and `−` are one step each.

### Focal set library

Lists WORKING DRAFTS, not database rows, so a set the user started is visible
before it is saved. Rows show `[n=x]` from the working membership. An inactive
dirty draft carries the pink asterisk; the active one does not, because the Save
control beside its title already says so.

* selecting a row never saves the row being left
* the pencil selects the draft and focuses its title; the rename is local until
  Save, like any other working-copy edit
* locking requires a SAVED, CLEAN draft — locking a dirty one would either lock
  the stale saved copy or auto-save a version the user never approved, so it is
  refused with an explanation instead
* a locked row keeps its lock and loses rename/delete, and stays selectable and
  runnable
* deleting is confirmed; a saved set goes through `project.deleteFocalSet`, a
  never-saved draft is dropped locally, and a dirty saved set says explicitly
  that the unsaved edits go too
* after a deletion a neighbour is selected; if nothing remains, one blank LOCAL
  draft appears and no database row is created for it

**Affordance gap.** The design shows no dedicated "create another focal set"
control. Rather than invent a panel, the capability is a modest `+ New focal
set` action in the library footer. If the intended affordance is something else,
this is the one place to change.

### Sequence visualizer — shell only

No transport, no renderer; the interior is the existing placeholder. What is
implemented is the layout contract:

* the pane OVERLAYS the analysis column (absolute, not a flex sibling), so
  dragging it never reflows or squeezes the controls underneath — verified: the
  analysis column stays 1049px wide at every snap
* the left edge is a drag handle moving between four snaps (34%, 54.6% default,
  70%, 85%) rather than to an arbitrary pixel, so it cannot be parked half
  across something; the handle is fully on screen at both extremes
* double-clicking the handle resets to the design's default; arrow keys and Home
  do the same from the keyboard
* showing the library takes VERTICAL space only — the frame's top edge drops
  from y187 to y349 and its width stays 850px

Tooltips portal to `<body>` at `z-index: 100` against the viewer's `2`, so they
render over it; and the analysis column takes `pointer-events: none` while a
drag is in progress, so sweeping the pointer across it cannot pop tooltips.

### Still drawn but not wired

The Search row under the focal box is present and disabled — a field that looks
live and does nothing is worse than one that says it is not ready. "Estimated
run time" shows an em dash for the same reason.


## 19. Interaction corrections (pass 10)

No redesign. Seven specific behaviours were wrong or inconsistent against
`docs/new/all_pages-20Aug2026.svg` and against each other; this pass fixed them
and left everything else alone.

### The project page behaves the same before and after Create

There was already one `ProjectScreen` and one `SourceTable`. What was still
different was the ROW: a candidate could not be locked, because a candidate has
no `fastaFileId` and the lock is a backend call.

**Pending source lock.** `pending.lockedPaths` holds the paths the user locked
before the project existed - by path, because that is the only identity a
candidate has. The row reads its lock from there and calls
`setPendingCandidateLocked`; a linked row calls `project.setFastaFileLocked`
exactly as before. One `SourceRow`, one lock control, one set of states:

* locking a pending row keeps the lock visible,
* the pencil and the delete disappear while it is locked,
* unlocking restores them,
* hover behaves identically to a linked row.

**Carry-over.** `linkCandidates(candidates, lockedPaths)` links each file and
then, for a file that was locked, immediately calls `project.setFastaFileLocked`
with the id the link produced. So the local intent becomes the persisted lock as
soon as there is a row to lock, and the row transitions pending to linked without
changing what it offers. Tested in `projectScreen.test.tsx`, including the
ordering (`linkFasta` before `setFastaFileLocked`) and the case where nothing was
locked.

A candidate lock is deliberately NOT persisted anywhere else. If the user never
creates the project, it was a decision about a row that never existed.

### Source rows: the pencil, and the help dots

**The pencil belongs to the name.** It used to occupy its own grid column, which
put it hundreds of pixels right of a short file name. Name and pencil are now one
inline group (`.source-row__name-group`) sharing the `1fr` track: the name is a
flex item sized to its content, so the pencil starts where the rendered text
ends, and a long name ellipsises with the pencil still against it. The group's
track ends before the seq column, so neither can collide with the metrics.

**The help dots on a row are a different treatment.** The standard dot is a
`--c-surface` disc, which is exactly the FASTA row's own background - so on a row
it vanished. Page 3 draws those dots the other way round, sampled from the SVG:
a `#525f72` (`--c-bg`) disc with a `#7f92ae` (`--c-muted`) question mark.
`HelpButton` takes `variant="on-row"` for that, rather than every dot in the app
being changed to suit one context.

**Row control rules**, identical for pending and linked rows:

| State | Appearance |
|---|---|
| unlocked, not hovered | no controls |
| row hovered | controls appear muted grey |
| pencil / lock hovered | white |
| delete hovered | red |
| locked | lock only, grey, visible without hover |
| locked, lock hovered | white |

The rename pencil takes the white hover even though renaming is not implemented:
it is part of the row's hover language. It keeps `cursor: default` and a tooltip
that says renaming is not available yet, so it never claims to work.

### Focal library rows: one rule for every row

The active row used to carry `.focal-row.is-active .focal-row__icon { opacity: 1;
color: muted }`. Three classes beat `.focal-row__icon:hover`'s two, so the active
row's controls were always visible and could never light up under the pointer.
The fix is the deletion of that rule, not another override on top of it: every
row now follows one set - hidden at rest, grey on row hover, white for
pencil/lock, red for delete, and a locked row keeps its lock (grey at rest, white
on hover) whether or not the row is hovered.

The pencil moved into a name group here too, for the same reason and with the
same shape. The select button and the pencil are SIBLINGS inside that group,
never nested: a button inside a button is invalid and unpredictable for assistive
tech. Row selection, the active arrow, the dirty asterisk, `[n=x]` and the
accessible names are unchanged.

### The upper focal controls

* **Title to STRING spacing.** Measured on screens 4 and 6: both draw a 33px
  title row whose top is 53px above the STRING row's top - a 20px gap. The saved
  title state was leaving 6px (a `min-height` plus a 6px margin), which is why it
  read as cramped while the input state did not. Both states now take the gap
  from one declaration, so they cannot drift apart again. The parameter block
  below is untouched.
* **The title help dot is gone.** "FOCAL SET TITLE" is a name for a focal set;
  a tooltip explaining that is noise. `HELP_TEXT.focalSetTitle` was removed with
  it rather than left as dead copy.
* **The STRING help dot** sat on top of the Enter chip. The chip overhangs the
  field wrapper it belongs to, so the flex line measured short and the dot landed
  on it - and widening the gap only squeezed Enter instead. The dot is now placed
  at a measured x (`1011px + 6px`), out of the flex line: Enter keeps its measured
  914..1011 and the dot ends at 1047, one pixel before the viewer's rail at 1048.
  The analysis column's right inset became `padding` instead of a transparent
  border so the dot is not clipped, and its scrollbar is inset by a transparent
  border on the scrollbar parts instead. STRING behaviour is unchanged.

### The library decides where the viewer starts

The viewer's top used to be one of two constants (73 without the library, 302
with it), so a library of two rows pushed it as far down as a library of eight
would and the band between stayed permanently empty.

The selectors, the library and the viewer are now ONE flex column
(`.workspace-right`) overlaying the analysis controls. The viewer is the item
after the pane, with a 16px top margin - the design's own gap between the
selectors (ending y=167) and the frame (starting y=183). So:

* one or two focal sets and the viewer starts high;
* more sets and it moves down as the list grows;
* at `--h-library-max` the list scrolls internally and the viewer stops
  descending;
* nothing reserves space for a library larger than the one on screen.

No ResizeObserver and no measured pixel value: the layout does it, which also
means it cannot disagree with the zoom in §19's scaling section. Horizontal
independence is unchanged - the pane is right-anchored at its own width and the
viewer's left edge is a margin percentage.

### One left rail

The workspace had two left-edge elements: a `viewer__edge` drag handle and the
placeholder's own `alignment__rail`, drawn side by side, so there were two thin
edges where the reference has one. The rail now belongs to the viewer
(`viewer__rail`), is the drag target over its whole height, keeps the darker
thumb, and the placeholder draws only what is inside the frame. The stage begins
at the rail's right edge with no gap.

Everything the handle did, the rail does: drag to resize, release snaps,
double-click resets to the design default, arrow keys and Home from the keyboard,
the full-expansion snap reaching the workspace's left boundary, tooltips
suppressed while dragging, and no reflow of the analysis column.

### Undo grouping: a 700ms idle timer

A burst used to stay open until a focus or action boundary, so a user who typed
for two minutes without leaving the field lost all of it to one Undo.

`app/state/typingBurst.ts` closes the burst after ~700ms without a document
change. Every change restarts the countdown, so a burst never expires mid-flow.
Boundaries still close it explicitly and immediately: blur, `+`, `-`, Save,
switching focal sets, undo/redo, and any programmatic canonical replacement.

**The idle timer does not canonicalise.** It closes the history group and
nothing else. A user who pauses halfway through typing a header must not have
their unfinished text rewritten; blur and Save keep their existing
canonicalisation boundaries. Timers are cleared on unmount and when the active
draft changes, and the payload travels through the timer so the burst that ends
is the one that was typed in.

Fake-timer tests in `state/typingBurst.test.tsx` cover rapid typing as one step,
a pause starting a new step, one Undo removing only the most recent burst, and
the absence of one-step-per-keystroke; `projectFlow.test.tsx` proves the wiring
through the real CodeMirror editor.

### 125% Windows scaling

**The authoritative appearance is still 1920x1080 at 100%.** Nothing about the
reference layout was changed to accommodate a smaller viewport.

At 125% Windows scaling the SAME physical screen reports a 1536x864 CSS
viewport - exactly 0.8 of the design in both axes. The window did not get
narrower; the pixels got bigger. Reflowing for that would produce a second design
nobody drew, and the 1920 frame would then have to be maintained against it.

So `app/uiScale.ts` keeps `--ui-zoom` at `min(1, innerWidth / 1920)` and
`global.css` applies it to `#root` with a compensating `width`/`height`, using
`zoom` rather than `transform: scale()` because zoom participates in layout -
viewport units, fixed positioning, scrolling and hit testing all follow it.

* At 1920 the factor is exactly 1 and nothing is touched.
* At 1536 it is 0.8, so the renderer lays out at 1920x1080 and paints into
  1536x864 - which on a 125% display is the same physical size as the reference
  at 100%. Text included: 27px at 0.8 zoom on a 1.25 device ratio is 27
  device-independent pixels again.
* It never magnifies a wider window, and it stops at a 0.55 floor.
* **The launcher is excluded.** It has its own 540x340 window sized to its own
  540x289 design, so measuring it against the workspace frame would read a small
  window as a scaled-down large one and shrink it to a third of its size.
  `useWorkspaceZoom(state.screen !== 'launcher')`.
* Height is deliberately not part of the factor. A short window is a different
  problem, and the analysis column already scrolls internally; scaling for height
  would shrink a layout with no horizontal problem and leave a wide empty margin
  down the right.

Two consequences worth knowing:

* Pointer deltas are in painted pixels while CSS lengths are in layout pixels.
  The focal box's resize grip divides by `currentZoom()` for that reason. The
  viewer's drag works in ratios, so it needs no correction.
* Tooltips portal to `<body>`, outside the zoomed root, so they render at 100%
  over a scaled app. They are overlays and stay readable; positioning is
  unaffected because it comes from `getBoundingClientRect`, which is already in
  painted coordinates.

### How the visual states were checked

The preview pane used for screenshots does not deliver pointer hover to the
page, so hover states cannot be photographed by moving a cursor. They were
verified by rewriting every `:hover` selector in the loaded stylesheets IN PLACE
to a class and applying that class: a pseudo-class and a class have the same
specificity and the rule keeps its position in the sheet, so the cascade that
decides a hover colour is exactly the app's own. That is how the active-row
specificity bug was confirmed fixed rather than assumed fixed.

### Still deferred, deliberately

The real sequence renderer; screen 5 (progress / pause / cancel); PIS
computation; run-time estimation; FASTA rename; source reordering; the Search
row under the focal box. All are drawn inert rather than faked.


## 20. Run observability (pass 11)

A run could sit on "Running..." with nothing to say for itself. On one machine
that was fine; on another it looked like a hang. This pass makes a run report
what it is doing, without changing what it does.

### What the user sees

The existing Run section's status slot, which used to hold one fixed sentence,
now holds one live line:

    Running — verifying FASTA 1/3 — 00:07
    Running — DMC search, size 3 — 12,000 / 266,916 — 00:19
    Running — 5-site optimisation 428,000/2,118,760 — 00:41
    Running — writing comparison workbook — 01:03

Same typography, same muted colour, same place. Three parts: the state
(`Running` / `Continuing`), the stage with whatever it counts, and an elapsed
clock. The clock is driven by a local one-second interval rather than by
arriving progress, because a stage that reports nothing for a minute is exactly
when the user needs to see that the app is alive.

**The renderer owns every word.** Python sends a stage NAME and numbers;
`STAGE_LABELS` in `DiagnosisRunPanel.tsx` is an exhaustive `Record` over the
`DiagnosisProgressStage` union, so a stage added in `progress.py` and to the
union cannot compile without wording, and the line can never go blank because
the two sides drifted.

Not built, deliberately: the deferred progress screen, pause, cancel, and
runtime estimation. This is a status line, not a progress UI.

### How progress travels

    core/pipeline (RunObserver)
      -> RunDiagnostics (service/diagnostics.py)
      -> ProgressChannel -> one JSON line on stdout, tagged with the request id
      -> PythonBridge routes it by envelope kind, leaving the request PENDING
      -> main sends 'analysis:progress' to the window
      -> preload's project.onDiagnosisProgress
      -> ProjectProvider dispatches 'diagnosisProgress'
      -> the reducer applies it, IF the token matches the running run

Typed the whole way: `DiagnosisProgress` in `backendContract.ts`. Nothing parses
stderr, and stderr never reaches React.

**Framing is unchanged.** A progress line carries `type: "progress"` and no
`ok`; a response carries `ok`. The bridge, and both Python test clients,
distinguish them by shape rather than by position, so any number of progress
lines may precede a response.

**Correlation is a token the renderer mints.** `runToken` goes out with the
request and comes back on every notification. The reducer applies progress only
when the run is still `running` AND the token matches, which drops both a
notification that arrives after the run finished and one from a run the user has
already replaced. The backend also stamps its own short `runId`, which is what
appears in the stderr log.

### When the backend dies

The bridge already resolved pending requests on child exit; what matters here is
that the UI acts on it. A `BACKEND_UNAVAILABLE` failure clears the running state
through the normal `diagnosisFailed` path, shows the message and the exit detail,
and re-enables the button. There is still NO timeout on a run: a legitimate
search can take hours, and inventing a deadline would abort real work.

### Developer diagnostics

Structured stderr lines, one JSON object each behind `[mdx]`, carrying the run
id, elapsed time, stage, durations and counts — plus an environment banner at
service start and a traceback with the failing stage on a crash. Electron writes
them to `<userData>/logs/backend-YYYY-MM-DD.log` as well as the console;
`MOLECULAR_TOOL_LOG_DIR` overrides the directory.

A parent that pipes stderr **must drain it**: the pipe fills at 64 KB and the
child then blocks mid-write, which looks precisely like the hang this pass
exists to diagnose. The bridge does; so do the test clients.

### Where the time goes

Measured, not guessed — see REPO_MAP section 14 for the numbers. Both expensive
stages are `C(n, k)` searches whose per-combination cost scales with the contrast
set, and both now announce their total before they start. Neither algorithm was
touched.

### All files with overlapping FASTAs

Investigated, not redesigned: REPO_MAP section 15 records exactly what the
backend does with a full alignment plus its focal subset (refused, one problem
per shared header, no sequence comparison, first-file-wins pooling underneath).
The only change made was capping the refusal's detail, which used to concatenate
one sentence per duplicate into the message this UI displays.
