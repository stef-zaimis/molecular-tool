# UI_NOTES.md

Notes for the Electron/React desktop UI in `desktop/`.

This file is where the uncertainty lives. If something in the code looks like a
decision, the reasoning is here; if something is unresolved, it is in §10 rather
than hidden behind a plausible-looking implementation.

Scope of this phase: **UI only**. No Python was modified, no backend transport
exists, and no analysis can be executed. `python main.py` still launches the
existing Tkinter application, unchanged.

Reference material: `docs/01-launcher.png`, `docs/02-project-creation.png`,
`docs/03-molecular-diagnosis.png`. Behavioural reference for the existing
Python: `REPO_MAP.md` at the repo root.

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
- The workspace "home" button returns to Project Creation, not to the launcher.
  The launcher is reachable from Project Creation's home button. This is an
  assumption — the mockups show a home icon on both screens but not where it goes.
- `setActiveAnalysis` is guarded **in the reducer**, not just on the button, so a
  disabled analysis cannot be activated by any code path.
- If you untick the analysis you are currently viewing, the active tab falls back
  to the first still-enabled analysis rather than showing an empty workspace.

---

## 3. Figma controls whose behaviour was unclear

Implemented visually; behaviour is either minimal or deferred. None of these
invent semantics.

| Control | Where | What it does now |
|---|---|---|
| Hamburger (title bar) | Screens 2, 3 | Renders and has hover/focus states. **No menu** — there is no menu content in the mockups. See Q7. |
| `Save` next to FOCAL SET TITLE | Screen 3 | Shows a "not available yet" notice. The title is already held in session state; there is nothing to save to. |
| Import icon (left of the STRING box) | Screen 3 | Shows a "not available yet" notice. Assumed to mean "load focal strings from a file"; no such format exists. See Q3. |
| `+` / `−` (STRING row) | Screen 3 | Assumed to add/remove **focal sets**. `+` adds an empty set and makes it active; `−` removes the active one and is disabled when only one remains. There is no visible focal-set switcher in the mockup, so added sets are currently reachable only via these buttons. See Q5. |
| Undo / redo arrows | Screen 3 | Fully wired to the editor's own history stack, shared with Ctrl+Z / Ctrl+Shift+Z / Ctrl+Y. Redo is disabled until something is undone — matching the mockup, where redo is drawn dimmed. |
| Resize grip (STRING box) | Screen 3 | Real: the box is vertically resizable. |
| Scrollbar rail left of the alignment panel | Screen 3 | **Presentational only.** The real scrollbar belongs to the future viewer. |

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

4. **Icons are hand-drawn, not exported assets.** `docs/` had no SVGs. Each glyph
   in `src/components/icons/Icons.tsx` was traced from the reference: the
   hamburger (35×26, three 6px bars), window controls (26px ink), and the home
   glyph (40×36, overhanging roof with a door notch) match measured geometry.
   Known differences: the reference roof appears to be drawn as an outlined
   triangle with an inner filled triangle and small 2px gaps near the eaves; ours
   is a single filled triangle. At 100% this is not visible. The import, undo,
   redo and plus/minus glyphs are close but not pixel-traced.

5. **The tooltip is 317px wide with a left-pointing notch**, matching the one
   open tooltip in screen 2. Only one tooltip style exists, so this is applied to
   every help dot.

6. **The alignment placeholder shows the text "Alignment Viewer"** plus the
   selected file name. The reference panel is empty. Per the brief the area
   should "look intentional, not like unfinished white space", so it carries a
   recessive (42% opacity) label. Remove it when the viewer lands.

7. **Numbers in the spinners are visible in our build.** The mockup's spinners
   appear empty. Ours show the actual values (1 and 2, the Python defaults),
   because an empty numeric control with no value would be misleading.

8. **A notice bar appears for unfinished actions.** Not in the design. It only
   renders after the user activates such a control, so the resting UI stays clean.

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
    offscreen Electron.

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
