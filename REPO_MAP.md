# REPO_MAP.md

Map of `molecular-tool`, branch `feature/ui-revamp`.
Last updated 2026-08-22, after the run-observability pass (section 14).

The repository is three layers, not one:

1. **The scientific core** — the original Python package (`molecular_diagnosis/`),
   unchanged in intent since the first map and still the only place DMC science
   happens.
2. **The project layer** — `molecular_diagnosis/project/` and
   `molecular_diagnosis/service/`: a SQLite-backed project format, linked FASTA
   sources with a live status model, persistent focal sets, and a
   newline-delimited-JSON service that exposes all of it.
3. **The desktop UI** — `desktop/`: an Electron + React renderer that never
   imports Python and never touches the database; it calls named service
   methods through a preload bridge.

Conventions used below:
- Line references are `path:line` against the working tree.
- `(UNVERIFIED)` marks a claim not confirmed from the code alone; each one says
  what would settle it.
- Line counts are `wc -l` and may be off by one where a file lacks a trailing
  newline.
- Sections 4, 6, 7, 10 and 12 analyse the scientific core. They were written
  against the pre-revamp tree and still describe it: the modules they cover
  (`core`, `punishments`, `consensus`, `sequence_subsets`, `excel`, `reports`,
  `pipeline`) have not been restructured since, and `tests/test_parity.py`
  exists to prove exactly that.

---

## 1. INVENTORY

### The scientific core (Python)

| File | Lines | Purpose | Class |
|---|---|---|---|
| `main.py` | 3 | Entry point for the LEGACY Tkinter app; imports `launch_gui` and calls it. | **LIVE (legacy UI)** |
| `molecular_diagnosis/__init__.py` | 6 | Package docstring and `__version__`. | **LIVE** |
| `molecular_diagnosis/constants.py` | 81 | Output filenames, base to hex colour table, IUPAC tables, punishment thresholds/weights. | **LIVE** |
| `molecular_diagnosis/models.py` | 81 | Frozen dataclasses: `DMCResult`, `FiveSiteResult`, `PunishmentEvent`, `PunishmentResult`, `PunishmentPipelineResult`, `PipelineResult`. | **LIVE** |
| `molecular_diagnosis/utils.py` | 26 | `next_available_filename` - non-clobbering output naming (`name(2).txt`). | **LIVE** |
| `molecular_diagnosis/fasta_io.py` | 68 | FASTA parsing, aligned-length validation, focal/non-focal header split. | **LIVE** |
| `molecular_diagnosis/progress.py` | 233 | **Observation hooks.** `RunObserver` (a no-op base class), the stage-name vocabulary that crosses the process boundary, `ProgressTicker` (throttled reporting for million-iteration loops), and `describe_path`. Imports nothing, so the science can report without knowing what a service is. | **LIVE** |
| `molecular_diagnosis/focal.py` | 177 | **The single focal matcher.** Case-sensitive substring containment, OR across selectors, selectors trimmed and deduplicated, empty selector refused. The four places that once wrote `target_string in header` by hand all route through it. | **LIVE** |
| `molecular_diagnosis/core.py` | 519 | DMC search: site scoring, focal consensus per column, candidate filtering, n-site combination search, 5-site optimisation, formatting. | **LIVE** |
| `molecular_diagnosis/punishments.py` | 508 | Focal-only punishment/anomaly scoring (POLY / BAL / PRL / INS / EW / BD). | **LIVE** |
| `molecular_diagnosis/consensus.py` | 347 | Focal consensus sequence (untrimmed and PRL/INS-trimmed) and its text report. | **LIVE** |
| `molecular_diagnosis/sequence_subsets.py` | 217 | Ungapped-substring grouping of focal sequences plus its Excel writer. | **LIVE** |
| `molecular_diagnosis/excel.py` | 239 | `comparison_output.xlsx` (Full/Gap5/Avg5 sheets) and `punishment_output.xlsx`. | **LIVE** |
| `molecular_diagnosis/reports.py` | 138 | `DMCs_output.txt` text report writer. | **LIVE** |
| `molecular_diagnosis/pipeline.py` | 300 | Orchestration: `load_inputs`, `run_pipeline_core`, `run_punishment_core`. The only UI-free composition layer. | **LIVE** |
| `molecular_diagnosis/gui.py` | 338 | **Legacy** Tkinter main window. Still runnable; superseded by `desktop/`. | **LIVE (legacy UI)** |
| `molecular_diagnosis/viewer.py` | 750 | **Legacy** bitmap-backed FASTA/alignment viewer, imported by `gui.py:7`. | **LIVE (legacy UI)** |

### The project layer (Python)

Everything durable lives here. The renderer holds no SQL and the Electron main
process holds no project state; this package owns the database exclusively.

| File | Lines | Purpose |
|---|---|---|
| `molecular_diagnosis/project/__init__.py` | 25 | Package exports. |
| `molecular_diagnosis/project/paths.py` | 36 | Where a project's `project.sqlite` and `outputs/` live. |
| `molecular_diagnosis/project/db.py` | 208 | Connection handling and the migration ladder (`001_initial.sql`, `002_fasta_file_locked.sql`). |
| `molecular_diagnosis/project/repository.py` | 535 | Every SQL statement in the application. Row dataclasses for projects, FASTA files, the header index, focal sets and entries. |
| `molecular_diagnosis/project/sources.py` | 304 | The live SOURCE STATUS model: `current`, `unverified`, `stale`, `never_indexed`, `missing`, `unreadable`, plus the persisted `locked` flag. Nothing about availability is stored. |
| `molecular_diagnosis/project/indexing.py` | 229 | Header indexing and fingerprinting (`indexed_sha256`, size, mtime). |
| `molecular_diagnosis/project/alignments.py` | 95 | Reading a linked alignment for a run, and vetting a candidate before it is linked. |
| `molecular_diagnosis/project/locations.py` | 193 | Cached header-to-file locations behind presence answers. |
| `molecular_diagnosis/project/search.py` | 152 | Header search over the index: the capped preview and the uncapped resolve. |
| `molecular_diagnosis/project/service.py` | 1421 | `ProjectService` - the operations the UI calls. Owns the rules: a locked focal set refuses every mutation, a locked FASTA refuses unlinking, `+` refuses a partial expansion, a run refuses an out-of-scope focal entry. |

### The service layer (Python)

| File | Lines | Purpose |
|---|---|---|
| `molecular_diagnosis/service/__init__.py` | 15 | Package docstring and exports. |
| `molecular_diagnosis/service/__main__.py` | 71 | `python -m molecular_diagnosis.service`: the stdio loop Electron spawns. Redirects `sys.stdout` to stderr so a stray `print` cannot corrupt the protocol. |
| `molecular_diagnosis/service/protocol.py` | 199 | Newline-delimited JSON framing: requests, responses, and progress notifications. |
| `molecular_diagnosis/service/diagnostics.py` | 330 | Run diagnostics (structured stderr lines, environment snapshot, run ids) and `ProgressChannel`, which turns an observer's progress into a protocol envelope for the request in flight. |
| `molecular_diagnosis/service/errors.py` | 200 | `ServiceError` and the error codes the UI writes sentences for. |
| `molecular_diagnosis/service/handlers.py` | 435 | `METHODS` and `dispatch`. Non-project methods: `ping`, `loadFasta`, `validateFocalStrings`, `runMolecularDiagnosis`, `runSequencePunishment`. |
| `molecular_diagnosis/service/projects.py` | 508 | `PROJECT_METHODS` - the 28 `project.*` methods (listed in section 3d). |

### Database schema

| File | Lines | Purpose |
|---|---|---|
| `db/001_initial.sql` | 190 | Projects, FASTA files, header index, focal sets and entries. |
| `db/002_fasta_file_locked.sql` | 18 | Adds `fasta_file.locked`, additive and defaulted so an existing project opens unchanged. Sets `PRAGMA user_version = 2`. |
| `db/optional_fts5_trigram.sql` | 25 | Optional accelerated search, applied only where SQLite supports it. |

### The desktop UI (`desktop/`)

Electron main plus preload plus React renderer. `src/main.ts` spawns the Python
service, `src/preload.ts` exposes `window.desktop.*`, and the renderer calls
that and nothing else.

| Area | Files (lines) | Notes |
|---|---|---|
| Shell | `src/main.ts` (303), `src/preload.ts` (339), `src/renderer.tsx` (25), `src/backend/pythonBridge.ts` (262) | Process spawn, IPC forwarding, interpreter resolution. |
| Contracts | `src/backendContract.ts` (458), `src/contract.ts` (674) | `backendContract.ts` is what actually crosses the boundary; `contract.ts` describes the wider intended architecture and carries the `// BACKEND GAP:` notes. |
| State | `src/app/state/projectState.ts` (968), `ProjectContext.tsx` (816), `focalDrafts.ts` (271), `sourceStatus.ts` (198), `runGate.ts` (135), `typingBurst.ts` (79) | One reducer plus context. `focalDrafts.ts` holds the saved-versus-working distinction; `typingBurst.ts` is the undo-grouping idle timer. |
| Scaling | `src/app/uiScale.ts` (80) | `--ui-zoom`: the 1920x1080 design frame, scaled for a high-DPI viewport (section 5b). |
| Screens | `src/screens/LauncherScreen.tsx` (104), `ProjectScreen.tsx` (258), `MolecularDiagnosisScreen.tsx` (480) | `ProjectScreen` is ONE screen for both project states, before and after Create. |
| Components | `src/components/sources/SourceTable.tsx` (347), `focal/FocalSetEditor.tsx` (396), `focal/FocalSetLibrary.tsx` (225), `alignment/WorkspaceRightPane.tsx` (211), `alignment/AlignmentPlaceholder.tsx` (37), `analysis/DiagnosisRunPanel.tsx` (159), `controls/*`, `chrome/*`, `icons/Icons.tsx` (340) | `SourceTable` serves candidates and linked sources from ONE `SourceRow`. |
| Styles | `src/styles/tokens.css` (241), `src/styles/global.css` (220), plus one CSS file per component | Every value in `tokens.css` was measured from the design SVG. |

### Design references

| Path | Status |
|---|---|
| `docs/new/all_pages-20Aug2026.svg` | **Authoritative.** Six pages; page 3 is the project/FASTA list, page 4 the new focal set, page 6 the saved focal set with the library. |
| `docs/new/loading.svg` | Authoritative; not implemented yet. |
| `docs/all_pages-20Aug2026.svg` | The same drawing at the old path, kept while references migrate. |
| `docs/legacy/` | **History only.** `01-launcher.png`, `02-project-creation.png`, `02-project-creation.svg`, `03-molecular-diagnosis.png`, `all_pages.svg`. Must not drive new work. |

### Not live

| File | Lines | Purpose | Class | Evidence |
|---|---|---|---|---|
| `molecular_diagnosis/viewer_old.py` | 473 | Earlier canvas-item-based virtualised viewer. | **LEGACY** | Nothing imports it; a grep for `viewer_old` across the project returns only its own file. |
| `pipeline_verbose.py` | 746 | Self-contained pre-refactor monolith with its own parser, Excel writer, Tk GUI and viewer. | **LEGACY** | Nothing imports it. Untouched since the first commit. |
| `consensus_dmc_pipeline.py` | 358 | Standalone consensus, DMC-site-mapping and Word-document script. | **LEGACY, UNRUNNABLE** | Imports `docx`, which is neither installed nor in `requirements.txt`, and hardcodes an absolute Linux path. Its consensus rule was ported into `consensus.py`, which says so at `consensus.py:44`. |

### Config and data

Unchanged: `pyproject.toml` holds pytest settings only, `requirements.txt`
lists `openpyxl` and `pillow` unpinned, and `input/` and `output/` are
gitignored. `desktop/package.json` holds the Electron Forge, Vite, Vitest and
ESLint scripts; `@codemirror/state` and `@codemirror/view` are the only
non-React runtime dependencies.

---

## 2. ENTRY POINTS AND CONTROL FLOW

### The current entry point

**The desktop app.** `npm start` in `desktop/` runs Electron Forge, which loads
`src/main.ts`. Main spawns `python -m molecular_diagnosis.service` as a child
process and talks to it in newline-delimited JSON over stdin/stdout; the
renderer never spawns anything and never sees the scientific code. Walkthrough D
below follows a run from a click to a written file.

Interpreter resolution, in order: `MOLECULAR_TOOL_PYTHON`, then the project
`.venv`, then `python`/`python3` on PATH. `MOLECULAR_TOOL_ROOT` overrides the
repo root. Packaging a bundled interpreter is not addressed.

The Tkinter entry point below still works and is described as it stands, but it
is the LEGACY UI: it has no project database, no persistent focal sets and no
multi-file support.

### Is `main.py` the real entry point?

Yes, for the application. `main.py` is three lines: it imports `launch_gui` from `molecular_diagnosis.gui` and calls it. There is **no CLI** — no `argparse`, no `console_scripts` entry (`pyproject.toml` holds pytest config only). The GUI is the only supported way to run the tool.

Two **other** executable entry points exist, both legacy (§1):
- `pipeline_verbose.py:746` → its own separate `launch_gui()`.
- `consensus_dmc_pipeline.py:357` → `main()` (batch script, currently unrunnable — no `docx`).

`tests/test_pipeline.py:3` is a fourth way in: it calls `run_pipeline_core` directly, bypassing the GUI. This is the existing proof that the core is callable headless today.

### Walkthrough A — "Run DMC Analysis"

1. `main.py:4` → `gui.launch_gui()` (`gui.py:10`) builds the Tk root and widgets, then blocks in `root.mainloop()` (`gui.py:339`).
2. User fills three text vars, two checkboxes and two spinboxes (§7), clicks **Run DMC Analysis** (`gui.py:327-331`) → callback `run()` (`gui.py:118`).
3. `run()` re-reads the Tk vars, does its *own* empty-string checks (`gui.py:123-133`) and its *own* min/max validation (`gui.py:139-152`) — duplicating checks that `pipeline.load_inputs` and `core.find_dmc_information` also perform.
4. `run()` enters a `while True:` loop (`gui.py:159`) and calls `pipeline.run_pipeline_core(...)` (`gui.py:160`).
5. `run_pipeline_core` (`pipeline.py:91`) → `load_inputs` (`pipeline.py:31`):
   - `fasta_io.parse_fasta` (`fasta_io.py:4`) → `dict[header, sequence]`, uppercased.
   - `fasta_io.validate_aligned_fasta` (`fasta_io.py:32`) → single alignment length, or `ValueError`.
   - `fasta_io.split_focal_headers` (`fasta_io.py:44`) → focal/non-focal header lists by substring match.
6. `pipeline.py:118` sets `ref_id = focal_headers[0]` — the first focal header **in FASTA file order**.
7. `consensus.build_focal_consensus_result` (`consensus.py:197`) builds untrimmed + trimmed focal consensus (§4).
8. `core.find_dmc_information` (`core.py:188`) runs the whole DMC search (§4) and returns a `DMCResult`.
9. `core.find_best_five_site_sets` (`core.py:441`) brute-forces every 5-subset of `dmc.unique` (§4, §10).
10. Three output paths are reserved via `utils.next_available_filename` (`pipeline.py:149-153`).
11. `reports.write_text_report` (`reports.py:14`) writes `DMCs_output.txt`. Note `punishment_result=None` is passed unconditionally at `pipeline.py:167`.
12. `consensus.write_consensus_text_report` (`consensus.py:233`) writes `focal_consensus_output.txt`.
13. `excel.write_excel_report` (`excel.py:97`) writes `comparison_output.xlsx` with sheet `Full`, plus `Gap5`/`Avg5` when the 5-site search produced results.
14. Back in the GUI: if `result.dmc.stop_reason != "reached_maximum_length"` the loop shows a summary and breaks (`gui.py:175-177`). Otherwise it asks whether to continue (`gui.py:179`), prompts for a new maximum via `simpledialog.askinteger` (`gui.py:192`), carries `diagnostic_combinations` and `combinations_tested_by_length` forward (`gui.py:205-206`), and **loops back to step 5**.
    - Each iteration re-parses the FASTA, re-runs the consensus, re-runs the 5-site search, and writes **a fresh numbered set of all three output files**. A three-round search leaves `DMCs_output.txt`, `DMCs_output(2).txt`, `DMCs_output(3).txt`, and so on.
15. Any exception in steps 5–14 is caught at `gui.py:208` and shown as `messagebox.showerror("Error", str(error))`; the traceback is discarded.

### Walkthrough B — "Run Punishment Analysis"

1. Button at `gui.py:333-337` → `run_punishments()` (`gui.py:211`).
2. `pipeline.run_punishment_core` (`pipeline.py:198`) → `load_inputs` again — a **second full parse** of the same FASTA if the user already ran A.
3. `punishments.find_focal_punishments` (`punishments.py:389`) — focal sequences only; non-focal sequences are never read.
4. `sequence_subsets.group_sequence_subsets` (`sequence_subsets.py:25`), also focal-only (`pipeline.py:222-229`).
5. `sequence_subsets.write_sequence_subset_excel_report` (`sequence_subsets.py:138`) → `sequence_subsets_output.xlsx`.
6. `excel.write_punishment_excel_report` (`excel.py:144`) → `punishment_output.xlsx`.
7. Output paths are shown in a messagebox (`gui.py:235`).

Note the ordering at `pipeline.py:231-249`: both output paths are reserved *before* either file is written, and the **subsets** workbook is written before the **punishment** workbook.

### Walkthrough C — "View FASTA"

`gui.py:56` `view_fasta()` → `parse_fasta` (a third, independent parse — nothing is cached) → `viewer.open_fasta_viewer(root, sequences)` (`viewer.py:10`). The viewer never receives any analysis result; it gets only the raw `dict[header, sequence]`.

### Walkthrough D - a project-backed Molecular Diagnosis run

1. Launcher: **Create New** or **Open Existing**. Both open a directory picker
   in main (`dialog.showOpenDialog`), then call `project.create` or
   `project.open`. They are separate methods on purpose: `create` refuses a
   folder that already holds a project (`PROJECT_ALREADY_EXISTS`), `open`
   refuses one that does not (`PROJECT_NOT_FOUND`).
2. The response carries the project metadata, its capabilities and a full
   source sweep, so a file that vanished while the app was closed is visible on
   the first paint.
3. **Project page** (`ProjectScreen`). Browse vets every chosen file through
   `project.validateFastaCandidate` BEFORE anything is linked; an unaligned or
   unreadable file is reported inline and the rest of the selection survives.
   With a project open, accepted files go straight to `project.linkFasta`;
   before Create they are held as candidates and linked by Create in order.
4. Analyses ticked here become the workspace tabs.
5. **Molecular Diagnosis workspace**. The FASTA pool selects the run scope;
   `+` resolves through `project.resolveFocalAddQuery` (uncapped, scope-checked)
   and `-` through `project.matchFocalHeaders`; entry colours come from
   `project.headerPresence`, debounced and generation-tagged in the renderer.
6. Editing writes nothing. `project.saveFocalSet` is the only write ordinary
   editing performs, and `runGate.canRun` refuses a new or dirty draft - checked
   again inside `runDiagnosis`, not just on the button.
7. **Run** calls `project.runMolecularDiagnosis` with the focal set id and the
   scope file ids. The service re-reads and verifies the linked files, then
   calls the same `pipeline.run_pipeline_core` the Tkinter app calls, writing
   into the project's `outputs/` directory.

### Call graph (live code)

```
main.py
└── gui.launch_gui
    ├── fasta_io.parse_fasta                      (View FASTA path)
    ├── viewer.open_fasta_viewer
    │   └── constants.COLORS
    ├── pipeline.run_pipeline_core
    │   ├── pipeline.load_inputs
    │   │   └── fasta_io.{parse_fasta, validate_aligned_fasta, split_focal_headers}
    │   ├── consensus.build_focal_consensus_result
    │   │   ├── consensus.build_untrimmed_consensus → assign_consensus_state
    │   │   ├── consensus.classify_prl_ins_sites
    │   │   └── consensus.build_trimmed_consensus  → assign_consensus_state
    │   ├── core.find_dmc_information
    │   ├── core.find_best_five_site_sets → core.compute_metrics → compute_similarity → score_state
    │   ├── utils.next_available_filename  (x3)
    │   ├── reports.write_text_report → core.format_diag_from_states
    │   ├── consensus.write_consensus_text_report
    │   └── excel.write_excel_report → excel.build_sheet → core.{extract_sites, compute_match_score, compute_similarity}
    └── pipeline.run_punishment_core
        ├── pipeline.load_inputs
        ├── punishments.find_focal_punishments
        ├── sequence_subsets.group_sequence_subsets
        ├── utils.next_available_filename  (x2)
        ├── sequence_subsets.write_sequence_subset_excel_report
        └── excel.write_punishment_excel_report
```

---

## 3. MODULE MAP

Sections 3a-3c below map the scientific core. 3d and 3e map the project/service
layer and the desktop renderer.

### 3a-3c. The scientific core

**No circular imports.** The live dependency graph is a DAG: `constants` and `models` are leaves; `pipeline` is the only node touching every branch; `gui` sits above `pipeline` and `viewer`. Verified by collecting every `import`/`from` line in the project (excluding `.venv`).

One structural note, stated without recommendation: `excel.py` (output layer) imports `core.py` (science layer) at `excel.py:8-12`, so Excel generation currently depends on the scoring functions.

### `constants.py`
Owns: output basenames, `COLORS` (16 entries incl. gap and `?`), `IUPAC` (16 entries; `-` and `?` map to `set()`), `IUPAC_FROM_BASES` (derived comprehension, `constants.py:62-66`), `STRICT_BASES`, punishment thresholds/weights.
Public names: `TXT_OUTPUT_BASENAME`, `XLSX_OUTPUT_BASENAME`, `PUNISHMENT_XLSX_OUTPUT_BASENAME`, `CONSENSUS_TXT_OUTPUT_BASENAME`, `SEQUENCE_SUBSETS_XLSX_OUTPUT_BASENAME`, `COLORS`, `IUPAC`, `IUPAC_FROM_BASES`, `STRICT_BASES`, `POLYMORPHISM_EMPTY_LIMIT`, `BALANCING_EMPTY_LIMIT`, `POLYMORPHISM_EMPTY_WEIGHT`, `BALANCING_EMPTY_WEIGHT`, `PROLONGATION_WEIGHT`, `BD_AMBIGUOUS_FRACTION_LIMIT`.
Imports: nothing. Imported by: `core`, `punishments`, `consensus`, `excel`, `pipeline`, `viewer`, `viewer_old`.

### `models.py`
Owns the six frozen dataclasses listed in §1. All are `@dataclass(frozen=True)`, but every one holds **mutable** containers (dicts/lists), so freezing prevents attribute rebinding only, not mutation of the contents.
Imports: `dataclasses`, `pathlib`. Imported by: `core`, `punishments`, `reports`, `excel`, `pipeline`.

### `utils.py`
- `next_available_filename(path: str | Path) -> Path` (`utils.py:4`)
Imports: `pathlib`. Imported by: `pipeline`.

### `fasta_io.py`
- `parse_fasta(path: str | Path) -> dict[str, str]` (`fasta_io.py:4`)
- `validate_aligned_fasta(sequences: dict[str, str]) -> int` (`fasta_io.py:32`)
- `split_focal_headers(sequences: dict[str, str], target_string: str) -> tuple[list[str], list[str]]` (`fasta_io.py:44`)
Imports: `pathlib`. Imported by: `pipeline`, `gui`, `tests/test_fasta_io.py`.

### `core.py`
Owns the DMC science. Public functions:
- `state_possibilities(state: str) -> set[str]` (`core.py:11`)
- `is_empty_state(state: str) -> bool` (`core.py:15`)
- `is_ambiguous_state(state: str) -> bool` (`core.py:19`)
- `column_empty_fraction(column: list[str]) -> float` (`core.py:24`)
- `column_has_empty_state(column: list[str]) -> bool` (`core.py:31`)
- `column_has_ambiguous_state(column: list[str]) -> bool` (`core.py:35`)
- `state_matches_dmc_base(state, base, *, include_ambiguous_dmc_bd) -> bool` (`core.py:39`)
- `consensus_base_for_column(column, *, include_ambiguous_dmc_bd, allow_gaps) -> str | None` (`core.py:56`)
- `column_has_non_empty_variation_from_base(column, base, *, include_ambiguous_dmc_bd) -> bool` (`core.py:108`)
- `extract_sites(seq: str, sites) -> list[str]` (`core.py:128`)
- `score_state(ref: str, query: str) -> float` (`core.py:132`)
- `compute_similarity(ref_states, query_states) -> float` (`core.py:144`)
- `compute_match_score(ref_states, query_states) -> float` (`core.py:152`)
- `compute_metrics(sequences, ref_id, sites, target_string, diagnostic_states=None) -> tuple[float, float, float]` (`core.py:156`)
- `find_dmc_information(sequences, target_string, *, include_ambiguous_dmc_bd=False, include_gappy_consensus_dmc_sites=False, min_combination_length=1, max_combination_length=2, start_combination_length=1, initial_diagnostic_combinations=None, initial_combinations_tested_by_length=None) -> DMCResult` (`core.py:188`)
- `find_best_five_site_sets(sequences, ref_id, sites, target_string, diagnostic_states=None) -> FiveSiteResult` (`core.py:441`)
- `format_diag(sites, seq) -> str` (`core.py:491`)
- `format_diag_from_states(sites, states) -> str` (`core.py:498`)
Imports: `itertools.combinations`, `constants`, `models`. Imported by: `pipeline`, `excel`, `reports`, `tests/test_core.py`.

### `punishments.py`
Owns punishment scoring. Public functions:
- `state_possibilities` / `is_empty_state` (`punishments.py:17`, `:21`) — **independent duplicates** of the `core.py` versions.
- `empty_fraction(column) -> float` (`punishments.py:25`) — duplicate of `core.column_empty_fraction`.
- `compatible_strict_base_for_bd(column) -> str | None` (`punishments.py:33`)
- `weighted_counts_for_column(column, *, focal_index=None, assumed_base=None) -> dict[str, float]` (`punishments.py:86`)
- `score_sequence_as_base(column, sequence_index, assumed_base) -> float` (`punishments.py:126`)
- `score_non_empty_state(column, sequence_index, state) -> float` (`punishments.py:156`)
- `add_score_event(...) -> None` (`punishments.py:185`)
- `process_polymorphism_or_balancing_column(...) -> None` (`punishments.py:219`)
- `process_gap_dominated_column(...) -> None` (`punishments.py:331`)
- `find_focal_punishments(sequences: Mapping[str, str], focal_headers: list[str]) -> PunishmentResult` (`punishments.py:389`)
Imports: `collections.defaultdict`, `collections.abc.Mapping`, `constants`, `models`. Imported by: `pipeline`.

### `consensus.py`
Owns consensus building **and** its own report writer **and** its own result dataclass.
- `ConsensusResult` dataclass (`consensus.py:13`) — lives here, not in `models.py`.
- `wrap_sequence(sequence, width=80) -> str` (`consensus.py:24`)
- `is_empty_state` / `empty_fraction` (`consensus.py:31`, `:35`) — a **third copy** of these predicates.
- `assign_consensus_state(non_gap_column) -> str` (`consensus.py:42`)
- `build_untrimmed_consensus(focal_sequences) -> tuple[str, int]` (`consensus.py:82`)
- `classify_prl_ins_sites(focal_sequences) -> tuple[set[int], set[int]]` (`consensus.py:110`)
- `build_trimmed_consensus(focal_sequences, removed_prl, removed_ins) -> tuple[str, tuple[int, ...], int]` (`consensus.py:160`)
- `build_focal_consensus_result(focal_sequences) -> ConsensusResult` (`consensus.py:197`)
- `write_consensus_text_report(output_path, *, target_string, focal_headers, alignment_length, consensus_result, dmc_sites=None) -> None` (`consensus.py:233`)
Imports: `dataclasses`, `pathlib`, `collections.abc.Sequence`, `constants`. Imported by: `pipeline`.

### `sequence_subsets.py`
Owns grouping **and** its Excel writer **and** its dataclass.
- `SequenceSubsetGroup` dataclass with a `count` property (`sequence_subsets.py:10`)
- `ungap(sequence) -> str` (`sequence_subsets.py:21`)
- `group_sequence_subsets(sequences: Mapping[str, str]) -> list[SequenceSubsetGroup]` (`sequence_subsets.py:25`)
- `autosize_columns(ws, *, max_width=100) -> None` (`sequence_subsets.py:112`) — **near-duplicate** of `excel.autosize_columns` (`excel.py:16`), differing only by the `max_width` cap.
- `excel_safe_sequence(sequence) -> str` (`sequence_subsets.py:124`)
- `write_sequence_subset_excel_report(output_path, sequence_subset_groups) -> None` (`sequence_subsets.py:138`)
Imports: `dataclasses`, `pathlib`, `collections.abc.Mapping`, `openpyxl`. Imported by: `pipeline`.

### `excel.py`
- `autosize_columns(ws) -> None` (`excel.py:16`)
- `build_sheet(workbook, name, sequences, ref_id, sites, target_string, diagnostic_states=None) -> None` (`excel.py:28`)
- `write_excel_report(output_path, sequences, ref_id, full_sites, target_string, best_gap_sites, best_avg_sites, diagnostic_states=None) -> None` (`excel.py:97`)
- `write_punishment_excel_report(output_path, focal_headers, alignment_length, punishment_result) -> None` (`excel.py:144`)
Imports: `pathlib`, `openpyxl`, `constants.COLORS`, `core.{compute_match_score, compute_similarity, extract_sites}`, `models.PunishmentResult`. Imported by: `pipeline`.

### `reports.py`
- `format_combo(combo, states) -> str` (`reports.py:7`)
- `write_text_report(output_path, fasta_path, output_dir, target_string, sequences, alignment_length, focal_headers, non_focal_headers, ref_id, dmc, five_site_result, punishment_result=None) -> None` (`reports.py:14`)
Imports: `pathlib`, `core.format_diag_from_states`, `models`. Imported by: `pipeline`.
The `punishment_result` parameter is **dead code**: `pipeline.py:167` always passes `None`, and the only branch consuming it (`reports.py:126-131`) writes a placeholder saying detailed punishment reporting "is not currently implemented".

### `pipeline.py`
- `load_inputs(fasta_path, target_string, output_dir) -> tuple[...]` (`pipeline.py:31`) — returns a positional 7-tuple.
- `run_pipeline_core(...) -> PipelineResult` (`pipeline.py:91`)
- `run_punishment_core(fasta_path, target_string, output_dir) -> PunishmentPipelineResult` (`pipeline.py:198`)
Imports: every package module except `gui`, `viewer`, `viewer_old`. Imported by: `gui`, `tests/test_pipeline.py`.

### `gui.py`
- `launch_gui() -> None` (`gui.py:10`) — the *only* module-level definition. Every widget, callback, and piece of parameter state is a local or closure inside this single function.
Imports: `pathlib`, `tkinter`, `fasta_io`, `pipeline`, `viewer`. Imported by: `main.py`.

### `viewer.py`
- `open_fasta_viewer(root: tk.Tk, sequences: Mapping[str, str], *, color_bases=True, show_letters=True, show_grid=False) -> None` (`viewer.py:10`) — again a single function, holding roughly two dozen nested closures.
Imports: `tkinter`, `tkinter.font`, `collections.abc.Mapping`, `PIL.{Image, ImageDraw, ImageTk}`, `constants.COLORS`. Imported by: `gui`.

### `viewer_old.py` (LEGACY)
- `open_fasta_viewer(root, sequences, *, color_bases=True, show_grid=False) -> None` (`viewer_old.py:8`) — same name, missing the `show_letters` keyword.
Imports: `tkinter`, `tkinter.font`, `constants.COLORS`. **Imported by: nothing.**

---

### 3d. The project and service layers

`service/handlers.py` merges `PROJECT_METHODS` into one `METHODS` table, so the
frontend has a single request channel. The 28 project methods:

| Group | Methods |
|---|---|
| Lifecycle | `project.create`, `project.open`, `project.close`, `project.setTitle` |
| Sources | `project.refreshSources`, `project.validateFastaCandidate`, `project.linkFasta`, `project.unlinkFasta`, `project.setFastaFileLocked`, `project.relinkFasta`, `project.reindexFasta` |
| Search | `project.searchHeaders` (capped preview), `project.resolveFocalAddQuery` (uncapped, what `+` uses), `project.matchFocalHeaders` (what `-` uses) |
| Focal sets | `project.listFocalSets`, `project.getFocalSet`, `project.createFocalSet`, `project.renameFocalSet`, `project.setFocalSetLocked`, `project.deleteFocalSet`, `project.replaceFocalEntries`, `project.saveFocalSet`, `project.addFocalEntries`, `project.removeFocalEntries` |
| Presence | `project.headerPresence`, `project.focalPresence` |
| Run | `project.runMolecularDiagnosis` |

Rules that live in Python rather than in disabled buttons:

- A locked focal set may be selected, presence-checked and analysed; every
  mutation is refused with `FOCAL_SET_LOCKED`.
- A locked FASTA source stays fully analysable; unlinking and renaming it are
  refused. Locking protects the LINK, not the file.
- `project.searchHeaders` is a capped preview and `+` must not truncate a focal
  set, so `+` resolves through `project.resolveFocalAddQuery` instead.
- A `+` whose scope contains a file that could not be searched is refused
  outright (`SEARCH_SCOPE_UNAVAILABLE`) rather than persisting half an
  expansion.
- `project.replaceFocalEntries` applies the big textbox as a DIFF, so surviving
  entries keep their ids and cached locations, and a header that matches no
  FASTA is KEPT so it can show red.

### 3e. The desktop renderer

    open/create project -> linked FASTA sources -> selected FASTA scope
      -> working focal-set draft -> saved persistent focal set
      -> project.runMolecularDiagnosis

- `projectState.ts` is the only reducer. Screens dispatch into it; no screen
  holds its own copy of project data.
- `focalDrafts.ts` holds `persistedId`, the saved title/headers, the working
  title/headers, `locked`, the undo history and `burstOpen`. `dirty` is DERIVED
  by comparing normalised working values against saved ones, never stored.
- `sourceStatus.ts` keeps `available`, `indexUsable` and `activity` apart -
  conflating them is what produces a UI that lies.
- `runGate.ts` answers one question, "may this run start", with a reason.
- `typingBurst.ts` closes a manual typing burst after ~700ms idle so one Undo
  cannot swallow an entire editing session. It closes the history group only and
  never rewrites the document.
- `uiScale.ts` keeps `--ui-zoom` in step with the window so the 1920x1080 design
  frame survives 125% Windows scaling.


## 4. THE SCIENTIFIC CORE

### 4.1 FASTA parsing and alignment representation

`fasta_io.parse_fasta` (`fasta_io.py:4-29`).

**The alignment in memory is a plain `dict[str, str]`**: header (the `>` line minus `>`, stripped) → the concatenated sequence, `.upper()`-ed at `fasta_io.py:19` and `:27`. There is no matrix, no per-column index, no NumPy array, no sequence object. Every column access anywhere in the codebase is an ad-hoc list comprehension of the form `[seq[site] for seq in sequences]`, rebuilt from scratch each time it is needed (`core.py:231-232`, `core.py:262-263`, `punishments.py:433`, `consensus.py:93`, `consensus.py:131`, `consensus.py:179`).

Consequences that are load-bearing:
- Blank lines are skipped (`fasta_io.py:14-15`); no other line is validated.
- **Duplicate headers silently collapse** — later records overwrite earlier ones, and the sequence count shrinks without warning. The real input file happens to have 2354 unique headers for 2354 `>` lines, so this is latent, not currently biting.
- Dict insertion order == FASTA file order, and that order is depended on for `ref_id` (`pipeline.py:118`) and for row ordering in the viewer.
- `validate_aligned_fasta` (`fasta_io.py:32-41`) is the only shape check: non-empty, and all sequences the same length. It returns that length.

### 4.2 Focal-set selection

`fasta_io.split_focal_headers` (`fasta_io.py:44-56`).

The focal set is chosen by **plain case-sensitive substring containment on the header string**: `target_string in header`. There is no regex, no field parsing, no delimiter awareness. `split_focal_headers` raises if either side is empty.

The same substring test is re-implemented, independently, in three more places:
- `core.py:215-217` inside `find_dmc_information` (which re-derives focal/non-focal from the full `sequences` dict rather than accepting the already-split lists).
- `core.py:171` inside `compute_metrics`, used to exclude focal sequences from similarity comparison.
- `excel.py:55` inside `build_sheet`, used to exclude non-reference focal sequences from the workbook.

`compute_metrics` at `core.py:171` tests `header == ref_id or target_string in header`. Since `ref_id` is itself a focal header, the `header == ref_id` clause is redundant.

### 4.3 IUPAC ambiguity handling

The single source of truth is `constants.IUPAC` (`constants.py:33-60`): every code maps to its set of possible strict bases; `-` and `?` map to the **empty set**.

The whole codebase encodes "gap or missing" as "`state_possibilities()` returned an empty set". Three separate copies of this predicate exist: `core.py:11-16`, `punishments.py:17-22`, `consensus.py:31-33`.

Two behaviours worth naming explicitly:
- **Any unrecognised character is silently treated as a gap.** `IUPAC.get(state.upper(), set())` returns the default `set()` for `X`, whitespace, digits, or anything else, and every caller then classifies it as empty. Nothing warns.
- `core.score_state` (`core.py:132-141`) is the one lookup that **does not uppercase**: it calls `IUPAC.get(query, set())` directly. In the running pipeline this never matters because `parse_fasta` uppercases everything, but it diverges from the neighbouring `state_possibilities` and would silently score lowercase input as 0.0 if the core were called as a library.

`is_ambiguous_state` (`core.py:19-21`) means "non-empty and not strict ACGT" — so gaps are *not* ambiguous.

`IUPAC_FROM_BASES` (`constants.py:62-66`) is the reverse map (frozenset of bases → code), built by comprehension, excluding empty states. It is used only by consensus generation.

### 4.4 Diagnostic site finding

All in `core.find_dmc_information` (`core.py:188-438`). Three phases.

**Phase 1 — focal consensus base per column (`core.py:230-255`).**
For each column:
1. `core.py:234` — if the **full column across every sequence in the file** (focal *and* non-focal) contains any ambiguous state and `include_ambiguous_dmc_bd` is off, the site is skipped and `skipped_non_acgt` is incremented. This is a whole-alignment veto: a single `N` in one of 1990 non-focal sequences disqualifies the site. On the real dataset this alone kills **243 of 736 sites** (measured).
2. `core.py:238` — `focal_gap_dominated` = focal-column empty fraction ≥ `BALANCING_EMPTY_LIMIT` (2/3).
3. `core.py:240-243` — gaps are tolerated only if `include_gappy_consensus_dmc_sites` is on **and** the column is not gap-dominated.
4. `core.py:245-249` → `consensus_base_for_column` (`core.py:56-105`), which returns a single strict base or `None`. With ambiguity enabled it intersects the possibility sets of all non-empty states and requires the intersection to be exactly one base (`core.py:92-98`); with ambiguity disabled it requires every state to be strict and identical (`core.py:100-105`).

With default options this means: **a candidate site must have a focal column that is entirely strict ACGT, gap-free, and invariant, and the same column must be free of ambiguity codes across the entire alignment.** On the real dataset 226 of 736 sites survive this phase.

**Phase 2 — global-conservation filter (`core.py:261-288`).**
A surviving site is kept as a *candidate* if either:
- the full column contains a non-empty state that does not match the focal base (`column_has_non_empty_variation_from_base`, `core.py:108-125`); or
- `gappy_consensus_signal` (`core.py:275-280`) — the gappy-sites option is on, the full column has gaps, the focal column is not gap-dominated, and there is *no* base variation. That is: a site with no base-level variation is still admitted purely because non-focal sequences have gaps there.
Otherwise the site is discarded as globally conserved.

**Phase 3 — combination search (`core.py:290-382`).**
Matching is defined by three nested closures:
- `sequence_matches_site` (`core.py:293`) → `state_matches_dmc_base` (`core.py:39-53`). With BD off this is exact equality; with BD on it is `base in possibilities`, so an `N` matches every base.
- **Any empty state returns `False`** (`core.py:47-48`). A non-focal sequence with a gap at a site therefore never "matches", so gappy non-focal sequences are automatically treated as discriminated-against. This is a major driver of what looks diagnostic.
- `combination_is_diagnostic` (`core.py:303`) — true iff **no non-focal sequence matches every site in the combination**. Focal sequences are never re-checked against the combination.

The search loop (`core.py:358-382`):
```
for combo_length in range(start_combination_length, max_combination_length + 1):
    for combo in combinations(candidate_sites, combo_length):
        if is_pruned_by_existing_dmc(combo): continue      # superset pruning
        combinations_tested_by_length[combo_length] += 1
        if combination_is_diagnostic(combo): found_this_length.append(combo)
    # results merged only AFTER the length completes
    if found_this_length and combo_length >= min_combination_length:
        stop_reason = "found_at_or_above_minimum_length"; break
```
Key semantics, as implemented:
- **`min_combination_length` is not a floor on which lengths get searched.** The loop always starts at `start_combination_length` (default 1). `min_combination_length` only gates *when the loop is allowed to stop early*. Setting min=3 does not skip lengths 1 and 2; it searches them, records any hits, and keeps going to 3.
- **Pruning is cross-length only.** `is_pruned_by_existing_dmc` (`core.py:337-344`) filters `found_combo` by `len(found_combo) < len(combo)`, and `found_this_length` is merged into `diagnostic_combinations` only after the inner loop finishes (`core.py:371-374`). So combinations found at length *k* never prune other length-*k* combinations, only length >*k* ones.
- **The search stops at the first productive length ≥ min.** With defaults (min=1, max=2), if any single-site DMC exists the search stops at length 1 and **no pairs are ever tested** — `pairs_tested` is 0 and `pairs` is empty. Verified by direct call.
- `stop_reason` is one of `"reached_maximum_length"` (initial value, `core.py:346`), `"no_candidate_sites"`, `"start_length_exceeds_maximum_length"`, `"found_at_or_above_minimum_length"`. Only the first triggers the GUI's continue-prompt.

Resumption: `initial_diagnostic_combinations` / `initial_combinations_tested_by_length` (`core.py:315-335`) let the GUI feed the previous round's results back in. Incoming combos are normalised (`tuple(sorted(set(combo)))`, `core.py:309`) and **dropped if they are not a subset of the current candidate set** (`core.py:321-322`). Because the candidate set is recomputed identically each round from unchanged inputs, this filter is a no-op in practice — but it would silently discard history if the options were changed mid-search.

**Result assembly (`core.py:384-438`):** `single` = length-1 combos flattened; `pairs` = length-2 combos; `unique` = the sorted set of all sites appearing in any recovered combination (`core.py:399-403`). Since the search usually stops at the first productive length, `unique` in practice contains only sites from that one length.

`ambiguous_bd_sites_included` (`core.py:405-409`) re-scans every candidate column against the full alignment — a repeat of work already done in phase 1.

### 4.5 Exhaustive vs. heuristic search paths, and how one is chosen

**There is no heuristic path. Every search in this codebase is exhaustive `itertools.combinations` enumeration.** The only things that bound the work are:
1. The candidate-site filter (phases 1–2) shrinking the input set.
2. Superset pruning across lengths (`core.py:337-344`).
3. The early `break` at the first productive length ≥ min (`core.py:380-382`).
4. The user-set `max_combination_length`, and the GUI's interactive "continue with a new maximum?" prompt (`gui.py:179-206`), which is the only mechanism resembling an adaptive search budget.

The 5-site search (`core.py:441-488`) is separately and unconditionally exhaustive over `C(len(dmc.unique), 5)`, with no pruning and no early exit, and it re-runs on **every** iteration of the GUI continue-loop. It short-circuits only when fewer than 5 unique sites exist (`core.py:454-461`). Each combination calls `compute_metrics`, which iterates the entire sequence dict — so the cost is `C(n,5) × total_sequences`. With the real file's 2354 sequences, 20 unique sites means ~36M state comparisons and 30 sites means ~340M. See §10.

Note the two objectives at `core.py:474-480`: `best_gap_sites` minimises `max_similarity` (the variable named `best_gap_score` actually stores a *max similarity*, not the gap `1 - max_similarity` that `compute_metrics` returns as its third element and which is discarded); `best_avg_sites` minimises mean similarity.

### 4.6 The punishment system

`punishments.find_focal_punishments` (`punishments.py:389-509`). It uses **focal sequences only** — non-focal sequences are never touched (the docstring at `punishments.py:393-398` says so and the code agrees).

**The left-to-right then right-to-left scan (`punishments.py:490-496`):**
```
left_stop       = (alignment_length + 1) // 2
left_positions  = range(0, left_stop)                    # 0 .. mid-1, ascending
right_positions = range(alignment_length - 1, left_stop - 1, -1)   # end .. mid, descending
```
The two ranges are **disjoint and together cover every column exactly once** — so this is not a double pass over the data. Its only purpose is to establish direction-dependent state: `process_scan_positions` (`punishments.py:435`) initialises `terminal_gap_run_open = True` at the start of each half and clears it (`punishments.py:461`) at the first column that is not gap-dominated. Each half therefore detects its own leading run of gap-dominated columns, inwards from that end.

Consequence worth flagging: if one entire half of the alignment is gap-dominated, `terminal_gap_run_open` never clears in that half, and gap-dominated columns adjacent to the **midpoint** are classified PRL rather than INS. The midpoint is a purely arithmetic boundary with no biological meaning, so the PRL/INS boundary in that scenario is an artefact of `(L+1)//2`.

**Category dispatch per column (`punishments.py:438-488`), by focal empty fraction:**

| empty fraction | category | handler |
|---|---|---|
| ≥ 2/3 (`BALANCING_EMPTY_LIMIT`) | `PRL` if the terminal run is still open, else `INS` | `process_gap_dominated_column` |
| < 1/3 (`POLYMORPHISM_EMPTY_LIMIT`) | `POLY` | `process_polymorphism_or_balancing_column`, weight 0.20 |
| 1/3 ≤ f < 2/3 | `BAL` | `process_polymorphism_or_balancing_column`, weight 0.10 |

Note the gap-dominated test runs **first**, so the ≥2/3 band never reaches the POLY/BAL handler.

**Gap-dominated scoring (`punishments.py:331-386`)** — penalises the sequences that *have* a base where most others have none:
- `PRL`: `score = (b_empty / b_any) * PROLONGATION_WEIGHT` (0.10) — `punishments.py:354`.
- `INS`: `score = float(b_any)` — `punishments.py:359`. **A raw count, with no weight applied at all**, unlike every other category. `b_any` here is the number of non-empty states in the column, so a lone insertion in a 364-sequence focal set scores 1.0 while a widely shared one scores proportionally more per sequence.
- The score is charged to every non-empty sequence in the column, and `prl_counts`/`ins_counts` are incremented per sequence.
- The two formulas are **reciprocal**: PRL uses empty/non-empty, INS uses non-empty (undivided). These are deliberately different scales.

**POLY/BAL scoring (`punishments.py:219-328`):**
- Empty-state charge: `empty_score = (b_any / b_empty) * empty_weight` (`punishments.py:240`), charged to each *empty* sequence. Note the ratio is **non-empty over empty**, so a column with a single gap among 364 sequences charges that one sequence `363 × weight`, while a column that is half gaps charges roughly `1 × weight` each. Penalty per gapped sequence therefore *decreases* as gaps become common.
- Base-disagreement charge: `score_non_empty_state` (`punishments.py:156-182`) averages `score_sequence_as_base` over each possible strict interpretation of the sequence's own state, so an `R` is scored as the mean of its A-interpretation and its G-interpretation.
- `score_sequence_as_base` (`punishments.py:126-153`) computes **`b_dif / b_same`**, where counts come from `weighted_counts_for_column` (`punishments.py:86-123`): each ambiguous state in *other* sequences contributes fractionally (`R` = 0.5 A + 0.5 G, `N` = 0.25 each), empty states contribute nothing, and the scored sequence itself contributes a full 1.0 to the assumed base (`punishments.py:111-116`). Because of that forced 1.0, `b_same` can never be zero, so the `if b_same == 0: return 0.0` guard at `punishments.py:149-150` is unreachable.
- Only positive scores generate an event (`punishments.py:315`).

**Benefit-of-the-doubt (`compatible_strict_base_for_bd`, `punishments.py:33-84`):**
BD applies to a POLY/BAL column when all three hold: exactly one distinct strict base appears among non-empty states; every ambiguous non-empty state includes that base in its possibility set; and ambiguous states are at most `BD_AMBIGUOUS_FRACTION_LIMIT` (0.20) of the non-empty states. Empty states are excluded from that fraction.
When BD fires (`punishments.py:244-285`), the function **returns early** (`punishments.py:285`). The effect is that **no base-disagreement score is charged for that column at all** — every sequence is forgiven, not just the ambiguous ones. Ambiguous sequences get a `BD` event with `score=0.0` and a `bd_counts` increment; empty sequences still get their `empty_score`. This is the single largest behavioural lever in the punishment system, because it converts a scored column into an unscored one.

**Score bookkeeping (`add_score_event`, `punishments.py:185-216`):** every charge is added to `total_scores` *and* to the passed-in `category_scores`, and appended to `events`.
Two things follow that matter for interpreting the Excel output:
- Empty-state charges are added to `empty_weight_scores` **and** (via `add_score_event`) to `total_scores` **and** to the POLY or BAL category bucket. So the `EW` column and the `POLY_PS`/`BAL_PS` columns overlap: EW is a subset of POLY_PS + BAL_PS, not a separate addend. `excel.py:197` subtracts EW from the total when computing `PS/bp`, but the `POLY_PS`/`BAL_PS` columns still contain it.
- BD events bypass `add_score_event` (they are appended directly at `punishments.py:273-283`), so they never affect any total.

`events` (`punishments.py:419`) accumulates one `PunishmentEvent` per charge — on the real dataset this is on the order of hundreds of thousands of frozen dataclass instances. **Nothing consumes it.** `write_punishment_excel_report` uses only the aggregate dicts, and `reports.write_text_report` never receives a `PunishmentResult` at all. The event log is built and discarded.

### 4.7 Consensus sequence generation

`consensus.build_focal_consensus_result` (`consensus.py:197-230`), focal sequences only, in `focal_headers` order.

`assign_consensus_state` (`consensus.py:42-79`) — the rule, exactly as coded:
1. Count only strict A/C/G/T among the already gap-filtered column; ties broken by base letter (`consensus.py:64`).
2. No strict bases at all → `"N"` (`consensus.py:60-61`).
3. Exactly one distinct base → that base (`consensus.py:66-67`).
4. Otherwise, if `freq(top) - freq(second) > 0.5` → the top base (`consensus.py:75-76`).
5. Otherwise → the IUPAC code for the **exact set of observed strict bases**, or `"N"` if that set has no code (`consensus.py:78-79`).

Two implementation details that diverge from the plain reading of the docstring at `consensus.py:44-51`:
- Ambiguity codes in the input are **not** decomposed. An `R` is neither counted as A nor as G; it is simply ignored by the `state in STRICT_BASES` filter at `consensus.py:57`. So a column of 100 `R`s and one `A` yields consensus `A`.
- The gap filter upstream (`consensus.py:94-98`, `:180-184`) uses `is_empty_state`, so `-`, `?` and any unrecognised character are dropped, but ambiguity codes reach `assign_consensus_state` and are then discarded there.

Trimming: `classify_prl_ins_sites` (`consensus.py:110-157`) **re-implements the punishment scan** — same `(L+1)//2` split, same `BALANCING_EMPTY_LIMIT`, same `terminal_gap_run_open` flag — to decide which columns are PRL or INS, and `build_trimmed_consensus` (`consensus.py:160-194`) drops all of them. This is a second, independent copy of the scan logic in `punishments.py:435-496`; the two are currently identical in behaviour but are not shared code, so they can drift.

`write_consensus_text_report` (`consensus.py:233-345`) emits the settings, the rule text, counts, the removed PRL/INS site lists (1-based), both consensus sequences wrapped at 80 columns, and a DMC-site → trimmed-position mapping. It is called with `dmc_sites=dmc.unique` (`pipeline.py:176`), so the mapping section reflects the DMC search that just ran.

### 4.8 Unique-sequence deduplication / clustering

`sequence_subsets.group_sequence_subsets` (`sequence_subsets.py:25-109`), focal sequences only (`pipeline.py:222-229`).

What it actually does:
1. Ungap every sequence — `sequence.replace("-", "")` (`sequence_subsets.py:21`). **Only the `-` character is removed**; `?` and `N` survive into the "ungapped" string and participate in containment testing.
2. Sort descending by ungapped length, ties broken by header (`sequence_subsets.py:40`).
3. Greedily assign each sequence to the first existing group whose **representative** string contains it as a substring (`sequence_subsets.py:53`); otherwise start a new group with this sequence as representative (`sequence_subsets.py:73-80`). Empty ungapped sequences group with other empty ones (`sequence_subsets.py:63-71`).

This is substring containment, not equality — exact duplicates are a special case. Because containment is not an equivalence relation and the assignment is first-match greedy, **the grouping depends on input order**, which is FASTA order after the length sort. The tie-break on header name makes it deterministic for a fixed input.

This differs from the legacy `cluster_sequences` (`consensus_dmc_pipeline.py:191-219`), which tested containment in *both* directions and promoted a longer member to representative mid-pass. The new version relies on the descending sort to make the first member the longest, and never promotes.

### 4.9 Report and Excel output generation

**`DMCs_output.txt`** — `reports.write_text_report` (`reports.py:14-131`). Settings, sequence summary, run diagnostics, combinations-tested-by-length, combinations grouped by length, the unique-site list (printed twice, at `reports.py:97` and `reports.py:100`, under different headings — "Unique diagnostic sites" and "Full diagnosis"), and the two 5-site selections. All site numbers are 1-based via `+ 1` in the formatters (`core.py:495`, `core.py:503`, `reports.py:11`). States printed are the **consensus** states from `dmc.states`, not any one sequence's bases.

**`focal_consensus_output.txt`** — `consensus.write_consensus_text_report` (§4.7).

**`comparison_output.xlsx`** — `excel.write_excel_report` (`excel.py:97-142`) → `build_sheet` (`excel.py:28-94`) per sheet (`Full`, then `Gap5`/`Avg5` if available).
- Row filter (`excel.py:54-56`): keep the row if it is the reference **or** not focal. So the sheet contains one focal sequence (the reference) plus every non-focal sequence.
- Baseline states (`excel.py:47-50`): when `diagnostic_states` is supplied — and `pipeline.py:187` always supplies it — the comparison baseline is the **consensus** state per site, *not* the reference sequence's own bases. But each row's displayed cells come from `extract_sites(sequence, sites)` (`excel.py:58`), i.e. that sequence's actual bases. **The reference row is therefore scored against the consensus rather than against itself, and can show less than 100% similarity to its own row.** This is a genuine divergence between what the row is labelled and what its score means.
- Sorting (`excel.py:69-70`): reference row pinned first, others by descending match score, then descending similarity, then header.
- The dataset-wide average similarity is written **only into row 2** (`excel.py:73`), leaving the rest of that column blank.
- Cell fills come from `constants.COLORS`, defaulting to white for unknown states (`excel.py:86`).

**`punishment_output.xlsx`** — `excel.write_punishment_excel_report` (`excel.py:144-233`). Columns `ID, PS, PS/bp, EW, BD, INS, PRL, POLY_PS, BAL_PS, INS_PS, PRL_PS`; rows sorted by descending total score then header (`excel.py:184-190`). `PS/bp = (PS - EW) / alignment_length` (`excel.py:197`) — see §4.6 on the EW/POLY_PS overlap. `alignment_length` here is the **full** alignment length from `load_inputs`, not the number of scored columns.

**`sequence_subsets_output.xlsx`** — `sequence_subsets.write_sequence_subset_excel_report` (`sequence_subsets.py:138-217`). One representative row per group (bold, green fill) followed by its subset rows; sequences longer than Excel's 32,767-character cell limit are truncated with a marker (`sequence_subsets.py:124-135`).

### 4.10 Where the README diverges from the code

`README.md` claims "Reports 1-site and 2-site diagnostic characters" and "Searches 5-site combinations for similarity optimisation". The 2-site claim is now conditional: with defaults, if any 1-site DMC exists the search stops at length 1 and no pairs are computed (§4.4). The README does not mention the punishment system, the consensus output, the sequence-subset output, or the configurable min/max combination lengths — all of which exist.

---

## 5. UI LAYER (1) - the legacy Tkinter app

This section describes `gui.py` and `viewer.py`, which are the LEGACY UI. They
still run, and the analysis below is still accurate for them. The current UI is
section 5b.

### Which files are tkinter

| File | Role |
|---|---|
| `molecular_diagnosis/gui.py` | Main window. `import tkinter as tk`, `filedialog`, `messagebox`, `simpledialog` (`gui.py:2-3`). |
| `molecular_diagnosis/viewer.py` | **The custom alignment viewer.** `tkinter`, `tkinter.font`, and PIL (`viewer.py:1-5`). |
| `molecular_diagnosis/viewer_old.py` | Legacy viewer, unimported (`viewer_old.py:1-2`). |
| `pipeline_verbose.py` | Legacy monolith with its own GUI + inline viewer (`pipeline_verbose.py:1-3`). |

**No tkinter import exists in any live non-UI module** — `core`, `punishments`, `consensus`, `sequence_subsets`, `excel`, `reports`, `pipeline`, `fasta_io`, `utils`, `models`, `constants` are all clean. That is the one clean boundary in the codebase.

---

### (a) Overall architecture of the viewer

`viewer.py:10` `open_fasta_viewer(root, sequences, *, color_bases=True, show_letters=True, show_grid=False)`.

One module-level function containing everything: a `tk.Toplevel`, a 4-canvas grid, and roughly two dozen nested closures sharing state through `nonlocal`. There is no class, no model object, no separation between view state and view logic. It is called exactly once, from `gui.py:83`, and returns immediately after wiring up bindings — the window then lives inside the parent's `mainloop`.

Layout (`viewer.py:68-107`) is a 3×3 `grid` inside a `Frame`:

```
            column 0            column 1                 column 2
row 0   corner_canvas       top_canvas (ruler)             —
row 1   name_canvas         seq_canvas (the alignment)   y_scroll
row 2       —               x_scroll                       —
```
Only row 1 / column 1 is given weight (`viewer.py:106-107`), so the sequence pane absorbs all resizing while the name column and ruler keep their fixed extents.

Fixed geometry, computed once at open time from font metrics (`viewer.py:47-66`):
- `seq_font` = Courier New 16; `seq_bold_font` = the same, bold.
- `base_width = max(seq_font.measure("W"), seq_bold_font.measure("W"))` → **13 px** on this machine (both fonts measure 13; the `max` is defensive).
- `row_height = max(36, linespace + 12)` → **36 px** (linespace 23).
- `name_width` = widest header, measured in the default UI font, + 16.
- `header_height` = 42, hardcoded.
- `total_width = seq_length * base_width`, `total_height = n_rows * row_height` — the virtual canvas extent, never materialised.

The three panes are kept in register purely by shared scalars: the ruler consumes `first_col`/`x_remainder`, the name column consumes `first_row`/`y_remainder`, and the sequence pane consumes all four. There is no Tk-level scroll linking between the canvases.

### (b) How alignment data is stored in memory

Two parallel Python lists, built once at `viewer.py:36-37`:
```python
headers = list(sequences.keys())
seqs    = [str(sequences[h]) for h in headers]
```
That is the entire model. No per-column index, no colour matrix, no tile cache, no precomputed image. Row *i* is `seqs[i]`, a plain `str`; column *j* of row *i* is `seqs[i][j]`.

`seq_length = max(len(seq) for seq in seqs)` (`viewer.py:40`) — the viewer takes the **maximum** length and pads shorter rows on the fly, so unlike the analysis path it does **not** require an aligned FASTA. `padded_slice` (`viewer.py:191-202`) right-pads any short row with `-` for display.

The viewer's entire interface to the rest of the system is the `Mapping[str, str]` it is handed, plus `constants.COLORS`. It has no knowledge of DMC sites, punishment scores, consensus, or the focal set, and it produces no output that flows back. It is a pure read-only sink.

Mutable view state (all `nonlocal` closure variables, `viewer.py:109-121`):
`x_offset`, `y_offset` (pixel scroll position), `hovered_row`, `hovered_col`, `redraw_pending`, `selection_anchor`, `selection_active`, plus `color_cache` and the `seq_canvas._viewport_image` reference.

### (c) How rendering works

This is the part that matters most for a Canvas port. The design is **"render the viewport into one PIL image, then blit it as a single canvas item."**

`redraw()` (`viewer.py:503-542`) per frame:
1. `clamp_offsets()` (`viewer.py:146`) — clamp scroll offsets to `[0, total - viewport]`.
2. `update_scrollbars()` (`viewer.py:158`) — push fractions into the scrollbars manually.
3. `visible_range()` (`viewer.py:176-189`) — integer-divide the offsets by cell size to get `first_row/last_row/first_col/last_col`, **plus two sub-cell remainders** `x_remainder`/`y_remainder`. The `+2` on each `last_*` (`viewer.py:181`, `:184`) is the overdraw margin that keeps partially-scrolled edge cells filled.
4. `draw_header(...)` and `draw_names(...)` — ordinary canvas items on the small panes.
5. `seq_canvas.delete("all")` — the sequence pane is fully cleared every frame.
6. `render_sequence_bitmap(...)` (`viewer.py:361-468`) — allocate a fresh `Image.new("RGB", (viewport_w, viewport_h), "white")`, draw filled rectangles into it with `ImageDraw`, wrap in `ImageTk.PhotoImage`, return.
7. Store the photo on the widget (`viewer.py:530`) and place it with a single `create_image(0, 0, anchor="nw")`.
8. `draw_sequence_letters(...)` (`viewer.py:470-501`) — **one `create_text` item per visible row**, each carrying the whole visible substring, positioned at `x = -x_remainder`.

So a viewport showing 12 rows × 90 columns costs **1 image item + 12 text items**, not 1080 rectangles + 1080 text items. The docstring at `viewer.py:18-27` describes exactly this and matches the implementation.

The colour path (`viewer.py:384-405`): when `color_bases` is false the row is filled white in one rectangle; otherwise it iterates the visible substring and draws one `draw.rectangle` per residue into the PIL image, with an optional 1px `(230,230,230)` outline when `show_grid` is on. Colours come from `parse_color` (`viewer.py:123-144`), which resolves `constants.COLORS` hex strings to RGB tuples and memoises them in `color_cache`.

Two rendering details that look accidental:
- `is_hovered` is computed at `viewer.py:382` inside the bitmap renderer and then **never used**. The row-hover highlight is not painted into the bitmap at all; the only visible hover feedback in the sequence pane is the bold font swap at `viewer.py:499`. (The name pane does draw its own hover fill, `viewer.py:329`.)
- The hovered-column guide lines (`viewer.py:444-460`) are nested **inside** the `if selected is not None:` block. Column guides therefore appear only while a selection rectangle exists, which is unlikely to be the intent given they track `hovered_col`, not the selection.

Text/cell alignment depends on Courier New being monospaced and on `base_width` equalling the font's advance width. Since the letters for a row are drawn as one string, Tk lays them out at the font's own advance; the coloured cells are laid out at `base_width`. These agree at 13px here. Note that the hovered row is drawn in the **bold** font while its cells keep the same `base_width` — this only stays aligned because bold Courier New has the same advance as regular on this machine (both measured at 13px). On a system where the bold advance differs, the hovered row's letters would drift out of their cells. (UNVERIFIED as a cross-platform claim — would need the same `font.measure("W")` comparison run on macOS/Linux with whatever Courier New substitute is resolved there.)

### (d) How scrolling works

Fully hand-rolled. The canvases have **no `scrollregion`**, no `xview`/`yview` delegation, no `xscrollcommand`. The viewer keeps `x_offset`/`y_offset` in pixels and implements the scrollbar protocol itself:

- `xview(*args)` / `yview(*args)` (`viewer.py:544-588`) are installed as the scrollbars' `command` (`viewer.py:590-591`). They parse the Tk protocol manually: `"moveto"` → `offset = fraction * total`; `"scroll"` with `"pages"` → one viewport minus one cell; `"scroll"` with units → `base_width * 5` horizontally, `row_height * 3` vertically.
- `update_scrollbars()` (`viewer.py:158-174`) computes `first`/`last` fractions and calls `scrollbar.set` directly.
- Mouse wheel: `on_mousewheel` (`viewer.py:648`) normalises `event.delta` through `wheel_units` (`viewer.py:637-646`, handling both the ±120 Windows convention and finer deltas), then moves 8 cells horizontally or 3 rows vertically. `Button-4`/`Button-5` handlers (`viewer.py:664`, `:675`) cover X11. Shift **or** Ctrl selects horizontal scrolling (`wants_horizontal_scroll`, `viewer.py:632-635`).
- Every handler returns `"break"` to stop Tk's default canvas scrolling from also firing.

**There is no zoom.** Font size is hardcoded at `viewer.py:47`, `base_width`/`row_height` are computed once, and no binding changes them. If the port is expected to have zoom, it is new functionality, not a port. (The `show_letters`/`color_bases`/`show_grid` keyword arguments are also fixed at call time — `gui.py:83` passes none of them, so the defaults `True/True/False` always apply and no UI exposes them.)

### (e) Performance optimisations actually present

1. **Viewport virtualisation.** Only `first_row..last_row` × `first_col..last_col` is ever touched. Cost is bounded by window size, not alignment size — a 2354 × 736 alignment draws the same amount as a 10 × 50 one.
2. **Single-bitmap compositing.** All residue backgrounds become one `PhotoImage` and one canvas item, avoiding thousands of canvas objects. This is the headline optimisation and the reason `viewer_old.py` was replaced.
3. **One text item per row**, not per residue (`viewer.py:494`), with the sub-cell offset absorbed by the `-x_remainder` origin.
4. **Coalesced redraws.** `request_redraw` (`viewer.py:204-217`) sets a `redraw_pending` flag and schedules via `viewer.after_idle`, so a burst of motion/scroll/configure events produces exactly one repaint per idle cycle.
5. **Colour memoisation.** `color_cache` (`viewer.py:121`) parses each hex string once.
6. **Sub-pixel-accurate integer scrolling.** Keeping the fractional remainder separately (`viewer.py:186-187`) means scrolling is smooth per-pixel while iteration stays per-cell.
7. **Overdraw margin of +2 cells** (`viewer.py:181`, `:184`) instead of clipping maths.
8. **Early-out per row and per cell** on off-screen geometry (`viewer.py:379-380`, `:393-394`).
9. **Hover diffing** — `on_motion` (`viewer.py:602`) only requests a redraw when the hovered cell actually changes.

Costs that remain: a full-viewport `Image.new` + `PhotoImage` allocation on **every** frame (no image reuse, no dirty-rect updates), and `seq_canvas.delete("all")` + re-creation of every text item per frame.

### (f) User-facing functionality the viewer currently supports

- Scrollable alignment grid with a frozen header row (position ruler) and frozen left column (sequence names).
- Per-base background colouring from the 16-entry `COLORS` table, covering strict bases, all 2- and 3-base ambiguity codes, `N`, `-` and `?`.
- Monospaced sequence letters overlaid on the colour cells.
- A position ruler with minor ticks every 5 columns and numbered major ticks every 10 (`viewer.py:255-300`), with edge-clipping suppression so half-visible numbers are not drawn (`viewer.py:286`).
- Row hover: the name cell fills grey and both the name and the sequence letters switch to bold (`viewer.py:329`, `:344`, `:499`).
- Rectangular cell-range selection by click-drag, drawn as a 2px black outline (`viewer.py:462-466`); anchor/active tracked at `viewer.py:713-741`.
- Column guide lines flanking the hovered column — but only while a selection exists, see (c).
- Wheel scrolling (vertical, and horizontal with Shift/Ctrl), scrollbar dragging, page scrolling; Windows and X11 wheel conventions both handled.
- Graceful empty state: a "No sequences to display." label (`viewer.py:29-34`).
- **Click emits `print(...)` to stdout** (`viewer.py:726`) — the only user-visible feedback for a click is a console line the user of a GUI app will never see. There is no status bar, no copy, no export, no search/goto, no zoom, no column/row sorting, no highlighting of DMC sites.

### (g) What ports to HTML Canvas, and what is tkinter-specific

**Ports essentially unchanged (the concepts, not the API):**
- The whole virtualisation scheme: `x_offset`/`y_offset` in pixels, `visible_range()`'s integer-divide plus remainder, the +2 overdraw margin. This is the same maths a Canvas implementation needs.
- The data model (`headers` + `seqs` lists, `padded_slice`) — plain strings, no Tk types.
- `parse_color` and the `COLORS` table — hex strings are more natural in a browser than they are here.
- The frozen-header/frozen-name-column layout as three synchronised viewports driven by shared scalars.
- The `request_redraw` idle-coalescing pattern → `requestAnimationFrame`, which does the same job natively and better.
- The wheel-delta normalisation intent (though the ±120 arithmetic in `wheel_units` is a Windows/Tk artefact; browsers expose `deltaMode`).
- Hit-testing (`cell_from_event`, `viewer.py:704-711`) — identical arithmetic.
- The ruler tick logic (every 5 / every 10, with edge suppression).
- Selection rectangle bookkeeping (anchor/active, min/max normalisation at `viewer.py:347-359`).

**Directly better on Canvas, i.e. the workaround disappears:**
- The PIL-image-per-frame trick exists *because* Tk canvas items are expensive. A Canvas 2D context draws rectangles and text into one buffer natively, so the `Image.new` / `ImageDraw` / `ImageTk.PhotoImage` / `create_image` chain collapses into direct `fillRect`/`fillText` calls, and the per-frame allocation goes away.
- The `seq_canvas._viewport_image` monkey-patched reference (`viewer.py:119`, `:530`) exists solely to defeat Tk's image garbage collection. It has no analogue.
- Drawing letters as one string per row was a workaround for per-item cost; on Canvas, per-row `fillText` is still the right call for a different reason (text shaping cost), so the technique survives but for a new reason.

**Tkinter-specific, needs replacing rather than porting:**
- The manual scrollbar protocol (`xview`/`yview` parsing `"moveto"`/`"scroll"`/`"pages"`, and `scrollbar.set(first, last)`) — replaced by native overflow scrolling or a custom scrollbar component.
- `tkfont.Font.measure()` / `.metrics("linespace")` for geometry → `ctx.measureText` / font metrics, with different rounding. The `max(normal, bold)` advance trick (see (c)) needs rethinking because web font stacks are far less likely to give bold and regular the same advance.
- The 3×3 `grid` with `grid_rowconfigure/columnconfigure` weights → CSS grid/flex.
- `Toplevel`, `Scrollbar`, `after_idle`, `bind`, the `"break"` return convention, `<Button-4>`/`<Button-5>`, `event.state & 0x0001|0x0004` modifier bit masks, `winfo_width()`.
- `print()` on click (`viewer.py:726`).

**Present in neither and would be new work:** zoom, any linkage between the viewer and analysis results, copy/export, virtualised rendering of more rows than fit in a browser's max canvas dimension.

---

## 5b. UI LAYER (2) - the Electron/React desktop app

The current UI. No Python is imported by the renderer, no SQL exists anywhere in
`desktop/`, and the Electron main process holds no project state - it forwards.

### Process boundaries

    renderer (React)  ->  window.desktop.*  (preload, contextIsolated)
                      ->  ipcRenderer/ipcMain  (main.ts)
                      ->  pythonBridge  (newline-delimited JSON over stdio)
                      ->  python -m molecular_diagnosis.service

stdout is protocol-only; the service redirects its own `sys.stdout` to stderr,
and the bridge attaches captured stderr to failures as diagnostics.

### Screens

| Screen | File | Notes |
|---|---|---|
| Launcher | `screens/LauncherScreen.tsx` | Create New / Open Existing, plus Recents and Browse. Recents needs a store nothing writes yet, and says so. |
| Project page | `screens/ProjectScreen.tsx` | ONE screen for both project states. Before Create the title field carries a Create button and the FASTA rows are vetted CANDIDATES; after Create the title becomes a heading with a rename pencil, SELECT ANALYSES appears, and the SAME table starts linking into the project. No navigation happens on Create. |
| Molecular Diagnosis | `screens/MolecularDiagnosisScreen.tsx` | Design pages 4 and 6. A never-saved set shows the labelled title input; a saved one shows the heading with SAVE. Same component. |

### The FASTA table

`components/sources/SourceTable.tsx` renders one `SourceRow` for both kinds of
row. A candidate has no `fastaFileId`, so its lock is held locally by path in
`pending.lockedPaths` and applied through `project.setFastaFileLocked` as soon
as Create links the file. Everything else about the row - the name and pencil
inline group, the seq/bp/PIS columns with first-row-only unit labels, the hover
reveal, the in-row removal confirmation - is identical either side of Create.

Deferred and drawn as inert rather than faked: source rename, source reordering,
and PIS (shown as an em dash because nothing computes it).

### The focal set

- `components/focal/FocalSetEditor.tsx` - CodeMirror 6 (`@codemirror/state` +
  `@codemirror/view`, no language modes, no history extension). Entries are
  complete exact headers separated by `;`. Manual editing never expands a
  substring; `+` is the thing that searches. Undo/redo belong to the DRAFT, so
  one stack covers typing, `+` and `-`.
- `components/focal/FocalSetLibrary.tsx` - lists WORKING DRAFTS, not database
  rows, so a set the user started is visible before it is saved. Every row
  follows one control rule: hidden at rest, grey on row hover, lit under the
  pointer, and a locked row keeps its lock and loses rename and delete.

### The right-hand column

`components/alignment/WorkspaceRightPane.tsx` renders the display selectors, the
optional focal-set library, and the visualizer shell as ONE flex column
overlaying the analysis controls. The viewer is the item after the pane, so its
top follows the library's real height; its left edge is a margin percentage
driven by the drag snaps. There is exactly one left rail, and that rail is the
drag handle. `components/alignment/AlignmentPlaceholder.tsx` draws only what is
inside the frame - there is no sequence renderer yet, deliberately.

### Scaling for Windows DPI

`app/uiScale.ts` sets `--ui-zoom` to `min(1, innerWidth / 1920)`, and
`global.css` applies it to `#root` with a compensating width/height. At 1920 CSS
px it is exactly 1 and the reference layout is untouched; at the 1536x864
viewport a 1920x1080 screen reports at 125% scaling it is 0.8, which reproduces
the reference layout at the same physical size instead of reflowing it into a
second design nobody drew.


## 6. UI/LOGIC ENTANGLEMENT

The list to work from when extracting the core. Ordered roughly by how much it will get in the way.

1. **Search-control flow lives in a GUI callback.** `gui.py:154-206` implements the entire iterative-deepening protocol: initialise `start_combination_length`, call the core, inspect `result.dmc.stop_reason`, ask the user via `messagebox.askyesno`, prompt for a new bound via `simpledialog.askinteger`, thread `initial_diagnostic_combinations` and `initial_combinations_tested_by_length` into the next call. **There is no headless equivalent** — a library user gets one pass and must re-implement the loop.

2. **Every parameter's only home is a Tk variable.** `gui.py:15-23` — `fasta_var`, `target_var`, `output_dir_var`, `include_ambiguous_dmc_bd_var`, `include_gappy_consensus_dmc_sites_var`, `min_combination_length_var`, `max_combination_length_var`. There is no config object, no dict, no dataclass. Defaults are encoded in `tk.BooleanVar(value=False)` / `tk.IntVar(value=1)` constructor arguments and are duplicated as Python defaults on `run_pipeline_core` (`pipeline.py:96-100`) and `find_dmc_information` (`core.py:192-196`) — three places, currently in agreement.

3. **Validation is triplicated across the UI/core boundary.** Empty-string checks at `gui.py:123-133` and again at `pipeline.py:48-55`; min/max combination checks at `gui.py:139-152` and again at `core.py:200-210`. Path existence is checked only in `pipeline.py:60-70`.

4. **`print()` inside the viewer.** `viewer.py:726` writes click details to stdout from inside an event handler. This is the only `print` in live code.

5. **The report writer is a side-effecting function, not a formatter.** `reports.write_text_report` (`reports.py:14`) takes eleven arguments and writes straight to a file handle; there is no way to obtain the report as a string. Same shape in `consensus.write_consensus_text_report` (`consensus.py:233`), `excel.write_excel_report` (`excel.py:97`), `excel.write_punishment_excel_report` (`excel.py:144`), `sequence_subsets.write_sequence_subset_excel_report` (`sequence_subsets.py:138`).

6. **Computation and serialisation are fused in the pipeline.** `run_pipeline_core` (`pipeline.py:91`) cannot compute without writing: file-path reservation (`pipeline.py:149-153`) and three writes (`pipeline.py:155-188`) are unconditional and interleaved with the analysis. There is no "compute only" entry point. `PipelineResult` returns **paths**, not data — except `dmc`, which is the sole computed object that escapes. `five_site_result` and `consensus_result` are computed, written to disk, and then dropped.

7. **Domain models scattered outside `models.py`.** `ConsensusResult` is defined in `consensus.py:13`; `SequenceSubsetGroup` in `sequence_subsets.py:10`. Both files also contain their own Excel/text writers, so model + algorithm + serialisation share a module.

8. **Output filenames are hardcoded constants**, not parameters — `constants.py:1-5`, consumed at `pipeline.py:149-153` and `pipeline.py:231-237`. A caller can choose the directory but never the filename.

9. **Hardcoded absolute path** (legacy only): `consensus_dmc_pipeline.py:5` embeds `/home/spinoulis/Desktop/Platygastroidea/...`. Also `DMC_SITES` hardcoded at `consensus_dmc_pipeline.py:8-9`, output names at `:12-15`, and the focal string `"Leptacis_tipulae"` at `:291`.

10. **Exception handling as UI.** `gui.py:208-209` and `gui.py:242-243` swallow every exception into `messagebox.showerror(..., str(error))`. All errors from the core are raised as bare `ValueError` with prose messages (`fasta_io.py:34`, `:39`, `:52`, `:55`; `pipeline.py:49-70`; `core.py:201-210`, `:220-223`, `:179`; `consensus.py:86`, `:123`, `:166`, `:201`, `:206`; `punishments.py:113`, `:364`, `:401`; `sequence_subsets.py:51`, `:91-94`), so callers can only distinguish failures by matching message text — which the tests already do (`tests/test_fasta_io.py:43`, `:53`, `:75`, `:85`).

11. **No progress reporting or cancellation.** The core is synchronous and silent; the GUI freezes for the duration of a run (see §12). Any extracted library will need a callback/generator seam that does not exist today.

12. **Redundant parses driven by UI structure.** `view_fasta` (`gui.py:64`), `run_pipeline_core` (`pipeline.py:72`) and `run_punishment_core` (`pipeline.py:72` via `load_inputs`) each parse the FASTA independently. Nothing caches. The GUI continue-loop re-parses on every iteration.

**Not a problem:** there is no global mutable state in the science modules, and no tkinter import in any non-UI live module. The entanglement is concentrated in `gui.py` and in the write-fused pipeline, not sprayed through the algorithms.

---

## 7. CONFIGURATION AND INPUTS

Every user-supplied parameter. There is no config file, no environment variable, no CLI flag.

| Parameter | Read at | Validated at | Consumed at | Default | Notes |
|---|---|---|---|---|---|
| FASTA file path | `gui.py:15` `fasta_var`; set by `browse_fasta` (`gui.py:26-36`) or typed into the entry (`gui.py:246`) | non-empty `gui.py:123`, `pipeline.py:48`; exists/is-file `pipeline.py:60-64` | `parse_fasta` (`pipeline.py:72`) | none | Also read independently by `view_fasta` (`gui.py:57`). |
| Focal identifier string | `gui.py:16` `target_var`, entry at `gui.py:264` | non-empty `gui.py:127`, `pipeline.py:51` | `split_focal_headers` (`fasta_io.py:48-49`), `find_dmc_information` (`core.py:215-217`), `compute_metrics` (`core.py:171`), `build_sheet` (`excel.py:55`), report headers, consensus FASTA header names (`consensus.py:316`, `:321`) | none | Case-sensitive substring. No validation that it is a sensible taxon string. Empty focal set or empty non-focal set raises (`fasta_io.py:51-55`). |
| Output directory | `gui.py:17` `output_dir_var`; `browse_output` (`gui.py:38`) or `use_fasta_location` (`gui.py:47`) | non-empty `gui.py:131`, `pipeline.py:54`; exists/is-dir `pipeline.py:66-70` | `next_available_filename` (`pipeline.py:149-153`, `:231-237`) | none | Not created if missing — it must already exist. |
| Benefit of doubt to ambiguous bases | `gui.py:19` `BooleanVar`, checkbox `gui.py:269-273` | none | `gui.py:164` → `pipeline.py:96` → `core.py:192`; used at `core.py:234`, `:247`, `:270`, `:297`, and via `state_matches_dmc_base` (`core.py:50-53`) and `consensus_base_for_column` (`core.py:84`, `:92`) | `False` | Label reads "Give benefit of doubt to ambiguous bases". Distinct from, and unrelated to, the punishment BD rule (`punishments.py:33`), which has no toggle and is always on. |
| Ignore gaps | `gui.py:20` `BooleanVar`, checkbox `gui.py:275-279` | none | `gui.py:165` → `pipeline.py:97` → `core.py:193`; used at `core.py:241`, `:276` | `False` | Variable name is `include_gappy_consensus_dmc_sites`; the UI label is "Ignore gaps". The two describe the same switch from opposite directions. |
| Minimum combination length | `gui.py:22` `IntVar`, Spinbox 1–20 (`gui.py:286-292`) | `gui.py:139`, `:147`; `core.py:200`, `:209` | `core.py:380` — the early-stop gate only | `1` | Does **not** set the starting length; see §4.4. |
| Maximum combination length | `gui.py:23` `IntVar`, Spinbox 1–20 (`gui.py:294-302`) | `gui.py:143`, `:147`; `core.py:203`, `:209` | `core.py:358` loop bound | `2` | Spinbox caps at 20; the mid-search dialog caps at 100 (`gui.py:197`). |
| New maximum (mid-search) | `simpledialog.askinteger` (`gui.py:192-198`) | `minvalue=next_start`, `maxvalue=100` | `gui.py:204` → next `run_pipeline_core` call | `stopped_at_length + 1` | Only reachable when `stop_reason == "reached_maximum_length"`. |
| `start_combination_length` | not exposed in the UI | `core.py:206`, `:353` | `core.py:358` | `1` | Set programmatically by the continue-loop (`gui.py:203`). |
| `initial_diagnostic_combinations` | not exposed | subset check `core.py:321` | `core.py:315-330` | `None` | Continue-loop only (`gui.py:205`). |
| `initial_combinations_tested_by_length` | not exposed | none | `core.py:334-335` | `None` | Continue-loop only (`gui.py:206`); affects reported counts only. |

**Compile-time constants that behave like parameters but are not exposed anywhere** (`constants.py:70-82`): `POLYMORPHISM_EMPTY_LIMIT` (1/3), `BALANCING_EMPTY_LIMIT` (2/3), `POLYMORPHISM_EMPTY_WEIGHT` (0.20), `BALANCING_EMPTY_WEIGHT` (0.10), `PROLONGATION_WEIGHT` (0.10), `BD_AMBIGUOUS_FRACTION_LIMIT` (0.20). `BALANCING_EMPTY_LIMIT` is doing triple duty: punishment PRL/INS classification (`punishments.py:442`), consensus trimming (`consensus.py:140`), and the DMC gap-dominated test (`core.py:238`, `:265`). Also unexposed: the consensus dominance margin `0.5` (hardcoded at `consensus.py:75`), the FASTA wrap width 80 (`consensus.py:24`), and the 5-site combination size (hardcoded `5` at `core.py:454`, `:463`).

The punishment run takes **no** tunables at all: `run_punishment_core` (`pipeline.py:198-202`) accepts only path/target/output and ignores both checkbox options.

---

## 8. DEPENDENCIES

### Third-party packages imported by LIVE code

| Package | Imported at | Used for | Installed in `.venv` |
|---|---|---|---|
| `openpyxl` | `excel.py:3-5`, `sequence_subsets.py:5-7` | All three Excel outputs | 3.1.5 |
| `Pillow` (`PIL`) | `viewer.py:5` | `Image`, `ImageDraw`, `ImageTk` — the viewer's bitmap renderer | 12.2.0 |
| `pytest` | `tests/test_fasta_io.py:3` | Test runner and `pytest.raises` | 9.0.3 |

`tkinter` is stdlib. Everything else in the live code is stdlib (`pathlib`, `itertools`, `collections`, `dataclasses`).

### Declared vs. used

- `requirements.txt` lists `openpyxl` and `pillow`. Both are used. **Neither is pinned** — no version specifiers at all.
- `pytest` is **used but not declared** in `requirements.txt` (it is installed in `.venv`, and `pyproject.toml` configures it). Anyone installing from `requirements.txt` alone cannot run the tests.
- `pyproject.toml` declares **no dependencies whatsoever** — it contains only `[tool.pytest.ini_options]`. The project is not installable as a package.

### Imported by LEGACY code, not installed

- `python-docx` (`consensus_dmc_pipeline.py:2-3`) — **not installed** in `.venv`, **not in `requirements.txt`**. Verified: `import docx` raises `ModuleNotFoundError` under both `.venv/Scripts/python.exe` and the system interpreter. `consensus_dmc_pipeline.py` cannot be executed as-is.

### Installed but not imported by any project code

`colorama`, `iniconfig`, `packaging`, `pluggy`, `pygments` — all transitive pytest dependencies. `et_xmlfile` is a transitive openpyxl dependency. None are directly imported.

### Environment notes

- `.venv` is Python 3.12.1. The code uses PEP 604 unions (`str | Path`) and PEP 585 generics (`dict[str, str]`) throughout, so **Python ≥ 3.10 is required**; nothing declares this.
- The `python` on PATH is `C:\Python312\python.exe`, **not** the venv interpreter, and it has its own (older, 10.2.0) Pillow. Anyone running `python main.py` without activating `.venv` gets a different dependency set than the one recorded here.

---

## 9. TESTS AND DATA

### Python: 315 tests, all passing

Collected under `pyproject.toml`'s `testpaths = ["tests"]`.

| File | Tests | Covers |
|---|---|---|
| `tests/test_core.py` | 10 | `extract_sites`, `score_state`, `compute_similarity`, `compute_match_score`, `find_dmc_information`, `format_diag`. |
| `tests/test_fasta_io.py` | 8 | Parsing, aligned-length validation, header splitting. |
| `tests/test_focal.py` | 29 | The single focal matcher: substring containment, OR across selectors, trimming, deduplication, refusal of an empty selector. |
| `tests/test_pipeline.py` | 2 | `run_pipeline_core` writes both outputs; `next_available_filename` produces `(2)` names. |
| `tests/test_parity.py` | 13 | Runs `tests/parity_driver.py` against the current tree and against a read-only worktree of the pre-integration commit, and diffs the results: focal selection, every `DMCResult` field, five-site optimisation, consensus, both text reports, the logical contents of `comparison_output.xlsx`, and continuation/resume. Two of them drive `service.dispatch` - the route Electron actually takes. |
| `tests/test_project_db.py` | 21 | Schema creation, the migration ladder, `user_version`. |
| `tests/test_project_sources.py` | 45 | Linking, the six source states, re-indexing, relinking, duplicate-header refusal. |
| `tests/test_project_sources_lock.py` | 19 | The persisted per-source lock: it is stored, it survives reopening, it refuses unlinking, and it does NOT make a file unanalysable. |
| `tests/test_project_focal_api.py` | 33 | Focal sets: create, rename, lock, delete, entry replacement as a diff. |
| `tests/test_project_draft_api.py` | 29 | `saveFocalSet`, presence, `resolveFocalAddQuery`, `matchFocalHeaders` - the calls the working-draft UI is built on. |
| `tests/test_project_rpc.py` | 41 | The `project.*` methods through `dispatch`, including refusals and their codes. |
| `tests/test_service.py` | 32 | Protocol framing, error mapping, and a real stdio subprocess proving stdout stays protocol-only. |
| `tests/test_progress.py` | 17 | Observability: progress does not disturb the request/response framing, several notifications arrive before the response, each is tagged with its run and request, an observed run produces identical science, the ticker throttles, and a refusal names the stage it happened in. |
| `tests/test_multi_fasta_scope.py` | 10 | What All-files ACTUALLY does with overlapping files today (section 15). Documentation, not endorsement. |

`tests/parity_driver.py` (233 lines) is a helper, not a test file.

**The formerly failing test is fixed.** Earlier revisions of this map recorded
`tests/test_core.py::test_find_dmc_information` failing because its assertion
predated the min/max search rework. The suite is green: 288 passed.

Still not covered by any Python test: `punishments.py`, `sequence_subsets.py`,
`excel.py` and `reports.py` beyond what parity exercises, and both Tkinter
viewers.

### Desktop: 238 tests, all passing

`npm test` in `desktop/` (Vitest, jsdom). `vitest.setup.ts` stubs the two
`Range` geometry methods CodeMirror calls, because jsdom has no layout engine.

| File | Tests | Covers |
|---|---|---|
| `src/app/state/focalEditing.test.ts` | 55 | The reducer: working copies, dirty derivation, undo/redo, lock refusals, scope reconciliation. |
| `src/components/focal/focalMatching.test.ts` | 29 | The frontend matcher modes. |
| `src/app/state/sourceStatus.test.ts` | 18 | The six source states, their tone, and what each offers. |
| `src/components/focal/FocalSetLibrary.test.tsx` | 18 | The library through the real workspace, plus row-control structure and the right-hand column's order and single rail. |
| `src/app/projectScreen.test.tsx` | 17 | One screen before and after Create, candidate vetting, the per-source lock including the PENDING lock and its carry-over through Create. |
| `src/app/projectFlow.test.tsx` | 17 | Launcher to run: no writes while editing, Run refused until saved, and the 700ms typing-burst wiring. |
| `src/app/focalScopes.test.tsx` | 17 | The three scopes kept distinct: run, `+` search, presence comparison. |
| `src/app/state/runGate.test.ts` | 17 | Every reason a run may be refused. |
| `src/components/focal/focalText.test.ts` | 15 | Tokenising and serialising the `;`-separated document. |
| `src/backend/backendContract.test.ts` | 12 | The shapes that cross the process boundary. |
| `src/app/state/typingBurst.test.tsx` | 9 | Fake-timer proof of undo granularity: rapid typing is one step, a pause starts another, one Undo takes back only the most recent burst, and the timer never rewrites the document. |
| `src/app/uiScale.test.ts` | 6 | 1 at 1920, 0.8 at 1536, never magnifying, and a floor. |
| `src/app/runProgress.test.tsx` | 7 | The live status line, run-token correlation (stale and post-run notifications ignored), the clock that advances without progress, and a dead backend clearing the running state. |

### Sample / test FASTA files

All under `input/`, all gitignored via `.gitignore` line 221 (`input/*`). Dimensions measured directly:

| File | Records | Alignment columns | Notes |
|---|---|---|---|
| `Leptacis_allSequences-BOLD-09March2026_aln.fasta` | **2354** (2354 `>` lines, 2354 unique headers) | **736**, uniform | Lowercase on disk; `parse_fasta` uppercases. Alphabet: `t` 614226, `a` 485835, `-` 267885, `c` 192977, `g` 171108, `n` 436, `w` 30, `r` 17, `y` 16, `m` 6, `k` 5, `s` 3. No `?`. Headers look like `AACTA5253-20|AU|Leptacis`. **364** headers contain `Leptacis_tipulae` (the focal set used by the legacy script), leaving 1990 non-focal. |
| `Leptacis_trial_alignment.fasta` | **14** | **not uniform — 13 rows of 22, 1 row of 23** | `10_Leptacis_phantasmatica` is 23 characters. **`validate_aligned_fasta` rejects this file** (`fasta_io.py:38-39`), so it cannot be run through the pipeline as-is. Uppercase; contains `N`, `R`, `B`, `H`. Headers: `1..10_Leptacis_phantasmatica` and `Leptacis_1..4`. |
| `Leptacis_trial_alignment_output.fasta` | **14** | **22**, uniform | Lowercase, same 14 headers, gaps redistributed — it looks like a corrected/realigned version of the previous file. This one passes validation. (UNVERIFIED: whether "_output" means "produced by this tool", "hand-fixed input", or something else. Nothing in the code writes a `.fasta` output — the only FASTA writer in the repo is `consensus_dmc_pipeline.write_fasta` at `consensus_dmc_pipeline.py:81`, which is legacy and writes different filenames. You would need to tell me.) |

Measured column statistics on the real file with the focal string `Leptacis_tipulae` and default options:
- 243 of 736 sites rejected by the whole-alignment ambiguity veto (`core.py:234`).
- 226 sites survive to become focal-consensus sites.
- 88 focal columns have empty fraction ≥ 2/3 and would be classified PRL or INS.

`output/` exists and is empty; its contents are gitignored (`.gitignore` line 222).

---

## 10. RISK REGISTER

The five things most likely to break during a refactor.

### 1. The stop-at-first-productive-length semantics, and `min_combination_length`'s counter-intuitive role
`core.py:358-382`. The loop starts at `start_combination_length` (not `min`), and `min` only gates the early `break`. It is very easy to "clean this up" into a loop that starts at `min` — which would silently change which combinations are found, empty `dmc.pairs` in different circumstances, change `dmc.unique`, and therefore change the 5-site results and the whole Excel workbook. A test asserting `pairs_tested == 6` was invalidated by exactly this behaviour changing once already; `tests/test_core.py:50` now pins BOTH sides of it, the early stop and the case that does enumerate pairs.
**Accidental-but-load-bearing:** `stop_reason` is initialised to `"reached_maximum_length"` at `core.py:346` and only overwritten on the early break. The GUI's entire continue-prompt (`gui.py:175`) hangs off that initial value. Any change to how `stop_reason` is assigned changes the interactive workflow, not just a report string.

### 2. The empty-state-never-matches rule
`core.py:47-48`: `state_matches_dmc_base` returns `False` for any state with no possibilities. A non-focal sequence that has a **gap** at a candidate site is therefore counted as not matching, which makes the site look diagnostic. With 267,885 gap characters across the real alignment (roughly 15% of all cells), this is not an edge case — it is a primary driver of the results. It looks like a null-guard and reads as one, so a refactor that "handles gaps properly" (skipping gapped sequences, or treating gaps as wildcards) will change the DMC set on real data without any test noticing.
Same shape, second instance: `core.py:234` vetoes a site if **any** sequence in the file — including the 1990 non-focal ones — carries an ambiguity code there. That single line discards 243 of 736 sites. Anyone who assumes the check is focal-only (as the neighbouring `consensus_base_for_column` call is) will roughly double the candidate set.

### 3. The `(L+1)//2` scan split in the punishment and consensus code
`punishments.py:490-496` and, duplicated, `consensus.py:150-152`. The two halves are disjoint and cover every column once, so the split is *not* about double-scanning — it exists solely so each end of the alignment can detect its own leading gap-dominated run and thereby distinguish PRL (terminal prolongation) from INS (internal insertion). The arithmetic midpoint has no biological meaning, and if one half is entirely gap-dominated, columns at the midpoint get classified PRL rather than INS. Refactoring this into a single forward pass, or into a "find terminal runs from both ends" formulation, will reclassify columns near the midpoint — which changes punishment scores **and** which columns the trimmed consensus drops.
Compounding: the logic is written twice, independently. They agree today. A fix applied to one and not the other is invisible until the consensus report and the punishment workbook disagree.

### 4. The 5-site exhaustive search is an unguarded combinatorial cliff
`core.py:441-488` enumerates `C(len(sites), 5)` and calls `compute_metrics` for each, which iterates every sequence in the file (`core.py:170-176`). With 2354 sequences: 20 unique sites ≈ 36M state comparisons, 30 sites ≈ 340M, 40 sites ≈ 1.5 billion. There is no cap, no sampling, no progress reporting, and it runs **again on every iteration of the GUI continue-loop**. It runs synchronously on the Tk main thread (§12), so the window stops responding for the duration.
This is currently masked by the aggressive candidate filtering (§4.4) keeping `dmc.unique` small. **Any refactor that loosens the site filters — including "fixing" either behaviour in risk 2 — makes this search explode.** (UNVERIFIED: the actual size of `dmc.unique` on the real dataset with default options; I did not run the full search. Running `find_dmc_information` on `input/Leptacis_allSequences-BOLD-09March2026_aln.fasta` with `target_string="Leptacis_tipulae"` and timing it would settle both the size and the runtime.)

### 5. Score aggregation aliasing in the punishment system
Two independent traps in the same area:
- **EW is counted inside POLY_PS/BAL_PS as well as in its own column.** `punishments.py:251-265` and `:291-305` increment `empty_weight_scores` *and* call `add_score_event`, which adds the same value to `total_scores` and to the category bucket. `excel.py:197` compensates by subtracting EW when computing `PS/bp`, but the category columns are not compensated. Any "cleanup" that stops the double-add will silently change `PS/bp` (by removing the compensation's counterpart), `POLY_PS`, and `BAL_PS` — three columns of scientific output — while every visible total still looks plausible.
- **The BD early return forgives the entire column, not just the ambiguous sequences.** `punishments.py:285`. It reads like a per-sequence exemption and is in fact a per-column one. A refactor that "only exempts the ambiguous state" would start charging base-disagreement scores on every BD column, changing every sequence's total.
**Accidental-but-load-bearing nearby:** `INS` uses `score = float(b_any)` (`punishments.py:359`) — a raw count with no weight constant, unlike PRL, POLY and BAL. It is plausibly deliberate (INS should dominate), but it is the one category with no tunable, and it is on a completely different scale from the others. Anyone normalising the weights will change the ranking that the whole punishment workbook is sorted by.

### Also worth knowing (below the top five)

- `ref_id = focal_headers[0]` (`pipeline.py:118`) — the reference is whichever focal sequence appears first in the file. Reordering the FASTA changes the reference, which changes the Excel `Full` sheet's pinned row and the `compute_metrics` exclusion. Nothing records which sequence was chosen except the text report.
- The Excel reference row is scored against the **consensus** states, not its own bases (`excel.py:47-50` vs `:58`), so it can show <100% self-similarity. Anyone "fixing" that changes published numbers.
- `parse_fasta` silently collapses duplicate headers (§4.1). Currently harmless on the real file; a different input makes sequences vanish without a message.
- `group_sequence_subsets` is order-dependent greedy substring containment (`sequence_subsets.py:44-80`); changing the sort key changes group membership.
- `next_available_filename` (`utils.py:4-27`) makes output non-deterministic with respect to the directory's prior contents — every test or comparison that re-runs into a dirty directory gets differently-named files.
- `models.py` dataclasses are `frozen=True` but hold mutable dicts/lists; `gui.py:205` passes `result.dmc.diagnostic_combinations` (a live list from the previous result) straight back into the next call, where `core.py:315` iterates it. Nothing mutates it today.

---

## 11. OPEN QUESTIONS FOR ME

1. **`min_combination_length`** — is the current behaviour (searches from length 1 regardless; `min` only permits early stopping) what you intended, or did you mean "do not bother testing combinations shorter than N"? These give different results whenever a shorter DMC exists. This is the single most consequential ambiguity in the map.

2. ~~`tests/test_core.py::test_find_dmc_information` asserts `pairs_tested == 6`, which the early-stop behaviour invalidated.~~ **Answered and fixed:** the test was stale. `tests/test_core.py:50` now asserts the early stop explicitly and a second case covers the pair enumeration.

3. **INS scoring** — `score = float(b_any)` with no weight constant (`punishments.py:359`), while PRL/POLY/BAL all have one. Deliberate (insertions should dominate), or an unfinished spot where an `INSERTION_WEIGHT` was meant to go?

4. **EW double-counting** — `EW` is included in `POLY_PS`/`BAL_PS` and in `PS`, and subtracted again only for `PS/bp` (§4.6, §10.5). Is `POLY_PS` supposed to include the empty-site charge, or should it be base-disagreement only?

5. **Benefit-of-doubt scope** — when the BD rule fires, the whole column is exempted from base-disagreement scoring, not just the ambiguous sequences (`punishments.py:285`). Intended?

6. **The `(L+1)//2` midpoint** — is the arithmetic midpoint the right boundary for the PRL/INS scan, or should the two passes meet where the terminal runs actually end? And should `consensus.classify_prl_ins_sites` always track `punishments.find_focal_punishments`, or are they allowed to diverge?

7. **The whole-alignment ambiguity veto** (`core.py:234`) — should a single ambiguity code in a *non-focal* sequence disqualify a site for everyone? It currently removes 243 of 736 sites on your real data. If it should be focal-only, the candidate set roughly doubles and the 5-site search cost explodes with it.

8. **The Excel reference row** — is it intended that the pinned reference row displays its own bases but is scored against the consensus (`excel.py:47-50` vs `:58`), so it can read below 100%?

9. **`input/Leptacis_trial_alignment.fasta` has one 23-character row among thirteen 22-character rows**, so `validate_aligned_fasta` rejects it. Is that file a deliberate "should be rejected" fixture, or just broken? And what is `Leptacis_trial_alignment_output.fasta` — tool output, hand-corrected input, or output from something else? Nothing in the live code writes `.fasta`.

10. **`viewer_old.py`, `pipeline_verbose.py`, `consensus_dmc_pipeline.py`** — I have classified all three LEGACY-SUSPECTED on the evidence in §1. Do you want any of them kept as reference during the port, or is any of them still something you actually run?

11. **Zoom** — the viewer has none (font size hardcoded at `viewer.py:47`). Is zoom expected in the Canvas version, i.e. new functionality rather than a port?

12. **Viewer click behaviour** — `viewer.py:726` prints cell details to stdout. Was that debugging left in, or is it standing in for a status-bar/inspector feature that the Canvas version should have?

13. **Should the viewer know about analysis results?** Today it is a pure sink taking only `dict[header, sequence]`. Highlighting DMC sites, PRL/INS regions, or punishment scores in the alignment would be the natural reason to port it — but nothing in the current interface supports that.

14. **Failure modes** — should the extracted library keep raising bare `ValueError` with prose messages (which the tests match on by text), or do you want typed exceptions? This affects `gui.py:208`, `gui.py:242` and every test in `tests/test_fasta_io.py`.

15. **Output naming** — is `next_available_filename`'s "never overwrite, add `(2)`" behaviour a requirement, or an artefact? It makes runs non-reproducible by filename and produces a numbered set per iteration of the continue-loop.

16. **`PunishmentResult.events`** is fully populated and never consumed by anything (§4.6). Was a per-event report planned — the `reports.py:126-131` placeholder suggests so — or should it go?

---

## 12. STATE AND CONCURRENCY

### Threads, processes, async
**None.** No `threading`, no `multiprocessing`, no `asyncio`, no `concurrent.futures`, no subprocess anywhere in the project (grep-verified across all live and legacy files). Everything is single-threaded and synchronous.

The direct consequence: **all analysis runs on the Tk main thread.** `gui.py:160` and `gui.py:229` call into the pipeline from inside a button callback, so `mainloop` is blocked for the entire run — the window is unresponsive and unpaintable, with no progress indication and no way to cancel. On the real dataset with a large `dmc.unique` this can be a very long freeze (§10.4).

### Callbacks and scheduled work
- **Tk event callbacks** — `gui.py`: `browse_fasta`, `browse_output`, `use_fasta_location`, `view_fasta`, `run`, `run_punishments`. `viewer.py`: `on_motion`, `on_mousewheel`, `on_linux_scroll_up/down`, `on_configure`, `on_press`, `on_drag`, `on_release`, `xview`, `yview`, `clear_hover` lambdas (bound at `viewer.py:689-748`).
- **`after_idle`** — the only scheduled work in the codebase: `viewer.py:217` and `viewer_old.py:136`, both inside `request_redraw`. Correctness depends on the `redraw_pending` flag being set before scheduling and cleared inside the callback (`viewer.py:210-215`); since everything runs on one thread this is a coalescing latch, not a lock.
- **Modal dialogs re-enter the event loop.** `messagebox.askyesno` (`gui.py:179`) and `simpledialog.askinteger` (`gui.py:192`) run a nested event loop inside the `while True:` search loop, so the user can interact with the main window mid-search — including clicking **Run DMC Analysis** or **Run Punishment Analysis** again. Nothing disables the buttons or guards re-entry. (UNVERIFIED whether Tk's grab actually blocks this on Windows for these specific dialogs — worth testing by clicking the main window's buttons while the "Continue DMC search?" dialog is open.)

### Module-level globals and singletons
No mutable module-level state in any live module. The module-level names are all constants in `constants.py`: `COLORS`, `IUPAC`, `IUPAC_FROM_BASES`, `STRICT_BASES`, the five output basenames, the six thresholds/weights. `COLORS`, `IUPAC` and `IUPAC_FROM_BASES` are **mutable dict objects** shared by every importer (`core`, `punishments`, `consensus`, `excel`, `viewer`, `viewer_old`), but nothing writes to them. There are no singletons and no module-level caches.

`__init__.py` defines only `__version__`.

### Shared mutable state

**In the GUI (`gui.py:10-243`), all closure-scoped inside `launch_gui`:**
| State | Written by | Read by |
|---|---|---|
| `fasta_var`, `target_var`, `output_dir_var` | `browse_fasta` (`:36`), `browse_output` (`:45`), `use_fasta_location` (`:54`), and the Entry widgets | `view_fasta` (`:57`), `run` (`:119-121`), `run_punishments` (`:212-214`) |
| `include_ambiguous_dmc_bd_var`, `include_gappy_consensus_dmc_sites_var` | Checkbutton widgets | `run` (`:164-165`) |
| `min_combination_length_var`, `max_combination_length_var` | Spinbox widgets | `run` (`:136-137`) |
| `start_combination_length`, `current_max_combination_length`, `initial_diagnostic_combinations`, `initial_combinations_tested_by_length` | `run`'s loop body (`:154-157`, `:203-206`) | the next `run_pipeline_core` call (`:166-170`) |

The last row is genuinely order-dependent loop-carried state: `gui.py:205` assigns `result.dmc.diagnostic_combinations` — **a live list owned by the previous `DMCResult`** — into the variable that `core.py:315` will iterate on the next pass. Nothing mutates that list, so it is safe today, but it is a cross-object aliasing edge waiting to be tripped.

**In the viewer (`viewer.py`), closure-scoped:**
| State | Written by | Read by | Order-dependent? |
|---|---|---|---|
| `x_offset`, `y_offset` | `xview` (`:554`,`:561-563`), `yview` (`:577`,`:584-586`), `on_mousewheel` (`:657-659`), `on_linux_scroll_up/down` (`:668-670`,`:679-681`), `clamp_offsets` (`:155-156`) | `visible_range` (`:180-187`), `update_scrollbars` (`:165-173`), `cell_from_motion_event` (`:594-595`), `cell_from_event` (`:705-706`) | Yes — `clamp_offsets` must run before `visible_range`, which `redraw` guarantees at `:504`/`:507`. Hit-testing reads the offsets *outside* that sequence, so between an event and the next idle repaint the offsets can be un-clamped. |
| `hovered_row`, `hovered_col` | `on_motion` (`:618-619`), `clear_hover` (`:628-629`) | `draw_names` (`:328`), `render_sequence_bitmap` (`:382` — computed and discarded, `:444`), `draw_sequence_letters` (`:499`) | No |
| `redraw_pending` | `request_redraw` (`:210`), the inner `_redraw` (`:214`) | `request_redraw` (`:207`) | Yes — the latch must be cleared before `redraw()` runs, which `:214-215` does in that order. |
| `selection_anchor`, `selection_active` | `on_press` (`:719-720`), `on_drag` (`:738`) | `selected_rectangle` (`:348-359`) → `render_sequence_bitmap` (`:407`) | No |
| `seq_canvas._viewport_image` | `viewer.py:530` (and initialised `:119`) | Nothing reads it — it exists purely to hold a reference so Tk does not garbage-collect the `PhotoImage` mid-frame | Yes, implicitly: the assignment must happen before/with `create_image` (`:531`) |

### Caches
- **`color_cache`** (`viewer.py:121`) — per-`open_fasta_viewer`-call dict, base→RGB, written and read only by `parse_color` (`viewer.py:123-144`). Unbounded but keyed by single characters, so effectively ≤ ~30 entries. Not shared between viewer windows.
- **`kept_index_to_trimmed_position`** (`consensus.py:244-247`) — a derived lookup built inside the report writer, local, read-only after construction.
- **On-disk:** `__pycache__/` and `.pytest_cache/` only. No application-level disk cache.
- **No `functools.lru_cache` or `@cache` anywhere.** Every column extraction, focal/non-focal split, and FASTA parse is recomputed from scratch at each call site (§4.1, §6.12).

### Temporary files
**None.** No `tempfile` import, no scratch files. The only files written are the four named outputs, all via `next_available_filename` into the user-chosen directory.

### Values mutated in one module and read in another

This is the short list, and it is short because the science modules are functional in style:

1. **`punishments.py` → `excel.py`.** `find_focal_punishments` builds nine `defaultdict`s and an `events` list (`punishments.py:408-419`), passes them by reference through `process_polymorphism_or_balancing_column` and `process_gap_dominated_column` (which mutate them in place — `punishments.py:202-203`, `:251`, `:271`, `:371`), converts them to plain dicts at `punishments.py:498-508`, and hands them to `excel.write_punishment_excel_report`, which reads them at `excel.py:193-214`. The mutation is confined to one module and the boundary object is a snapshot, so this is safe — but the *accumulation order* determines floating-point rounding in the totals. Reordering the scan (risk 3) changes the last bits of every score.
2. **`core.py` → `gui.py` → `core.py`.** `DMCResult.diagnostic_combinations` and `.combinations_tested_by_length` cross out of the core into the GUI (`gui.py:205-206`) and back in as `initial_*` arguments (`core.py:315`, `:334`). `combinations_tested_by_length` is `.update()`d into a fresh dict at `core.py:334-335` (safe copy); `diagnostic_combinations` is iterated but not mutated (`core.py:315-328`).
3. **`core.py` → `pipeline.py` → `excel.py`/`reports.py`.** `dmc.states` and `dmc.unique` are read by three separate writers (`pipeline.py:176`, `:183`, `:187`; `reports.py:29`, `:92`, `:97`). Read-only.
4. **`constants.COLORS`/`IUPAC`** — shared mutable dicts read by six modules, written by none.

Nothing is mutated in one module and read in another across a genuine ownership boundary.

---

## 13. MISSING CHARACTERIZATION TESTS

Priority order. These pin current behaviour before it moves; none of them are written here.

**1. Golden-file test of the full DMC text report on a fixed small alignment.**
Pin: the entire byte content of `DMCs_output.txt` for a checked-in ~15×40 alignment with a known focal string, run with default options.
Input needed: a small aligned FASTA committed to the repo containing, deliberately, at least one single-site DMC, at least one gap-only-difference site, at least one ambiguity code in a non-focal sequence, and one globally conserved site. `input/Leptacis_trial_alignment_output.fasta` (14 × 22, has `n`/`r`/`b`/`h` and gaps) is close to ideal but is gitignored and lowercase — it should be copied into `tests/data/` and its focal string decided.
Why high-risk untested: this one file transitively pins the candidate filter, the search loop, the stop reason, the counts, and the 1-based formatting — the four places risks 1 and 2 live. Right now **no test asserts on any output file's content at all**.

**2. Parametrised table test of `find_dmc_information` search semantics.**
Pin: `single`, `pairs`, `unique`, `stop_reason`, `stopped_at_length`, and `combinations_tested_by_length` across the grid of `(min, max, start)` combinations — specifically `(1,2,1)`, `(2,3,1)`, `(3,5,1)`, and a resume case with `start=3` plus `initial_diagnostic_combinations`.
Input needed: two tiny hand-built dicts — one where a 1-site DMC exists (so the early break fires) and one where the shortest DMC is 3 sites (so it does not).
Why high-risk untested: this is risk 1. `tests/test_core.py:50` now covers the early stop, and `tests/test_parity.py::test_parity_non_default_min_and_max` covers the min/max interaction against the pre-integration baseline. Include the "min=3 still searches lengths 1 and 2" case explicitly — that is the behaviour most likely to be silently "corrected".

**3. Gap and ambiguity truth table for the two option flags.**
Pin: for a purpose-built alignment, the exact `candidate_count`, `skipped_non_acgt`, `globally_conserved_removed`, `ambiguous_bd_sites_included` and `gappy_consensus_sites_included` under all four combinations of `include_ambiguous_dmc_bd` × `include_gappy_consensus_dmc_sites`.
Input needed: an alignment with one column per case — focal-invariant-with-a-gap; focal-invariant with an `N` in a *non-focal* sequence only; focal column with an `R` compatible with the focal base; a column ≥2/3 gaps in the focal set; a globally conserved column.
Why high-risk untested: risk 2, both halves. The whole-alignment ambiguity veto and the empty-never-matches rule are each one line, both read like defensive guards, and together they decide 243 of 736 sites on the real data. Neither flag has a single test today.

**4. Punishment score snapshot, per category, on a hand-built focal set.**
Pin: `total_scores`, `polymorphism_scores`, `balancing_scores`, `prolongation_scores`, `insertion_scores`, `empty_weight_scores`, `bd_counts`, `prl_counts`, `ins_counts` — all nine dicts, to a fixed float tolerance — plus `len(events)`.
Input needed: a ~6×20 focal-only alignment engineered so that each column exercises exactly one path: a clean POLY column, a POLY column that triggers BD (one strict base + ≤20% compatible ambiguity), a POLY column with one gap (to pin the `b_any/b_empty` formula), a BAL column, a terminal gap-dominated run at each end (PRL), and an internal gap-dominated column (INS).
Why high-risk untested: `punishments.py` is 508 lines with **zero** tests and produces a primary scientific output. It contains risk 5 (the EW aliasing and the BD early return) and half of risk 3 (the scan split). The per-category breakdown is what makes a regression legible — a total-only assertion would let EW/POLY reshuffling pass.

**5. PRL/INS classification equivalence between `punishments.py` and `consensus.py`.**
Pin: for the same focal alignment, `classify_prl_ins_sites` returns exactly the site sets that `find_focal_punishments` charges as PRL and INS (derivable from `PunishmentEvent.category`).
Input needed: the same alignment as test 4, plus one case where an entire half is gap-dominated so a midpoint column is classified PRL.
Why high-risk untested: the logic is duplicated in two modules (risk 3) with no shared code and no test tying them together. They agree today; a one-sided fix would be invisible until the consensus report and the punishment workbook silently disagreed about which columns are real.

**6. Consensus rule truth table.**
Pin: `assign_consensus_state` output for: all-gap column → `N`; single base → that base; 60/40 split → the IUPAC code (not the majority base, since 0.6−0.4 = 0.2 ≯ 0.5); 80/20 split → the majority base; a column of ambiguity codes with one strict base → that strict base (pinning that ambiguity codes are *ignored*, not decomposed); a base set with no IUPAC code → `N`. Then pin `untrimmed_sequence`, `trimmed_sequence` and `kept_indices` for a whole small alignment.
Input needed: hand-built columns; no file required for the first half.
Why high-risk untested: the `> 0.5` dominance margin (`consensus.py:75`) and the ignore-ambiguity behaviour are both unobvious, both undocumented in the report text the tool writes, and both easy to "correct" toward a conventional consensus rule.

**7. Excel workbook structural snapshot.**
Pin: for `comparison_output.xlsx` — sheet names present, header row, row order, the reference row's position and its computed `% similarity` (deliberately capturing the consensus-vs-own-bases scoring), and the fill colour of a couple of known cells. For `punishment_output.xlsx` — column order and the sort order of rows. Read back with `openpyxl` rather than comparing bytes.
Input needed: the same small alignment as test 1.
Why high-risk untested: `excel.py` has no tests, and it contains the reference-row scoring divergence (§4.9) and the `PS/bp` EW compensation (risk 5) — two places where a "fix" changes numbers a scientist will read as authoritative.

**8. `next_available_filename` and the non-clobbering contract.**
Pin: first call returns the bare name; with the bare name present returns `(2)`; with `(2)` present returns `(3)`; with `(7)` present returns `(8)`; a non-numeric `name(x).txt` sibling is ignored.
Input needed: `tmp_path` only.
Why high-risk untested: only indirectly covered by `test_pipeline`'s `(2)` case. Every output path in the tool flows through it, and the GUI continue-loop depends on it producing a fresh set per iteration.

**9. `group_sequence_subsets` ordering and containment.**
Pin: group count, representative choice, and member ordering for a set containing an exact duplicate pair, a strict substring, an unrelated sequence, an empty-after-ungapping sequence, and a length tie broken by header name.
Input needed: a small dict of ungapped-with-gaps strings.
Why high-risk untested: the grouping is order-dependent greedy containment (§4.8) that changed semantics from the legacy `cluster_sequences`, has no tests, and feeds a delivered workbook.

**10. A slow-marked end-to-end run on the real alignment.**
Pin: not exact output, but invariants — the run completes; `len(dmc.unique)`; `total_combinations_tested`; `stop_reason`; the number of PRL/INS columns (88 at present); the number of sites skipped by the ambiguity veto (243 at present); and a wall-clock ceiling.
Input needed: `input/Leptacis_allSequences-BOLD-09March2026_aln.fasta` with `target_string="Leptacis_tipulae"` — currently gitignored, so it needs either un-ignoring, a Git LFS entry, or a fixture that skips when absent.
Why high-risk untested: this is the only thing that would catch risk 4 — the 5-site combinatorial cliff. Nothing today measures how long a real run takes, so a filter change that turns 30 seconds into 30 minutes would ship unnoticed.

## 14. RUN OBSERVABILITY

Added because a run on another machine could sit in "running" indefinitely with
no way to tell a hang from an expensive search. Nothing about the science
changed; what changed is that it now says what it is doing.

### Two channels, two audiences

| Channel | Who reads it | What it carries |
|---|---|---|
| **stderr**, one JSON object per line behind `[mdx]` | developers | run id, elapsed ms, stage, durations, counts, environment, tracebacks |
| **the protocol channel**, `{"id", "type":"progress", "progress":{...}}` | the UI | stage name plus numbers, correlated with the pending request |

stdout remains protocol-only. A progress line carries `type` and no `ok`, so it
cannot be confused with a response, and it leaves its request PENDING — which
is the point. `tests/test_progress.py` pins that framing.

**A parent that pipes stderr MUST drain it.** The service writes diagnostics for
every request; an undrained pipe fills at 64 KB and the child then blocks
mid-write and never answers, which looks exactly like a hung analysis. The
Electron bridge reads it line by line, and the two test clients run a drain
thread for the same reason.

### What is instrumented

`project.runMolecularDiagnosis` end to end:

request received and its parameters; the environment (interpreter, version,
platform, PID, package versions, cwd, repo root); per-file source status with
path/size/mtime/permissions; the read-and-hash of each alignment; record counts
and alignment dimensions; focal-entry validation; the pooled scope's size and
its problems; the focal/non-focal split; the focal consensus; the DMC search,
with `C(n,k)` announced per size and periodic counts inside it; the number of
candidate and unique sites; the five-site search, with its `C(n,5)` total and
periodic counts; output-path allocation; each of the three writers; completion;
and any refusal or exception, with the stage it happened in.

### Cost of being observable

`ProgressTicker.advance()` is a counter increment and a comparison; it consults
the clock every few thousand iterations and emits at most once a second. The
default observer everywhere is `NULL_OBSERVER`, whose methods are empty, so a
direct `run_pipeline_core` call (the Tkinter app, every existing test) is
unchanged. `tests/test_progress.py::test_an_observed_run_produces_identical_science`
compares an observed run against an unobserved one field by field.

### Where the time actually goes

Measured on `input/Leptacis_allSequences-BOLD-09March2026_aln.fasta`
(2,354 records, 736 columns, focal `Leptacis_tipulae` = 364 sequences,
1,990 contrast sequences):

| Stage | Cost |
|---|---|
| parse + validate | 0.02 s |
| DMC search, max size 2 (the default) | 0.94 s — 118 candidate sites, **0 unique sites**, stops with `reached_maximum_length` |
| DMC search, max size 3 | 4.94 s — 273,937 combinations tested, 8 unique sites |
| five-site search | trivial here (`C(8,5)` = 56) |

The two cliffs, both `C(n, k)` with a per-combination cost proportional to the
contrast set:

* **the DMC combination search.** 18 µs per combination on this data, so
  `C(118,4)` ≈ 7.7 M is ~2 minutes, `C(118,5)` ≈ 175 M is ~53 minutes, and
  `C(118,6)` ≈ 3.3 G is ~16 hours. A user who raises the maximum candidate
  size, or who keeps accepting the continuation prompt, walks straight into it.
* **the five-site optimisation.** `C(n,5)` where n is the number of unique
  diagnostic sites, each iteration comparing the reference against every
  non-focal sequence. 30 unique sites is 142,506 combinations — 14 s against 45
  sequences, and roughly 40x that against 1,990.

Both now announce their totals before they start and report progress while they
run, so the UI can distinguish "working" from "stuck". Neither algorithm was
changed.

### Collecting comparable logs from two machines

1. `MOLECULAR_TOOL_LOG_DIR=/some/dir npm start` in `desktop/` — every backend
   line is appended to `backend-YYYY-MM-DD.log` there (default:
   `<userData>/logs`).
2. Or drive the service directly, which needs no Electron:

   ```
   echo '{"id":"1","method":"diagnostics.environment"}' | python -m molecular_diagnosis.service
   ```

   The `service.start` line and the response both carry the full environment.
3. Compare `service.start`, then `run.start`, then the `stage.end:*` durations.
   The first line that differs is the answer.

---

## 15. MULTI-FASTA SCOPE: WHAT "ALL FILES" DOES TODAY

Findings, not a design. The overlapping-FASTA rules are an open question; this
section records what the implementation currently does so the decision can be
made against the code. `tests/test_multi_fasta_scope.py` holds the same facts
as executable tests.

### The reported failure

File A is a full alignment; File B is the focal-only subset of it. Under
All files the run is **refused** with `DUPLICATE_HEADER_ACROSS_FILES`:

    The header 'Leptacis_tipulae_1' appears in more than one selected file.
    detail: full.fasta and subset.fasta; (+N more)

`build_scope` (`project/service.py`) records **one problem per shared header**,
and `run_molecular_diagnosis` refuses on the first. Since every header of a
subset is by definition also in the superset, the refusal is guaranteed for
this shape of input. It is not a crash and not a silent wrong answer — it is
the deliberate refusal added when multi-file scope was built, meeting a case it
was not designed around.

### How several files are actually combined

* Files are pooled **in memory** into ONE `header -> sequence` dictionary, in
  selection order. No temporary combined FASTA is written and nothing is
  re-parsed; the per-file alignment cache is merged directly.
* There are **no record ids**. The scientific core is keyed by header, so a
  header is the only identity a sequence has.
* Each file is validated **independently before pooling**: unreadable, empty,
  unaligned, or containing repeated headers within itself
  (`DUPLICATE_HEADER_IN_FILE`) removes that file from the run entirely, and the
  other files still pool.
* After pooling, the files must share ONE alignment length, or
  `INCOMPATIBLE_ALIGNMENT_LENGTHS`.

### Duplicate headers, precisely

| Case | Current behaviour |
|---|---|
| Same header twice **in one file** | The file is refused (`DUPLICATE_HEADER_IN_FILE`). The indexer keeps the LAST record for a repeated header, so the count of records and the count of distinct headers differ, and analysing it would silently discard sequences. |
| Same header in **two selected files**, identical sequence | Refused (`DUPLICATE_HEADER_ACROSS_FILES`). |
| Same header in **two selected files**, different sequence | Refused with the SAME code. |

**No sequence comparison happens.** "The same record twice" and "two different
records under one name" are not currently distinguished — the check is purely
on the header. That is the distinction the next pass has to decide about.

Underneath the refusal, the pooling itself is **first-wins**: the first selected
file's record is kept and the later one is skipped (the skip is what becomes the
problem). Selecting the files in the other order would keep the other record.
That order-dependence is exactly why the run is refused rather than allowed
through.

### What is NOT broken

The established All-files UI semantics still hold and are covered by tests:
presence is green when a header is in ANY selected file, orange cannot occur
under All files (there is no "elsewhere" when everything is selected), `+`
searches the union of the scope, and a multi-file run over files with **no**
shared headers works normally.

### The one narrow fix made here

The refusal's `detail` used to be every problem's detail joined together — with
a real alignment, hundreds of near-identical sentences in the message the UI
shows. It is now capped at eight with a `(+N more)` count
(`summarise_problems`). No policy changed.

---
