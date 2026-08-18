/**
 * contract.ts — types-only frontend/backend contract.
 *
 * THIS FILE CONTAINS NO IMPLEMENTATION. It exists so that a Python adapter can
 * later be built against a stable, explicit shape, and so that the places where
 * the current Python code cannot yet satisfy the frontend are written down
 * rather than discovered during integration.
 *
 * Every place the existing Python package cannot currently supply what the
 * frontend needs is marked:
 *
 *     // BACKEND GAP: ...
 *
 * Source of truth for the current Python behaviour is REPO_MAP.md at the repo
 * root; section references below (§4, §7, ...) point into it.
 *
 * Naming convention: this file uses frontend camelCase names. Where a field
 * maps onto a differently-named Python parameter, the Python name is given in
 * the doc comment. See desktop/UI_NOTES.md for the full mapping table.
 */

/* ------------------------------------------------------------------ */
/* Shared primitives                                                    */
/* ------------------------------------------------------------------ */

/** Absolute filesystem path. Produced by the native picker; never parsed in the renderer. */
export type FilePath = string;

/** Opaque identifier minted by the frontend for a session-local entity. */
export type Id = string;

/** ISO-8601 timestamp. */
export type IsoTimestamp = string;

/** Zero-based alignment column index, matching Python's internal indexing (§4.9). */
export type SiteIndex = number;

/**
 * One-based alignment column, as shown to users and written into every Python
 * report (`site + 1`, REPO_MAP §4.9). Kept as a distinct alias so conversions
 * are visible at call sites.
 */
export type SiteNumber = number;

/* ------------------------------------------------------------------ */
/* Analyses                                                            */
/* ------------------------------------------------------------------ */

export type AnalysisKind =
  | 'molecularDiagnosis'
  | 'sequencePunishmentTest'
  | 'consensusSequenceGeneration';

/**
 * Which analyses a project has enabled.
 *
 * ASSUMPTION (see UI_NOTES §7): selection here drives which workspace tabs are
 * enabled. The Python side has no concept of a "project" with enabled analyses;
 * it exposes two independent entry points (`run_pipeline_core`,
 * `run_punishment_core`) that are called on demand.
 *
 * BACKEND GAP: consensus generation has no standalone entry point at all. It
 * runs as a side effect of `run_pipeline_core` (REPO_MAP §2 walkthrough A step
 * 7) and always writes `focal_consensus_output.txt`. Selecting consensus
 * without molecular diagnosis is therefore not currently executable.
 */
export type AnalysisSelection = Record<AnalysisKind, boolean>;

/**
 * BACKEND GAP: no execution ordering exists for multiple selected analyses.
 * The Tkinter UI has two separate buttons the user presses independently, and
 * each re-parses the FASTA from scratch (REPO_MAP §2, §6.12). Ordering,
 * shared parsing, and whether analyses may run concurrently are all undecided.
 */
export type AnalysisRunOrder = readonly AnalysisKind[];

/* ------------------------------------------------------------------ */
/* Project                                                             */
/* ------------------------------------------------------------------ */

/**
 * The in-session project draft. This is what the creation screen builds and
 * what the workspace reads. It is intentionally NOT the persisted shape.
 */
export interface ProjectDraft {
  readonly id: Id;
  name: string;
  fastaPath: FilePath | null;
  analyses: AnalysisSelection;

  /** Per-analysis settings, filled in progressively by the workspace screens. */
  molecularDiagnosis: MolecularDiagnosisConfig;

  /** Room for future settings without reshaping the draft. */
  punishment?: PunishmentConfig;
  consensus?: ConsensusConfig;

  readonly createdAt: IsoTimestamp;
}

/**
 * What a saved project file would need to contain.
 *
 * BACKEND GAP: no project format exists anywhere in the repository. Nothing in
 * the Python package reads or writes project metadata — it takes three loose
 * arguments (fasta path, focal string, output directory) per invocation
 * (REPO_MAP §7). The shape below is a proposal, not an implemented format.
 */
export interface PersistedProject {
  readonly formatVersion: 1;
  readonly id: Id;
  readonly name: string;
  readonly createdAt: IsoTimestamp;
  readonly updatedAt: IsoTimestamp;

  /**
   * BACKEND GAP: undecided whether the FASTA is referenced by path or copied
   * into the project. A bare path breaks when the project moves between
   * machines; copying duplicates a 1.8 MB alignment per project.
   */
  readonly fastaPath: FilePath;

  /** BACKEND GAP: no integrity check exists today. Needed to detect the FASTA changing under a saved project. */
  readonly fastaChecksum?: string;

  readonly analyses: AnalysisSelection;
  readonly focalSets: readonly FocalSet[];
  readonly molecularDiagnosis: MolecularDiagnosisConfig;
  readonly punishment?: PunishmentConfig;
  readonly consensus?: ConsensusConfig;

  /**
   * BACKEND GAP: results are currently only files on disk with
   * collision-avoiding names (`DMCs_output(2).txt`, REPO_MAP §10). There is no
   * record of which run produced which file, so a reopened project cannot
   * currently restore its results.
   */
  readonly runs?: readonly AnalysisRunRecord[];
}

/** What "Open Existing Project" will need in order to work. */
export interface OpenProjectRequest {
  readonly path: FilePath;
}

export interface OpenProjectResponse {
  readonly project: PersistedProject;
  /** Non-fatal problems, e.g. the referenced FASTA no longer exists at its recorded path. */
  readonly warnings: readonly string[];
}

/* ------------------------------------------------------------------ */
/* Focal selection                                                     */
/* ------------------------------------------------------------------ */

/**
 * A named group of focal search strings.
 *
 * BACKEND GAP: the Python API accepts EXACTLY ONE focal string. It is a single
 * `target_string` threaded through `split_focal_headers`,
 * `find_dmc_information`, `compute_metrics` and `build_sheet` (REPO_MAP §4.2),
 * and it is also used to name consensus FASTA records. The frontend models a
 * list because the design shows multiple tokens; reconciling the two is
 * explicitly out of scope for the UI phase and MUST NOT be resolved by
 * changing Python behaviour unilaterally.
 */
export interface FocalSet {
  readonly id: Id;
  /** The "FOCAL SET TITLE" field. Frontend-only today — Python has no equivalent. */
  title: string;
  /** BACKEND GAP: multiple strings; Python consumes one. */
  strings: readonly string[];
}

/**
 * How a focal string is matched against FASTA headers.
 *
 * `caseSensitiveSubstring` is the ONLY mode the Python code implements today:
 * `target_string in header`, case-sensitive, no field awareness (REPO_MAP §4.2).
 * The other modes are named so the frontend matcher can be swapped without
 * reshaping the editor — none of them exist in Python.
 */
export type FocalMatchMode =
  | 'caseSensitiveSubstring'
  // BACKEND GAP: not implemented in Python.
  | 'caseInsensitiveSubstring'
  // BACKEND GAP: not implemented in Python.
  | 'exactId'
  // BACKEND GAP: not implemented in Python. Would need header field parsing (`ID|COUNTRY|TAXON`).
  | 'fieldAware'
  // BACKEND GAP: not implemented in Python.
  | 'regex';

export interface FocalValidationRequest {
  readonly fastaPath: FilePath;
  readonly tokens: readonly string[];
  readonly mode: FocalMatchMode;
}

export interface FocalTokenValidation {
  readonly token: string;
  readonly matched: boolean;
  readonly matchCount: number;
  /**
   * Sample of matching headers for UI display. Capped by the backend — a token
   * like "Leptacis" matches 2354 headers in the real dataset (REPO_MAP §9).
   */
  readonly sampleHeaderIds: readonly string[];
  readonly truncated: boolean;
}

export interface FocalValidationResponse {
  readonly results: readonly FocalTokenValidation[];
  readonly totalHeaders: number;
  /**
   * Headers matched by at least one token — the union, i.e. the prospective
   * focal set.
   *
   * BACKEND GAP: Python has no union concept; with one string the focal set is
   * simply that string's matches. Union vs. intersection semantics for multiple
   * tokens is an open question (UI_NOTES Q4).
   */
  readonly unionMatchCount: number;
}

/* ------------------------------------------------------------------ */
/* Analysis configuration                                              */
/* ------------------------------------------------------------------ */

/**
 * Molecular Diagnosis parameters.
 *
 * Field-by-field mapping onto the Python keyword arguments of
 * `find_dmc_information` / `run_pipeline_core` (REPO_MAP §7).
 */
export interface MolecularDiagnosisConfig {
  /**
   * UI label: "Ignore gaps".
   * Python: `include_gappy_consensus_dmc_sites` (SAME polarity — not inverted).
   * true  => gap-containing focal columns may still yield a consensus base, and
   *          gap-only differences can admit a site as a candidate.
   */
  ignoreGaps: boolean;

  /**
   * UI label: "Give BoTD to ambiguous bases".
   * Python: `include_ambiguous_dmc_bd` (SAME polarity).
   */
  giveBenefitOfDoubtToAmbiguousBases: boolean;

  /**
   * UI label: "Min. candidate DNC size".
   * Python: `min_combination_length` (default 1).
   *
   * WARNING — the name is misleading and this is load-bearing. It is NOT a
   * floor on which combination lengths are searched. The Python loop always
   * starts at `start_combination_length` (default 1); `min` only gates whether
   * the search is allowed to STOP early once a productive length is found
   * (REPO_MAP §4.4, §10.1). Do not "fix" this in the frontend.
   */
  minCandidateSize: number;

  /**
   * UI label: "Max. candidate DNC size".
   * Python: `max_combination_length` (default 2). Upper bound of the search loop.
   */
  maxCandidateSize: number;

  /** The focal set this run uses. */
  focalSetId: Id | null;

  /**
   * BACKEND GAP: `output_dir` is a required Python argument with no UI in the
   * mockups. `load_inputs` raises if it is missing or does not exist
   * (REPO_MAP §7). Where outputs go is an unresolved product question.
   */
  outputDirectory?: FilePath;
}

/**
 * BACKEND GAP: `run_punishment_core` accepts NO tunables at all — only path,
 * focal string and output directory (REPO_MAP §7). The thresholds and weights
 * it uses are module-level constants in `constants.py`
 * (POLYMORPHISM_EMPTY_LIMIT, BALANCING_EMPTY_LIMIT, POLYMORPHISM_EMPTY_WEIGHT,
 * BALANCING_EMPTY_WEIGHT, PROLONGATION_WEIGHT, BD_AMBIGUOUS_FRACTION_LIMIT).
 * Exposing any of them requires a Python signature change.
 */
export interface PunishmentConfig {
  readonly focalSetId: Id | null;
}

/**
 * BACKEND GAP: consensus generation likewise has no parameters and no entry
 * point of its own. The dominance margin (0.5) and FASTA wrap width (80) are
 * hardcoded in `consensus.py` (REPO_MAP §7).
 */
export interface ConsensusConfig {
  readonly focalSetId: Id | null;
}

/* ------------------------------------------------------------------ */
/* Alignment transfer                                                  */
/* ------------------------------------------------------------------ */

/**
 * Metadata about a loaded alignment. Cheap, structured, safe to send as JSON.
 *
 * BACKEND GAP: nothing in Python returns this. `parse_fasta` produces a
 * `dict[str, str]` consumed in-process (REPO_MAP §4.1); no API exposes it.
 */
export interface AlignmentMetadata {
  readonly sequenceCount: number;
  readonly alignmentLength: number;
  readonly headers: readonly string[];

  /**
   * True when every sequence is the same length. Python's
   * `validate_aligned_fasta` raises rather than reporting, so a false value
   * here has no Python counterpart yet.
   */
  readonly isAligned: boolean;

  /** Distinct characters observed, uppercased. Drives the residue colour legend. */
  readonly alphabet: readonly string[];

  /**
   * BACKEND GAP: `parse_fasta` silently collapses duplicate headers because it
   * builds a dict (REPO_MAP §4.1). The viewer needs to know if this happened.
   */
  readonly duplicateHeaderCount: number;
}

/** How the residue matrix is encoded for transfer. */
export type AlignmentEncoding =
  /** One byte per residue, ASCII uppercase, row-major. */
  | 'ascii-uint8'
  /** One byte per residue, index into `AlignmentMatrixDescriptor.symbols`. */
  | 'symbol-index-uint8';

/**
 * Describes the binary payload WITHOUT carrying it. The bytes themselves are
 * expected to arrive out-of-band (shared memory, a file handle, or a binary
 * IPC channel) — deliberately not modelled as a JS array, because the real
 * dataset is 2354 x 736 = ~1.7M residues (REPO_MAP §9) and JSON-encoding that
 * per redraw is not viable.
 *
 * BACKEND GAP: no transport exists. Choosing one is part of the alignment
 * viewer task, not this one.
 */
export interface AlignmentMatrixDescriptor {
  readonly encoding: AlignmentEncoding;
  readonly rows: number;
  readonly columns: number;
  /** Row-major: residue (r, c) lives at byteOffset + r * columns + c. */
  readonly rowStride: number;
  readonly byteLength: number;
  /** Present for 'symbol-index-uint8'. Index -> residue character. */
  readonly symbols?: readonly string[];
}

/**
 * The shape the future alignment viewer will consume.
 * Deliberately separates cheap metadata from the bulk payload.
 */
export interface AlignmentData {
  readonly metadata: AlignmentMetadata;
  readonly matrix: AlignmentMatrixDescriptor;
  /**
   * Populated only once a transport exists. Left optional so the viewer can be
   * built against metadata first.
   */
  readonly bytes?: Uint8Array;
}

/** A viewport request, for a viewer that will virtualise like `viewer.py` does (REPO_MAP §5). */
export interface AlignmentSliceRequest {
  readonly firstRow: number;
  readonly rowCount: number;
  readonly firstColumn: number;
  readonly columnCount: number;
}

/* ------------------------------------------------------------------ */
/* Results — modelled as DATA, not file paths                          */
/* ------------------------------------------------------------------ */

/**
 * BACKEND GAP — THE BIG ONE. The current pipeline fuses computation and output
 * writing: `run_pipeline_core` reserves paths and writes three files inline and
 * unconditionally, and returns PATHS rather than data (REPO_MAP §6.6). Of
 * everything it computes, only `DMCResult` escapes the function; the
 * `FiveSiteResult` and the consensus result are written to disk and dropped.
 *
 * Serving the types below requires splitting computation from serialisation in
 * Python. That is a backend task and is explicitly NOT part of the UI work.
 */
export interface DiagnosticCombination {
  readonly sites: readonly SiteNumber[];
  /** Consensus state per site, e.g. "A". Python's `DMCResult.states`. */
  readonly states: readonly string[];
}

export interface MolecularDiagnosisResult {
  readonly kind: 'molecularDiagnosis';

  readonly combinationsByLength: Readonly<Record<number, readonly DiagnosticCombination[]>>;
  readonly uniqueSites: readonly SiteNumber[];

  readonly stoppedAtLength: number;
  /** Mirrors Python's `DMCResult.stop_reason` string union (REPO_MAP §4.4). */
  readonly stopReason:
    | 'reached_maximum_length'
    | 'found_at_or_above_minimum_length'
    | 'no_candidate_sites'
    | 'start_length_exceeds_maximum_length';

  readonly diagnostics: {
    readonly totalSequences: number;
    readonly alignmentLength: number;
    readonly focalCount: number;
    readonly nonFocalCount: number;
    readonly skippedSites: number;
    readonly candidateSites: number;
    readonly globallyConservedRemoved: number;
    readonly combinationsTestedByLength: Readonly<Record<number, number>>;
  };

  /**
   * The reference sequence Python picked.
   * NOTE: this is `focal_headers[0]` — whichever focal sequence appears first
   * in the file (REPO_MAP §10). Surfacing it matters because it silently
   * changes if the FASTA is reordered.
   */
  readonly referenceHeader: string;

  readonly fiveSite?: {
    readonly combinationsTested: number;
    readonly bestGapSites: readonly SiteNumber[] | null;
    readonly bestAverageSites: readonly SiteNumber[] | null;
  };

  /** Files the run wrote, retained during the transition away from path-only results. */
  readonly writtenFiles: readonly FilePath[];
}

export interface PunishmentScoreRow {
  readonly sequenceId: string;
  readonly total: number;
  readonly perBase: number;
  readonly emptyWeight: number;
  readonly benefitOfDoubtCount: number;
  readonly insertionCount: number;
  readonly prolongationCount: number;
  readonly polymorphismScore: number;
  readonly balancingScore: number;
  readonly insertionScore: number;
  readonly prolongationScore: number;
}

export interface PunishmentResultData {
  readonly kind: 'sequencePunishmentTest';
  readonly rows: readonly PunishmentScoreRow[];
  /**
   * BACKEND GAP: `PunishmentResult.events` is fully populated (one entry per
   * charge, hundreds of thousands on the real dataset) and then discarded —
   * nothing consumes it (REPO_MAP §4.6). If the UI is to show a per-site
   * explanation, the backend must expose and probably paginate it.
   */
  readonly events?: readonly unknown[];
  readonly writtenFiles: readonly FilePath[];
}

export interface ConsensusResultData {
  readonly kind: 'consensusSequenceGeneration';
  readonly untrimmedSequence: string;
  readonly trimmedSequence: string;
  readonly keptSites: readonly SiteNumber[];
  readonly removedProlongationSites: readonly SiteNumber[];
  readonly removedInsertionSites: readonly SiteNumber[];
  readonly writtenFiles: readonly FilePath[];
}

export type AnalysisResult =
  | MolecularDiagnosisResult
  | PunishmentResultData
  | ConsensusResultData;

export interface AnalysisRunRecord {
  readonly runId: Id;
  readonly kind: AnalysisKind;
  readonly startedAt: IsoTimestamp;
  readonly finishedAt: IsoTimestamp | null;
  readonly status: AnalysisRunStatus;
  readonly result?: AnalysisResult;
}

/* ------------------------------------------------------------------ */
/* Long-running work                                                   */
/* ------------------------------------------------------------------ */

export type AnalysisRunStatus =
  | 'queued'
  | 'running'
  | 'succeeded'
  | 'failed'
  | 'cancelled';

/**
 * BACKEND GAP: the Python core has NO progress seam. Every analysis is a single
 * synchronous call with no callbacks, no generators and no logging hooks, and
 * it runs on the Tk main thread today (REPO_MAP §12). Emitting these events
 * requires threading the work off the UI thread and adding reporting points
 * inside the search loops.
 */
export interface AnalysisProgressEvent {
  readonly runId: Id;
  readonly kind: AnalysisKind;
  /** 0..1, or null when the backend cannot estimate. */
  readonly fraction: number | null;
  /** Short human-readable stage, e.g. "Testing 3-site combinations". */
  readonly stage: string;
  /**
   * Optional detail for the combinatorial phases, where cost is unbounded
   * (REPO_MAP §10.4 — the 5-site search is C(n,5) x sequence count).
   */
  readonly combinationsTested?: number;
  readonly at: IsoTimestamp;
}

/**
 * BACKEND GAP: there is no cancellation mechanism. A running search cannot be
 * interrupted; the Tkinter UI simply freezes for the duration (REPO_MAP §12).
 * Honouring this request requires cooperative cancellation checks inside the
 * Python loops.
 */
export interface AnalysisCancelRequest {
  readonly runId: Id;
  readonly reason?: string;
}

export interface AnalysisCompletedEvent {
  readonly runId: Id;
  readonly result: AnalysisResult;
  readonly durationMs: number;
}

export interface AnalysisFailedEvent {
  readonly runId: Id;
  readonly error: AppError;
}

export type AnalysisEvent =
  | { readonly type: 'progress'; readonly payload: AnalysisProgressEvent }
  | { readonly type: 'completed'; readonly payload: AnalysisCompletedEvent }
  | { readonly type: 'failed'; readonly payload: AnalysisFailedEvent }
  | { readonly type: 'cancelled'; readonly payload: { readonly runId: Id } };

/* ------------------------------------------------------------------ */
/* Errors                                                              */
/* ------------------------------------------------------------------ */

/**
 * Typed error codes the frontend can branch on.
 *
 * BACKEND GAP: Python raises bare `ValueError`s carrying prose messages — e.g.
 * "Sequences are not all the same length. Input must be an aligned FASTA."
 * (REPO_MAP §6.10). Callers can currently only distinguish failures by matching
 * message text, which the existing Python tests already do. Producing these
 * codes requires typed exceptions on the Python side.
 */
export type AppErrorCode =
  | 'FASTA_NOT_FOUND'
  | 'FASTA_NOT_A_FILE'
  | 'FASTA_EMPTY'
  | 'FASTA_NOT_ALIGNED'
  | 'FOCAL_NO_MATCH'
  | 'FOCAL_MATCHES_EVERYTHING'
  | 'OUTPUT_DIR_MISSING'
  | 'OUTPUT_DIR_NOT_A_DIRECTORY'
  | 'INVALID_PARAMETER'
  | 'ANALYSIS_CANCELLED'
  | 'BACKEND_UNAVAILABLE'
  | 'UNKNOWN';

export interface AppError {
  readonly code: AppErrorCode;
  /** Safe to show to the user. */
  readonly message: string;
  /** Which control caused it, for inline field errors. */
  readonly field?: string;
  /**
   * BACKEND GAP: the original Python message/traceback. Retained so nothing is
   * lost while `code` is still being inferred from prose.
   */
  readonly rawMessage?: string;
  readonly detail?: string;
}

/* ------------------------------------------------------------------ */
/* The adapter surface a Python backend will eventually implement       */
/* ------------------------------------------------------------------ */

/**
 * The complete set of operations the renderer expects from a backend.
 *
 * NOT IMPLEMENTED. No transport is chosen; per the task brief there is no
 * FastAPI, no localhost HTTP and no Python bridge in this phase. This interface
 * exists so the eventual adapter has a target, and so the preload API can be
 * widened deliberately rather than ad hoc.
 */
export interface BackendApi {
  loadAlignmentMetadata(path: FilePath): Promise<AlignmentMetadata>;
  validateFocalStrings(request: FocalValidationRequest): Promise<FocalValidationResponse>;

  // BACKEND GAP: requires computation/serialisation split (see MolecularDiagnosisResult).
  runMolecularDiagnosis(
    project: ProjectDraft,
    config: MolecularDiagnosisConfig,
  ): Promise<{ readonly runId: Id }>;

  // BACKEND GAP: no progress seam exists (see AnalysisProgressEvent).
  subscribeToRun(runId: Id, listener: (event: AnalysisEvent) => void): () => void;

  // BACKEND GAP: no cancellation exists (see AnalysisCancelRequest).
  cancelRun(request: AnalysisCancelRequest): Promise<void>;

  // BACKEND GAP: no project format exists (see PersistedProject).
  openProject(request: OpenProjectRequest): Promise<OpenProjectResponse>;
  saveProject(project: PersistedProject): Promise<void>;
}
