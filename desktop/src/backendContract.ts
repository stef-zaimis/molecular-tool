/**
 * backendContract.ts — the shapes that actually cross the Python boundary.
 *
 * These mirror, field for field, what `molecular_diagnosis.service.handlers`
 * returns. Unlike `contract.ts` (which describes the eventual architecture,
 * including parts that do not exist yet), everything in this file is
 * IMPLEMENTED and exercised end to end.
 *
 * Keep it free of Electron and React types: it is shared by the preload, the
 * main process and the renderer.
 */

/** Stable error codes produced by `service/errors.py`. */
export type BackendErrorCode =
  | 'INVALID_REQUEST'
  | 'UNKNOWN_METHOD'
  | 'FASTA_NOT_FOUND'
  | 'FASTA_NOT_A_FILE'
  | 'FASTA_UNREADABLE'
  | 'FASTA_EMPTY'
  | 'FASTA_NOT_ALIGNED'
  | 'FOCAL_EMPTY'
  | 'FOCAL_NO_MATCH'
  | 'FOCAL_MATCHES_EVERYTHING'
  | 'OUTPUT_DIR_MISSING'
  | 'OUTPUT_DIR_NOT_A_DIRECTORY'
  | 'OUTPUT_WRITE_FAILED'
  | 'INVALID_PARAMETER'
  | 'BACKEND_UNAVAILABLE'
  | 'UNKNOWN'
  /* Project-level codes. Raised by `project/service.py` and passed through the
     boundary with their own code, so the renderer can branch on them. */
  | 'NO_PROJECT_OPEN'
  | 'PROJECT_NOT_FOUND'
  | 'PROJECT_ALREADY_EXISTS'
  | 'UNKNOWN_FILE'
  | 'PATH_ALREADY_LINKED'
  | 'SOURCE_UNAVAILABLE'
  | 'SOURCE_CHANGED_DURING_READ'
  | 'SEARCH_SCOPE_UNAVAILABLE'
  | 'UNKNOWN_FOCAL_SET'
  | 'FOCAL_SET_LOCKED'
  | 'FOCAL_QUERY_EMPTY'
  | 'FOCAL_ENTRIES_NOT_IN_FILE'
  | 'FOCAL_ENTRIES_NOT_IN_SCOPE'
  | 'FOCAL_PRESENCE_UNKNOWN'
  | 'NO_FILES_SELECTED'
  | 'SCOPE_MISMATCH'
  | 'DUPLICATE_HEADER_IN_FILE'
  | 'DUPLICATE_HEADER_ACROSS_FILES'
  | 'INCOMPATIBLE_ALIGNMENT_LENGTHS'
  | 'FASTA_FILE_LOCKED';

export interface BackendError {
  readonly code: BackendErrorCode | string;
  /** Safe to show a user. */
  readonly message: string;
  /** Original Python wording. For logs and debugging, not for the UI. */
  readonly detail?: string;
  /** Python traceback, when one was available. */
  readonly traceback?: string;
}

/**
 * Backend failures are DATA, not exceptions. Every call resolves; callers
 * branch on `ok` rather than catching.
 */
export type BackendResult<T> =
  | { readonly ok: true; readonly result: T }
  | { readonly ok: false; readonly error: BackendError };

/* ------------------------------------------------------------------ */
/* loadFasta                                                          */
/* ------------------------------------------------------------------ */

export interface FastaLoadResult {
  readonly path: string;
  /** Distinct headers. Repeats collapse in the parser, as they do in analysis. */
  readonly sequenceCount: number;
  readonly alignmentLength: number;
  readonly headers: readonly string[];
  /** '>' lines in the file, which may exceed `sequenceCount`. */
  readonly headerLineCount: number;
  readonly duplicateHeaderCount: number;
  readonly alphabet: readonly string[];
}

/* ------------------------------------------------------------------ */
/* validateFocalStrings                                               */
/* ------------------------------------------------------------------ */

export interface FocalStringValidation {
  readonly focalString: string;
  /** Green when true, red when false. */
  readonly matched: boolean;
  readonly matchCount: number;
  readonly sampleHeaders: readonly string[];
  readonly truncated: boolean;
  /** True for an entry that cannot be matched at all, e.g. an empty string. */
  readonly invalid: boolean;
}

export interface FocalValidationResult {
  readonly path: string;
  readonly totalHeaders: number;
  /** One entry per supplied string, in the order supplied. */
  readonly results: readonly FocalStringValidation[];
  /** Headers matched by AT LEAST ONE string — the prospective focal set. */
  readonly unionMatchCount: number;
  readonly unionSampleHeaders: readonly string[];
  readonly matchedBy: Readonly<Record<string, readonly string[]>>;
}

/* ------------------------------------------------------------------ */
/* runMolecularDiagnosis                                              */
/* ------------------------------------------------------------------ */

/**
 * Continuation state.
 *
 * Opaque to the frontend: it is produced by the backend and handed straight
 * back to continue a search. Passing it back is what makes a continuation
 * resume rather than restart.
 */
export interface DiagnosisResumeState {
  readonly startCombinationLength: number;
  readonly diagnosticCombinations: readonly (readonly number[])[];
  readonly combinationsTestedByLength: Readonly<Record<string, number>>;
}

export interface MolecularDiagnosisRequest {
  readonly fastaPath: string;
  readonly focalStrings: readonly string[];
  readonly outputDirectory: string;
  /** Python: `include_gappy_consensus_dmc_sites` (same polarity). */
  readonly ignoreGaps: boolean;
  /** Python: `include_ambiguous_dmc_bd` (same polarity). */
  readonly giveBenefitOfDoubtToAmbiguousBases: boolean;
  /** Python: `min_combination_length`. */
  readonly minCandidateSize: number;
  /** Python: `max_combination_length`. */
  readonly maxCandidateSize: number;
  /** Omit to start a fresh search; pass a previous result's `resume` to continue. */
  readonly resume?: DiagnosisResumeState | null;
}

export interface DiagnosisDiagnostics {
  readonly fixedCount: number;
  readonly skippedSites: number;
  readonly globallyConservedRemoved: number;
  readonly candidateCount: number;
  readonly pairsTested: number;
  readonly totalCombinationsTested: number;
  readonly combinationsTestedByLength: Readonly<Record<string, number>>;
  readonly ambiguousBdSitesIncluded: readonly number[];
  readonly gappyConsensusSitesIncluded: readonly number[];
}

export type DiagnosisStopReason =
  | 'reached_maximum_length'
  | 'found_at_or_above_minimum_length'
  | 'no_candidate_sites'
  | 'start_length_exceeds_maximum_length';

export interface DiagnosisDmcResult {
  /** Combination size (as a string key) -> the combinations found at that size. */
  readonly combinationsByLength: Readonly<Record<string, readonly (readonly number[])[]>>;
  readonly singleSites: readonly number[];
  readonly pairs: readonly (readonly number[])[];
  readonly uniqueSites: readonly number[];
  /** Zero-based site index (as a string key) -> consensus state. */
  readonly states: Readonly<Record<string, string>>;
  readonly stopReason: DiagnosisStopReason;
  readonly stoppedAtLength: number;
  readonly minCombinationLength: number;
  readonly maxCombinationLength: number;
  readonly startCombinationLength: number;
  readonly diagnostics: DiagnosisDiagnostics;
}

export interface MolecularDiagnosisResult {
  /** The selectors actually used, after trimming and deduplication. */
  readonly focalStrings: readonly string[];
  readonly outputs: {
    readonly reportTxt: string;
    readonly workbookXlsx: string;
    readonly consensusTxt: string | null;
  };
  readonly dmc: DiagnosisDmcResult;
  /** True when the search stopped only because it hit the configured maximum. */
  readonly canContinue: boolean;
  readonly resume: DiagnosisResumeState | null;
}

/* ------------------------------------------------------------------ */
/* Projects: persistent, SQLite-backed                                */
/* ------------------------------------------------------------------ */

/**
 * Runtime state of a linked FASTA. Mirrors `project/sources.py::SourceState`.
 *
 * None of this is persisted: it is recomputed from the filesystem every time.
 * The renderer must treat it as a live reading, not as a stored property of
 * the project.
 */
export type SourceState =
  | 'missing'
  | 'unreadable'
  | 'never_indexed'
  | 'unverified'
  | 'current'
  | 'stale';

/** Mirrors `SourceStatus.to_payload()`. */
export interface SourceStatusPayload {
  readonly fastaFileId: string;
  readonly sourcePath: string;
  readonly displayName: string;
  readonly state: SourceState;
  /** Readable right now. Says nothing about whether the index still matches. */
  readonly available: boolean;
  /** May the stored header index be trusted for arbitrary searching? */
  readonly indexUsable: boolean;
  readonly exists: boolean;
  readonly currentSizeBytes: number | null;
  readonly currentMtimeNs: number | null;
  readonly indexedSizeBytes: number | null;
  readonly indexedMtimeNs: number | null;
  readonly indexRevision: number;
  readonly sequenceCount: number | null;
  readonly alignmentLength: number | null;
  readonly duplicateHeaderCount: number;
  /**
   * A locked source may be analysed but not unlinked or renamed.
   *
   * Locking protects the LINK, not the data — the backend refuses the unlink,
   * so this is a property of the project rather than a hidden button.
   */
  readonly locked: boolean;
  readonly message: string | null;
}

/**
 * A FASTA vetted BEFORE it is linked.
 *
 * Produced by `project.validateFastaCandidate`, which needs no open project —
 * the new-project screen has to check files before any database exists. It runs
 * the same scan the indexer does, so a file accepted here cannot be rejected at
 * link time for a reason the user was never shown.
 */
export interface FastaCandidate {
  readonly path: string;
  readonly displayName: string;
  readonly sequenceCount: number;
  readonly alignmentLength: number;
  readonly duplicateHeaderCount: number;
}

export interface ProjectCapabilities {
  readonly sqliteVersion: string;
  readonly fts5: boolean;
  readonly trigram: boolean;
  /** True only when the FTS accelerator is actually installed and usable. */
  readonly acceleratedSearch: boolean;
}

export interface ProjectMetadata {
  readonly projectUuid: string;
  readonly title: string;
}

export interface OpenProjectResult {
  readonly projectDir: string;
  readonly outputsDir: string;
  readonly metadata: ProjectMetadata;
  readonly capabilities: ProjectCapabilities;
  /** Every linked file, already status-checked. Missing files appear here. */
  readonly sources: readonly SourceStatusPayload[];
}

export interface HeaderHitPayload {
  readonly fastaFileId: string;
  readonly ordinal: number;
  readonly header: string;
  readonly recordStartByte: number | null;
}

/**
 * Search is non-mutating, so a file it could not read is reported rather than
 * fatal. `unavailable` is what the UI must show before the user trusts the
 * hit list as complete.
 *
 * `+` does NOT behave this way: it refuses instead of persisting a partial
 * expansion. See `addFocalEntries`.
 */
export interface SearchHeadersResult {
  readonly query: string;
  readonly hits: readonly HeaderHitPayload[];
  readonly unavailable: readonly SourceStatusPayload[];
}

/**
 * Where a focal entry currently is. Mirrors `locations.py::PresenceState`.
 *
 * `present_other` only ever appears when a specific file is selected: it means
 * "not in the selected file, but present elsewhere in the project".
 */
export type FocalPresenceState = 'present_current' | 'present_other' | 'missing' | 'unknown';

export interface FocalPresencePayload {
  readonly entryId: string;
  readonly header: string;
  readonly state: FocalPresenceState;
  /** fastaFileId -> how many records in that file carry this exact header. */
  readonly occurrences: Readonly<Record<string, number>>;
}

export interface FocalEntryPayload {
  readonly id: string;
  readonly header: string;
}

/**
 * A focal set and its complete explicit membership.
 *
 * `locked` is enforced in Python, not by disabling controls: a locked set may
 * be selected, presence-checked and analysed, but every mutation is refused
 * with `FOCAL_SET_LOCKED`. Grey the buttons out for clarity, never for safety.
 */
export interface FocalSetPayload {
  readonly id: string;
  readonly title: string;
  readonly locked: boolean;
  readonly entries: readonly FocalEntryPayload[];
}

/**
 * The result of replacing the large textbox's contents.
 *
 * Applied as a diff: `kept` headers keep their entry ids and their cached
 * locations, so a debounced editor may resend freely.
 */
export interface FocalReplacementResult {
  readonly focalSetId: string;
  /** Exactly what is now stored, in textbox order, trimmed and deduplicated. */
  readonly headers: readonly string[];
  readonly added: readonly string[];
  readonly removed: readonly string[];
  readonly kept: readonly string[];
}

export interface FocalQueryResult {
  readonly query: string;
  /** Every exact header the query resolved to. Never the query itself. */
  readonly matched: readonly string[];
  readonly added: readonly string[];
}

/**
 * The explicit Save the workspace performs.
 *
 * `focalSetId: null` creates the set and its entries in one transaction; an id
 * updates title + membership in place, preserving entry ids (and their cached
 * locations) for headers that survive the diff.
 *
 * Nothing else in focal editing writes: the renderer holds a working copy, and
 * this is the only call that commits it.
 */
export interface SaveFocalSetRequest {
  readonly focalSetId?: string | null;
  readonly title: string;
  /** Exact headers. Ones matching no FASTA are still valid members. */
  readonly headers: readonly string[];
}

/**
 * Presence for arbitrary headers, including ones in an UNSAVED draft.
 *
 * Answered from the header index with one batched lookup: no FASTA is opened,
 * hashed or reindexed, and no `focal_entry_location` row is created. Safe to
 * call on a debounce while the user types.
 */
export interface HeaderPresencePayload {
  readonly header: string;
  readonly state: FocalPresenceState;
  /** fastaFileId -> how many records in that file carry this exact header. */
  readonly occurrences: Readonly<Record<string, number>>;
}

/**
 * What `+` WOULD add, resolved but not added.
 *
 * Uncapped and non-mutating. `SearchHeadersResult` is a capped PREVIEW and must
 * never be used to expand `+`: a query matching 300 headers has to yield all
 * 300, or the focal set silently ends up smaller than the user asked for.
 */
export interface ResolveFocalAddQueryResult {
  readonly query: string;
  /** Every matching complete header, in project-file then record order. */
  readonly headers: readonly string[];
}

/**
 * `-` over a working copy.
 *
 * Runs in Python so the renderer does not need a second implementation of
 * `str.casefold()`; `toLowerCase()` is not the same function, and the two
 * would eventually disagree about which member a `-` removes.
 */
export interface MatchFocalHeadersResult {
  readonly query: string;
  readonly matched: readonly string[];
}

export interface RelinkResult {
  /** True when the new path holds byte-identical content, so nothing was rebuilt. */
  readonly identical: boolean;
  readonly reindexed: boolean;
  readonly indexRevision: number;
  readonly source: SourceStatusPayload;
}

export interface ProjectDiagnosisRequest {
  readonly focalSetId: string;
  readonly fastaFileIds: readonly string[];
  /**
   * Any run is refused unless every focal entry is present somewhere in the
   * files being analysed, because an absent entry would silently shrink the
   * focal group. `singleFile` only changes which refusal code comes back
   * (`FOCAL_ENTRIES_NOT_IN_FILE` vs `FOCAL_ENTRIES_NOT_IN_SCOPE`), and asserts
   * that exactly one file was selected.
   */
  readonly singleFile: boolean;
  readonly options: {
    readonly ignoreGaps: boolean;
    readonly giveBenefitOfDoubtToAmbiguousBases: boolean;
    readonly minCandidateSize: number;
    readonly maxCandidateSize: number;
  };
  readonly resume?: DiagnosisResumeState | null;
}

export interface ProjectDiagnosisResult {
  readonly focalSetId: string;
  /** The exact headers used as the focal group. Not substrings. */
  readonly focalHeaders: readonly string[];
  readonly fastaFileIds: readonly string[];
  readonly sequenceCount: number;
  readonly alignmentLength: number | null;
  readonly outputs: {
    readonly reportTxt: string;
    readonly workbookXlsx: string;
    readonly consensusTxt: string | null;
  };
  readonly dmc: DiagnosisDmcResult;
  readonly canContinue: boolean;
  readonly resume: DiagnosisResumeState | null;
}
