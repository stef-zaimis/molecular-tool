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
  | 'UNKNOWN';

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
