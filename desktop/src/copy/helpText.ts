/**
 * Every tooltip string in the application, in one place.
 *
 * WHY CENTRAL: the final scientific wording has not been signed off. Keeping
 * the copy here means it can be reviewed and corrected in one pass, instead of
 * invented ad hoc inside components and then quietly diverging.
 *
 * STATUS KEY
 *   FROM DESIGN — reproduced verbatim from the reference mockups.
 *   PROVISIONAL — written from the observed behaviour of the Python code
 *                 (REPO_MAP.md); factual about what the code does, but NOT
 *                 reviewed as user-facing scientific wording. Needs sign-off.
 *
 * Nothing here should assert a scientific claim that the code does not
 * implement. Where the implementation is surprising, the copy says so plainly
 * rather than describing the intent.
 */

export const HELP_TEXT = {
  /* ---- Analyses (project creation) ---- */

  /** PROVISIONAL. */
  molecularDiagnosis:
    'Molecular Diagnosis searches an aligned FASTA for diagnostic molecular characters: ' +
    'alignment positions, or combinations of positions, whose states are shared by the ' +
    'focal group and by no sequence outside it. It reports the combinations found, the ' +
    'sites involved, and a similarity comparison against the non-focal sequences.',

  /** PROVISIONAL. */
  sequencePunishmentTest:
    'The Sequence Punishment Test scores each sequence in the focal group for internal ' +
    'disagreement. It examines the focal group only, weighing base-state polymorphism, ' +
    'gap-heavy regions, and terminal versus internal gaps, and ranks sequences by the ' +
    'total score so that anomalous members can be reviewed.',

  /** FROM DESIGN — verbatim from docs/02-project-creation.png. */
  consensusSequenceGeneration:
    'A consensus sequence generator combines aligned sequences from a selected focal ' +
    'group into a single representative sequence. It identifies the most common bases at ' +
    'each position, uses IUPAC ambiguity codes where variation exists, optionally removes ' +
    'gap-heavy columns, and produces consensus FASTA files plus a report comparing unique ' +
    'sequence groups with the consensus.',

  /* ---- Molecular Diagnosis parameters ---- */

  /**
   * PROVISIONAL.
   * Frontend name `ignoreGaps` maps 1:1 onto the Python keyword
   * `include_gappy_consensus_dmc_sites` — same polarity, no inversion.
   */
  ignoreGaps:
    'When enabled, alignment columns where part of the focal group has a gap can still ' +
    'contribute a diagnostic character, and a column may qualify on gap differences ' +
    'alone. When disabled, any gap in the focal column disqualifies the site. Columns ' +
    'that are mostly gaps within the focal group are excluded either way.',

  /**
   * PROVISIONAL.
   * Frontend name `giveBenefitOfDoubtToAmbiguousBases` maps 1:1 onto
   * `include_ambiguous_dmc_bd`.
   */
  benefitOfDoubt:
    'When enabled, IUPAC ambiguity codes are treated as compatible with any base they ' +
    'could represent, so an ambiguous state does not by itself rule a site in or out. ' +
    'When disabled, a site is skipped if an ambiguity code appears anywhere in that ' +
    'column, including in non-focal sequences.',

  /**
   * PROVISIONAL — and the note about the search floor is load-bearing, not a
   * caveat. See REPO_MAP.md §4.4 and §10.1, and UI_NOTES Q6.
   */
  minCandidateSize:
    'The smallest combination size the search is allowed to stop at. Note that the ' +
    'search still begins at size 1 regardless of this value: raising it does not skip ' +
    'smaller combinations, it only prevents the search from finishing before this size ' +
    'is reached.',

  /** PROVISIONAL. */
  maxCandidateSize:
    'The largest combination size the search will test. The search stops at the first ' +
    'size that yields results at or above the minimum; if it reaches this limit without ' +
    'stopping, it can be continued from the next size up.',

  /** PROVISIONAL. */
  focalSetTitle:
    'A name for this focal set, so it can be recognised later. The name is used for ' +
    'labelling only and does not affect which sequences are selected.',

  /** PROVISIONAL. */
  focalString:
    'Each string is matched against the FASTA headers. A string matches a header if it ' +
    'appears anywhere within it, and matching is case sensitive. Green means the string ' +
    'matches at least one header, red means it matches none.',
} as const;

export type HelpTextKey = keyof typeof HELP_TEXT;

/**
 * Copy shown when an action exists in the design but is not implemented yet.
 * Kept here so unfinished states read consistently instead of leaking
 * developer placeholders into the UI.
 */
export const UNFINISHED_TEXT = {
  openExistingProject:
    'Opening a saved project is not available yet. This build has no project file ' +
    'format, so there is nothing to open. Create a new project to continue.',
  loadFocalStrings:
    'Loading focal strings from a file is not available yet. There is no focal-set ' +
    'file format in this build.',
  saveFocalSet:
    'Saving a focal set is not available yet. The set is kept for this session only ' +
    'and is not written to disk.',
} as const;
