/**
 * Mock FASTA headers for the UI-only phase.
 *
 * A small representative fixture — NOT a generated dataset. The real alignment
 * has 2354 headers (REPO_MAP.md §9); reproducing that volume here would prove
 * nothing about the interaction and would slow every render.
 *
 * Header shape is taken from the real file
 * (input/Leptacis_allSequences-BOLD-09March2026_aln.fasta):
 *
 *     <BOLD process id>|<country code>|<taxon>
 *     AACTA5253-20|AU|Leptacis
 *
 * with a subset carrying the species-level focal taxon `Leptacis_tipulae`,
 * which is the focal string the legacy analysis script uses.
 *
 * BACKEND GAP: real headers will come from the Python FASTA parser. Nothing in
 * this app parses FASTA — see contract.ts `AlignmentMetadata`.
 */
export const MOCK_FASTA_HEADERS: readonly string[] = [
  'AACTA5253-20|AU|Leptacis',
  'ABOTH6382-22|CA|Leptacis',
  'ABWYT17523-24|CA|Leptacis',
  'ACGAT8891-21|US|Leptacis',
  'ADFRT2210-19|FR|Leptacis',
  'AEKLM4417-23|DE|Leptacis',
  'BBHYQ1043-18|GB|Leptacis_tipulae',
  'BCJKR7756-20|NL|Leptacis_tipulae',
  'BDMNP3391-21|SE|Leptacis_tipulae',
  'BEQRS9028-22|NO|Leptacis_tipulae',
  'BFTUV1164-23|FI|Leptacis_tipulae',
  'CGHIJ5502-19|ES|Leptacis_phantasmatica',
  'CHKLM6613-20|PT|Leptacis_phantasmatica',
  'CIMNO7724-21|IT|Leptacis_phantasmatica',
  'DJPQR8835-22|GR|Platygaster',
  'DKSTU9946-23|TR|Platygaster',
  'ELVWX1057-24|PL|Synopeas',
  'EMYZA2168-24|CZ|Synopeas',
];

/** Headers a user is likely to try in the STRING editor, for tooltip examples. */
export const MOCK_HEADER_SOURCE_NOTE =
  'Validated against a representative fixture of FASTA headers, not a loaded file.';
