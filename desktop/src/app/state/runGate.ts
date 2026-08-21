/**
 * When may Molecular Diagnosis start?
 *
 * The backend refuses an unsound run on its own — that is the guarantee, and
 * this file does not replace it. What this adds is a REASON, before the click:
 * the user should see "one focal entry is not in the selected file" while they
 * can still fix it, not as a failed run afterwards.
 *
 * Every rule here mirrors one the service enforces:
 *
 *   dirty draft        -> nothing to run; the saved set is what a run names
 *   empty membership   -> FOCAL_EMPTY
 *   presence not green -> FOCAL_ENTRIES_NOT_IN_FILE / _NOT_IN_SCOPE
 *   presence unknown   -> FOCAL_PRESENCE_UNKNOWN
 *   unusable source    -> SOURCE_UNAVAILABLE
 *
 * If the two ever disagree, the backend wins and the run fails with its own
 * message. This is a courtesy, never a substitute.
 */

import type { HeaderPresencePayload, SourceStatusPayload } from '../../backendContract';
import type { FocalDraft } from './focalDrafts';
import { isDirty, isNewDraft, normaliseHeaders } from './focalDrafts';
import type { FocalPresenceState } from './projectState';

export interface RunGate {
  readonly canRun: boolean;
  /** One sentence naming what to fix, or null when the run may start. */
  readonly reason: string | null;
}

interface RunGateInput {
  readonly draft: FocalDraft;
  /**
   * Presence measured against the RUN SCOPE, not the display comparison scope.
   *
   * Passing the display answer here was a correctness bug: it is computed
   * against every linked file, so an unrelated unavailable FASTA made entries
   * unknown and the gate blamed the selected file for it.
   */
  readonly presence: FocalPresenceState;
  /** True when exactly one FASTA is selected, which is what makes orange real. */
  readonly singleFileScope: boolean;
  /** The linked files this run would actually read. */
  readonly scopeSources: readonly SourceStatusPayload[];
  readonly running: boolean;
}

function unusable(sources: readonly SourceStatusPayload[]): readonly SourceStatusPayload[] {
  // `available` is "readable right now". A stale index is fine for a run: the
  // run re-reads the file anyway and refreshes the index from what it read.
  return sources.filter((source) => !source.available);
}

export function evaluateRunGate({
  draft,
  presence,
  singleFileScope,
  scopeSources,
  running,
}: RunGateInput): RunGate {
  if (running) return { canRun: false, reason: 'A run is already in progress.' };

  if (scopeSources.length === 0) {
    return { canRun: false, reason: 'Link at least one FASTA file to this project first.' };
  }

  const blocked = unusable(scopeSources);
  if (blocked.length > 0) {
    return {
      canRun: false,
      reason: `${blocked[0].displayName} is not available. Relink it before running.`,
    };
  }

  const headers = normaliseHeaders(draft.headers);
  if (headers.length === 0) {
    return { canRun: false, reason: 'Add at least one focal entry before running.' };
  }

  // Checked before dirtiness so a brand-new draft with entries says "save it",
  // not "add entries".
  if (isNewDraft(draft)) {
    return { canRun: false, reason: 'Save this focal set before running an analysis with it.' };
  }
  if (isDirty(draft)) {
    return {
      canRun: false,
      reason: 'This focal set has unsaved changes. Save it before running.',
    };
  }

  if (presence.status === 'failed') {
    return { canRun: false, reason: 'The focal entries could not be checked against the FASTA files.' };
  }
  if (presence.status !== 'loaded') {
    return { canRun: false, reason: 'Checking where the focal entries are...' };
  }

  const verdicts: (HeaderPresencePayload | undefined)[] = headers.map(
    (header) => presence.byHeader[header],
  );

  if (verdicts.some((verdict) => verdict === undefined)) {
    return { canRun: false, reason: 'Checking where the focal entries are...' };
  }

  /*
   * `unknown` here means the RUN scope could not answer — i.e. a file this run
   * would actually read is unavailable. A file elsewhere in the project being
   * unavailable cannot reach this, because it is not in the run scope.
   */
  const unknown = verdicts.filter((verdict) => verdict?.state === 'unknown');
  if (unknown.length > 0) {
    return {
      canRun: false,
      reason: singleFileScope
        ? `${unknown[0]?.header} cannot be located: the selected FASTA is unavailable.`
        : `${unknown[0]?.header} cannot be located: a FASTA in this run is unavailable.`,
    };
  }

  const absent = verdicts.filter((verdict) => verdict?.state !== 'present_current');
  if (absent.length > 0) {
    const first = absent[0];
    return {
      canRun: false,
      reason: singleFileScope
        ? `${first?.header} is not in the selected FASTA. Choose All files, or a file that contains it.`
        : `${first?.header} is not in any of the selected FASTA files.`,
    };
  }

  return { canRun: true, reason: null };
}
