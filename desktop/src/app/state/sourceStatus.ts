/**
 * Turning backend source states into what the workspace actually shows.
 *
 * The rule this file exists to enforce: a linked FASTA's state is a LIVE
 * READING of the filesystem, never a stored property of the project. The
 * renderer therefore has no opinion of its own about whether a file is
 * current — it only decides how to phrase, colour and gate what Python
 * reported.
 *
 * Three ideas are kept deliberately separate, because conflating them is what
 * produces a UI that lies:
 *
 *  - `available`   — the file can be read right now.
 *  - `indexUsable` — the stored header index provably describes those bytes.
 *  - `activity`    — what the app is doing about it at this instant.
 *
 * A file can be available with an unusable index (it was edited), and a file
 * can have a usable index while a check is in flight. Neither implies the
 * other.
 */

import type { SourceState, SourceStatusPayload } from '../../backendContract';

/** What the app is currently doing to a source, as opposed to what it is. */
export type SourceActivity =
  | 'idle'
  /** A cheap stat check, or a strong hash, is in flight. */
  | 'verifying'
  /** The header index is being rebuilt from the file. */
  | 'reindexing'
  /** The user is pointing the link at a different path. */
  | 'relinking';

/** Severity, so the panel can order and colour without re-deriving meaning. */
export type SourceTone = 'ok' | 'pending' | 'warning' | 'error';

export interface SourceView {
  readonly payload: SourceStatusPayload;
  readonly activity: SourceActivity;
  /** Short status word shown next to the file name. */
  readonly label: string;
  /** One sentence explaining the consequence, not the mechanism. */
  readonly detail: string;
  readonly tone: SourceTone;
  /** True when the user must do something before this file can be analysed. */
  readonly needsAttention: boolean;
  /** True when this file may take part in an analysis right now. */
  readonly analysable: boolean;
  readonly canRelink: boolean;
  readonly canReindex: boolean;
}

/**
 * Activity wins over state for the LABEL, because a stale-looking file that is
 * actively being rebuilt should not be announced as broken. It never wins over
 * `analysable`: work in progress is not a reason to let a run start.
 */
const ACTIVITY_LABELS: Readonly<Record<Exclude<SourceActivity, 'idle'>, string>> = {
  verifying: 'Checking',
  reindexing: 'Re-indexing',
  relinking: 'Relinking',
};

const ACTIVITY_DETAILS: Readonly<Record<Exclude<SourceActivity, 'idle'>, string>> = {
  verifying: 'Confirming this file still matches what was indexed.',
  reindexing: 'Reading the file and rebuilding its header index.',
  relinking: 'Pointing this project at the file you chose.',
};

interface StateDescription {
  readonly label: string;
  readonly detail: string;
  readonly tone: SourceTone;
}

const STATE_DESCRIPTIONS: Readonly<Record<SourceState, StateDescription>> = {
  missing: {
    label: 'Missing',
    detail: 'This file is not at the path it was linked from. Relink it to continue.',
    tone: 'error',
  },
  unreadable: {
    label: 'Unreadable',
    detail: 'The path exists but cannot be read as a file.',
    tone: 'error',
  },
  never_indexed: {
    label: 'Not indexed',
    detail: 'This file has not been read yet.',
    tone: 'warning',
  },
  unverified: {
    // Not an error: the size and timestamp still match. It simply has not been
    // proven byte for byte during this session yet.
    label: 'Unverified',
    detail: 'Looks unchanged, but has not been fully checked yet this session.',
    tone: 'pending',
  },
  current: {
    label: 'Current',
    detail: 'Verified against the file on disk.',
    tone: 'ok',
  },
  stale: {
    label: 'Changed on disk',
    detail: 'The file has been edited since it was indexed. Re-index it to use it.',
    tone: 'warning',
  },
};

export function describeSource(
  payload: SourceStatusPayload,
  activity: SourceActivity = 'idle',
): SourceView {
  const base = STATE_DESCRIPTIONS[payload.state];
  const busy = activity !== 'idle';

  return {
    payload,
    activity,
    label: busy ? ACTIVITY_LABELS[activity] : base.label,
    // Python's own message is more specific when it bothered to send one.
    detail: busy ? ACTIVITY_DETAILS[activity] : payload.message ?? base.detail,
    tone: busy ? 'pending' : base.tone,
    needsAttention: !busy && (payload.state === 'missing' || payload.state === 'unreadable' || payload.state === 'stale' || payload.state === 'never_indexed'),
    // `unverified` is analysable: the run itself reads and verifies the file,
    // so blocking on it would just mean hashing everything twice.
    analysable: !busy && payload.available && payload.state !== 'stale' && payload.state !== 'never_indexed',
    canRelink: !busy,
    canReindex: !busy && payload.available,
  };
}

export function describeSources(
  payloads: readonly SourceStatusPayload[],
  activities: Readonly<Record<string, SourceActivity>> = {},
): readonly SourceView[] {
  return payloads.map((payload) =>
    describeSource(payload, activities[payload.fastaFileId] ?? 'idle'),
  );
}

/** Sources the user has to deal with, worst first. */
const TONE_ORDER: Readonly<Record<SourceTone, number>> = {
  error: 0,
  warning: 1,
  pending: 2,
  ok: 3,
};

export function sourcesNeedingAttention(views: readonly SourceView[]): readonly SourceView[] {
  return views
    .filter((view) => view.needsAttention)
    .slice()
    .sort((a, b) => TONE_ORDER[a.tone] - TONE_ORDER[b.tone]);
}

/**
 * The single banner line shown on entering a project, or null when everything
 * is fine. Missing files are called out by name because "some files are
 * missing" gives the user nothing to act on.
 */
export function sourceBannerMessage(views: readonly SourceView[]): string | null {
  const missing = views.filter(
    (view) => view.payload.state === 'missing' || view.payload.state === 'unreadable',
  );
  const stale = views.filter((view) => view.payload.state === 'stale');

  const parts: string[] = [];
  if (missing.length > 0) {
    parts.push(
      missing.length === 1
        ? `${missing[0].payload.displayName} is missing from disk.`
        : `${missing.length} linked files are missing from disk.`,
    );
  }
  if (stale.length > 0) {
    parts.push(
      stale.length === 1
        ? `${stale[0].payload.displayName} has changed since it was indexed.`
        : `${stale.length} linked files have changed since they were indexed.`,
    );
  }

  return parts.length > 0 ? parts.join(' ') : null;
}

/** Whether a run over this exact scope may start. */
export function scopeIsAnalysable(
  views: readonly SourceView[],
  fastaFileIds: readonly string[],
): boolean {
  if (fastaFileIds.length === 0) return false;
  return fastaFileIds.every((id) => {
    const view = views.find((candidate) => candidate.payload.fastaFileId === id);
    return view !== undefined && view.analysable;
  });
}
