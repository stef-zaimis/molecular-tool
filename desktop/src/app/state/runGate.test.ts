import { describe, expect, it } from 'vitest';
import { evaluateRunGate } from './runGate';
import { blankDraft, draftFromPayload, withHeaders } from './focalDrafts';
import type { FocalDraft } from './focalDrafts';
import type { FocalPresenceState } from './projectState';
import type {
  FocalPresenceState as PresenceVerdict,
  SourceStatusPayload,
} from '../../backendContract';

/**
 * Every rule here mirrors one the Python service enforces. These tests are
 * about the REASON arriving before the click — the backend refusal is still
 * the guarantee.
 */

const SAVED: FocalDraft = draftFromPayload({
  id: 'set-1',
  title: 'Targets',
  locked: false,
  entries: [
    { id: 'e1', header: 'focal_1' },
    { id: 'e2', header: 'focal_2' },
  ],
});

function source(overrides: Partial<SourceStatusPayload> = {}): SourceStatusPayload {
  return {
    fastaFileId: 'f1',
    sourcePath: '/data/a.fasta',
    displayName: 'a.fasta',
    state: 'current',
    available: true,
    indexUsable: true,
    exists: true,
    currentSizeBytes: 10,
    currentMtimeNs: 1,
    indexedSizeBytes: 10,
    indexedMtimeNs: 1,
    indexRevision: 1,
    sequenceCount: 4,
    alignmentLength: 8,
    duplicateHeaderCount: 0,
    locked: false,
    message: null,
    ...overrides,
  };
}

function presence(states: Record<string, PresenceVerdict>): FocalPresenceState {
  return {
    status: 'loaded',
    generation: 1,
    error: null,
    byHeader: Object.fromEntries(
      Object.entries(states).map(([header, state]) => {
        const occurrences: Record<string, number> = state === 'missing' ? {} : { f1: 1 };
        return [header, { header, state, occurrences }];
      }),
    ),
  };
}

const ALL_PRESENT = presence({ focal_1: 'present_current', focal_2: 'present_current' });

function gate(overrides: Partial<Parameters<typeof evaluateRunGate>[0]> = {}) {
  return evaluateRunGate({
    draft: SAVED,
    presence: ALL_PRESENT,
    singleFileScope: true,
    scopeSources: [source()],
    running: false,
    ...overrides,
  });
}

describe('run gating', () => {
  it('allows a run when the set is saved and every entry is present', () => {
    expect(gate()).toEqual({ canRun: true, reason: null });
  });

  it('blocks while a run is already going', () => {
    expect(gate({ running: true }).canRun).toBe(false);
  });

  it('blocks when no FASTA is linked', () => {
    const result = gate({ scopeSources: [] });
    expect(result.canRun).toBe(false);
    expect(result.reason).toMatch(/link at least one fasta/i);
  });

  it('blocks when a file in the scope is unavailable, naming it', () => {
    const result = gate({
      scopeSources: [source({ displayName: 'gone.fasta', available: false, state: 'missing' })],
    });
    expect(result.canRun).toBe(false);
    expect(result.reason).toContain('gone.fasta');
  });

  it('allows a stale-but-readable file: the run re-reads it anyway', () => {
    expect(gate({ scopeSources: [source({ state: 'stale', indexUsable: false })] }).canRun).toBe(
      true,
    );
  });

  it('blocks an empty focal set', () => {
    const result = gate({ draft: withHeaders(SAVED, []) });
    expect(result.canRun).toBe(false);
    expect(result.reason).toMatch(/add at least one focal entry/i);
  });

  it('blocks a set that has never been saved, even with entries', () => {
    const draft = withHeaders(blankDraft('New'), ['focal_1', 'focal_2']);
    const result = gate({ draft });
    expect(result.canRun).toBe(false);
    expect(result.reason).toMatch(/save this focal set/i);
  });

  it('blocks a saved set with unsaved changes', () => {
    const result = gate({ draft: withHeaders(SAVED, ['focal_1']) });
    expect(result.canRun).toBe(false);
    expect(result.reason).toMatch(/unsaved changes/i);
  });

  it('blocks until presence has been answered', () => {
    const result = gate({
      presence: { status: 'checking', byHeader: {}, generation: 1, error: null },
    });
    expect(result.canRun).toBe(false);
    expect(result.reason).toMatch(/checking/i);
  });

  it('blocks when an entry is absent from every file in an All-files run', () => {
    const result = gate({
      singleFileScope: false,
      presence: presence({ focal_1: 'present_current', focal_2: 'missing' }),
    });
    expect(result.canRun).toBe(false);
    expect(result.reason).toContain('focal_2');
    expect(result.reason).toMatch(/not in any/i);
  });

  it('blocks when an entry is absent from the one selected file', () => {
    const result = gate({
      presence: presence({ focal_1: 'present_current', focal_2: 'missing' }),
    });
    expect(result.canRun).toBe(false);
    expect(result.reason).toMatch(/not in the selected fasta/i);
  });

  it('blames the RUN scope for an unknown, never an unrelated file', () => {
    /*
     * The gate only ever sees run-scope presence, so an unknown here can only
     * have come from a file this run would read.
     */
    const single = gate({ presence: presence({ focal_1: 'unknown', focal_2: 'present_current' }) });
    expect(single.reason).toMatch(/the selected fasta is unavailable/i);

    const all = gate({
      singleFileScope: false,
      presence: presence({ focal_1: 'unknown', focal_2: 'present_current' }),
    });
    expect(all.reason).toMatch(/a fasta in this run is unavailable/i);
  });

  it('blocks when an entry is only in a file outside the selection, and says so', () => {
    const result = gate({
      presence: presence({ focal_1: 'present_current', focal_2: 'present_other' }),
    });
    expect(result.canRun).toBe(false);
    expect(result.reason).toContain('focal_2');
    // The fix is different from "it does not exist", so the wording is too.
    expect(result.reason).toMatch(/all files/i);
  });

  it('separates "cannot be checked" from "is not there"', () => {
    const result = gate({
      presence: presence({ focal_1: 'present_current', focal_2: 'unknown' }),
    });
    expect(result.canRun).toBe(false);
    expect(result.reason).toMatch(/unavailable/i);
  });

  it('blocks when presence could not be fetched at all', () => {
    const result = gate({
      presence: {
        status: 'failed',
        byHeader: {},
        generation: 1,
        error: { code: 'UNKNOWN', message: 'boom' },
      },
    });
    expect(result.canRun).toBe(false);
  });

  it('blocks when an entry has no verdict yet, rather than assuming it is fine', () => {
    const result = gate({ presence: presence({ focal_1: 'present_current' }) });
    expect(result.canRun).toBe(false);
  });

  it('a locked set that is saved and present may still run', () => {
    const locked = draftFromPayload({
      id: 'set-1',
      title: 'Targets',
      locked: true,
      entries: [
        { id: 'e1', header: 'focal_1' },
        { id: 'e2', header: 'focal_2' },
      ],
    });
    expect(gate({ draft: locked }).canRun).toBe(true);
  });
});
