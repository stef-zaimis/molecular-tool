import { describe, expect, it } from 'vitest';
import {
  describeSource,
  describeSources,
  scopeIsAnalysable,
  sourceBannerMessage,
  sourcesNeedingAttention,
} from './sourceStatus';
import type { SourceState, SourceStatusPayload } from '../../backendContract';

function payload(
  state: SourceState,
  overrides: Partial<SourceStatusPayload> = {},
): SourceStatusPayload {
  const available = state !== 'missing' && state !== 'unreadable';
  return {
    fastaFileId: overrides.fastaFileId ?? `file-${state}`,
    sourcePath: `/data/${state}.fasta`,
    displayName: `${state}.fasta`,
    state,
    available,
    indexUsable: state === 'current',
    exists: available,
    currentSizeBytes: available ? 100 : null,
    currentMtimeNs: available ? 1 : null,
    indexedSizeBytes: 100,
    indexedMtimeNs: 1,
    indexRevision: 1,
    sequenceCount: 4,
    alignmentLength: 12,
    duplicateHeaderCount: 0,
    message: null,
    ...overrides,
  };
}

describe('describeSource', () => {
  it('reports a verified file as current and analysable', () => {
    const view = describeSource(payload('current'));
    expect(view.label).toBe('Current');
    expect(view.tone).toBe('ok');
    expect(view.needsAttention).toBe(false);
    expect(view.analysable).toBe(true);
  });

  it('treats an unverified file as pending, not as a problem', () => {
    // Size and mtime still match; it simply has not been hashed this session.
    // Showing it as broken would train users to ignore the warning that means
    // something.
    const view = describeSource(payload('unverified'));
    expect(view.tone).toBe('pending');
    expect(view.needsAttention).toBe(false);
    // The run reads and verifies the file anyway, so blocking here would only
    // mean hashing it twice.
    expect(view.analysable).toBe(true);
  });

  it('blocks a changed file until it is re-indexed', () => {
    const view = describeSource(payload('stale'));
    expect(view.tone).toBe('warning');
    expect(view.needsAttention).toBe(true);
    expect(view.analysable).toBe(false);
    expect(view.canReindex).toBe(true);
  });

  it('blocks a missing file and offers only relinking', () => {
    const view = describeSource(payload('missing'));
    expect(view.tone).toBe('error');
    expect(view.analysable).toBe(false);
    expect(view.canRelink).toBe(true);
    // Nothing to read, so re-indexing is not on the table.
    expect(view.canReindex).toBe(false);
  });

  it('prefers the backend message over the generic wording', () => {
    const view = describeSource(payload('missing', { message: 'D:/old/x.fasta is gone.' }));
    expect(view.detail).toBe('D:/old/x.fasta is gone.');
  });
});

describe('activity', () => {
  it('announces work in progress instead of the underlying state', () => {
    // A file being rebuilt should not still read as "Changed on disk".
    const view = describeSource(payload('stale'), 'reindexing');
    expect(view.label).toBe('Re-indexing');
    expect(view.tone).toBe('pending');
    expect(view.needsAttention).toBe(false);
  });

  it('does not let work in progress make a file analysable', () => {
    expect(describeSource(payload('stale'), 'reindexing').analysable).toBe(false);
    expect(describeSource(payload('current'), 'verifying').analysable).toBe(false);
  });

  it('disables the recovery actions while they are running', () => {
    const view = describeSource(payload('missing'), 'relinking');
    expect(view.label).toBe('Relinking');
    expect(view.canRelink).toBe(false);
  });

  it('distinguishes checking from re-indexing', () => {
    expect(describeSource(payload('unverified'), 'verifying').label).toBe('Checking');
    expect(describeSource(payload('current'), 'reindexing').label).toBe('Re-indexing');
  });
});

describe('describeSources', () => {
  it('applies per-file activity by id', () => {
    const views = describeSources(
      [payload('current', { fastaFileId: 'a' }), payload('stale', { fastaFileId: 'b' })],
      { b: 'reindexing' },
    );
    expect(views[0].activity).toBe('idle');
    expect(views[1].activity).toBe('reindexing');
  });
});

describe('sourcesNeedingAttention', () => {
  it('lists only actionable files, worst first', () => {
    const views = describeSources([
      payload('current', { fastaFileId: 'a' }),
      payload('stale', { fastaFileId: 'b' }),
      payload('missing', { fastaFileId: 'c' }),
      payload('unverified', { fastaFileId: 'd' }),
    ]);
    expect(sourcesNeedingAttention(views).map((view) => view.payload.fastaFileId)).toEqual([
      'c',
      'b',
    ]);
  });
});

describe('sourceBannerMessage', () => {
  it('is silent when every file is fine', () => {
    expect(sourceBannerMessage(describeSources([payload('current')]))).toBeNull();
    expect(sourceBannerMessage(describeSources([payload('unverified')]))).toBeNull();
  });

  it('names a single missing file so the user can act on it', () => {
    const message = sourceBannerMessage(
      describeSources([payload('missing', { displayName: 'beetles.fasta' })]),
    );
    expect(message).toContain('beetles.fasta');
  });

  it('counts rather than lists when several are missing', () => {
    const message = sourceBannerMessage(
      describeSources([
        payload('missing', { fastaFileId: 'a' }),
        payload('missing', { fastaFileId: 'b' }),
      ]),
    );
    expect(message).toContain('2 linked files are missing');
  });

  it('reports missing and changed files together', () => {
    const message = sourceBannerMessage(
      describeSources([
        payload('missing', { fastaFileId: 'a', displayName: 'gone.fasta' }),
        payload('stale', { fastaFileId: 'b', displayName: 'edited.fasta' }),
      ]),
    );
    expect(message).toContain('gone.fasta');
    expect(message).toContain('edited.fasta');
  });
});

describe('scopeIsAnalysable', () => {
  const views = describeSources([
    payload('current', { fastaFileId: 'ok' }),
    payload('stale', { fastaFileId: 'changed' }),
    payload('missing', { fastaFileId: 'gone' }),
  ]);

  it('allows a scope of healthy files', () => {
    expect(scopeIsAnalysable(views, ['ok'])).toBe(true);
  });

  it('refuses a scope containing any unusable file', () => {
    expect(scopeIsAnalysable(views, ['ok', 'changed'])).toBe(false);
    expect(scopeIsAnalysable(views, ['ok', 'gone'])).toBe(false);
  });

  it('refuses an empty scope and unknown ids', () => {
    expect(scopeIsAnalysable(views, [])).toBe(false);
    expect(scopeIsAnalysable(views, ['not-in-project'])).toBe(false);
  });
});
