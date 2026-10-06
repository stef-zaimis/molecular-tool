import { describe, expect, it } from 'vitest';
import { parseFocalText, scanFocalText, serialiseFocalText, textMatchesHeaders } from './focalText';

/**
 * The rules the editable field obeys, mirroring `parse_focal_text` in
 * `project/service.py`.
 */

describe('parsing the focal document', () => {
  it('splits on semicolons and trims around them', () => {
    expect(parseFocalText('a ; b;c')).toEqual(['a', 'b', 'c']);
  });

  it('keeps whitespace INSIDE an entry, which FASTA headers legitimately have', () => {
    expect(parseFocalText('gi|123 Homo sapiens; other')).toEqual([
      'gi|123 Homo sapiens',
      'other',
    ]);
  });

  it('drops blank segments rather than storing empty entries', () => {
    expect(parseFocalText(' ; ;;a;; ')).toEqual(['a']);
    expect(parseFocalText('')).toEqual([]);
    expect(parseFocalText('   ')).toEqual([]);
  });

  it('collapses duplicates onto the first occurrence', () => {
    expect(parseFocalText('a;b;a')).toEqual(['a', 'b']);
  });

  it('leaves a trailing separator alone, so typing the next entry is uninterrupted', () => {
    expect(parseFocalText('a; ')).toEqual(['a']);
  });
});

describe('scanning for decoration offsets', () => {
  it('reports the offsets of the trimmed text, not of the raw segment', () => {
    const spans = scanFocalText('  abc ; def');
    expect(spans[0]).toMatchObject({ text: 'abc', from: 2, to: 5, isMember: true });
    expect(spans[1]).toMatchObject({ text: 'def', from: 8, to: 11, isMember: true });
  });

  it('marks a repeat as present but not a member', () => {
    const spans = scanFocalText('a;a');
    expect(spans[0].isMember).toBe(true);
    // Shown — the user typed it — but it is not a second member.
    expect(spans[1]).toMatchObject({ text: 'a', isMember: false });
  });

  it('offsets a decorator can use land on the right substring', () => {
    const text = 'first ; second';
    for (const span of scanFocalText(text)) {
      expect(text.slice(span.from, span.to)).toBe(span.text);
    }
  });
});

describe('serialising back', () => {
  it('uses the design’s "a; b" form', () => {
    expect(serialiseFocalText(['a', 'b'])).toBe('a; b');
  });

  it('round-trips', () => {
    const headers = ['focal_1|AU|Target', 'gi|9 Homo sapiens'];
    expect(parseFocalText(serialiseFocalText(headers))).toEqual(headers);
  });
});

describe('deciding whether the document needs replacing', () => {
  it('accepts a document that already expresses the membership', () => {
    // Including one still being typed around: replacing here would move the caret.
    expect(textMatchesHeaders('a; b; ', ['a', 'b'])).toBe(true);
    expect(textMatchesHeaders('a;b', ['a', 'b'])).toBe(true);
  });

  it('rejects a different membership or a different order', () => {
    expect(textMatchesHeaders('a; b', ['b', 'a'])).toBe(false);
    expect(textMatchesHeaders('a', ['a', 'b'])).toBe(false);
  });
});

describe('canonicalising at a boundary', () => {
  /*
   * The promise the editor makes: between boundaries the document is left
   * exactly as typed, and at a boundary it is squared up with the stored
   * membership — so a struck-out duplicate cannot survive a Save for an entry
   * the database does not hold.
   */
  it('a duplicate is shown while typing but is not a member', () => {
    const typed = 'a; b; a';
    const spans = scanFocalText(typed);

    expect(spans.map((span) => span.text)).toEqual(['a', 'b', 'a']);
    expect(spans.map((span) => span.isMember)).toEqual([true, true, false]);
    // The membership the caller stores has no duplicate in it.
    expect(parseFocalText(typed)).toEqual(['a', 'b']);
  });

  it('the document is left alone while the membership still matches', () => {
    // Mid-word, trailing separator, a duplicate: none of these change what is
    // stored, so none of them justify rewriting the text under the caret.
    expect(textMatchesHeaders('a; b; a', ['a', 'b'])).toBe(true);
    expect(textMatchesHeaders('a; b; ', ['a', 'b'])).toBe(true);
  });

  it('canonical text for a saved set carries no duplicate at all', () => {
    const saved = parseFocalText('a; b; a');
    const canonical = serialiseFocalText(saved);

    expect(canonical).toBe('a; b');
    expect(scanFocalText(canonical).every((span) => span.isMember)).toBe(true);
  });
});
