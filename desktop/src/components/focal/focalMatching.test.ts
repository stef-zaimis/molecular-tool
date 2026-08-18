import { describe, expect, it } from 'vitest';
import {
  countMatchingHeaders,
  evaluateFocalTokens,
  headerMatchesToken,
  matchingHeaders,
  serializeFocalTokens,
  tokenizeFocalInput,
  unionMatchCount,
} from './focalMatching';
import { MOCK_FASTA_HEADERS } from '../../fixtures/mockHeaders';

describe('tokenizeFocalInput', () => {
  it('splits on the token separator and trims each token', () => {
    expect(tokenizeFocalInput('Leptacis_tipulae; Platygaster')).toEqual([
      'Leptacis_tipulae',
      'Platygaster',
    ]);
  });

  it('returns a single token when there is no separator', () => {
    expect(tokenizeFocalInput('Leptacis_tipulae')).toEqual(['Leptacis_tipulae']);
  });

  it('drops empty tokens from trailing, leading and doubled separators', () => {
    expect(tokenizeFocalInput(';Leptacis;;Platygaster;')).toEqual(['Leptacis', 'Platygaster']);
  });

  it('drops whitespace-only tokens', () => {
    // An empty token would substring-match every header, so it must never survive.
    expect(tokenizeFocalInput('   ;  \t ;Leptacis')).toEqual(['Leptacis']);
  });

  it('returns an empty list for empty input', () => {
    expect(tokenizeFocalInput('')).toEqual([]);
    expect(tokenizeFocalInput('   ')).toEqual([]);
    expect(tokenizeFocalInput(';;;')).toEqual([]);
  });

  it('preserves internal spaces and pipe characters inside a token', () => {
    expect(tokenizeFocalInput('AACTA5253-20|AU|Leptacis')).toEqual(['AACTA5253-20|AU|Leptacis']);
  });

  it('round-trips through serializeFocalTokens', () => {
    const tokens = ['Leptacis_tipulae', 'Platygaster'];
    expect(tokenizeFocalInput(serializeFocalTokens(tokens))).toEqual(tokens);
  });
});

describe('headerMatchesToken — substring semantics', () => {
  it('matches on substring containment, not exact identity', () => {
    // This is the behaviour of the current Python code: `target_string in header`.
    expect(headerMatchesToken('BBHYQ1043-18|GB|Leptacis_tipulae', 'Leptacis_tipulae')).toBe(true);
    expect(headerMatchesToken('BBHYQ1043-18|GB|Leptacis_tipulae', 'tipulae')).toBe(true);
    expect(headerMatchesToken('BBHYQ1043-18|GB|Leptacis_tipulae', 'GB')).toBe(true);
  });

  it('is case sensitive by default', () => {
    expect(headerMatchesToken('AACTA5253-20|AU|Leptacis', 'leptacis')).toBe(false);
    expect(headerMatchesToken('AACTA5253-20|AU|Leptacis', 'Leptacis')).toBe(true);
  });

  it('never matches on an empty token', () => {
    expect(headerMatchesToken('AACTA5253-20|AU|Leptacis', '')).toBe(false);
  });

  it('supports the alternative modes without affecting the default', () => {
    const header = 'AACTA5253-20|AU|Leptacis';
    expect(headerMatchesToken(header, 'leptacis', 'caseInsensitiveSubstring')).toBe(true);
    expect(headerMatchesToken(header, 'Leptacis', 'exactId')).toBe(false);
    expect(headerMatchesToken(header, header, 'exactId')).toBe(true);
    expect(headerMatchesToken(header, 'Leptacis', 'fieldAware')).toBe(true);
    expect(headerMatchesToken(header, 'AU', 'fieldAware')).toBe(true);
    expect(headerMatchesToken(header, 'Lept', 'fieldAware')).toBe(false);
  });
});

describe('countMatchingHeaders / matchingHeaders', () => {
  it('counts every matching header for a positive match', () => {
    expect(countMatchingHeaders(MOCK_FASTA_HEADERS, 'Leptacis_tipulae')).toBe(5);
  });

  it('counts zero for a token that matches nothing', () => {
    expect(countMatchingHeaders(MOCK_FASTA_HEADERS, 'Lptacis')).toBe(0);
  });

  it('counts broad substrings across genus and species headers', () => {
    // "Leptacis" is a substring of "Leptacis_tipulae" too, so this is a superset.
    expect(countMatchingHeaders(MOCK_FASTA_HEADERS, 'Leptacis')).toBe(14);
  });

  it('returns the matching headers themselves', () => {
    expect(matchingHeaders(MOCK_FASTA_HEADERS, 'Synopeas')).toEqual([
      'ELVWX1057-24|PL|Synopeas',
      'EMYZA2168-24|CZ|Synopeas',
    ]);
  });
});

describe('evaluateFocalTokens', () => {
  it('marks a token that matches at least one header as a match', () => {
    const [token] = evaluateFocalTokens(['Leptacis_tipulae'], MOCK_FASTA_HEADERS);
    expect(token.state).toBe('match');
    expect(token.matchCount).toBe(5);
  });

  it('marks a token that matches zero headers as noMatch', () => {
    const [token] = evaluateFocalTokens(['Lptacis'], MOCK_FASTA_HEADERS);
    expect(token.state).toBe('noMatch');
    expect(token.matchCount).toBe(0);
  });

  it('marks every token neutral when there are no headers to validate against', () => {
    const tokens = evaluateFocalTokens(['Leptacis_tipulae', 'Lptacis'], null);
    expect(tokens.map((token) => token.state)).toEqual(['neutral', 'neutral']);
    expect(tokens.every((token) => token.matchCount === 0)).toBe(true);
  });

  it('resolves mixed states independently across multiple tokens', () => {
    // Mirrors the reference mockup: one green token, one red token.
    const tokens = evaluateFocalTokens(
      ['AACTA5253-20|AU|Leptacis', 'ALGFA5279-26|FR|Lptacis', 'Synopeas'],
      MOCK_FASTA_HEADERS,
    );
    expect(tokens.map((token) => token.state)).toEqual(['match', 'noMatch', 'match']);
    expect(tokens.map((token) => token.matchCount)).toEqual([1, 0, 2]);
  });

  it('returns an empty result for no tokens', () => {
    expect(evaluateFocalTokens([], MOCK_FASTA_HEADERS)).toEqual([]);
  });

  it('is case sensitive, so a case-mismatched token reads as noMatch', () => {
    const [token] = evaluateFocalTokens(['leptacis_tipulae'], MOCK_FASTA_HEADERS);
    expect(token.state).toBe('noMatch');
  });
});

describe('unionMatchCount', () => {
  it('counts each header once even when several tokens match it', () => {
    // "Leptacis" and "Leptacis_tipulae" overlap on the 5 tipulae headers.
    expect(unionMatchCount(['Leptacis', 'Leptacis_tipulae'], MOCK_FASTA_HEADERS)).toBe(14);
  });

  it('sums disjoint token matches', () => {
    expect(unionMatchCount(['Synopeas', 'Platygaster'], MOCK_FASTA_HEADERS)).toBe(4);
  });

  it('is zero when no tokens are supplied', () => {
    expect(unionMatchCount([], MOCK_FASTA_HEADERS)).toBe(0);
  });
});
