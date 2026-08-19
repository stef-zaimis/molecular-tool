import type { FocalMatchMode } from '../../contract';

/**
 * Focal-string tokenisation and validation.
 *
 * This module is the ONLY place that decides what "a token matches a header"
 * means. FocalStringEditor renders state; it never compares strings itself.
 * That separation is deliberate: the matching rule is expected to change
 * (case-insensitive, exact id, field-aware, regex — see contract.ts
 * `FocalMatchMode`) and the editor must not need rewriting when it does.
 *
 * The default mode reproduces the current Python behaviour exactly:
 * case-sensitive substring containment, `target_string in header`
 * (REPO_MAP.md §4.2). It is NOT exact-id equality.
 */

export type FocalTokenState = 'match' | 'noMatch' | 'neutral';

export interface FocalToken {
  readonly text: string;
  readonly state: FocalTokenState;
  /** Number of known headers matched. 0 when unmatched or unvalidatable. */
  readonly matchCount: number;
}

/**
 * Focal entries are DISPLAYED separated by ';', mirroring the design, but the
 * separator is presentation only — the focal set is stored as an array of
 * strings and never round-tripped through a delimited string. That is why a
 * ';' inside a single entry is rejected outright rather than silently split.
 */
export const TOKEN_SEPARATOR = ';';

/** Why a candidate focal string cannot be added, or null when it is fine. */
export type FocalStringProblem =
  | { readonly kind: 'empty' }
  | { readonly kind: 'duplicate' };

export function describeFocalStringProblem(problem: FocalStringProblem): string {
  switch (problem.kind) {
    case 'empty':
      return 'Enter a search string first.';
    case 'duplicate':
      return 'That string is already in the focal set.';
  }
}

/**
 * Validates one candidate entry against the current set.
 *
 * `existing` is only consulted for the duplicate check, so this same function
 * serves both the add and remove paths (remove passes an empty list).
 */
export function validateFocalString(
  value: string,
  existing: readonly string[] = [],
): FocalStringProblem | null {
  const trimmed = value.trim();
  if (trimmed.length === 0) return { kind: 'empty' };
  // A ';' is allowed inside an entry: the separator shown between entries is a
  // rendered element, never a stored character, so there is nothing to confuse.
  if (existing.includes(trimmed)) return { kind: 'duplicate' };
  return null;
}

/**
 * Split raw editor text into trimmed, non-empty tokens.
 * Empty segments (from a trailing ';', a doubled ';;', or whitespace-only
 * input) are dropped rather than becoming empty tokens — an empty token would
 * match every header under substring semantics, which is never intended.
 */
export function tokenizeFocalInput(raw: string): string[] {
  return raw
    .split(TOKEN_SEPARATOR)
    .map((part) => part.trim())
    .filter((part) => part.length > 0);
}

/** Render tokens back to editor text. */
export function serializeFocalTokens(tokens: readonly string[]): string {
  return tokens.join(`${TOKEN_SEPARATOR} `);
}

/**
 * Does one header match one token?
 *
 * Only `caseSensitiveSubstring` reflects real backend behaviour. The other
 * modes are frontend-only conveniences and are marked as backend gaps in
 * contract.ts — do not assume the Python side honours them.
 */
export function headerMatchesToken(
  header: string,
  token: string,
  mode: FocalMatchMode = 'caseSensitiveSubstring',
): boolean {
  if (token.length === 0) return false;

  switch (mode) {
    case 'caseSensitiveSubstring':
      return header.includes(token);
    case 'caseInsensitiveSubstring':
      return header.toLowerCase().includes(token.toLowerCase());
    case 'exactId':
      return header === token;
    case 'fieldAware':
      return header.split('|').some((field) => field === token);
    case 'regex':
      try {
        return new RegExp(token).test(header);
      } catch {
        // An in-progress regex is not an error state for the user.
        return false;
      }
    default:
      return header.includes(token);
  }
}

export function countMatchingHeaders(
  headers: readonly string[],
  token: string,
  mode: FocalMatchMode = 'caseSensitiveSubstring',
): number {
  if (token.length === 0) return 0;
  return headers.reduce(
    (total, header) => (headerMatchesToken(header, token, mode) ? total + 1 : total),
    0,
  );
}

export function matchingHeaders(
  headers: readonly string[],
  token: string,
  mode: FocalMatchMode = 'caseSensitiveSubstring',
): string[] {
  if (token.length === 0) return [];
  return headers.filter((header) => headerMatchesToken(header, token, mode));
}

/**
 * Evaluate every token against the known headers.
 *
 * `headers === null` means "no alignment is loaded, so nothing can be
 * validated" and yields NEUTRAL for every token — distinct from "validated and
 * found nothing", which is NO MATCH. Keeping those apart matters: showing red
 * for an unloaded file would tell the user their string is wrong when it may
 * be fine.
 */
export function evaluateFocalTokens(
  tokens: readonly string[],
  headers: readonly string[] | null,
  mode: FocalMatchMode = 'caseSensitiveSubstring',
): FocalToken[] {
  return tokens.map((text) => {
    if (headers === null) return { text, state: 'neutral' as const, matchCount: 0 };

    const matchCount = countMatchingHeaders(headers, text, mode);
    return {
      text,
      state: matchCount > 0 ? ('match' as const) : ('noMatch' as const),
      matchCount,
    };
  });
}

/** Headers matched by at least one token — the prospective focal set. */
export function unionMatchCount(
  tokens: readonly string[],
  headers: readonly string[],
  mode: FocalMatchMode = 'caseSensitiveSubstring',
): number {
  if (tokens.length === 0) return 0;
  return headers.reduce(
    (total, header) =>
      tokens.some((token) => headerMatchesToken(header, token, mode)) ? total + 1 : total,
    0,
  );
}
