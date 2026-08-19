import { describe, expect, it } from 'vitest';
import { defaultOutputDirectory } from '../app/state/ProjectContext';
import type {
  MolecularDiagnosisRequest,
  MolecularDiagnosisResult,
} from '../backendContract';

/**
 * Serialization tests for the Electron -> Python boundary.
 *
 * Everything crossing it travels as JSON over a pipe, so anything that does not
 * survive `JSON.stringify` / `JSON.parse` — undefined, Map, Set, tuples with
 * identity — is a bug that would only show up at runtime.
 */

const REQUEST: MolecularDiagnosisRequest = {
  fastaPath: '/data/aligned.fasta',
  focalStrings: ['Leptacis_tipulae', 'Synopeas'],
  outputDirectory: '/data',
  ignoreGaps: true,
  giveBenefitOfDoubtToAmbiguousBases: false,
  minCandidateSize: 1,
  maxCandidateSize: 3,
  resume: null,
};

const RESULT: MolecularDiagnosisResult = {
  focalStrings: ['Leptacis_tipulae'],
  outputs: {
    reportTxt: '/data/DMCs_output.txt',
    workbookXlsx: '/data/comparison_output.xlsx',
    consensusTxt: '/data/focal_consensus_output.txt',
  },
  dmc: {
    combinationsByLength: { '1': [[1], [3]] },
    singleSites: [1, 3],
    pairs: [],
    uniqueSites: [1, 3],
    states: { '1': 'A', '3': 'C' },
    stopReason: 'found_at_or_above_minimum_length',
    stoppedAtLength: 1,
    minCombinationLength: 1,
    maxCombinationLength: 2,
    startCombinationLength: 1,
    diagnostics: {
      fixedCount: 4,
      skippedSites: 0,
      globallyConservedRemoved: 0,
      candidateCount: 4,
      pairsTested: 0,
      totalCombinationsTested: 4,
      combinationsTestedByLength: { '1': 4 },
      ambiguousBdSitesIncluded: [],
      gappyConsensusSitesIncluded: [],
    },
  },
  canContinue: false,
  resume: null,
};

const roundTrip = <T,>(value: T): T => JSON.parse(JSON.stringify(value)) as T;

describe('Molecular Diagnosis request serialization', () => {
  it('survives a JSON round trip unchanged', () => {
    expect(roundTrip(REQUEST)).toEqual(REQUEST);
  });

  it('keeps focal strings as an array, never a joined string', () => {
    const wire = JSON.parse(JSON.stringify(REQUEST)) as { focalStrings: unknown };
    expect(Array.isArray(wire.focalStrings)).toBe(true);
    expect(wire.focalStrings).toEqual(['Leptacis_tipulae', 'Synopeas']);
  });

  it('carries a focal string containing a semicolon without splitting it', () => {
    const request = { ...REQUEST, focalStrings: ['wei;rd', 'plain'] };
    expect(roundTrip(request).focalStrings).toEqual(['wei;rd', 'plain']);
  });

  it('serializes a continuation request', () => {
    const request: MolecularDiagnosisRequest = {
      ...REQUEST,
      maxCandidateSize: 5,
      resume: {
        startCombinationLength: 3,
        diagnosticCombinations: [[1, 3], [0, 2, 5]],
        combinationsTestedByLength: { '1': 4, '2': 6 },
      },
    };

    const revived = roundTrip(request);
    expect(revived.resume?.startCombinationLength).toBe(3);
    expect(revived.resume?.diagnosticCombinations).toEqual([[1, 3], [0, 2, 5]]);
    expect(revived.resume?.combinationsTestedByLength).toEqual({ '1': 4, '2': 6 });
  });
});

describe('Molecular Diagnosis result serialization', () => {
  it('survives a JSON round trip unchanged', () => {
    expect(roundTrip(RESULT)).toEqual(RESULT);
  });

  it('keeps a null consensus path as null rather than dropping the key', () => {
    const result = {
      ...RESULT,
      outputs: { ...RESULT.outputs, consensusTxt: null },
    };
    const revived = roundTrip(result);
    expect('consensusTxt' in revived.outputs).toBe(true);
    expect(revived.outputs.consensusTxt).toBeNull();
  });

  it('round-trips a continuable result together with its resume state', () => {
    const result: MolecularDiagnosisResult = {
      ...RESULT,
      dmc: { ...RESULT.dmc, stopReason: 'reached_maximum_length', stoppedAtLength: 2 },
      canContinue: true,
      resume: {
        startCombinationLength: 3,
        diagnosticCombinations: [],
        combinationsTestedByLength: { '1': 4, '2': 6 },
      },
    };

    const revived = roundTrip(result);
    expect(revived.canContinue).toBe(true);
    expect(revived.resume?.startCombinationLength).toBe(3);
    // The blob is handed straight back to continue, so it must be identical.
    expect(revived.resume).toEqual(result.resume);
  });

  it('keeps integer-keyed maps addressable after the round trip', () => {
    // JSON turns integer keys into strings; the contract already says string.
    const revived = roundTrip(RESULT);
    expect(revived.dmc.combinationsByLength['1']).toEqual([[1], [3]]);
    expect(revived.dmc.states['3']).toBe('C');
  });
});

describe('defaultOutputDirectory', () => {
  it('uses the folder holding the FASTA on Windows paths', () => {
    expect(defaultOutputDirectory('C:\\data\\runs\\aligned.fasta')).toBe('C:\\data\\runs');
  });

  it('uses the folder holding the FASTA on POSIX paths', () => {
    expect(defaultOutputDirectory('/home/me/data/aligned.fasta')).toBe('/home/me/data');
  });

  it('falls back to the input when there is no directory part', () => {
    expect(defaultOutputDirectory('aligned.fasta')).toBe('aligned.fasta');
  });
});
