"""
Service operations.

These are the only functions the desktop frontend can reach. They take plain
dictionaries and return plain dictionaries; nothing here imports Electron, and
nothing here imports `gui.py` or any Tkinter machinery.

This module owns the orchestration that used to live inside the Tkinter run
button — notably the DMC search-continuation loop — but as a pure function of
its inputs, with continuation state travelling in and out as data.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from molecular_diagnosis.fasta_io import parse_fasta, validate_aligned_fasta
from molecular_diagnosis.focal import (
    header_matches_focal,
    matching_focal_strings,
    normalise_focal_strings,
)
from molecular_diagnosis.pipeline import run_pipeline_core, run_punishment_core
from molecular_diagnosis.service.errors import ErrorCode, ServiceError, to_service_error
from molecular_diagnosis.service.protocol import (
    resume_state_from_payload,
    resume_state_to_payload,
)

__all__ = ["METHODS", "dispatch"]

#: Sample of matching headers returned per focal string. A broad selector can
#: match thousands of headers; the UI only needs enough to show the user what
#: it hit.
HEADER_SAMPLE_LIMIT = 20


# ---------------------------------------------------------------------------
# Parameter helpers
# ---------------------------------------------------------------------------


def _require_str(params: dict[str, Any], key: str) -> str:
    value = params.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER,
            f"The request is missing '{key}'.",
            detail=f"{key}={value!r}",
        )
    return value


def _require_focal_strings(params: dict[str, Any], key: str = "focalStrings") -> list[str]:
    value = params.get(key)

    if isinstance(value, str):
        raw: Sequence[str] = [value]
    elif isinstance(value, list) and all(isinstance(item, str) for item in value):
        raw = value
    else:
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER,
            "The focal strings must be a list of strings.",
            detail=f"{key}={value!r}",
        )

    # normalise_focal_strings enforces "no empty selector" and deduplicates.
    return normalise_focal_strings(raw)


def _optional_bool(params: dict[str, Any], key: str, default: bool) -> bool:
    value = params.get(key, default)
    if not isinstance(value, bool):
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER,
            f"'{key}' must be true or false.",
            detail=f"{key}={value!r}",
        )
    return value


def _optional_int(params: dict[str, Any], key: str, default: int) -> int:
    value = params.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int):
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER,
            f"'{key}' must be a whole number.",
            detail=f"{key}={value!r}",
        )
    return value


def _count_header_lines(path: Path) -> int:
    """
    How many '>' lines the file has.

    `parse_fasta` returns a dict, so repeated headers collapse into one entry
    and the sequence count silently drops. Reporting both numbers lets the UI
    say so instead of the user wondering where sequences went.
    """
    count = 0
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


# ---------------------------------------------------------------------------
# loadFasta
# ---------------------------------------------------------------------------


def load_fasta(params: dict[str, Any]) -> dict[str, Any]:
    """
    Parse and validate a FASTA through the existing backend, and return the
    metadata the frontend needs.

    Validation is `fasta_io.validate_aligned_fasta` — the same check the
    analysis pipeline performs — so a file that loads here is a file the
    analysis will accept.
    """
    path_text = _require_str(params, "path")
    path = Path(path_text)

    if not path.exists():
        raise ServiceError(
            ErrorCode.FASTA_NOT_FOUND,
            "The FASTA file could not be found.",
            detail=path_text,
        )
    if not path.is_file():
        raise ServiceError(
            ErrorCode.FASTA_NOT_A_FILE,
            "The FASTA path is not a file.",
            detail=path_text,
        )

    sequences = parse_fasta(path)
    alignment_length = validate_aligned_fasta(sequences)

    headers = list(sequences)
    header_line_count = _count_header_lines(path)

    alphabet = sorted({residue for sequence in sequences.values() for residue in sequence})

    return {
        "path": str(path),
        "sequenceCount": len(sequences),
        "alignmentLength": alignment_length,
        "headers": headers,
        # Repeated headers collapse in the parser; report the discrepancy.
        "headerLineCount": header_line_count,
        "duplicateHeaderCount": max(0, header_line_count - len(sequences)),
        "alphabet": alphabet,
    }


# ---------------------------------------------------------------------------
# validateFocalStrings
# ---------------------------------------------------------------------------


def validate_focal_strings(params: dict[str, Any]) -> dict[str, Any]:
    """
    Report, per focal string, whether it matches any header of the given FASTA.

    This is the authority for the green/red state in the UI. It runs the exact
    matcher the analysis uses (`focal.header_matches_focal`), so the colours
    cannot disagree with which sequences the analysis will actually select.
    """
    path_text = _require_str(params, "path")
    path = Path(path_text)

    if not path.is_file():
        raise ServiceError(
            ErrorCode.FASTA_NOT_FOUND,
            "The FASTA file could not be found.",
            detail=path_text,
        )

    raw = params.get("focalStrings")
    if not isinstance(raw, list) or not all(isinstance(item, str) for item in raw):
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER,
            "The focal strings must be a list of strings.",
            detail=repr(raw),
        )

    sequences = parse_fasta(path)
    headers = list(sequences)

    # Deliberately NOT normalised: the UI asks about the entries it is holding,
    # including any that are empty or duplicated, and needs an answer for each
    # one in the order it supplied them.
    results = []
    for selector in raw:
        if not selector.strip():
            results.append(
                {
                    "focalString": selector,
                    "matched": False,
                    "matchCount": 0,
                    "sampleHeaders": [],
                    "truncated": False,
                    "invalid": True,
                }
            )
            continue

        matches = [header for header in headers if selector in header]
        results.append(
            {
                "focalString": selector,
                "matched": len(matches) > 0,
                "matchCount": len(matches),
                "sampleHeaders": matches[:HEADER_SAMPLE_LIMIT],
                "truncated": len(matches) > HEADER_SAMPLE_LIMIT,
                "invalid": False,
            }
        )

    usable = [selector for selector in raw if selector.strip()]
    union = [header for header in headers if header_matches_focal(header, usable)]

    return {
        "path": str(path),
        "totalHeaders": len(headers),
        "results": results,
        # The prospective focal set: headers matched by AT LEAST ONE selector.
        "unionMatchCount": len(union),
        "unionSampleHeaders": union[:HEADER_SAMPLE_LIMIT],
        "matchedBy": {
            header: matching_focal_strings(header, usable)
            for header in union[:HEADER_SAMPLE_LIMIT]
        },
    }


# ---------------------------------------------------------------------------
# runMolecularDiagnosis
# ---------------------------------------------------------------------------


def _dmc_to_payload(dmc: Any) -> dict[str, Any]:
    """Flatten a DMCResult into plain JSON-safe data."""
    return {
        "combinationsByLength": {
            str(length): [list(combo) for combo in combos]
            for length, combos in dmc.diagnostic_combinations_by_length.items()
        },
        "singleSites": list(dmc.single),
        "pairs": [list(pair) for pair in dmc.pairs],
        "uniqueSites": list(dmc.unique),
        # Site -> consensus state, 0-based keys as strings for JSON.
        "states": {str(site): state for site, state in dmc.states.items()},
        "stopReason": dmc.stop_reason,
        "stoppedAtLength": dmc.stopped_at_length,
        "minCombinationLength": dmc.min_combination_length,
        "maxCombinationLength": dmc.max_combination_length,
        "startCombinationLength": dmc.start_combination_length,
        "diagnostics": {
            "fixedCount": dmc.fixed_count,
            "skippedSites": dmc.skipped_non_acgt,
            "globallyConservedRemoved": dmc.globally_conserved_removed,
            "candidateCount": dmc.candidate_count,
            "pairsTested": dmc.pairs_tested,
            "totalCombinationsTested": dmc.total_combinations_tested,
            "combinationsTestedByLength": {
                str(length): count
                for length, count in dmc.combinations_tested_by_length.items()
            },
            "ambiguousBdSitesIncluded": list(dmc.ambiguous_bd_sites_included),
            "gappyConsensusSitesIncluded": list(dmc.gappy_consensus_sites_included),
        },
    }


def run_molecular_diagnosis(params: dict[str, Any]) -> dict[str, Any]:
    """
    Run one pass of the Molecular Diagnosis pipeline.

    Continuation: when the search stops only because it reached the configured
    maximum combination length, `canContinue` is true and `resume` carries the
    state needed to pick up from the next length. Passing that blob back as
    `resume` on a later call continues the search instead of restarting it —
    the already-found combinations and already-counted tests are threaded into
    `find_dmc_information` exactly as the Tkinter loop did.

    The decision to continue belongs to the caller; this function never prompts.
    """
    fasta_path = _require_str(params, "fastaPath")
    focal_strings = _require_focal_strings(params)
    output_directory = _require_str(params, "outputDirectory")

    include_ambiguous = _optional_bool(params, "giveBenefitOfDoubtToAmbiguousBases", False)
    ignore_gaps = _optional_bool(params, "ignoreGaps", False)
    min_length = _optional_int(params, "minCandidateSize", 1)
    max_length = _optional_int(params, "maxCandidateSize", 2)

    resume = resume_state_from_payload(params.get("resume"))

    if resume is None:
        start_length = 1
        initial_combinations: list[tuple[int, ...]] | None = None
        initial_tested: dict[int, int] | None = None
    else:
        start_length, initial_combinations, initial_tested = resume
        # A continuation must extend the search, never repeat it.
        if max_length < start_length:
            raise ServiceError(
                ErrorCode.INVALID_PARAMETER,
                "The new maximum combination size must be at least the size the search resumes from.",
                detail=f"maxCandidateSize={max_length}, resumes at {start_length}",
            )

    result = run_pipeline_core(
        fasta_path=fasta_path,
        focal_strings=focal_strings,
        output_dir=output_directory,
        include_ambiguous_dmc_bd=include_ambiguous,
        include_gappy_consensus_dmc_sites=ignore_gaps,
        min_combination_length=min_length,
        max_combination_length=max_length,
        start_combination_length=start_length,
        initial_diagnostic_combinations=initial_combinations,
        initial_combinations_tested_by_length=initial_tested,
    )

    dmc = result.dmc
    can_continue = dmc.stop_reason == "reached_maximum_length"

    payload: dict[str, Any] = {
        "focalStrings": focal_strings,
        "outputs": {
            "reportTxt": str(result.txt_output_path),
            "workbookXlsx": str(result.xlsx_output_path),
            "consensusTxt": (
                str(result.consensus_txt_output_path)
                if result.consensus_txt_output_path is not None
                else None
            ),
        },
        "dmc": _dmc_to_payload(dmc),
        "canContinue": can_continue,
        "resume": None,
    }

    if can_continue:
        payload["resume"] = resume_state_to_payload(
            stopped_at_length=dmc.stopped_at_length,
            diagnostic_combinations=dmc.diagnostic_combinations,
            combinations_tested_by_length=dmc.combinations_tested_by_length,
        )

    return payload


# ---------------------------------------------------------------------------
# runSequencePunishment (wired for completeness; no UI reaches it yet)
# ---------------------------------------------------------------------------


def run_sequence_punishment(params: dict[str, Any]) -> dict[str, Any]:
    fasta_path = _require_str(params, "fastaPath")
    focal_strings = _require_focal_strings(params)
    output_directory = _require_str(params, "outputDirectory")

    result = run_punishment_core(
        fasta_path=fasta_path,
        focal_strings=focal_strings,
        output_dir=output_directory,
    )

    return {
        "focalStrings": focal_strings,
        "outputs": {
            "punishmentXlsx": str(result.xlsx_output_path),
            "sequenceSubsetsXlsx": (
                str(result.sequence_subsets_xlsx_output_path)
                if result.sequence_subsets_xlsx_output_path is not None
                else None
            ),
        },
    }


def ping(_params: dict[str, Any]) -> dict[str, Any]:
    """Liveness probe used by the frontend to confirm the backend started."""
    import sys

    return {
        "ok": True,
        "python": sys.version.split()[0],
    }


METHODS = {
    "ping": ping,
    "loadFasta": load_fasta,
    "validateFocalStrings": validate_focal_strings,
    "runMolecularDiagnosis": run_molecular_diagnosis,
    "runSequencePunishment": run_sequence_punishment,
}


def dispatch(method: str, params: dict[str, Any]) -> dict[str, Any]:
    """
    Run one operation.

    This is the boundary: it always raises `ServiceError`, never a bare
    exception from the scientific modules. The original exception is kept as the
    cause so no context is lost on the way out.
    """
    handler = METHODS.get(method)
    if handler is None:
        raise ServiceError(
            ErrorCode.UNKNOWN_METHOD,
            f"The backend does not support '{method}'.",
        )

    try:
        return handler(params)
    except ServiceError:
        raise
    except Exception as error:  # noqa: BLE001 - translated, not swallowed
        raise to_service_error(error) from error
