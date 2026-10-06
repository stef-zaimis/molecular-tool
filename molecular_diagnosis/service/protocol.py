"""
Wire format for the Electron <-> Python boundary.

One JSON object per line, in both directions:

    request   {"id": "7", "method": "runMolecularDiagnosis", "params": {...}}
    success   {"id": "7", "ok": true,  "result": {...}}
    failure   {"id": "7", "ok": false, "error": {"code": "...", "message": "..."}}
    progress  {"id": "7", "type": "progress", "progress": {...}}

A PROGRESS line may appear any number of times between a request and its
response, and never in place of one. It carries no `ok`, so a reader that only
knows about requests and responses cannot mistake it for either: the framing is
unchanged, one JSON object per line, and the pending request stays pending
until its `ok` arrives.

Everything crossing this boundary is plain JSON-serializable data. No Python
objects, no pickling, no knowledge of Electron on this side and no knowledge of
the scientific modules on the other.

Newline-delimited JSON is deliberate: it needs no framing library, is trivially
readable in a log, and lets the parent read responses incrementally.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from molecular_diagnosis.service.errors import ErrorCode, ServiceError

__all__ = [
    "Request",
    "decode_request",
    "encode_failure",
    "encode_progress",
    "encode_success",
    "resume_state_from_payload",
    "resume_state_to_payload",
]


@dataclass(frozen=True)
class Request:
    id: str
    method: str
    params: dict[str, Any]


def decode_request(line: str) -> Request:
    """Parse one request line. Raises `ServiceError` on anything malformed."""
    try:
        payload = json.loads(line)
    except json.JSONDecodeError as error:
        raise ServiceError(
            ErrorCode.INVALID_REQUEST,
            "The backend received a malformed request.",
            detail=str(error),
            cause=error,
        ) from error

    if not isinstance(payload, dict):
        raise ServiceError(
            ErrorCode.INVALID_REQUEST,
            "The backend received a request that was not an object.",
        )

    request_id = payload.get("id")
    method = payload.get("method")
    params = payload.get("params", {})

    if not isinstance(request_id, str) or not request_id:
        raise ServiceError(ErrorCode.INVALID_REQUEST, "The request has no id.")

    if not isinstance(method, str) or not method:
        raise ServiceError(ErrorCode.INVALID_REQUEST, "The request has no method.")

    if params is None:
        params = {}

    if not isinstance(params, dict):
        raise ServiceError(
            ErrorCode.INVALID_REQUEST,
            "The request parameters were not an object.",
        )

    return Request(id=request_id, method=method, params=params)


def encode_success(request_id: str, result: Any) -> str:
    return json.dumps({"id": request_id, "ok": True, "result": result}, ensure_ascii=False)


def encode_progress(request_id: str, progress: dict[str, Any]) -> str:
    """
    One progress notification for a request that is still running.

    `type` marks it and there is no `ok` field, so this can never be read as a
    response. The payload is plain numbers and stage names; the wording of what
    the user sees belongs to the renderer.
    """
    return json.dumps(
        {"id": request_id, "type": "progress", "progress": progress},
        ensure_ascii=False,
        default=str,
    )


def encode_failure(request_id: str, error: ServiceError) -> str:
    return json.dumps(
        {"id": request_id, "ok": False, "error": error.to_payload()},
        ensure_ascii=False,
    )


# ---------------------------------------------------------------------------
# Continuation state
# ---------------------------------------------------------------------------

# The DMC search can stop because it hit the configured maximum combination
# length rather than because it found anything. The caller may then continue
# from the next length. Carrying the prior search state forward is what stops
# the continuation from re-testing combinations that were already ruled out.
#
# The state is passed through the boundary rather than held in this process, so
# a continuation survives a backend restart and there is no hidden session.
#
# JSON has no tuples and no integer object keys, so the conversion is explicit
# in both directions.


def resume_state_to_payload(
    *,
    stopped_at_length: int,
    diagnostic_combinations: Sequence[Sequence[int]],
    combinations_tested_by_length: dict[int, int],
) -> dict[str, Any]:
    """Build the blob the caller sends back to continue this search."""
    return {
        "startCombinationLength": stopped_at_length + 1,
        "diagnosticCombinations": [list(combo) for combo in diagnostic_combinations],
        "combinationsTestedByLength": {
            str(length): count for length, count in combinations_tested_by_length.items()
        },
    }


def resume_state_from_payload(
    payload: Any,
) -> tuple[int, list[tuple[int, ...]], dict[int, int]] | None:
    """
    Read a continuation blob back into the shapes `find_dmc_information` wants.

    Returns None when no continuation was supplied. Raises `ServiceError` when
    one was supplied but is not usable, rather than silently restarting the
    search from scratch.
    """
    if payload is None:
        return None

    if not isinstance(payload, dict):
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER,
            "The continuation state was not an object.",
        )

    start = payload.get("startCombinationLength")
    if not isinstance(start, int) or isinstance(start, bool) or start < 1:
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER,
            "The continuation state has no usable start length.",
            detail=f"startCombinationLength={start!r}",
        )

    raw_combos = payload.get("diagnosticCombinations", [])
    if not isinstance(raw_combos, list):
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER,
            "The continuation state has an unusable combination list.",
        )

    combinations: list[tuple[int, ...]] = []
    for combo in raw_combos:
        if not isinstance(combo, list) or not all(
            isinstance(site, int) and not isinstance(site, bool) for site in combo
        ):
            raise ServiceError(
                ErrorCode.INVALID_PARAMETER,
                "The continuation state contains a malformed combination.",
                detail=repr(combo),
            )
        combinations.append(tuple(combo))

    raw_tested = payload.get("combinationsTestedByLength", {})
    if not isinstance(raw_tested, dict):
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER,
            "The continuation state has an unusable tested-count map.",
        )

    tested: dict[int, int] = {}
    for length, count in raw_tested.items():
        try:
            tested[int(length)] = int(count)
        except (TypeError, ValueError) as error:
            raise ServiceError(
                ErrorCode.INVALID_PARAMETER,
                "The continuation state has an unusable tested-count entry.",
                detail=f"{length!r}: {count!r}",
                cause=error,
            ) from error

    return start, combinations, tested
