"""
Structured errors for the service boundary.

The scientific modules raise bare `ValueError`s carrying prose messages. That is
fine inside Python but useless to a UI, which needs to branch on *what went
wrong*, not on message text. This module is the single translation point: prose
in, a stable code out, with the original message and traceback preserved so
nothing is lost for logging.

Nothing here changes when or why the science raises — only how the failure is
described at the boundary.
"""

from __future__ import annotations

import traceback

__all__ = ["ErrorCode", "ServiceError", "to_service_error"]


class ErrorCode:
    """Stable codes the frontend may branch on."""

    INVALID_REQUEST = "INVALID_REQUEST"
    UNKNOWN_METHOD = "UNKNOWN_METHOD"

    FASTA_NOT_FOUND = "FASTA_NOT_FOUND"
    FASTA_NOT_A_FILE = "FASTA_NOT_A_FILE"
    FASTA_UNREADABLE = "FASTA_UNREADABLE"
    FASTA_EMPTY = "FASTA_EMPTY"
    FASTA_NOT_ALIGNED = "FASTA_NOT_ALIGNED"

    FOCAL_EMPTY = "FOCAL_EMPTY"
    FOCAL_NO_MATCH = "FOCAL_NO_MATCH"
    FOCAL_MATCHES_EVERYTHING = "FOCAL_MATCHES_EVERYTHING"

    OUTPUT_DIR_MISSING = "OUTPUT_DIR_MISSING"
    OUTPUT_DIR_NOT_A_DIRECTORY = "OUTPUT_DIR_NOT_A_DIRECTORY"
    OUTPUT_WRITE_FAILED = "OUTPUT_WRITE_FAILED"

    INVALID_PARAMETER = "INVALID_PARAMETER"
    UNKNOWN = "UNKNOWN"


class ServiceError(Exception):
    """
    An error with a code the frontend can act on.

    `message` is safe to show a user. `detail` and `traceback_text` carry the
    technical context for logs and are never required for the UI to behave
    correctly.
    """

    def __init__(
        self,
        code: str,
        message: str,
        *,
        detail: str | None = None,
        cause: BaseException | None = None,
    ) -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.detail = detail
        self.cause = cause
        self.traceback_text = (
            "".join(traceback.format_exception(type(cause), cause, cause.__traceback__))
            if cause is not None
            else None
        )

    def to_payload(self) -> dict[str, object]:
        payload: dict[str, object] = {"code": self.code, "message": self.message}
        if self.detail:
            payload["detail"] = self.detail
        if self.traceback_text:
            payload["traceback"] = self.traceback_text
        return payload


# Prose emitted by the scientific modules -> stable code. Substring matched
# because the messages are written for humans and may gain punctuation.
_VALUE_ERROR_CODES: tuple[tuple[str, str, str], ...] = (
    (
        "No sequences were read",
        ErrorCode.FASTA_EMPTY,
        "The FASTA file contains no sequences.",
    ),
    (
        "not all the same length",
        ErrorCode.FASTA_NOT_ALIGNED,
        "The sequences are not all the same length, so this is not an aligned FASTA.",
    ),
    (
        "No sequences matched the identifier",
        ErrorCode.FOCAL_NO_MATCH,
        "No FASTA header contains any of the focal strings.",
    ),
    (
        "All sequences match the identifier",
        ErrorCode.FOCAL_MATCHES_EVERYTHING,
        "Every FASTA header matches the focal strings, so there is no contrast set to compare against.",
    ),
    (
        "Focal strings cannot be empty",
        ErrorCode.FOCAL_EMPTY,
        "A focal string is empty.",
    ),
    (
        "No identifier string entered",
        ErrorCode.FOCAL_EMPTY,
        "No focal strings were supplied.",
    ),
    (
        "FASTA file does not exist",
        ErrorCode.FASTA_NOT_FOUND,
        "The FASTA file could not be found.",
    ),
    (
        "FASTA path is not a file",
        ErrorCode.FASTA_NOT_A_FILE,
        "The FASTA path is not a file.",
    ),
    (
        "No FASTA file selected",
        ErrorCode.FASTA_NOT_FOUND,
        "No FASTA file was supplied.",
    ),
    (
        "Output directory does not exist",
        ErrorCode.OUTPUT_DIR_MISSING,
        "The output directory does not exist.",
    ),
    (
        "Output path is not a directory",
        ErrorCode.OUTPUT_DIR_NOT_A_DIRECTORY,
        "The output path is not a directory.",
    ),
    (
        "No output directory selected",
        ErrorCode.OUTPUT_DIR_MISSING,
        "No output directory was supplied.",
    ),
    (
        "combination length",
        ErrorCode.INVALID_PARAMETER,
        "A combination-length setting is out of range.",
    ),
)


def to_service_error(error: BaseException) -> ServiceError:
    """
    Wrap any exception as a `ServiceError`, preserving its cause.

    A `ServiceError` passes through unchanged. A `ValueError` from the science
    is matched against the table above so the UI gets a code; its original
    wording is kept as `detail` rather than discarded.
    """
    if isinstance(error, ServiceError):
        return error

    if isinstance(error, ValueError):
        text = str(error)
        for needle, code, message in _VALUE_ERROR_CODES:
            if needle in text:
                return ServiceError(code, message, detail=text, cause=error)
        return ServiceError(ErrorCode.INVALID_PARAMETER, text, detail=text, cause=error)

    if isinstance(error, FileNotFoundError):
        return ServiceError(
            ErrorCode.FASTA_NOT_FOUND,
            "The file could not be found.",
            detail=str(error),
            cause=error,
        )

    if isinstance(error, PermissionError):
        return ServiceError(
            ErrorCode.OUTPUT_WRITE_FAILED,
            "Permission was denied while reading or writing a file.",
            detail=str(error),
            cause=error,
        )

    if isinstance(error, OSError):
        return ServiceError(
            ErrorCode.FASTA_UNREADABLE,
            "The file could not be read.",
            detail=str(error),
            cause=error,
        )

    return ServiceError(
        ErrorCode.UNKNOWN,
        "The analysis backend failed unexpectedly.",
        detail=f"{type(error).__name__}: {error}",
        cause=error,
    )
