"""
Project RPC methods.

Thin adapters: they validate the incoming plain data, call `ProjectService`,
and convert its results and errors into the service protocol. No SQL, no
filesystem logic and no scientific work happens here.

One project is open at a time per service process, which matches the desktop
app: one window, one project.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from molecular_diagnosis.project.service import (
    PROJECT_DB_NAME,
    ProjectError,
    ProjectService,
    validate_fasta_candidate,
)
from molecular_diagnosis.progress import STAGE_STARTING
from molecular_diagnosis.service.diagnostics import (
    RunDiagnostics,
    environment_snapshot,
    log_line,
    new_run_id,
)
from molecular_diagnosis.service.errors import ErrorCode, ServiceError

__all__ = ["PROJECT_METHODS", "close_open_project", "current_project"]

_open: dict[str, ProjectService] = {}

#: Project failures that map onto an existing boundary code; anything else
#: keeps its own project-level code, which the frontend can still branch on.
_CODE_ALIASES = {
    "FASTA_NOT_FOUND": ErrorCode.FASTA_NOT_FOUND,
    "FASTA_UNREADABLE": ErrorCode.FASTA_UNREADABLE,
    "FASTA_EMPTY": ErrorCode.FASTA_EMPTY,
    "FASTA_NOT_ALIGNED": ErrorCode.FASTA_NOT_ALIGNED,
    "FOCAL_EMPTY": ErrorCode.FOCAL_EMPTY,
    "FOCAL_NO_MATCH": ErrorCode.FOCAL_NO_MATCH,
    "FOCAL_MATCHES_EVERYTHING": ErrorCode.FOCAL_MATCHES_EVERYTHING,
}


def _as_service_error(error: ProjectError) -> ServiceError:
    # No `cause`, deliberately. A ProjectError is a decision the service made
    # on purpose — "this file is missing", "that entry is not in this file" —
    # and it already carries a code and a sentence written for the user. A
    # traceback would describe where we chose to refuse, which tells nobody
    # anything. Unexpected exceptions still travel with theirs.
    return ServiceError(
        _CODE_ALIASES.get(error.code, error.code),
        error.message,
        detail=error.detail,
    )


def current_project() -> ProjectService:
    project = _open.get("current")
    if project is None:
        raise ServiceError("NO_PROJECT_OPEN", "No project is open.")
    return project


def close_open_project() -> None:
    project = _open.pop("current", None)
    if project is not None:
        project.close()


def _require_str(params: dict[str, Any], key: str) -> str:
    value = params.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER, f"The request is missing '{key}'.", detail=repr(value)
        )
    return value


def _require_ids(params: dict[str, Any], key: str) -> list[str]:
    value = params.get(key)
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER, f"'{key}' must be a list of ids.", detail=repr(value)
        )
    return value


def _optional_scope(params: dict[str, Any], key: str = "fastaFileIds") -> list[str] | None:
    """
    Absent/null means the whole project; a list means exactly those files.

    The two are different requests, so a missing key is never folded into an
    empty list (or vice versa) anywhere on this boundary.
    """
    value = params.get(key)
    if value is None:
        return None
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER,
            f"'{key}' must be a list of ids, or omitted for the whole project.",
            detail=repr(value),
        )
    return value


def _guard(call):
    try:
        return call()
    except ProjectError as error:
        raise _as_service_error(error) from error


# ---------------------------------------------------------------------------
# Lifecycle
# ---------------------------------------------------------------------------


def create_project(params: dict[str, Any]) -> dict[str, Any]:
    """
    Deliberately initialise a new project in a directory.

    Refuses a directory that already holds one: overwriting or adopting an
    existing project.sqlite under the banner of "create" is how someone loses
    the project they meant to open.
    """
    directory = _require_str(params, "projectDir")
    title = params.get("title")
    database = Path(directory) / PROJECT_DB_NAME
    if database.exists():
        raise ServiceError(
            "PROJECT_ALREADY_EXISTS",
            "That folder already contains a project. Open it instead of creating a new one.",
            detail=str(database),
        )

    close_open_project()
    project = ProjectService(Path(directory), create_if_missing=True)
    _open["current"] = project
    return _guard(
        lambda: project.open_project(
            title if isinstance(title, str) and title.strip() else "Untitled project"
        )
    )


def open_project(params: dict[str, Any]) -> dict[str, Any]:
    """
    Open an EXISTING project. Never creates one.

    A folder without a project.sqlite is a mistake to report — the user picked
    the wrong folder — not a reason to silently create an empty project that
    then looks like their work went missing.
    """
    directory = _require_str(params, "projectDir")
    if not (Path(directory) / PROJECT_DB_NAME).is_file():
        # Checked BEFORE closing: picking the wrong folder must not also throw
        # away the project the user already had open.
        raise ServiceError(
            "PROJECT_NOT_FOUND",
            "That folder does not contain a project.",
            detail=str(Path(directory) / PROJECT_DB_NAME),
        )

    close_open_project()
    project = _guard(lambda: ProjectService(Path(directory), create_if_missing=False))
    _open["current"] = project
    return _guard(lambda: project.open_project())


def set_project_title(params: dict[str, Any]) -> dict[str, Any]:
    project = current_project()
    title = _require_str(params, "title")
    return {"metadata": _guard(lambda: project.set_title(title))}


def close_project(_params: dict[str, Any]) -> dict[str, Any]:
    close_open_project()
    return {"closed": True}


def refresh_sources(params: dict[str, Any]) -> dict[str, Any]:
    """Cheap stat refresh for every linked file. Safe to call often."""
    project = current_project()
    strong = bool(params.get("strong", False))
    statuses = _guard(lambda: project.refresh_sources(strong=strong))
    return {"sources": [status.to_payload() for status in statuses]}


# ---------------------------------------------------------------------------
# FASTA files
# ---------------------------------------------------------------------------


def link_fasta(params: dict[str, Any]) -> dict[str, Any]:
    project = current_project()
    path = _require_str(params, "path")
    row, status = _guard(lambda: project.link_fasta(path))
    return {"fastaFileId": row.id, "source": status.to_payload()}


def validate_fasta_candidate_method(params: dict[str, Any]) -> dict[str, Any]:
    """
    Vet a file before it is linked. Writes nothing, and needs no open project.

    The new-project screen runs before any database exists, so this deliberately
    does NOT call `current_project()`. It is the same scan the indexer performs,
    which is the point: a file accepted here cannot be rejected at link time for
    a reason the user was never shown.
    """
    path = _require_str(params, "path")
    return {"candidate": _guard(lambda: validate_fasta_candidate(path))}


def set_fasta_file_locked(params: dict[str, Any]) -> dict[str, Any]:
    """Lock or unlock a linked source. A locked source stays analysable."""
    project = current_project()
    file_id = _require_str(params, "fastaFileId")
    locked = params.get("locked")
    if not isinstance(locked, bool):
        raise ServiceError(ErrorCode.INVALID_PARAMETER, "'locked' must be true or false.")
    status = _guard(lambda: project.set_fasta_file_locked(file_id, locked))
    return {"source": status.to_payload()}


def unlink_fasta(params: dict[str, Any]) -> dict[str, Any]:
    project = current_project()
    file_id = _require_str(params, "fastaFileId")
    _guard(lambda: project.unlink_fasta(file_id))
    return {"removed": file_id}


def relink_fasta(params: dict[str, Any]) -> dict[str, Any]:
    project = current_project()
    file_id = _require_str(params, "fastaFileId")
    path = _require_str(params, "path")
    return _guard(lambda: project.relink_fasta(file_id, path))


def reindex_fasta(params: dict[str, Any]) -> dict[str, Any]:
    project = current_project()
    file_id = _require_str(params, "fastaFileId")
    status = _guard(lambda: project.reindex(file_id))
    return {"source": status.to_payload()}


# ---------------------------------------------------------------------------
# Search and focal sets
# ---------------------------------------------------------------------------


def search_headers(params: dict[str, Any]) -> dict[str, Any]:
    """
    Plain, non-mutating search.

    Unlike `+`, this reports the files it could not search rather than
    refusing: the caller is looking, not committing anything to the project.
    """
    project = current_project()
    query = params.get("query")
    if not isinstance(query, str):
        raise ServiceError(ErrorCode.INVALID_PARAMETER, "'query' must be a string.")
    scope_ids = _optional_scope(params)
    limit = params.get("limit")
    return _guard(
        lambda: project.search_headers(
            query,
            fasta_file_ids=scope_ids,
            limit=limit if isinstance(limit, int) else 200,
        )
    )


def list_focal_sets(_params: dict[str, Any]) -> dict[str, Any]:
    project = current_project()
    return {"focalSets": _guard(project.list_focal_sets)}


def get_focal_set(params: dict[str, Any]) -> dict[str, Any]:
    project = current_project()
    set_id = _require_str(params, "focalSetId")
    return {"focalSet": _guard(lambda: project.get_focal_set(set_id))}


def create_focal_set(params: dict[str, Any]) -> dict[str, Any]:
    project = current_project()
    title = _require_str(params, "title")
    focal_set = _guard(lambda: project.create_focal_set(title))
    # focalSetId is kept alongside the full set so existing callers that only
    # read the id keep working.
    return {"focalSetId": focal_set["id"], "title": focal_set["title"], "focalSet": focal_set}


def rename_focal_set(params: dict[str, Any]) -> dict[str, Any]:
    project = current_project()
    set_id = _require_str(params, "focalSetId")
    title = _require_str(params, "title")
    return {"focalSet": _guard(lambda: project.rename_focal_set(set_id, title))}


def set_focal_set_locked(params: dict[str, Any]) -> dict[str, Any]:
    project = current_project()
    set_id = _require_str(params, "focalSetId")
    locked = params.get("locked")
    if not isinstance(locked, bool):
        raise ServiceError(ErrorCode.INVALID_PARAMETER, "'locked' must be true or false.")
    return {"focalSet": _guard(lambda: project.set_focal_set_locked(set_id, locked))}


def delete_focal_set(params: dict[str, Any]) -> dict[str, Any]:
    project = current_project()
    set_id = _require_str(params, "focalSetId")
    return _guard(lambda: project.delete_focal_set(set_id))


def replace_focal_entries(params: dict[str, Any]) -> dict[str, Any]:
    """
    The large textbox: replace the set's explicit contents wholesale.

    `text` is the raw semicolon-separated box; `headers` is the same thing
    already split. Exactly one is supplied. Neither is a query — these are
    complete headers, and one that matches no FASTA is kept so the UI can
    show it red.
    """
    project = current_project()
    set_id = _require_str(params, "focalSetId")
    text = params.get("text")
    headers = params.get("headers")

    if isinstance(text, str) and headers is None:
        return _guard(lambda: project.replace_focal_entries(set_id, text=text))
    if headers is not None and text is None:
        if not isinstance(headers, list) or not all(isinstance(item, str) for item in headers):
            raise ServiceError(ErrorCode.INVALID_PARAMETER, "'headers' must be a list of strings.")
        return _guard(lambda: project.replace_focal_entries(set_id, headers=headers))

    raise ServiceError(
        ErrorCode.INVALID_PARAMETER, "Supply exactly one of 'text' or 'headers'."
    )


def _require_headers(params: dict[str, Any], key: str = "headers") -> list[str]:
    value = params.get(key)
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER, f"'{key}' must be a list of strings.", detail=repr(value)
        )
    return value


def save_focal_set(params: dict[str, Any]) -> dict[str, Any]:
    """
    The explicit Save the workspace calls. Nothing else in editing writes.

    `focalSetId` absent or null creates; an id updates in place.
    """
    project = current_project()
    set_id = params.get("focalSetId")
    if set_id is not None and not isinstance(set_id, str):
        raise ServiceError(
            ErrorCode.INVALID_PARAMETER, "'focalSetId' must be an id or null.", detail=repr(set_id)
        )
    title = _require_str(params, "title")
    headers = _require_headers(params)

    return {
        "focalSet": _guard(
            lambda: project.save_focal_set(focal_set_id=set_id, title=title, headers=headers)
        )
    }


def header_presence(params: dict[str, Any]) -> dict[str, Any]:
    """
    Presence for arbitrary headers, including ones not saved to any focal set.

    Index-only and non-mutating: safe to call on every debounced keystroke.
    """
    project = current_project()
    headers = _require_headers(params)
    selected = params.get("selectedFastaFileId")
    scope_ids = _optional_scope(params)

    return {
        "entries": _guard(
            lambda: project.header_presence(
                headers,
                selected_file_id=selected if isinstance(selected, str) else None,
                fasta_file_ids=scope_ids,
            )
        )
    }


def match_focal_headers(params: dict[str, Any]) -> dict[str, Any]:
    """`-` against a working copy, using Python's casefold rather than JS's."""
    project = current_project()
    query = params.get("query")
    if not isinstance(query, str):
        raise ServiceError(ErrorCode.INVALID_PARAMETER, "'query' must be a string.")
    headers = _require_headers(params)
    return _guard(lambda: project.match_focal_headers(query, headers))


def resolve_focal_add_query(params: dict[str, Any]) -> dict[str, Any]:
    """
    What `+` would add, resolved but NOT added.

    Non-mutating, uncapped, and scope-verified exactly as the mutating `+` is:
    the renderer appends the result to a working draft and saves separately.
    """
    project = current_project()
    query = params.get("query")
    if not isinstance(query, str):
        raise ServiceError(ErrorCode.INVALID_PARAMETER, "'query' must be a string.")
    scope_ids = _optional_scope(params)
    return _guard(lambda: project.resolve_focal_add_query(query, fasta_file_ids=scope_ids))


def add_focal_entries(params: dict[str, Any]) -> dict[str, Any]:
    """`+` mode: a query expands to complete headers, which become members."""
    project = current_project()
    set_id = _require_str(params, "focalSetId")
    query = params.get("query")
    if not isinstance(query, str):
        raise ServiceError(ErrorCode.INVALID_PARAMETER, "'query' must be a string.")
    scope_ids = _optional_scope(params)
    return _guard(
        lambda: project.add_focal_entries_from_query(set_id, query, fasta_file_ids=scope_ids)
    )


def remove_focal_entries(params: dict[str, Any]) -> dict[str, Any]:
    """`-` mode: drop members whose exact stored header contains the query."""
    project = current_project()
    set_id = _require_str(params, "focalSetId")
    query = params.get("query")
    if not isinstance(query, str):
        raise ServiceError(ErrorCode.INVALID_PARAMETER, "'query' must be a string.")
    return _guard(lambda: project.remove_focal_entries_by_query(set_id, query))


def focal_presence(params: dict[str, Any]) -> dict[str, Any]:
    project = current_project()
    set_id = _require_str(params, "focalSetId")
    selected = params.get("selectedFastaFileId")
    scope_ids = _optional_scope(params)
    presence = _guard(
        lambda: project.focal_presence(
            set_id,
            selected_file_id=selected if isinstance(selected, str) else None,
            fasta_file_ids=scope_ids,
        )
    )
    return {"entries": [item.to_payload() for item in presence]}


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------


def run_project_molecular_diagnosis(params: dict[str, Any]) -> dict[str, Any]:
    """
    The one long-running method, and therefore the instrumented one.

    Every run gets a short id that appears in each of its diagnostic lines, so
    a log holding several runs can still be read. `runToken` is the CALLER's
    correlation id: the renderer mints it, it comes back on every progress
    notification, and that is how the UI knows a late notification belongs to a
    run it has already replaced.
    """
    project = current_project()
    set_id = _require_str(params, "focalSetId")
    file_ids = _require_ids(params, "fastaFileIds")
    single = bool(params.get("singleFile", len(file_ids) == 1))
    options = params.get("options")
    resume = params.get("resume")
    token = params.get("runToken")

    run = RunDiagnostics(new_run_id(), run_token=token if isinstance(token, str) else None)
    log_line(
        "run.start",
        run=run.run_id,
        runToken=run.run_token,
        focalSetId=set_id,
        fastaFileIds=file_ids,
        singleFile=single,
        environment=environment_snapshot(),
    )
    run.progress(STAGE_STARTING)

    try:
        result = _guard(
            lambda: project.run_molecular_diagnosis(
                focal_set_id=set_id,
                fasta_file_ids=file_ids,
                single_file=single,
                options=options if isinstance(options, dict) else {},
                resume=resume if isinstance(resume, dict) else None,
                observer=run,
            )
        )
    except ServiceError as error:
        # An expected refusal: the code and the stage it happened in, no
        # traceback. The stage is the part that is hard to guess afterwards.
        log_line(
            "run.refused",
            run=run.run_id,
            stage=run.last_stage,
            code=str(error.code),
            message=error.message,
            elapsedMs=run.elapsed_ms,
        )
        raise
    except BaseException as error:
        run.failure(error)
        raise

    log_line("run.end", run=run.run_id, ok=True, elapsedMs=run.elapsed_ms)
    return result


def diagnostics_environment(_params: dict[str, Any]) -> dict[str, Any]:
    """
    The environment this service is actually running in.

    Exists so two machines can be compared without having to reproduce a
    failure first: run it from a terminal (see REPO_MAP section 14) and diff
    the two answers.
    """
    return {"environment": environment_snapshot()}


PROJECT_METHODS = {
    "diagnostics.environment": diagnostics_environment,
    "project.create": create_project,
    "project.open": open_project,
    "project.close": close_project,
    "project.setTitle": set_project_title,
    "project.refreshSources": refresh_sources,
    "project.validateFastaCandidate": validate_fasta_candidate_method,
    "project.linkFasta": link_fasta,
    "project.unlinkFasta": unlink_fasta,
    "project.setFastaFileLocked": set_fasta_file_locked,
    "project.relinkFasta": relink_fasta,
    "project.reindexFasta": reindex_fasta,
    "project.searchHeaders": search_headers,
    "project.listFocalSets": list_focal_sets,
    "project.getFocalSet": get_focal_set,
    "project.createFocalSet": create_focal_set,
    "project.renameFocalSet": rename_focal_set,
    "project.setFocalSetLocked": set_focal_set_locked,
    "project.deleteFocalSet": delete_focal_set,
    "project.replaceFocalEntries": replace_focal_entries,
    "project.saveFocalSet": save_focal_set,
    "project.headerPresence": header_presence,
    "project.matchFocalHeaders": match_focal_headers,
    "project.resolveFocalAddQuery": resolve_focal_add_query,
    "project.addFocalEntries": add_focal_entries,
    "project.removeFocalEntries": remove_focal_entries,
    "project.focalPresence": focal_presence,
    "project.runMolecularDiagnosis": run_project_molecular_diagnosis,
}
