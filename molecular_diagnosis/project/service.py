"""
ProjectService — the façade the RPC layer talks to.

It owns the database connection, the source verifier, and the alignment cache
for one open project, and it is the only place that composes them. Handlers
call methods here; they never touch SQL, the filesystem, or the cache.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from molecular_diagnosis.focal import ExactHeaders
from molecular_diagnosis.project.alignments import AlignmentCache, CachedAlignment
from molecular_diagnosis.project.db import Capabilities, fts_available, open_project_db, optimize
from molecular_diagnosis.project.indexing import (
    FastaScan,
    SourceChangedDuringRead,
    scan_fasta,
)
from molecular_diagnosis.project.locations import (
    EntryPresence,
    PresenceState,
    presence_for_entries,
    verify_entry_locations,
)
from molecular_diagnosis.project.paths import display_name_for, normalise_path_key
from molecular_diagnosis.project.repository import FastaFileRow, FocalSetRow, Repository
from molecular_diagnosis.project.search import resolve_headers, search_headers
from molecular_diagnosis.project.sources import SourceState, SourceStatus, SourceVerifier, hash_file

__all__ = [
    "AnalysisScope",
    "ProjectError",
    "ProjectService",
    "ScopeProblem",
    "dedupe_headers",
    "parse_focal_text",
]

PROJECT_DB_NAME = "project.sqlite"
OUTPUTS_DIR = "outputs"

#: How many times a racing writer may defeat a read before we give up.
SCAN_ATTEMPTS = 3


class ProjectError(RuntimeError):
    """A project-level failure with a stable code for the service boundary."""

    def __init__(self, code: str, message: str, *, detail: str | None = None) -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.detail = detail


def validate_fasta_candidate(path: str | Path) -> dict[str, object]:
    """
    Decide whether a file may be linked at all — WITHOUT linking it.

    Deliberately module-level: the new-project screen has to vet files before
    any project exists, so this cannot depend on an open database. It also runs
    no differently from the real thing, because it uses the same `scan_fasta`
    the indexer does. Judging a FASTA by its extension would accept a renamed
    spreadsheet and reject a valid `.txt`.

    Returns the metadata the pending row shows. Raises a coded `ProjectError`
    for every reason a file cannot be used, so the caller can say which one.
    """
    candidate = Path(path)
    if not candidate.is_file():
        raise ProjectError(
            "FASTA_NOT_FOUND", "That file could not be found.", detail=str(path)
        )

    try:
        scan = scan_fasta(candidate)
    except SourceChangedDuringRead as error:
        raise ProjectError(
            "SOURCE_CHANGED_DURING_READ",
            "The file kept changing while it was being read.",
            detail=str(error),
        ) from error
    except OSError as error:
        raise ProjectError(
            "FASTA_UNREADABLE", "That file could not be read.", detail=str(error)
        ) from error

    # No '>' line anywhere means this is not a FASTA at all, which is the same
    # outcome as an empty one: there is nothing to analyse.
    if not scan.records:
        raise ProjectError(
            "FASTA_EMPTY",
            "That file contains no FASTA sequences.",
            detail=candidate.name,
        )
    if not scan.sequences:
        raise ProjectError(
            "FASTA_EMPTY", "That FASTA file contains no sequences.", detail=candidate.name
        )
    if scan.alignment_length is None:
        raise ProjectError(
            "FASTA_NOT_ALIGNED",
            "The FASTA file needs to be aligned.",
            detail=(
                f"{candidate.name}: the sequences are not all the same length, so this is "
                "not an aligned FASTA."
            ),
        )

    return {
        "path": str(candidate),
        "displayName": display_name_for(candidate),
        "sequenceCount": scan.sequence_count,
        "alignmentLength": scan.alignment_length,
        "duplicateHeaderCount": scan.duplicate_header_count,
    }


def dedupe_headers(headers: Sequence[str]) -> list[str]:
    """
    Set semantics, first occurrence wins.

    `focal_set_entry` is UNIQUE on (set, header), so a repeated header cannot
    be stored twice anyway; collapsing here means the caller gets back exactly
    what was stored instead of a silent constraint failure.
    """
    seen: set[str] = set()
    ordered: list[str] = []
    for raw in headers:
        header = str(raw).strip()
        if not header or header in seen:
            continue
        seen.add(header)
        ordered.append(header)
    return ordered


def parse_focal_text(text: str) -> list[str]:
    """
    Parse the large focal textbox into exact headers.

    `a ; b;c` is three headers. Whitespace immediately around a separator is
    trimmed, empty entries are ignored, and duplicates collapse to their first
    occurrence. Whitespace INSIDE a header is preserved, because FASTA headers
    routinely contain spaces.
    """
    return dedupe_headers(text.split(";"))


@dataclass(frozen=True)
class ScopeProblem:
    code: str
    message: str
    detail: str | None = None

    def to_payload(self) -> dict[str, object]:
        return {"code": self.code, "message": self.message, "detail": self.detail}


@dataclass
class AnalysisScope:
    """A validated, in-memory combination of one or more alignments."""

    fasta_file_ids: list[str]
    sequences: dict[str, str]
    alignment_length: int | None
    problems: list[ScopeProblem]

    @property
    def ok(self) -> bool:
        return not self.problems


class ProjectService:
    def __init__(self, project_dir: str | Path, *, create_if_missing: bool = True) -> None:
        self.project_dir = Path(project_dir)

        # "Create a project" and "open an existing project" are different user
        # intentions, so they must be different operations here. Opening a
        # folder that holds no project.sqlite is a mistake to report, not a
        # silent invitation to initialise a second, empty project on top of
        # whatever the user actually meant to open.
        if not create_if_missing and not (self.project_dir / PROJECT_DB_NAME).is_file():
            raise ProjectError(
                "PROJECT_NOT_FOUND",
                "That folder does not contain a project.",
                detail=str(self.project_dir / PROJECT_DB_NAME),
            )

        self.project_dir.mkdir(parents=True, exist_ok=True)
        self.outputs_dir = self.project_dir / OUTPUTS_DIR
        self.outputs_dir.mkdir(parents=True, exist_ok=True)

        con, capabilities, fts = open_project_db(self.project_dir / PROJECT_DB_NAME)
        self.connection = con
        self.capabilities: Capabilities = capabilities
        self.repository = Repository(con, fts_enabled=fts)
        self.verifier = SourceVerifier(self.repository)
        self.alignments = AlignmentCache()

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self) -> None:
        optimize(self.connection)
        self.connection.close()

    def open_project(self, title: str = "Untitled project") -> dict[str, object]:
        """
        Open (or initialise) the project and immediately report source status.

        The cheap stat check runs for every linked FASTA here, so a missing
        file is visible the moment the project opens rather than when an
        analysis fails.

        `title` is only used when the project metadata row does not exist yet:
        opening an existing project never renames it. Use `set_title` for that.
        """
        metadata = self.repository.ensure_project(title)
        statuses = self.refresh_sources(strong=False)
        return {
            "projectDir": str(self.project_dir),
            "outputsDir": str(self.outputs_dir),
            "metadata": {
                "projectUuid": metadata["project_uuid"],
                "title": metadata["title"],
            },
            "capabilities": {
                "sqliteVersion": self.capabilities.sqlite_version,
                "fts5": self.capabilities.fts5,
                "trigram": self.capabilities.trigram,
                "acceleratedSearch": self.repository.fts_enabled,
            },
            "sources": [status.to_payload() for status in statuses],
        }

    def set_title(self, title: str) -> dict[str, object]:
        """Rename the project. The only mutable piece of project metadata."""
        cleaned = title.strip()
        if not cleaned:
            raise ProjectError("INVALID_PARAMETER", "A project title cannot be empty.")

        self.repository.begin()
        try:
            self.repository.set_project_title(cleaned)
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise

        metadata = self.repository.ensure_project(cleaned)
        return {"projectUuid": metadata["project_uuid"], "title": metadata["title"]}

    # ------------------------------------------------------------------
    # Sources
    # ------------------------------------------------------------------

    def refresh_sources(self, *, strong: bool = False) -> list[SourceStatus]:
        """Cheap status for every linked file. Never hashes unless asked."""
        return [
            self.verifier.ensure_current(row, strong=strong)
            for row in self.repository.list_fasta_files()
        ]

    def link_fasta(self, path: str) -> tuple[FastaFileRow, SourceStatus]:
        """
        Register a FASTA and build its first complete header index.

        Linking is ALL OR NOTHING for a new file. It used to insert the row
        first and index afterwards, so a file that turned out to be ragged or
        unreadable left a permanent `fasta_file` row behind that the user then
        had to clean up. Now the candidate is vetted before anything is
        written, and if indexing still fails — the file can change between the
        two reads — the row this call created is removed again.

        Re-linking a path that is already in the project is left alone on
        failure: that row is not ours to delete.
        """
        validate_fasta_candidate(path)

        existing = self.repository.find_by_path(path)
        self.repository.begin()
        try:
            row = self.repository.add_fasta_file(path)
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise

        try:
            if existing is None or not row.has_index:
                status = self.reindex(row.id)
            else:
                status = self.verifier.ensure_current(row, strong=True)
        except Exception:
            if existing is None:
                # We created this row moments ago; nothing else can depend on
                # it yet, so removing it restores the project exactly.
                self.repository.begin()
                try:
                    self.repository.remove_fasta_file(row.id)
                    self.repository.commit()
                except Exception:
                    self.repository.rollback()
                    raise
                self.verifier.forget(row.id)
                self.alignments.invalidate(row.id)
            raise

        row = self.repository.get_fasta_file(row.id)
        assert row is not None
        return row, status

    def set_fasta_file_locked(self, file_id: str, locked: bool) -> SourceStatus:
        """
        Lock or unlock a linked source.

        A locked source stays analysable — locking protects the LINK, not the
        data. It is what stops a settled input being unlinked by a mis-click.
        """
        self._require_file(file_id)
        self.repository.begin()
        try:
            self.repository.set_fasta_file_locked(file_id, locked)
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise
        return self.verifier.cheap_status(self._require_file(file_id))

    def unlink_fasta(self, file_id: str) -> None:
        row = self._require_file(file_id)
        if row.locked:
            # Enforced here, not only by hiding the button: a lock that only
            # exists in React is not a property of the project.
            raise ProjectError(
                "FASTA_FILE_LOCKED",
                f"{row.display_name} is locked. Unlock it before removing it.",
                detail=file_id,
            )

        self.repository.begin()
        try:
            self.repository.remove_fasta_file(file_id)
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise
        self.verifier.forget(file_id)
        self.alignments.invalidate(file_id)

    # ------------------------------------------------------------------
    # Indexing
    # ------------------------------------------------------------------

    def _stable_scan(self, row: FastaFileRow) -> FastaScan:
        """
        Read the file completely, retrying if it moved underneath us.

        A snapshot assembled while the file was being written is discarded, not
        committed: the index must represent exactly one stable source version.
        """
        last: Exception | None = None
        for _attempt in range(SCAN_ATTEMPTS):
            try:
                return scan_fasta(row.source_path)
            except SourceChangedDuringRead as error:
                last = error
            except OSError as error:
                raise ProjectError(
                    "FASTA_UNREADABLE", "The FASTA file could not be read.", detail=str(error)
                ) from error
        raise ProjectError(
            "SOURCE_CHANGED_DURING_READ",
            "The FASTA file kept changing while it was being read, so it was not indexed.",
            detail=str(last),
        )

    def reindex(self, file_id: str, *, force: bool = False) -> SourceStatus:
        """
        Rebuild one file's complete header index from its current bytes.

        Parse and validate first; only a complete, stable, valid scan is
        allowed to replace the previous snapshot. A failure leaves the old
        snapshot intact but the source stale, i.e. runtime-unusable.
        """
        row = self._require_file(file_id)
        scan = self._stable_scan(row)

        if not scan.sequences:
            raise ProjectError("FASTA_EMPTY", "The FASTA file contains no sequences.")
        if scan.alignment_length is None:
            raise ProjectError(
                "FASTA_NOT_ALIGNED",
                "The sequences are not all the same length, so this is not an aligned FASTA.",
            )

        content_changed = scan.sha256 != row.indexed_sha256

        self.repository.begin()
        try:
            self.repository.replace_index(file_id, scan, bump_revision=content_changed or force)
            # Focal locations are NOT wiped here. They are a separate cache and
            # are repaired individually below, from this same scan.
            entries = self._all_focal_entries()
            if entries:
                verify_entry_locations(self.repository, row, entries, scan=scan)
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise

        if content_changed:
            self.alignments.invalidate(file_id)
        self.verifier.note_scan(file_id, scan, matched_index=True)
        self.alignments.put(file_id, self._require_file(file_id).index_revision, scan)

        return self.verifier.ensure_current(self._require_file(file_id), strong=False)

    def ensure_index_usable(self, file_id: str) -> SourceStatus:
        """
        Guarantee the complete header index describes the current file.

        Used before any arbitrary header search, which is the one operation
        that genuinely needs the whole snapshot to be trustworthy.
        """
        row = self._require_file(file_id)
        status = self.verifier.ensure_current(row, strong=True)

        if status.state in (SourceState.MISSING, SourceState.UNREADABLE):
            return status
        if status.state is SourceState.CURRENT:
            return status
        return self.reindex(file_id)

    # ------------------------------------------------------------------
    # Relinking
    # ------------------------------------------------------------------

    def relink_fasta(self, file_id: str, new_path: str) -> dict[str, object]:
        """
        Point a linked FASTA at a replacement file.

        A byte-identical replacement is proved by SHA-256 and keeps everything:
        the header index, the focal locations and index_revision. A different
        but valid FASTA is reindexed transactionally. An invalid candidate is
        rejected without touching project state.
        """
        row = self._require_file(file_id)
        candidate = Path(new_path)
        if not candidate.is_file():
            raise ProjectError(
                "FASTA_NOT_FOUND", "The replacement file could not be found.", detail=new_path
            )

        clash = self.repository.find_by_path(new_path)
        if clash is not None and clash.id != file_id:
            raise ProjectError(
                "PATH_ALREADY_LINKED",
                "That file is already linked to this project as a different entry.",
                detail=clash.display_name,
            )

        try:
            digest, _size, _mtime = hash_file(candidate)
        except (OSError, SourceChangedDuringRead) as error:
            raise ProjectError(
                "FASTA_UNREADABLE", "The replacement file could not be read.", detail=str(error)
            ) from error

        if row.indexed_sha256 is not None and digest == row.indexed_sha256:
            # Same bytes at a new path: nothing derived can have changed.
            self.repository.begin()
            try:
                self.repository.update_source_path(file_id, str(candidate))
                stat = os.stat(candidate)
                self.repository.touch_indexed_stat(file_id, stat.st_size, stat.st_mtime_ns)
                self.repository.commit()
            except Exception:
                self.repository.rollback()
                raise
            self.verifier.forget(file_id)
            updated = self._require_file(file_id)
            status = self.verifier.ensure_current(updated, strong=True)
            return {
                "identical": True,
                "reindexed": False,
                "indexRevision": updated.index_revision,
                "source": status.to_payload(),
            }

        # Different content: validate before committing anything.
        probe = scan_fasta(candidate)
        if not probe.sequences:
            raise ProjectError("FASTA_EMPTY", "The replacement file contains no sequences.")
        if probe.alignment_length is None:
            raise ProjectError(
                "FASTA_NOT_ALIGNED",
                "The replacement file is not an aligned FASTA, so it was not linked.",
            )

        self.repository.begin()
        try:
            self.repository.update_source_path(file_id, str(candidate))
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise

        self.verifier.forget(file_id)
        self.alignments.invalidate(file_id)
        status = self.reindex(file_id)
        return {
            "identical": False,
            "reindexed": True,
            "indexRevision": self._require_file(file_id).index_revision,
            "source": status.to_payload(),
        }

    # ------------------------------------------------------------------
    # Search and focal sets
    # ------------------------------------------------------------------

    def _resolve_scope(self, fasta_file_ids: list[str] | None) -> list[str]:
        """
        Turn a requested scope into concrete file ids.

        `None` means the whole project. `[]` means an empty scope. Those are
        different requests, and the `ids or all_files` idiom this replaces
        could not express the second one: it quietly promoted "nothing
        selected" to "everything".
        """
        if fasta_file_ids is None:
            return [row.id for row in self.repository.list_fasta_files()]
        return list(fasta_file_ids)

    def _partition_scope(self, scope: list[str]) -> tuple[list[str], list[SourceStatus]]:
        """
        Split a scope into files whose header index may be trusted and files
        whose may not, verifying (and refreshing) each one on the way.

        Unknown ids raise rather than being dropped: a scope naming a file that
        is not in this project is a caller bug, not a partial result.
        """
        usable: list[str] = []
        unusable: list[SourceStatus] = []
        for file_id in scope:
            self._require_file(file_id)
            status = self.ensure_index_usable(file_id)
            if status.index_usable:
                usable.append(file_id)
            else:
                unusable.append(status)
        return usable, unusable

    def search_headers(
        self, query: str, *, fasta_file_ids: list[str] | None = None, limit: int | None = 200
    ) -> dict[str, object]:
        """
        Substring search over headers.

        The complete header index has to be trustworthy for an arbitrary
        search, so every in-scope file is verified (and refreshed if stale)
        first.

        This is the non-mutating endpoint, so a file it could not use is
        REPORTED rather than fatal: the user gets the matches that do exist
        plus an explicit list of what could not be searched. `+` cannot afford
        that, and refuses instead — see `add_focal_entries_from_query`.
        """
        usable, unusable = self._partition_scope(self._resolve_scope(fasta_file_ids))
        hits = search_headers(self.repository, query, fasta_file_ids=usable, limit=limit)
        return {
            "query": query,
            "hits": [hit.to_payload() for hit in hits],
            "unavailable": [status.to_payload() for status in unusable],
        }

    # ------------------------------------------------------------------
    # Focal sets
    # ------------------------------------------------------------------

    def _require_focal_set(self, set_id: str) -> FocalSetRow:
        """
        Every focal operation starts here.

        Without it an unknown id reaches SQLite and comes back as a foreign-key
        or integrity error — a database detail leaking through the service
        boundary, with no stable code for the UI to branch on.
        """
        row = self.repository.get_focal_set(set_id)
        if row is None:
            raise ProjectError(
                "UNKNOWN_FOCAL_SET", "That focal set is not in this project.", detail=set_id
            )
        return row

    def _require_unlocked(self, set_id: str) -> FocalSetRow:
        """
        A locked focal set may be selected, checked and analysed — never edited.

        Enforced here rather than by disabling buttons: the renderer is one
        client of this service, and a rule that only exists in React is not a
        rule about the data.
        """
        row = self._require_focal_set(set_id)
        if row.locked:
            raise ProjectError(
                "FOCAL_SET_LOCKED",
                f"The focal set {row.title!r} is locked, so its contents cannot be changed.",
                detail=set_id,
            )
        return row

    @staticmethod
    def _require_query(query: str) -> str:
        """
        A mutation query must have content.

        An empty `-` query is the dangerous one: every stored header contains
        the empty string, so it would match — and delete — the entire set.
        """
        cleaned = query.strip()
        if not cleaned:
            raise ProjectError(
                "FOCAL_QUERY_EMPTY", "Type something to search for before adding or removing."
            )
        return cleaned

    def focal_set_payload(self, row: FocalSetRow) -> dict[str, object]:
        return {
            "id": row.id,
            "title": row.title,
            "locked": row.locked,
            "entries": [
                {"id": entry.id, "header": entry.header}
                for entry in self.repository.list_focal_entries(row.id)
            ],
        }

    def list_focal_sets(self) -> list[dict[str, object]]:
        return [self.focal_set_payload(row) for row in self.repository.list_focal_sets()]

    def get_focal_set(self, set_id: str) -> dict[str, object]:
        return self.focal_set_payload(self._require_focal_set(set_id))

    def create_focal_set(self, title: str) -> dict[str, object]:
        cleaned = title.strip()
        if not cleaned:
            raise ProjectError("INVALID_PARAMETER", "A focal set needs a title.")
        self.repository.begin()
        try:
            row = self.repository.create_focal_set(cleaned)
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise
        return self.focal_set_payload(row)

    def rename_focal_set(self, set_id: str, title: str) -> dict[str, object]:
        self._require_unlocked(set_id)
        cleaned = title.strip()
        if not cleaned:
            raise ProjectError("INVALID_PARAMETER", "A focal set needs a title.")
        self.repository.begin()
        try:
            self.repository.rename_focal_set(set_id, cleaned)
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise
        return self.focal_set_payload(self._require_focal_set(set_id))

    def set_focal_set_locked(self, set_id: str, locked: bool) -> dict[str, object]:
        """
        Lock or unlock. Unlocking a locked set is the one mutation it allows,
        because otherwise a lock could never be undone.
        """
        self._require_focal_set(set_id)
        self.repository.begin()
        try:
            self.repository.set_focal_set_locked(set_id, locked)
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise
        return self.focal_set_payload(self._require_focal_set(set_id))

    def delete_focal_set(self, set_id: str) -> dict[str, object]:
        self._require_unlocked(set_id)
        self.repository.begin()
        try:
            self.repository.delete_focal_set(set_id)
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise
        return {"deleted": set_id}

    def replace_focal_entries(
        self, set_id: str, *, text: str | None = None, headers: list[str] | None = None
    ) -> dict[str, object]:
        """
        Make the focal set's explicit contents exactly what the textbox says.

        The large textbox is the set's actual membership, not a query field, so
        this takes complete headers. A manually typed header does NOT have to
        exist in any FASTA: it is kept so presence can report it as missing
        (red) instead of the edit silently discarding it.

        Applied as one transactional diff — see
        `Repository.replace_focal_entries` — so the caller may debounce and
        resend freely without churning entry ids or the location cache.
        """
        self._require_unlocked(set_id)

        if (text is None) == (headers is None):
            raise ProjectError(
                "INVALID_PARAMETER", "Supply exactly one of 'text' or 'headers'."
            )

        wanted = parse_focal_text(text) if text is not None else dedupe_headers(headers or [])

        self.repository.begin()
        try:
            diff = self.repository.replace_focal_entries(set_id, wanted)
            if diff["added"]:
                added = set(diff["added"])
                new_entries = [
                    entry
                    for entry in self.repository.list_focal_entries(set_id)
                    if entry.header in added
                ]
                # Seed only from indexes already proven current this session.
                # A debounced keystroke must never trigger a rehash or a
                # reindex; an entry with no seeded location is not wrong, it
                # just verifies itself the next time presence is asked for.
                self._seed_locations_from_index(new_entries, self._trusted_file_ids())
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise

        return {
            "focalSetId": set_id,
            "headers": wanted,
            "added": diff["added"],
            "removed": diff["removed"],
            "kept": diff["kept"],
        }

    def save_focal_set(
        self, *, focal_set_id: str | None, title: str, headers: Sequence[str]
    ) -> dict[str, object]:
        """
        The UI's explicit Save boundary: title + membership, in one transaction.

        Editing a focal set in the workspace writes NOTHING; the renderer holds
        a working copy and calls this when the user saves. That is what makes
        "unsaved changes" a real state rather than a label on data already
        committed.

        `focal_set_id=None` creates the set and its entries together — a create
        that half-succeeded would leave a titled, empty set behind. An existing
        id updates in place, as a diff, so surviving headers keep their entry
        ids and their `focal_entry_location` rows.

        Headers are trimmed and deduplicated, first occurrence winning. A header
        that matches no FASTA is still stored: it is a member the user typed,
        and presence will report it as missing rather than the save discarding
        it.
        """
        cleaned_title = title.strip()
        if not cleaned_title:
            raise ProjectError("INVALID_PARAMETER", "A focal set needs a title.")

        wanted = dedupe_headers(headers)

        if focal_set_id is not None:
            self._require_unlocked(focal_set_id)

        self.repository.begin()
        try:
            if focal_set_id is None:
                row = self.repository.create_focal_set(cleaned_title)
                set_id = row.id
            else:
                set_id = focal_set_id
                self.repository.rename_focal_set(set_id, cleaned_title)

            diff = self.repository.replace_focal_entries(set_id, wanted)
            if diff["added"]:
                added = set(diff["added"])
                new_entries = [
                    entry
                    for entry in self.repository.list_focal_entries(set_id)
                    if entry.header in added
                ]
                self._seed_locations_from_index(new_entries, self._trusted_file_ids())
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise

        return self.focal_set_payload(self._require_focal_set(set_id))

    def header_presence(
        self,
        headers: Sequence[str],
        *,
        selected_file_id: str | None = None,
        fasta_file_ids: list[str] | None = None,
    ) -> list[dict[str, object]]:
        """
        Presence of arbitrary exact headers, including UNSAVED draft headers.

        Read-only in every sense: it answers from the header index with a
        batched lookup, creates no `focal_entry_location` rows, opens no FASTA,
        and — critically — never hashes or reindexes. It is called on a debounce
        while the user types, so anything that touched the filesystem per header
        would make typing cost I/O.

        A file counts as answerable when its stat still matches the fingerprint
        the index was built from. That is evidence, not proof, and it is the
        right bar HERE: these colours are editing feedback. The run gate is
        `focal_presence` + `require_focal_completeness`, which verify locations
        against the real bytes before any analysis starts.
        """
        wanted = dedupe_headers(headers)
        if not wanted:
            return []

        scope = self._resolve_scope(fasta_file_ids)
        answerable: list[str] = []
        unanswerable: set[str] = set()
        for row in self.repository.list_fasta_files():
            if row.id not in scope:
                continue
            status = self.verifier.cheap_status(row)
            if status.state in (SourceState.CURRENT, SourceState.UNVERIFIED):
                answerable.append(row.id)
            else:
                unanswerable.add(row.id)

        counts = self.repository.count_headers_in_files(wanted, answerable)

        results: list[dict[str, object]] = []
        for header in wanted:
            occurrences = {
                file_id: counts[(file_id, header)]
                for file_id in answerable
                if (file_id, header) in counts
            }

            if selected_file_id is None:
                # All files: presence anywhere in scope is success, so there is
                # no "elsewhere" to distinguish and orange is not meaningful.
                if occurrences:
                    state = PresenceState.PRESENT_CURRENT
                elif unanswerable:
                    state = PresenceState.UNKNOWN
                else:
                    state = PresenceState.MISSING
            elif selected_file_id in unanswerable:
                state = PresenceState.UNKNOWN
            elif occurrences.get(selected_file_id):
                state = PresenceState.PRESENT_CURRENT
            elif occurrences:
                state = PresenceState.PRESENT_OTHER
            elif unanswerable:
                state = PresenceState.UNKNOWN
            else:
                state = PresenceState.MISSING

            results.append(
                {
                    "header": header,
                    "state": state.value,
                    "occurrences": occurrences,
                }
            )

        return results

    def match_focal_headers(self, query: str, headers: Sequence[str]) -> dict[str, object]:
        """
        Which of these headers contain `query`, case-insensitively.

        This exists so `-` can operate on an unsaved working copy without
        reimplementing Python's `str.casefold()` in JavaScript. Unicode case
        folding is not the same as `toLowerCase()`, and two implementations of
        "the same" rule would eventually disagree about which member a `-`
        removes.
        """
        needle = self._require_query(query).casefold()
        return {
            "query": query,
            "matched": [header for header in headers if needle in header.casefold()],
        }

    def _trusted_file_ids(self) -> list[str]:
        """Files whose header index is already known-good, without any I/O beyond stat."""
        return [
            row.id
            for row in self.repository.list_fasta_files()
            if self.verifier.cheap_status(row).index_usable
        ]

    def _expand_add_query(
        self, query: str, fasta_file_ids: list[str] | None
    ) -> tuple[str, list[str], list[str]]:
        """
        The shared body of `+`: verify the scope, then resolve the query.

        Returns (cleaned query, complete matching headers, usable file ids).

        `+` promises "every header in this scope that matches". If a requested
        file cannot be searched the whole expansion is refused rather than
        answered from the subset that happened to be readable — that subset
        looks complete and is not reproducible.

        There is deliberately NO result limit here. A cap belongs to a preview
        listing; applying one to `+` would silently drop members from the focal
        set the user asked for.
        """
        cleaned = self._require_query(query)

        usable, unusable = self._partition_scope(self._resolve_scope(fasta_file_ids))
        if unusable:
            raise ProjectError(
                "SEARCH_SCOPE_UNAVAILABLE",
                "Some of the selected files cannot be searched right now, so nothing was "
                "added.",
                detail="; ".join(
                    f"{status.display_name} ({status.state.value})" for status in unusable
                ),
            )

        return cleaned, resolve_headers(self.repository, cleaned, fasta_file_ids=usable), usable

    def resolve_focal_add_query(
        self, query: str, *, fasta_file_ids: list[str] | None = None
    ) -> dict[str, object]:
        """
        What `+` WOULD add, without adding it or writing anything.

        Exists because the renderer edits a working copy: it needs the complete
        expansion to append to the draft, and the capped `search_headers`
        preview cannot supply that — a query matching 300 headers must yield
        all 300, not the first 200.

        Same substring rule, same scope verification, same deterministic
        project-file/record ordering as `add_focal_entries_from_query`; they
        share `_expand_add_query`, so the two cannot drift.
        """
        cleaned, headers, _usable = self._expand_add_query(query, fasta_file_ids)
        return {"query": cleaned, "headers": headers}

    def add_focal_entries_from_query(
        self, set_id: str, query: str, *, fasta_file_ids: list[str] | None = None
    ) -> dict[str, object]:
        """
        `+` : query -> substring search -> complete headers -> dedupe -> add.

        The query itself never becomes a member. This is the MUTATING form,
        kept for the direct-to-database path; the renderer resolves through
        `resolve_focal_add_query` and saves explicitly.
        """
        self._require_unlocked(set_id)
        query, headers, usable = self._expand_add_query(query, fasta_file_ids)

        self.repository.begin()
        try:
            added = self.repository.add_focal_entries(set_id, headers)
            # Seed the persistent location cache straight from the header
            # index, which was just verified above. Without this a brand-new
            # entry would have no saved location to verify later, and the
            # first presence check would have to rescan every file.
            if added:
                new_entries = [
                    entry
                    for entry in self.repository.list_focal_entries(set_id)
                    if entry.header in set(added)
                ]
                self._seed_locations_from_index(new_entries, usable)
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise

        return {"query": query, "matched": headers, "added": added}

    def _seed_locations_from_index(self, entries, fasta_file_ids: list[str]) -> None:
        """Record where each entry's exact header sits, per the current index."""
        for entry in entries:
            for file_id in fasta_file_ids:
                rows = self.repository.locate_header(file_id, entry.header)
                self.repository.replace_entry_locations_in_file(
                    entry.id,
                    file_id,
                    [
                        (row["ordinal"], row["record_start_byte"], row["record_end_byte"])
                        for row in rows
                    ],
                )

    def remove_focal_entries_by_query(self, set_id: str, query: str) -> dict[str, object]:
        """
        `-` : drop members whose EXACT stored header contains the query.

        The comparison is against the stored headers, not against the FASTA:
        this removes from the set, it does not re-search the file.
        """
        self._require_unlocked(set_id)
        needle = self._require_query(query).casefold()
        entries = self.repository.list_focal_entries(set_id)
        doomed = [entry.header for entry in entries if needle in entry.header.casefold()]

        self.repository.begin()
        try:
            removed = self.repository.remove_focal_entries(set_id, doomed)
            self.repository.commit()
        except Exception:
            self.repository.rollback()
            raise

        return {"query": query, "removed": removed}

    def focal_presence(
        self,
        set_id: str,
        *,
        selected_file_id: str | None,
        fasta_file_ids: list[str] | None = None,
    ) -> list[EntryPresence]:
        """
        Presence of every focal entry, verified location-first.

        A saved location is checked directly against the current file by byte
        offset; only entries that fail that fall back to a scan. The whole-file
        hash is not a prerequisite for keeping a location that still resolves.

        `fasta_file_ids` bounds which files may confer presence: `None` is the
        whole project, a list is exactly those files. Run gating passes the
        analysis scope, so an entry cannot be counted as present because it
        happens to exist in some project file the run will not read.
        """
        self._require_focal_set(set_id)
        entries = self.repository.list_focal_entries(set_id)
        if not entries:
            return []

        scope = set(self._resolve_scope(fasta_file_ids))
        per_file: dict[str, dict[str, list[int]]] = {}
        unavailable: set[str] = set()

        for row in self.repository.list_fasta_files():
            if row.id not in scope:
                continue
            status = self.verifier.cheap_status(row)
            if not status.available:
                unavailable.add(row.id)
                continue
            self.repository.begin()
            try:
                per_file[row.id] = verify_entry_locations(self.repository, row, entries)
                self.repository.commit()
            except Exception:
                self.repository.rollback()
                raise

        return presence_for_entries(
            entries, per_file, selected_file_id=selected_file_id, unavailable_files=unavailable
        )

    # ------------------------------------------------------------------
    # Alignments and analysis scope
    # ------------------------------------------------------------------

    def load_alignment(self, file_id: str) -> CachedAlignment:
        """
        A verified, parsed alignment, reusing the session cache when possible.

        The cheap current-file check runs first every time; a missing or
        changed source can never be served from cache.
        """
        row = self._require_file(file_id)
        status = self.verifier.cheap_status(row)

        if status.state in (SourceState.MISSING, SourceState.UNREADABLE):
            raise ProjectError(
                "SOURCE_UNAVAILABLE",
                f"{row.display_name} is not available, so it cannot be analysed.",
                detail=row.source_path,
            )

        if status.state is SourceState.CURRENT:
            cached = self.alignments.get(file_id, row.index_revision, row.indexed_sha256)
            if cached is not None:
                return cached

        # The bytes are needed anyway, so this read also verifies the file.
        scan = self._stable_scan(row)
        if not scan.sequences:
            raise ProjectError("FASTA_EMPTY", f"{row.display_name} contains no sequences.")
        if scan.alignment_length is None:
            raise ProjectError(
                "FASTA_NOT_ALIGNED",
                f"{row.display_name} is not an aligned FASTA.",
            )

        if scan.sha256 != row.indexed_sha256:
            # Content genuinely changed: refresh the snapshot from this scan
            # rather than reading the file a second time.
            self.repository.begin()
            try:
                self.repository.replace_index(file_id, scan, bump_revision=True)
                entries = self._all_focal_entries()
                if entries:
                    verify_entry_locations(self.repository, row, entries, scan=scan)
                self.repository.commit()
            except Exception:
                self.repository.rollback()
                raise
            self.alignments.invalidate(file_id)
            row = self._require_file(file_id)

        self.verifier.note_scan(file_id, scan, matched_index=True)
        return self.alignments.put(file_id, row.index_revision, scan)

    def build_scope(self, fasta_file_ids: list[str]) -> AnalysisScope:
        """
        Combine verified alignments in memory for a multi-file run.

        No temporary combined FASTA is written and nothing is re-parsed: the
        cached dictionaries are merged directly.
        """
        problems: list[ScopeProblem] = []
        combined: dict[str, str] = {}
        owner: dict[str, str] = {}
        lengths: dict[str, int] = {}

        for file_id in fasta_file_ids:
            row = self.repository.get_fasta_file(file_id)
            if row is None:
                problems.append(ScopeProblem("UNKNOWN_FILE", "A selected file is not in this project."))
                continue
            try:
                alignment = self.load_alignment(file_id)
            except ProjectError as error:
                problems.append(ScopeProblem(error.code, error.message, error.detail))
                continue

            if alignment.duplicate_header_count > 0:
                # scan_fasta keeps every record, but `sequences` is a dict, so
                # a repeated header collapses to one entry — and the scientific
                # core is header-keyed, so it would analyse whichever record
                # won and never mention the others. Choosing one silently,
                # renaming, or dropping records would all invent data, so the
                # run is refused instead. The file stays linked and indexable;
                # it just cannot take part until a record-id core exists.
                kept = len(alignment.sequences)
                problems.append(
                    ScopeProblem(
                        "DUPLICATE_HEADER_IN_FILE",
                        f"{row.display_name} contains repeated FASTA headers, so it cannot be "
                        "analysed without silently discarding sequences.",
                        detail=(
                            f"{row.display_name}: {alignment.duplicate_header_count} duplicate "
                            f"record(s) across {kept} distinct header(s)"
                        ),
                    )
                )
                continue

            if alignment.alignment_length is not None:
                lengths[row.display_name] = alignment.alignment_length

            for header, sequence in alignment.sequences.items():
                if header in combined:
                    # The scientific core is keyed by header; silently letting
                    # one file overwrite another's record would change results
                    # invisibly. Renaming would be inventing data.
                    problems.append(
                        ScopeProblem(
                            "DUPLICATE_HEADER_ACROSS_FILES",
                            f"The header {header!r} appears in more than one selected file.",
                            detail=f"{owner[header]} and {row.display_name}",
                        )
                    )
                    continue
                combined[header] = sequence
                owner[header] = row.display_name

        distinct = set(lengths.values())
        if len(distinct) > 1:
            summary = ", ".join(f"{name}: {length}" for name, length in sorted(lengths.items()))
            problems.append(
                ScopeProblem(
                    "INCOMPATIBLE_ALIGNMENT_LENGTHS",
                    "The selected files do not share one alignment length.",
                    detail=summary,
                )
            )

        return AnalysisScope(
            fasta_file_ids=list(fasta_file_ids),
            sequences=combined,
            alignment_length=next(iter(distinct)) if len(distinct) == 1 else None,
            problems=problems,
        )

    def require_focal_completeness(
        self, set_id: str, fasta_file_ids: list[str], *, single_file: bool
    ) -> list[EntryPresence]:
        """
        Every explicit focal entry must exist inside the ACTUAL analysis scope.

        This is a scientific safety rule, not a UI convenience: an entry that
        resolves to nothing does not fail loudly downstream, it just shrinks
        the focal group, and a run on a smaller focal group is a different
        experiment reported under the same name.

        Presence is asked for with the included file ids, so "present" always
        means "present in something this run will actually read".
        """
        selected = fasta_file_ids[0] if single_file and len(fasta_file_ids) == 1 else None
        presence = self.focal_presence(
            set_id, selected_file_id=selected, fasta_file_ids=fasta_file_ids
        )

        unresolved = [
            entry for entry in presence if entry.state is not PresenceState.PRESENT_CURRENT
        ]
        if not unresolved:
            return presence

        detail = "; ".join(entry.header for entry in unresolved[:10])

        # "Cannot be established" is checked first and reported separately from
        # "is not there". Calling an unreadable file's contents absent would be
        # asserting something we did not observe.
        if any(entry.state is PresenceState.UNKNOWN for entry in unresolved):
            raise ProjectError(
                "FOCAL_PRESENCE_UNKNOWN",
                "Some focal entries cannot be located because a selected file is unavailable, "
                "so the run was not started.",
                detail=detail,
            )

        if single_file:
            raise ProjectError(
                "FOCAL_ENTRIES_NOT_IN_FILE",
                "Every focal entry must be present in the selected file before it can be "
                "analysed on its own.",
                detail=detail,
            )

        raise ProjectError(
            "FOCAL_ENTRIES_NOT_IN_SCOPE",
            "Every focal entry must be present in at least one of the selected files.",
            detail=detail,
        )

    def focal_selector(self, set_id: str) -> ExactHeaders:
        """The focal set as EXACT headers, for the scientific core."""
        self._require_focal_set(set_id)
        entries = self.repository.list_focal_entries(set_id)
        if not entries:
            raise ProjectError("FOCAL_EMPTY", "The focal set is empty.")
        return ExactHeaders([entry.header for entry in entries])

    # ------------------------------------------------------------------

    def _require_file(self, file_id: str) -> FastaFileRow:
        row = self.repository.get_fasta_file(file_id)
        if row is None:
            raise ProjectError("UNKNOWN_FILE", "That FASTA is not linked to this project.")
        return row

    def _all_focal_entries(self):
        entries = []
        for focal_set in self.repository.list_focal_sets():
            entries.extend(self.repository.list_focal_entries(focal_set.id))
        return entries

    # ------------------------------------------------------------------
    # Molecular Diagnosis over a project scope
    # ------------------------------------------------------------------

    def run_molecular_diagnosis(
        self,
        *,
        focal_set_id: str,
        fasta_file_ids: list[str],
        single_file: bool,
        options: dict[str, object],
        resume: dict[str, object] | None = None,
    ) -> dict[str, object]:
        """
        Run the pipeline over a project scope with EXACT focal membership.

        Any run — one file or twenty — is blocked unless every focal entry is
        present somewhere in the scope being analysed, because an absent entry
        silently shrinks the focal group and changes the science without
        telling anyone.
        """
        from molecular_diagnosis.focal import partition_headers
        from molecular_diagnosis.pipeline import run_pipeline_on_sequences
        from molecular_diagnosis.service.handlers import _dmc_to_payload
        from molecular_diagnosis.service.protocol import (
            resume_state_from_payload,
            resume_state_to_payload,
        )

        if not fasta_file_ids:
            raise ProjectError("NO_FILES_SELECTED", "Select at least one FASTA file to analyse.")

        # Resolve the scope before anything is judged against it, so an id that
        # is not in this project is reported as such rather than as a focal
        # entry mysteriously missing from the scope.
        for file_id in fasta_file_ids:
            self._require_file(file_id)

        selector = self.focal_selector(focal_set_id)

        if single_file and len(fasta_file_ids) != 1:
            raise ProjectError(
                "SCOPE_MISMATCH", "A single-file run must select exactly one file."
            )

        self.require_focal_completeness(
            focal_set_id, list(fasta_file_ids), single_file=single_file
        )

        scope = self.build_scope(fasta_file_ids)
        if not scope.ok:
            first = scope.problems[0]
            raise ProjectError(
                first.code,
                first.message,
                detail="; ".join(
                    filter(None, [problem.detail for problem in scope.problems])
                ) or None,
            )

        focal_headers, non_focal_headers = partition_headers(list(scope.sequences), selector)
        if not focal_headers:
            raise ProjectError(
                "FOCAL_NO_MATCH", "No sequence in the selected scope matches the focal set."
            )
        if not non_focal_headers:
            raise ProjectError(
                "FOCAL_MATCHES_EVERYTHING",
                "Every sequence in the selected scope is focal, so there is no contrast set.",
            )

        parsed_resume = resume_state_from_payload(resume)
        if parsed_resume is None:
            start_length, initial_combos, initial_tested = 1, None, None
        else:
            start_length, initial_combos, initial_tested = parsed_resume

        names = [
            row.display_name
            for row in (self.repository.get_fasta_file(fid) for fid in fasta_file_ids)
            if row is not None
        ]

        result = run_pipeline_on_sequences(
            sequences=scope.sequences,
            selectors=selector,
            focal_headers=focal_headers,
            non_focal_headers=non_focal_headers,
            alignment_length=scope.alignment_length or 0,
            source_label=" + ".join(names),
            output_dir=self.outputs_dir,
            include_ambiguous_dmc_bd=bool(options.get("giveBenefitOfDoubtToAmbiguousBases", False)),
            include_gappy_consensus_dmc_sites=bool(options.get("ignoreGaps", False)),
            min_combination_length=int(options.get("minCandidateSize", 1)),
            max_combination_length=int(options.get("maxCandidateSize", 2)),
            start_combination_length=start_length,
            initial_diagnostic_combinations=initial_combos,
            initial_combinations_tested_by_length=initial_tested,
        )

        dmc = result.dmc
        can_continue = dmc.stop_reason == "reached_maximum_length"
        return {
            "focalSetId": focal_set_id,
            "focalHeaders": list(selector),
            "fastaFileIds": list(fasta_file_ids),
            "sequenceCount": len(scope.sequences),
            "alignmentLength": scope.alignment_length,
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
            "resume": (
                resume_state_to_payload(
                    stopped_at_length=dmc.stopped_at_length,
                    diagnostic_combinations=dmc.diagnostic_combinations,
                    combinations_tested_by_length=dmc.combinations_tested_by_length,
                )
                if can_continue
                else None
            ),
        }
