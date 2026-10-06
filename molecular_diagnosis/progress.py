"""
Observation hooks for long-running work.

The scientific modules do not know what a service, a socket or a UI is, and
must not learn. What they do know is when they start an expensive stage and how
far through it they are — so they accept an OBSERVER and tell it, and something
further out decides whether that becomes a log line, a progress event, or
nothing at all.

Three rules keep this honest:

1. **The default observes nothing.** `NULL_OBSERVER` is the parameter default
   everywhere, so every existing caller and every test keeps its exact previous
   behaviour and cost.
2. **Observation never changes results.** Nothing here is read back by the
   science; an observer that raises would be a bug in the observer, but the
   loops that report are written so a callback cannot alter their state.
3. **Reporting is throttled, not per-iteration.** A `C(n,5)` search runs
   millions of iterations. `ProgressTicker` reports on a wall-clock interval,
   and only checks the clock every few thousand iterations, so the cost of
   being observable stays in the noise.
"""

from __future__ import annotations

import os
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

__all__ = [
    "NULL_OBSERVER",
    "describe_path",
    "PROGRESS_STAGES",
    "ProgressTicker",
    "RunObserver",
    "STAGE_BUILDING_SCOPE",
    "STAGE_CONSENSUS",
    "STAGE_DMC_SEARCH",
    "STAGE_FINISHING",
    "STAGE_FIVE_SITE",
    "STAGE_LOADING_ALIGNMENT",
    "STAGE_STARTING",
    "STAGE_VALIDATING_FOCAL",
    "STAGE_VERIFYING_SOURCES",
    "STAGE_WRITING_CONSENSUS",
    "STAGE_WRITING_REPORT",
    "STAGE_WRITING_WORKBOOK",
]

# ---------------------------------------------------------------------------
# Stage names
# ---------------------------------------------------------------------------
#
# These strings cross the process boundary and are matched by name in the
# renderer (`DiagnosisProgressStage` in desktop/src/backendContract.ts). The
# TypeScript side maps every one of them to a label through an exhaustive
# record, so adding a stage here without adding it there fails to compile.
# `tests/test_progress.py` pins the list so a rename is deliberate.

STAGE_STARTING = "starting"
STAGE_VERIFYING_SOURCES = "verifying_sources"
STAGE_LOADING_ALIGNMENT = "loading_alignment"
STAGE_BUILDING_SCOPE = "building_scope"
STAGE_VALIDATING_FOCAL = "validating_focal"
STAGE_CONSENSUS = "consensus"
STAGE_DMC_SEARCH = "dmc_search"
STAGE_FIVE_SITE = "five_site"
STAGE_WRITING_REPORT = "writing_report"
STAGE_WRITING_CONSENSUS = "writing_consensus"
STAGE_WRITING_WORKBOOK = "writing_workbook"
STAGE_FINISHING = "finishing"

PROGRESS_STAGES: tuple[str, ...] = (
    STAGE_STARTING,
    STAGE_VERIFYING_SOURCES,
    STAGE_LOADING_ALIGNMENT,
    STAGE_BUILDING_SCOPE,
    STAGE_VALIDATING_FOCAL,
    STAGE_CONSENSUS,
    STAGE_DMC_SEARCH,
    STAGE_FIVE_SITE,
    STAGE_WRITING_REPORT,
    STAGE_WRITING_CONSENSUS,
    STAGE_WRITING_WORKBOOK,
    STAGE_FINISHING,
)


class RunObserver:
    """
    Somewhere for a long operation to say what it is doing.

    The base implementation does nothing, which is also the default parameter
    everywhere it is accepted. Subclasses turn these calls into diagnostics
    (developer-facing, stderr) and progress (user-facing, protocol channel);
    the two are deliberately separate methods because they have different
    audiences and different budgets.
    """

    # -- developer diagnostics ------------------------------------------

    def event(self, name: str, /, **fields: Any) -> None:
        """A point-in-time fact worth a log line."""

    @contextmanager
    def stage(self, name: str, /, **fields: Any) -> Iterator[None]:
        """
        Bracket a phase of work with a start and an end (and its duration).

        The base implementation still yields, so wrapping code in it costs
        nothing measurable when nobody is observing.
        """
        yield

    # -- user-facing progress -------------------------------------------

    def progress(
        self,
        stage: str,
        /,
        *,
        current: int | None = None,
        total: int | None = None,
        detail: str | None = None,
    ) -> None:
        """One live progress update. Must be cheap; may be dropped."""


NULL_OBSERVER = RunObserver()


class ProgressTicker:
    """
    Throttled progress for a loop with millions of iterations.

    `advance()` is called once per iteration and is nearly free: it increments a
    counter and only looks at the clock every `check_every` iterations. A report
    is emitted at most every `min_interval_s`, so a five-minute search produces
    a few hundred updates rather than a few million.
    """

    __slots__ = (
        "_observer",
        "_stage",
        "_total",
        "_detail",
        "_min_interval",
        "_check_every",
        "_count",
        "_next_check",
        "_last_report",
        "_started",
    )

    def __init__(
        self,
        observer: RunObserver,
        stage: str,
        *,
        total: int | None = None,
        detail: str | None = None,
        min_interval_s: float = 1.0,
        check_every: int = 2048,
    ) -> None:
        self._observer = observer
        self._stage = stage
        self._total = total
        self._detail = detail
        self._min_interval = min_interval_s
        self._check_every = max(1, check_every)
        self._count = 0
        self._next_check = self._check_every
        self._started = time.monotonic()
        self._last_report = self._started

    @property
    def count(self) -> int:
        return self._count

    def advance(self, step: int = 1) -> None:
        self._count += step
        if self._count < self._next_check:
            return
        self._next_check = self._count + self._check_every
        now = time.monotonic()
        if now - self._last_report < self._min_interval:
            return
        self._last_report = now
        self._emit()

    def finish(self) -> None:
        """Report the final count once, whatever the interval says."""
        self._emit()

    def _emit(self) -> None:
        self._observer.progress(
            self._stage,
            current=self._count,
            total=self._total,
            detail=self._detail,
        )


def describe_path(path: str | Path) -> dict[str, Any]:
    """
    What the filesystem says about one path, without reading it.

    Used for FASTA sources, the project directory and the outputs directory:
    a permission or a mount difference shows up here long before it shows up
    as a mysterious failure two stages later.
    """
    target = Path(path)
    info: dict[str, Any] = {"path": str(target)}
    try:
        stat = target.stat()
    except OSError as error:
        info["exists"] = False
        info["error"] = f"{type(error).__name__}: {error}"
        return info

    info.update(
        {
            "exists": True,
            "sizeBytes": stat.st_size,
            "mtimeNs": stat.st_mtime_ns,
            "mode": oct(stat.st_mode & 0o7777),
            "readable": os.access(target, os.R_OK),
            "writable": os.access(target, os.W_OK),
        }
    )
    return info
