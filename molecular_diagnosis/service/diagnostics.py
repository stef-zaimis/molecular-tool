"""
Run diagnostics and live progress for the service.

Two audiences, two channels, deliberately not mixed:

* **Developers** get structured lines on **stderr**. One JSON object per line
  behind a `[mdx]` tag, carrying the run id, the elapsed time and the stage, so
  a log from a Windows machine and a log from a Linux machine can be diffed
  against each other. stderr is captured by the Electron parent and written to
  its log file; nothing here ever touches stdout directly.

* **Users** get progress envelopes on the **protocol channel**, correlated with
  the pending request id, so the UI can say what a long run is doing while it
  is still running. These are deliberately sparse and carry numbers, not
  prose — the renderer owns the wording.

The rule that makes this safe: stdout is the protocol. Diagnostics go to
stderr; progress goes through `ProgressChannel`, which writes a properly framed
JSON envelope through the same writer the responses use, and is a no-op when
no request is in flight (so a direct `ProjectService` call in a test emits
nothing).
"""

from __future__ import annotations

import json
import os
import platform
import sqlite3
import sys
import time
import traceback
import uuid
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from molecular_diagnosis.progress import RunObserver, describe_path

__all__ = [
    "CHANNEL",
    "ProgressChannel",
    "RunDiagnostics",
    # Re-exported: it lives in the observation layer so modules below the
    # service boundary can use it without importing the service package.
    "describe_path",
    "environment_snapshot",
    "log_line",
    "new_run_id",
]

#: Prefix every diagnostic line carries, so a log can be grepped for ours.
TAG = "[mdx]"


def log_line(event: str, **fields: Any) -> None:
    """
    One structured diagnostic line on stderr.

    Flushed immediately: a log that is still in a buffer when the process is
    killed describes nothing, and "where did it stop" is the whole point.
    """
    payload = {"ts": time.strftime("%Y-%m-%dT%H:%M:%S"), "event": event, **fields}
    try:
        line = json.dumps(payload, default=str, ensure_ascii=False)
    except (TypeError, ValueError):  # pragma: no cover - defensive
        line = json.dumps({"event": event, "unserialisable": True})
    print(f"{TAG} {line}", file=sys.stderr, flush=True)


def new_run_id() -> str:
    """Short, human-quotable, unique enough to correlate one run's lines."""
    return uuid.uuid4().hex[:8]


# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------


def _package_version(name: str) -> str | None:
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:  # pragma: no cover - Python < 3.8 only
        return None
    try:
        return version(name)
    except PackageNotFoundError:
        return None
    except Exception:  # pragma: no cover - defensive
        return None


def environment_snapshot() -> dict[str, Any]:
    """
    Everything needed to compare two machines that behave differently.

    Deliberately includes the boring things: which interpreter actually got
    picked, where it thinks the repo is, and what it is running as. Most
    "works on my machine" differences in this app are one of those three, not
    the science.
    """
    return {
        "python": {
            "executable": sys.executable,
            "version": sys.version.split()[0],
            "fullVersion": sys.version.replace("\n", " "),
            "implementation": platform.python_implementation(),
            "pid": os.getpid(),
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "platform": platform.platform(),
            "filesystemEncoding": sys.getfilesystemencoding(),
            "defaultEncoding": sys.getdefaultencoding(),
        },
        "packages": {
            "openpyxl": _package_version("openpyxl"),
            "pillow": _package_version("pillow"),
            # The SQLite LIBRARY version. `sqlite3.version` is the module's
            # own and is deprecated in 3.12, so it is deliberately not here.
            "sqlite3": sqlite3.sqlite_version,
        },
        "paths": {
            "cwd": os.getcwd(),
            "packageRoot": str(Path(__file__).resolve().parents[1]),
            "repoRoot": str(Path(__file__).resolve().parents[2]),
            "sysPath0": sys.path[0] if sys.path else None,
            "tempDir": os.environ.get("TMPDIR") or os.environ.get("TEMP") or "/tmp",
        },
        "env": {
            "MOLECULAR_TOOL_PYTHON": os.environ.get("MOLECULAR_TOOL_PYTHON"),
            "MOLECULAR_TOOL_ROOT": os.environ.get("MOLECULAR_TOOL_ROOT"),
            "PYTHONIOENCODING": os.environ.get("PYTHONIOENCODING"),
        },
    }


# ---------------------------------------------------------------------------
# Progress channel
# ---------------------------------------------------------------------------


class ProgressChannel:
    """
    Where progress envelopes go while a request is in flight.

    The stdio loop installs a writer and sets the current request id before
    dispatching; everything else in the process asks this object. With no
    writer installed (a test calling `ProjectService` directly, or the Tkinter
    app) every send is a no-op, so nothing outside the service pays for it.
    """

    def __init__(self) -> None:
        self._write: Callable[[str], None] | None = None
        self._request_id: str | None = None

    def install(self, write: Callable[[str], None] | None) -> None:
        self._write = write

    @contextmanager
    def request(self, request_id: str) -> Iterator[None]:
        previous = self._request_id
        self._request_id = request_id
        try:
            yield
        finally:
            self._request_id = previous

    @property
    def active(self) -> bool:
        return self._write is not None and self._request_id is not None

    def send(self, payload: dict[str, Any]) -> None:
        """
        Emit one progress envelope, correlated with the pending request.

        Never raises: progress is diagnostics, and a broken pipe or an
        unserialisable field must not take down a scientific run that is
        otherwise fine.
        """
        write = self._write
        request_id = self._request_id
        if write is None or request_id is None:
            return
        try:
            from molecular_diagnosis.service.protocol import encode_progress

            write(encode_progress(request_id, payload))
        except Exception:  # noqa: BLE001 - progress must never break a run
            pass


#: Process-wide channel. Installed by `service.__main__`.
CHANNEL = ProgressChannel()


# ---------------------------------------------------------------------------
# The observer the run path uses
# ---------------------------------------------------------------------------


class RunDiagnostics(RunObserver):
    """
    One analysis run's diagnostics and progress.

    Stage brackets and events become stderr lines; `progress` becomes a
    protocol envelope AND (throttled) a stderr line, because the same numbers
    answer both "what should the user see" and "where did it stop".

    `last_stage` is kept so an unexpected exception can be reported with the
    stage it happened in, which is the single most useful fact in a hang
    report.
    """

    def __init__(
        self,
        run_id: str,
        *,
        run_token: str | None = None,
        channel: ProgressChannel | None = None,
        log: Callable[..., None] = log_line,
    ) -> None:
        self.run_id = run_id
        self.run_token = run_token
        self.started = time.monotonic()
        self.last_stage: str = "starting"
        self._channel = channel if channel is not None else CHANNEL
        self._log = log
        self._stack: list[str] = []
        self._last_progress_log = 0.0

    # -- helpers ---------------------------------------------------------

    @property
    def elapsed_ms(self) -> int:
        return int((time.monotonic() - self.started) * 1000)

    def _emit(self, event: str, **fields: Any) -> None:
        self._log(
            event,
            run=self.run_id,
            elapsedMs=self.elapsed_ms,
            stage=self.last_stage,
            **fields,
        )

    # -- RunObserver -----------------------------------------------------

    def event(self, name: str, /, **fields: Any) -> None:
        self._emit(name, **fields)

    @contextmanager
    def stage(self, name: str, /, **fields: Any) -> Iterator[None]:
        previous = self.last_stage
        self.last_stage = name
        self._stack.append(name)
        started = time.monotonic()
        self._emit(f"stage.start:{name}", **fields)
        failed = False
        try:
            yield
        except Exception as error:  # noqa: BLE001 - re-raised after reporting
            failed = True
            self._emit(
                f"stage.error:{name}",
                durationMs=int((time.monotonic() - started) * 1000),
                error=f"{type(error).__name__}: {error}",
            )
            raise
        else:
            self._emit(
                f"stage.end:{name}",
                durationMs=int((time.monotonic() - started) * 1000),
            )
        finally:
            self._stack.pop()
            # On the way OUT of a failure, `last_stage` keeps the stage that
            # failed. Unwinding it back to the caller's stage would throw away
            # the one fact a hang or crash report is for.
            if not failed:
                self.last_stage = self._stack[-1] if self._stack else previous

    def progress(
        self,
        stage: str,
        /,
        *,
        current: int | None = None,
        total: int | None = None,
        detail: str | None = None,
    ) -> None:
        elapsed = self.elapsed_ms
        self._channel.send(
            {
                "kind": "diagnosisProgress",
                "runId": self.run_id,
                "runToken": self.run_token,
                "stage": stage,
                "current": current,
                "total": total,
                "detail": detail,
                "elapsedMs": elapsed,
            }
        )

        # The same numbers on stderr, but far more sparsely: the log is for
        # comparing two machines, not for animating a bar.
        now = time.monotonic()
        if now - self._last_progress_log >= 5.0:
            self._last_progress_log = now
            self._emit("progress", progressStage=stage, current=current, total=total)

    # -- failures --------------------------------------------------------

    def failure(self, error: BaseException) -> None:
        """Report an unexpected exception with the stage it happened in."""
        self._emit(
            "run.failed",
            failedStage=self.last_stage,
            error=f"{type(error).__name__}: {error}",
            traceback=traceback.format_exc(),
        )
