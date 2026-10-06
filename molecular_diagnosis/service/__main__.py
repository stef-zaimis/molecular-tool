"""
Stdio entry point for the service.

    python -m molecular_diagnosis.service

Reads one JSON request per line on stdin, writes one JSON response per line on
stdout, and runs until stdin closes.

**stdout is reserved for the protocol.** Anything a handler or a third-party
library prints would otherwise corrupt the stream, so `sys.stdout` is swapped
for `sys.stderr` while handlers run and responses are written through a private
reference to the real stream. Diagnostics belong on stderr, which the parent
captures for logs.

Requests are served one at a time, in order: the request loop is deliberately
single-threaded. Concurrency belongs to the parent, which can run more than one
child if it ever needs to.

The ONE exception is control traffic. A dedicated thread reads stdin and hands
ordinary requests to the loop through a queue; a request named in
`OUT_OF_BAND_METHODS` (today only `control.cancelRun`) is answered by that
reader thread immediately instead. Without this a Stop could only be read after
the run it is meant to stop had already finished. Out-of-band handlers only
flip a flag the running request polls at its own checkpoints — they never touch
the project, the database or the scientific state, so there is still exactly
one thread doing real work.
"""

from __future__ import annotations

import io
import json
import queue
import sys
import threading
import time
import traceback
from collections.abc import Callable, Iterable

from molecular_diagnosis.service.diagnostics import (
    CHANNEL,
    environment_snapshot,
    log_line,
)
from molecular_diagnosis.service.errors import ServiceError, to_service_error
from molecular_diagnosis.service.handlers import dispatch
from molecular_diagnosis.service.projects import OUT_OF_BAND_METHODS
from molecular_diagnosis.service.protocol import (
    decode_request,
    encode_failure,
    encode_success,
)


def answer_out_of_band(line: str, emit: Callable[[str], None]) -> bool:
    """
    If `line` is an out-of-band control request, answer it now and return True.

    Anything else — including a line that is not JSON at all — returns False and
    is left for the ordinary loop, which reports malformed input exactly as it
    always has.
    """
    try:
        payload = json.loads(line)
    except json.JSONDecodeError:
        return False
    if not isinstance(payload, dict):
        return False
    method = payload.get("method")
    handler = OUT_OF_BAND_METHODS.get(method) if isinstance(method, str) else None
    if handler is None:
        return False

    request_id = payload.get("id") if isinstance(payload.get("id"), str) else "unknown"
    try:
        request = decode_request(line)
        emit(encode_success(request.id, handler(request.params)))
    except ServiceError as error:
        emit(encode_failure(request_id, error))
    except Exception as error:  # noqa: BLE001 - the reader thread must not die
        emit(encode_failure(request_id, to_service_error(error)))
    return True


def read_requests(
    lines: Iterable[str],
    pending: "queue.Queue[str | None]",
    emit: Callable[[str], None],
) -> None:
    """Reader thread body: answer control traffic, queue everything else."""
    try:
        for raw_line in lines:
            line = raw_line.strip()
            if not line:
                continue
            if answer_out_of_band(line, emit):
                continue
            pending.put(line)
    finally:
        # stdin closed (or the read failed): let the loop finish what it has.
        pending.put(None)


def main() -> int:
    # Private handle on the real stdout, taken before anything can replace it.
    channel = sys.stdout
    if isinstance(channel, io.TextIOWrapper):
        # Line buffering so the parent sees each response as it is produced.
        channel.reconfigure(line_buffering=True)

    # Two threads write now (the loop, and the reader answering control
    # requests), so each line goes out whole under a lock.
    write_lock = threading.Lock()

    def emit(line: str) -> None:
        with write_lock:
            channel.write(line + "\n")
            channel.flush()

    # From here on, a stray print() lands in the log rather than the protocol.
    sys.stdout = sys.stderr

    # Progress notifications are written through the SAME private handle as
    # responses, so they cannot interleave with a half-written response line
    # and cannot be affected by the stdout swap above.
    CHANNEL.install(emit)

    # One environment banner per process. This is the line to compare when two
    # machines behave differently: interpreter, versions, platform, paths.
    log_line("service.start", environment=environment_snapshot())

    emit(encode_success("ready", {"ready": True, "python": sys.version.split()[0]}))

    pending: queue.Queue[str | None] = queue.Queue()
    reader = threading.Thread(
        target=read_requests,
        args=(sys.stdin, pending, emit),
        name="stdin-reader",
        daemon=True,
    )
    reader.start()

    while True:
        line = pending.get()
        if line is None:
            break

        request_id = "unknown"
        started = time.monotonic()
        method = "unknown"
        try:
            request = decode_request(line)
            request_id = request.id
            method = request.method
            log_line("request.start", id=request_id, method=method)
            # Progress emitted anywhere inside this dispatch is tagged with
            # this request's id, which is what lets the parent route it.
            with CHANNEL.request(request_id):
                result = dispatch(request.method, request.params)
            log_line(
                "request.end",
                id=request_id,
                method=method,
                ok=True,
                durationMs=int((time.monotonic() - started) * 1000),
            )
            emit(encode_success(request_id, result))
        except ServiceError as error:
            log_line(
                "request.end",
                id=request_id,
                method=method,
                ok=False,
                code=str(error.code),
                message=error.message,
                durationMs=int((time.monotonic() - started) * 1000),
            )
            emit(encode_failure(request_id, error))
        except Exception as error:  # noqa: BLE001 - boundary must not die
            log_line(
                "request.crashed",
                id=request_id,
                method=method,
                error=f"{type(error).__name__}: {error}",
                traceback=traceback.format_exc(),
                durationMs=int((time.monotonic() - started) * 1000),
            )
            # The cause and traceback are preserved inside the ServiceError.
            emit(encode_failure(request_id, to_service_error(error)))

    log_line("service.stop", reason="stdin closed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
