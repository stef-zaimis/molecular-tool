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

The loop is deliberately single-threaded and sequential: one request is served
at a time. Concurrency belongs to the parent, which can run more than one child
if it ever needs to.
"""

from __future__ import annotations

import io
import sys

from molecular_diagnosis.service.errors import ServiceError, to_service_error
from molecular_diagnosis.service.handlers import dispatch
from molecular_diagnosis.service.protocol import (
    decode_request,
    encode_failure,
    encode_success,
)


def main() -> int:
    # Private handle on the real stdout, taken before anything can replace it.
    channel = sys.stdout
    if isinstance(channel, io.TextIOWrapper):
        # Line buffering so the parent sees each response as it is produced.
        channel.reconfigure(line_buffering=True)

    def emit(line: str) -> None:
        channel.write(line + "\n")
        channel.flush()

    # From here on, a stray print() lands in the log rather than the protocol.
    sys.stdout = sys.stderr

    emit(encode_success("ready", {"ready": True, "python": sys.version.split()[0]}))

    for raw_line in sys.stdin:
        line = raw_line.strip()
        if not line:
            continue

        request_id = "unknown"
        try:
            request = decode_request(line)
            request_id = request.id
            result = dispatch(request.method, request.params)
            emit(encode_success(request_id, result))
        except ServiceError as error:
            emit(encode_failure(request_id, error))
        except Exception as error:  # noqa: BLE001 - boundary must not die
            # The cause and traceback are preserved inside the ServiceError.
            emit(encode_failure(request_id, to_service_error(error)))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
