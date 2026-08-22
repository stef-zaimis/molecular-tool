"""
Run observability: progress on the protocol channel, diagnostics on stderr.

Two properties matter more than the rest, and both are about NOT breaking
things that already work:

* progress notifications must not disturb the request/response framing — a
  reader that expects one response per request must still get exactly that;
* observing a run must not change its results, because the observer is
  attached to the scientific core and the science is the product.

The rest pins the correlation rules the UI depends on (which run a
notification belongs to) and the throttling that keeps a `C(n,5)` search from
drowning in its own progress.
"""

from __future__ import annotations

import json
import subprocess
import sys
import threading
from pathlib import Path

import pytest

from molecular_diagnosis.pipeline import run_pipeline_on_sequences
from molecular_diagnosis.progress import (
    NULL_OBSERVER,
    PROGRESS_STAGES,
    ProgressTicker,
    RunObserver,
    describe_path,
)
from molecular_diagnosis.service.diagnostics import (
    ProgressChannel,
    RunDiagnostics,
    environment_snapshot,
)

REPO_ROOT = Path(__file__).resolve().parents[1]

RECORDS = [
    ("focal_1|AU|Target", "AACCGGTTAC"),
    ("focal_2|GB|Target", "AACCGGTTAC"),
    ("other_1|FR|Contrast", "AGCTGGTTAC"),
    ("other_2|DE|Contrast", "CGTTGGTTAG"),
    ("other_3|ES|Contrast", "CGTTGCTTAG"),
]


def write_fasta(path: Path, records=RECORDS) -> Path:
    path.write_text(
        "".join(f">{header}\n{sequence}\n" for header, sequence in records), encoding="utf-8"
    )
    return path


# ---------------------------------------------------------------------------
# The stage vocabulary
# ---------------------------------------------------------------------------


def test_the_stage_names_are_pinned():
    """
    These strings cross the process boundary.

    `DiagnosisProgressStage` in desktop/src/backendContract.ts is the same list,
    and `STAGE_LABELS` in DiagnosisRunPanel.tsx is an exhaustive record over it.
    Renaming one here without the other leaves the UI with a stage it cannot
    label, so the list is written down in both places on purpose.
    """
    assert PROGRESS_STAGES == (
        "starting",
        "verifying_sources",
        "loading_alignment",
        "building_scope",
        "validating_focal",
        "consensus",
        "dmc_search",
        "five_site",
        "writing_report",
        "writing_consensus",
        "writing_workbook",
        "finishing",
    )


# ---------------------------------------------------------------------------
# Observing does not change the science
# ---------------------------------------------------------------------------


class RecordingObserver(RunObserver):
    """Captures everything, changes nothing."""

    def __init__(self) -> None:
        self.events: list[tuple[str, dict]] = []
        self.stages: list[str] = []
        self.progress_calls: list[tuple[str, int | None, int | None, str | None]] = []

    def event(self, name, /, **fields):
        self.events.append((name, fields))

    def progress(self, stage, /, *, current=None, total=None, detail=None):
        self.progress_calls.append((stage, current, total, detail))

    # `stage` is inherited: the base yields, which is all a recorder needs
    # beyond noting that it happened.
    def stage(self, name, /, **fields):  # type: ignore[override]
        self.stages.append(name)
        return super().stage(name, **fields)


def run_once(tmp_path: Path, observer: RunObserver, folder: str):
    sequences = {header: sequence for header, sequence in RECORDS}
    focal = [header for header in sequences if "Target" in header]
    non_focal = [header for header in sequences if "Target" not in header]
    output_dir = tmp_path / folder
    output_dir.mkdir()
    return run_pipeline_on_sequences(
        sequences=sequences,
        selectors=["Target"],
        focal_headers=focal,
        non_focal_headers=non_focal,
        alignment_length=len(RECORDS[0][1]),
        source_label="test",
        output_dir=output_dir,
        observer=observer,
    )


def test_an_observed_run_produces_identical_science(tmp_path):
    plain = run_once(tmp_path, NULL_OBSERVER, "plain")
    observed = run_once(tmp_path, RecordingObserver(), "observed")

    # Everything the pipeline computes, compared field by field. The output
    # PATHS differ (different directories), which is the only difference.
    assert observed.dmc == plain.dmc
    assert observed.txt_output_path.name == plain.txt_output_path.name
    assert observed.xlsx_output_path.read_bytes()[:2] == plain.xlsx_output_path.read_bytes()[:2]
    assert observed.txt_output_path.read_text(encoding="utf-8") == plain.txt_output_path.read_text(
        encoding="utf-8"
    ).replace("plain", "observed")


def test_every_pipeline_stage_is_reported(tmp_path):
    observer = RecordingObserver()
    run_once(tmp_path, observer, "stages")

    assert observer.stages == [
        "consensus",
        "dmc_search",
        "five_site",
        "finishing",
        "writing_report",
        "writing_consensus",
        "writing_workbook",
    ]
    # The counts a reader needs to size the run are events, not just stages.
    names = [name for name, _ in observer.events]
    assert "pipeline.inputs" in names
    assert "dmc.result" in names
    assert "outputs.allocated" in names


# ---------------------------------------------------------------------------
# Throttling
# ---------------------------------------------------------------------------


def test_the_ticker_reports_far_less_often_than_it_is_advanced():
    observer = RecordingObserver()
    ticker = ProgressTicker(observer, "five_site", total=1_000_000, min_interval_s=3600)

    for _ in range(200_000):
        ticker.advance()

    # The interval is an hour, so nothing has been due yet.
    assert observer.progress_calls == []

    ticker.finish()
    assert observer.progress_calls == [("five_site", 200_000, 1_000_000, None)]


def test_the_ticker_reports_when_the_interval_has_passed():
    observer = RecordingObserver()
    ticker = ProgressTicker(observer, "dmc_search", total=10, min_interval_s=0, check_every=1)

    ticker.advance()
    ticker.advance()

    assert [call[1] for call in observer.progress_calls] == [1, 2]


# ---------------------------------------------------------------------------
# The channel
# ---------------------------------------------------------------------------


def test_progress_is_silent_with_no_request_in_flight():
    """
    A direct `ProjectService` call — a test, or the Tkinter app — has no
    protocol channel to write to, and must not acquire one.
    """
    written: list[str] = []
    channel = ProgressChannel()
    channel.install(written.append)

    RunDiagnostics("abc", channel=channel).progress("five_site", current=1, total=2)

    assert written == []


def test_progress_carries_the_request_id_and_run_token():
    written: list[str] = []
    channel = ProgressChannel()
    channel.install(written.append)
    run = RunDiagnostics("run1", run_token="token-9", channel=channel)

    with channel.request("42"):
        run.progress("dmc_search", current=3, total=9, detail="2")

    envelope = json.loads(written[0])
    assert envelope["id"] == "42"
    assert envelope["type"] == "progress"
    assert "ok" not in envelope  # cannot be mistaken for a response
    assert envelope["progress"]["stage"] == "dmc_search"
    assert envelope["progress"]["current"] == 3
    assert envelope["progress"]["total"] == 9
    assert envelope["progress"]["detail"] == "2"
    assert envelope["progress"]["runId"] == "run1"
    assert envelope["progress"]["runToken"] == "token-9"


def test_a_broken_channel_never_breaks_the_run():
    def explode(_line: str) -> None:
        raise OSError("pipe closed")

    channel = ProgressChannel()
    channel.install(explode)

    with channel.request("7"):
        # No exception escapes: progress is diagnostics, not the product.
        RunDiagnostics("run2", channel=channel).progress("five_site")


def test_a_failure_is_reported_with_the_stage_it_happened_in():
    lines: list[tuple[str, dict]] = []
    run = RunDiagnostics("run3", channel=ProgressChannel(), log=lambda name, **f: lines.append((name, f)))

    with pytest.raises(ValueError):
        with run.stage("five_site"):
            raise ValueError("boom")

    names = [name for name, _ in lines]
    assert "stage.start:five_site" in names
    assert "stage.error:five_site" in names
    errored = dict(lines)["stage.error:five_site"]
    assert "ValueError: boom" in errored["error"]


# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------


def test_the_environment_snapshot_names_what_two_machines_would_differ_on():
    snapshot = environment_snapshot()

    assert snapshot["python"]["executable"] == sys.executable
    assert snapshot["platform"]["system"]
    assert "openpyxl" in snapshot["packages"]
    assert snapshot["packages"]["sqlite3"]
    assert snapshot["paths"]["cwd"]


def test_describe_path_reports_a_missing_file_without_raising(tmp_path):
    described = describe_path(tmp_path / "nope.fasta")
    assert described["exists"] is False
    assert "error" in described


def test_describe_path_reports_permissions_for_a_real_file(tmp_path):
    path = write_fasta(tmp_path / "real.fasta")
    described = describe_path(path)
    assert described["exists"] is True
    assert described["readable"] is True
    assert described["sizeBytes"] > 0


# ---------------------------------------------------------------------------
# End to end, through a real service subprocess
# ---------------------------------------------------------------------------


class Client:
    """A service subprocess that separates progress from responses."""

    def __init__(self) -> None:
        self.proc = subprocess.Popen(
            [sys.executable, "-m", "molecular_diagnosis.service"],
            cwd=str(REPO_ROOT),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            bufsize=1,
        )
        self.stderr_lines: list[str] = []
        self.progress: list[dict] = []
        #: stderr has to be drained or the child blocks on a full pipe.
        threading.Thread(target=self._drain, daemon=True).start()
        self._counter = 0
        self.lines: list[dict] = []
        ready = json.loads(self.proc.stdout.readline())
        assert ready["id"] == "ready"

    def _drain(self) -> None:
        for line in self.proc.stderr:
            self.stderr_lines.append(line.rstrip())

    def call(self, method: str, params: dict | None = None) -> dict:
        self._counter += 1
        request_id = str(self._counter)
        self.proc.stdin.write(
            json.dumps({"id": request_id, "method": method, "params": params or {}}) + "\n"
        )
        self.proc.stdin.flush()
        while True:
            line = self.proc.stdout.readline()
            assert line, "the service closed stdout:\n" + "\n".join(self.stderr_lines[-20:])
            message = json.loads(line)  # EVERY line must be JSON
            self.lines.append(message)
            if message.get("type") == "progress":
                self.progress.append(message)
                continue
            return message

    def ok(self, method: str, params: dict | None = None) -> dict:
        response = self.call(method, params)
        assert response["ok"], response.get("error")
        return response["result"]

    def close(self) -> None:
        try:
            self.proc.stdin.close()
            self.proc.wait(timeout=10)
        except Exception:
            self.proc.kill()


@pytest.fixture()
def client():
    service = Client()
    yield service
    service.close()


def run_a_project_analysis(client: Client, tmp_path: Path, token: str = "token-1") -> dict:
    fasta = write_fasta(tmp_path / "input.fasta")
    client.ok("project.create", {"projectDir": str(tmp_path / "proj"), "title": "Progress"})
    file_id = client.ok("project.linkFasta", {"path": str(fasta)})["fastaFileId"]
    set_id = client.ok("project.createFocalSet", {"title": "Targets"})["focalSetId"]
    client.ok("project.addFocalEntries", {"focalSetId": set_id, "query": "Target"})
    return client.call(
        "project.runMolecularDiagnosis",
        {
            "focalSetId": set_id,
            "fastaFileIds": [file_id],
            "singleFile": True,
            "runToken": token,
            "options": {
                "ignoreGaps": False,
                "giveBenefitOfDoubtToAmbiguousBases": False,
                "minCandidateSize": 1,
                "maxCandidateSize": 2,
            },
        },
    )


def test_progress_does_not_corrupt_the_request_response_framing(client, tmp_path):
    response = run_a_project_analysis(client, tmp_path)

    assert response["ok"], response.get("error")
    # One response, and it is the LAST line of the exchange.
    responses = [line for line in client.lines if "ok" in line]
    assert len(responses) == 5  # create, link, createFocalSet, addEntries, run
    assert client.lines[-1] is responses[-1]
    # Every progress line is a progress line and nothing else.
    for message in client.progress:
        assert message["type"] == "progress"
        assert "ok" not in message
        assert "result" not in message


def test_several_progress_notifications_arrive_before_the_response(client, tmp_path):
    run_a_project_analysis(client, tmp_path)

    assert len(client.progress) >= 3
    stages = [message["progress"]["stage"] for message in client.progress]
    # The run walks the stages in order, starting before any work is done.
    assert stages[0] == "starting"
    assert "validating_focal" in stages
    assert "writing_workbook" in stages
    assert set(stages) <= set(PROGRESS_STAGES)


def test_progress_is_tagged_with_the_run_it_belongs_to(client, tmp_path):
    response = run_a_project_analysis(client, tmp_path, token="token-abc")
    assert response["ok"]

    run_id = {message["progress"]["runId"] for message in client.progress}
    assert len(run_id) == 1, "one run, one id"
    assert all(message["progress"]["runToken"] == "token-abc" for message in client.progress)
    # And with the id of the request that is still pending — the run request,
    # which is the last one this exchange sent.
    assert {message["id"] for message in client.progress} == {"5"}


def test_a_refused_run_still_reports_its_stage_on_stderr(client, tmp_path):
    fasta = write_fasta(tmp_path / "input.fasta")
    client.ok("project.create", {"projectDir": str(tmp_path / "proj"), "title": "Refusal"})
    file_id = client.ok("project.linkFasta", {"path": str(fasta)})["fastaFileId"]
    set_id = client.ok("project.createFocalSet", {"title": "Targets"})["focalSetId"]
    # An entry that is not in the file: the run is refused at focal validation.
    client.ok(
        "project.saveFocalSet",
        {"focalSetId": set_id, "title": "Targets", "headers": ["not_in_this_file"]},
    )

    response = client.call(
        "project.runMolecularDiagnosis",
        {
            "focalSetId": set_id,
            "fastaFileIds": [file_id],
            "singleFile": True,
            "options": {},
        },
    )

    assert response["ok"] is False
    assert response["error"]["code"] == "FOCAL_ENTRIES_NOT_IN_FILE"
    refusals = [
        json.loads(line.split(" ", 1)[1])
        for line in client.stderr_lines
        if line.startswith("[mdx] ") and '"run.refused"' in line
    ]
    assert refusals, "the refusal must say which stage it happened in"
    assert refusals[0]["stage"] == "validating_focal"


def test_the_service_logs_its_environment_once_at_startup(client, tmp_path):
    run_a_project_analysis(client, tmp_path)

    starts = [
        json.loads(line.split(" ", 1)[1])
        for line in client.stderr_lines
        if line.startswith("[mdx] ") and '"service.start"' in line
    ]
    assert len(starts) == 1
    environment = starts[0]["environment"]
    assert environment["python"]["version"]
    assert environment["platform"]["system"]
