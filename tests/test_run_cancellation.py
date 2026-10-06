"""
Stopping a Molecular Diagnosis run, and continuing one safely.

Cancellation is cooperative: `control.cancelRun` marks a run token, and the
run's observer raises `RunCancelled` at the next checkpoint (a stage boundary
or a `ProgressTicker` clock check). Output writing is never interrupted; a stop
that lands there deletes the run's outputs afterwards. A cancelled run yields
no result, and the service stays usable.

Continuation: a resume state is bound to the inputs that produced it and is
refused against any other inputs, while "raise the maximum and continue" works.
"""

from __future__ import annotations

import json
import queue
import random
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

from molecular_diagnosis.progress import (
    STAGE_DMC_SEARCH,
    STAGE_FIVE_SITE,
    STAGE_WRITING_REPORT,
    STAGE_WRITING_WORKBOOK,
    ProgressTicker,
    RunCancelled,
    RunObserver,
)
from molecular_diagnosis.project.service import ProjectError, ProjectService
from molecular_diagnosis.service import projects as rpc
from molecular_diagnosis.service.__main__ import answer_out_of_band
from molecular_diagnosis.service.diagnostics import CancellationRegistry, RunDiagnostics
from molecular_diagnosis.service.errors import ServiceError

REPO = Path(__file__).resolve().parents[1]

# Every single site and every pair is matched by some contrast sequence; only
# the full 4-site combination is diagnostic. So max=2 stops at the maximum
# with nothing found, and continuing to 4 finds [0, 1, 2, 3].
CONTINUABLE = [
    ("F1_target", "AAAA"),
    ("F2_target", "AAAA"),
    ("C1_other", "CAAA"),
    ("C2_other", "ACAA"),
    ("C3_other", "AACA"),
    ("C4_other", "AAAC"),
]

# Enough candidate sites for a five-site search to run.
SEVEN_SITES = [
    ("F1_target", "AAAAAAA"),
    ("F2_target", "AAAAAAA"),
    ("C1_other", "CCCCCCC"),
    ("C2_other", "GGGGGGG"),
]

FOCAL = ["F1_target", "F2_target"]


def write_fasta(path, records):
    path.write_text(
        "".join(f">{header}\n{sequence}\n" for header, sequence in records), encoding="utf-8"
    )
    return path


@pytest.fixture()
def project(tmp_path):
    service = ProjectService(tmp_path / "proj")
    service.open_project("Cancel")
    yield service
    service.close()


def link(project, tmp_path, name, records):
    row, _status = project.link_fasta(str(write_fasta(tmp_path / name, records)))
    return row.id


def outputs(project) -> list[str]:
    return sorted(path.name for path in project.outputs_dir.iterdir())


class StopAt(RunDiagnostics):
    """Requests a stop the moment a given stage is entered or reported."""

    def __init__(self, stage: str, registry: CancellationRegistry) -> None:
        super().__init__("test-run", run_token="tok", cancellations=registry, log=lambda *a, **k: None)
        self._stop_stage = stage
        self._registry = registry

    def _maybe_stop(self, stage: str) -> None:
        if stage == self._stop_stage:
            self._registry.request("tok")

    def stage(self, name, /, **fields):
        self._maybe_stop(name)
        return super().stage(name, **fields)

    def progress(self, stage, /, **kwargs):
        self._maybe_stop(stage)
        super().progress(stage, **kwargs)


def run(project, set_id, file_ids, *, observer=None, resume=None, **options):
    opts = {"minCandidateSize": 1, "maxCandidateSize": 2, **options}
    kwargs = {"observer": observer} if observer is not None else {}
    return project.run_molecular_diagnosis(
        focal_set_id=set_id,
        fasta_file_ids=list(file_ids),
        single_file=len(file_ids) == 1,
        options=opts,
        resume=resume,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# The checkpoint mechanism
# ---------------------------------------------------------------------------


def test_the_ticker_checkpoints_only_at_its_clock_checks():
    calls = []

    class Counting(RunObserver):
        def checkpoint(self):
            calls.append(1)

    ticker = ProgressTicker(Counting(), STAGE_DMC_SEARCH, check_every=100)
    for _ in range(1000):
        ticker.advance()
    assert len(calls) == 10


def test_the_default_observer_never_stops():
    observer = RunObserver()
    observer.checkpoint()
    assert observer.cancel_requested is False


def test_run_diagnostics_stops_only_its_own_token():
    registry = CancellationRegistry()
    mine = RunDiagnostics("r1", run_token="mine", cancellations=registry, log=lambda *a, **k: None)
    other = RunDiagnostics("r2", run_token="other", cancellations=registry, log=lambda *a, **k: None)
    registry.request("mine")
    with pytest.raises(RunCancelled):
        mine.checkpoint()
    other.checkpoint()  # unaffected


def test_stops_are_deferred_once_outputs_are_being_written():
    registry = CancellationRegistry()
    observer = RunDiagnostics("r", run_token="t", cancellations=registry, log=lambda *a, **k: None)
    with observer.stage("finishing"):
        registry.request("t")
        observer.checkpoint()  # deferred, not raised
        with observer.stage(STAGE_WRITING_REPORT):
            pass
    assert observer.cancel_requested is True


def test_the_registry_is_bounded_and_discardable():
    registry = CancellationRegistry()
    for index in range(CancellationRegistry.LIMIT + 50):
        registry.request(f"t{index}")
    assert not registry.is_requested("t0")
    assert registry.is_requested(f"t{CancellationRegistry.LIMIT + 49}")
    registry.discard(f"t{CancellationRegistry.LIMIT + 49}")
    assert not registry.is_requested(f"t{CancellationRegistry.LIMIT + 49}")


# ---------------------------------------------------------------------------
# Cancelling a project run
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("stage", [STAGE_DMC_SEARCH, STAGE_FIVE_SITE, "validating_focal"])
def test_a_stop_before_outputs_leaves_no_files(project, tmp_path, stage):
    file_id = link(project, tmp_path, "a.fasta", SEVEN_SITES)
    set_id = project.save_focal_set(focal_set_id=None, title="t", headers=FOCAL)["id"]

    with pytest.raises(RunCancelled):
        run(project, set_id, [file_id], observer=StopAt(stage, CancellationRegistry()))

    assert outputs(project) == []


@pytest.mark.parametrize("stage", [STAGE_WRITING_REPORT, STAGE_WRITING_WORKBOOK])
def test_a_stop_while_writing_finishes_writing_then_removes_the_outputs(project, tmp_path, stage):
    file_id = link(project, tmp_path, "a.fasta", SEVEN_SITES)
    set_id = project.save_focal_set(focal_set_id=None, title="t", headers=FOCAL)["id"]

    with pytest.raises(RunCancelled):
        run(project, set_id, [file_id], observer=StopAt(stage, CancellationRegistry()))

    assert outputs(project) == []


def test_a_cancelled_run_does_not_remove_earlier_outputs(project, tmp_path):
    file_id = link(project, tmp_path, "a.fasta", SEVEN_SITES)
    set_id = project.save_focal_set(focal_set_id=None, title="t", headers=FOCAL)["id"]
    run(project, set_id, [file_id])
    before = outputs(project)
    assert before

    with pytest.raises(RunCancelled):
        run(project, set_id, [file_id], observer=StopAt(STAGE_WRITING_WORKBOOK, CancellationRegistry()))

    assert outputs(project) == before


def test_the_project_runs_normally_after_a_cancel(project, tmp_path):
    file_id = link(project, tmp_path, "a.fasta", SEVEN_SITES)
    set_id = project.save_focal_set(focal_set_id=None, title="t", headers=FOCAL)["id"]
    baseline = run(project, set_id, [file_id])

    with pytest.raises(RunCancelled):
        run(project, set_id, [file_id], observer=StopAt(STAGE_DMC_SEARCH, CancellationRegistry()))

    again = run(project, set_id, [file_id])
    assert again["dmc"] == baseline["dmc"]
    # Saving and reading still work: no transaction was left open.
    project.save_focal_set(focal_set_id=set_id, title="renamed", headers=FOCAL)
    assert project.get_focal_set(set_id)["title"] == "renamed"


def test_an_unstopped_observer_changes_nothing(project, tmp_path):
    """A run that is never stopped gives the same science with or without a registry."""
    file_id = link(project, tmp_path, "a.fasta", SEVEN_SITES)
    set_id = project.save_focal_set(focal_set_id=None, title="t", headers=FOCAL)["id"]
    plain = run(project, set_id, [file_id])
    observed = run(
        project,
        set_id,
        [file_id],
        observer=RunDiagnostics(
            "r", run_token="never", cancellations=CancellationRegistry(), log=lambda *a, **k: None
        ),
    )
    assert observed["dmc"] == plain["dmc"]
    first = Path(plain["outputs"]["reportTxt"]).read_text(encoding="utf-8")
    second = Path(observed["outputs"]["reportTxt"]).read_text(encoding="utf-8")
    assert first == second


# ---------------------------------------------------------------------------
# The RPC boundary
# ---------------------------------------------------------------------------


@pytest.fixture()
def rpc_project(project, monkeypatch):
    monkeypatch.setitem(rpc._open, "current", project)
    yield project
    monkeypatch.delitem(rpc._open, "current", raising=False)


def test_a_cancelled_rpc_run_is_its_own_outcome(rpc_project, tmp_path):
    file_id = link(rpc_project, tmp_path, "a.fasta", SEVEN_SITES)
    set_id = rpc_project.save_focal_set(focal_set_id=None, title="t", headers=FOCAL)["id"]

    # Stop requested before the run is even picked up: it still stops.
    assert rpc.cancel_run({"runToken": "early"}) == {"runToken": "early", "cancelRequested": True}
    with pytest.raises(ServiceError) as caught:
        rpc.run_project_molecular_diagnosis(
            {"focalSetId": set_id, "fastaFileIds": [file_id], "runToken": "early"}
        )
    assert caught.value.code == "RUN_CANCELLED"
    assert caught.value.message == "Analysis stopped."
    # The token is spent; a later run with a fresh token is unaffected.
    assert not rpc.CANCELLATIONS.is_requested("early")
    result = rpc.run_project_molecular_diagnosis(
        {"focalSetId": set_id, "fastaFileIds": [file_id], "runToken": "later"}
    )
    assert result["outputs"]["reportTxt"]


def test_repeated_cancels_are_harmless(rpc_project):
    for _ in range(3):
        assert rpc.cancel_run({"runToken": "same"})["cancelRequested"] is True
    rpc.CANCELLATIONS.discard("same")


def test_cancel_needs_a_token():
    with pytest.raises(ServiceError) as caught:
        rpc.cancel_run({})
    assert caught.value.code == "INVALID_PARAMETER"


def test_out_of_band_lines_are_answered_immediately_and_others_are_not():
    sent: list[str] = []
    line = json.dumps({"id": "9", "method": "control.cancelRun", "params": {"runToken": "x"}})
    assert answer_out_of_band(line, sent.append) is True
    assert json.loads(sent[0]) == {
        "id": "9",
        "ok": True,
        "result": {"runToken": "x", "cancelRequested": True},
    }
    rpc.CANCELLATIONS.discard("x")

    assert answer_out_of_band('{"id": "1", "method": "ping"}', sent.append) is False
    assert answer_out_of_band("not json", sent.append) is False
    assert len(sent) == 1


# ---------------------------------------------------------------------------
# The real process: Stop while a run is in flight
# ---------------------------------------------------------------------------


class Service:
    def __init__(self) -> None:
        self.proc = subprocess.Popen(
            [sys.executable, "-u", "-m", "molecular_diagnosis.service"],
            cwd=REPO,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            encoding="utf-8",
        )
        self.lines: queue.Queue[dict] = queue.Queue()
        threading.Thread(target=self._pump, daemon=True).start()
        self.wait(lambda m: m.get("id") == "ready")

    def _pump(self) -> None:
        for line in self.proc.stdout:
            self.lines.put(json.loads(line))

    def send(self, request_id: str, method: str, params: dict) -> None:
        self.proc.stdin.write(json.dumps({"id": request_id, "method": method, "params": params}) + "\n")
        self.proc.stdin.flush()

    def wait(self, predicate, timeout: float = 60.0) -> dict:
        deadline = time.monotonic() + timeout
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError("the service did not answer in time")
            message = self.lines.get(timeout=remaining)
            if predicate(message):
                return message

    def response(self, request_id: str, timeout: float = 60.0) -> dict:
        return self.wait(lambda m: m.get("id") == request_id and "ok" in m, timeout)

    def close(self) -> None:
        try:
            self.proc.stdin.close()
            self.proc.wait(timeout=10)
        finally:
            if self.proc.poll() is None:
                self.proc.kill()


def test_the_real_service_stops_a_running_search_and_keeps_serving(tmp_path):
    rng = random.Random(7)
    length = 300
    focal = "".join(rng.choice("ACGT") for _ in range(length))
    records = [("F1_target", focal), ("F2_target", focal)] + [
        (f"C{i}_other", "".join(rng.choice("ACGT") for _ in range(length))) for i in range(60)
    ]
    fasta = write_fasta(tmp_path / "big.fasta", records)

    service = Service()
    try:
        service.send("1", "project.create", {"projectDir": str(tmp_path / "proj"), "title": "p"})
        assert service.response("1")["ok"]
        service.send("2", "project.linkFasta", {"path": str(fasta)})
        file_id = service.response("2")["result"]["fastaFileId"]
        service.send("3", "project.saveFocalSet", {"focalSetId": None, "title": "t", "headers": FOCAL})
        set_id = service.response("3")["result"]["focalSet"]["id"]

        # A 4-site search over ~300 candidate sites: far longer than this test.
        service.send(
            "4",
            "project.runMolecularDiagnosis",
            {
                "focalSetId": set_id,
                "fastaFileIds": [file_id],
                "options": {"minCandidateSize": 4, "maxCandidateSize": 4},
                "runToken": "tok-A",
            },
        )
        service.wait(
            lambda m: m.get("id") == "4"
            and m.get("type") == "progress"
            and m["progress"]["stage"] == STAGE_DMC_SEARCH
        )

        started = time.monotonic()
        service.send("5", "control.cancelRun", {"runToken": "tok-A"})
        service.send("6", "control.cancelRun", {"runToken": "tok-A"})  # a second click
        assert service.response("5", timeout=5)["result"]["cancelRequested"] is True
        assert service.response("6", timeout=5)["ok"] is True

        stopped = service.response("4", timeout=30)
        assert stopped["ok"] is False
        assert stopped["error"]["code"] == "RUN_CANCELLED"
        assert time.monotonic() - started < 30
        assert list((tmp_path / "proj" / "outputs").iterdir()) == []

        # Same process, same project: an ordinary request and a fresh run work.
        service.send("7", "project.listFocalSets", {})
        assert service.response("7")["ok"] is True
        service.send(
            "8",
            "project.runMolecularDiagnosis",
            {
                "focalSetId": set_id,
                "fastaFileIds": [file_id],
                "options": {"minCandidateSize": 1, "maxCandidateSize": 1},
                "runToken": "tok-B",
            },
        )
        assert service.response("8", timeout=120)["ok"] is True
    finally:
        service.close()


# ---------------------------------------------------------------------------
# Continuation safety
# ---------------------------------------------------------------------------


@pytest.fixture()
def stopped_at_max(project, tmp_path):
    """A run that stopped at max=2 with a resume state, and how it was made."""
    file_id = link(project, tmp_path, "a.fasta", CONTINUABLE)
    set_id = project.save_focal_set(focal_set_id=None, title="t", headers=FOCAL)["id"]
    first = run(project, set_id, [file_id], maxCandidateSize=2)
    assert first["canContinue"] is True
    assert first["resume"]["inputsFingerprint"]
    return file_id, set_id, first["resume"]


def test_raising_the_maximum_and_continuing_still_works(project, tmp_path, stopped_at_max):
    file_id, set_id, resume = stopped_at_max
    continued = run(project, set_id, [file_id], maxCandidateSize=4, resume=resume)
    fresh = run(project, set_id, [file_id], maxCandidateSize=4)

    assert continued["dmc"]["combinationsByLength"] == fresh["dmc"]["combinationsByLength"]
    assert continued["dmc"]["combinationsByLength"]["4"] == [[0, 1, 2, 3]]
    assert continued["dmc"]["startCombinationLength"] == 3


@pytest.mark.parametrize(
    "change",
    [
        {"ignoreGaps": True},
        {"giveBenefitOfDoubtToAmbiguousBases": True},
        {"minCandidateSize": 2},
    ],
)
def test_changed_options_refuse_the_continuation(project, stopped_at_max, change):
    file_id, set_id, resume = stopped_at_max
    with pytest.raises(ProjectError) as caught:
        run(project, set_id, [file_id], maxCandidateSize=4, resume=resume, **change)
    assert caught.value.code == "RESUME_INPUTS_CHANGED"


def test_another_focal_set_refuses_the_continuation(project, stopped_at_max):
    file_id, _set_id, resume = stopped_at_max
    other = project.save_focal_set(focal_set_id=None, title="o", headers=["F1_target"])["id"]
    with pytest.raises(ProjectError) as caught:
        run(project, other, [file_id], maxCandidateSize=4, resume=resume)
    assert caught.value.code == "RESUME_INPUTS_CHANGED"


def test_edited_focal_membership_refuses_the_continuation(project, stopped_at_max):
    file_id, set_id, resume = stopped_at_max
    project.save_focal_set(focal_set_id=set_id, title="t", headers=["F1_target"])
    with pytest.raises(ProjectError) as caught:
        run(project, set_id, [file_id], maxCandidateSize=4, resume=resume)
    assert caught.value.code == "RESUME_INPUTS_CHANGED"


def test_a_changed_comparison_set_refuses_the_continuation(project, stopped_at_max):
    file_id, set_id, resume = stopped_at_max
    project.save_focal_set(
        focal_set_id=set_id, title="t", headers=FOCAL, comparison_headers=["C1_other", "C2_other"]
    )
    with pytest.raises(ProjectError) as caught:
        run(project, set_id, [file_id], maxCandidateSize=4, resume=resume)
    assert caught.value.code == "RESUME_INPUTS_CHANGED"


def test_a_different_scope_refuses_the_continuation(project, tmp_path, stopped_at_max):
    _file_id, set_id, resume = stopped_at_max
    copy = link(project, tmp_path, "copy.fasta", CONTINUABLE)
    with pytest.raises(ProjectError) as caught:
        run(project, set_id, [copy], maxCandidateSize=4, resume=resume)
    assert caught.value.code == "RESUME_INPUTS_CHANGED"


def test_an_edited_fasta_refuses_the_continuation(project, tmp_path, stopped_at_max):
    file_id, set_id, resume = stopped_at_max
    edited = [(h, s) for h, s in CONTINUABLE] + [("C5_other", "AAGG")]
    write_fasta(tmp_path / "a.fasta", edited)
    with pytest.raises(ProjectError) as caught:
        run(project, set_id, [file_id], maxCandidateSize=4, resume=resume)
    assert caught.value.code == "RESUME_INPUTS_CHANGED"


def test_a_resume_without_a_fingerprint_is_refused(project, stopped_at_max):
    file_id, set_id, resume = stopped_at_max
    stripped = {key: value for key, value in resume.items() if key != "inputsFingerprint"}
    with pytest.raises(ProjectError) as caught:
        run(project, set_id, [file_id], maxCandidateSize=4, resume=stripped)
    assert caught.value.code == "RESUME_INPUTS_CHANGED"
