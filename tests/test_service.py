"""
Tests for the Electron <-> Python service boundary.

These exercise the boundary as data in / data out: every input is a plain dict
and every output must be JSON-serializable, because that is the contract the
desktop frontend depends on.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from molecular_diagnosis.service.errors import ErrorCode, ServiceError, to_service_error
from molecular_diagnosis.service.handlers import dispatch
from molecular_diagnosis.service.protocol import (
    decode_request,
    encode_failure,
    encode_success,
    resume_state_from_payload,
    resume_state_to_payload,
)

ALIGNED = (
    ">alpha_1|AU|Leptacis\nAACC\n"
    ">alpha_2|GB|Leptacis\nAACC\n"
    ">beta_1|FR|Synopeas\nAGCT\n"
    ">beta_2|DE|Synopeas\nCGTT\n"
)


@pytest.fixture()
def fasta(tmp_path: Path) -> Path:
    path = tmp_path / "aligned.fasta"
    path.write_text(ALIGNED, encoding="utf-8")
    return path


def json_roundtrip(value):
    """Everything crossing the boundary must survive JSON."""
    return json.loads(json.dumps(value))


# ---------------------------------------------------------------------------
# loadFasta
# ---------------------------------------------------------------------------


def test_load_fasta_returns_headers_and_metadata(fasta: Path) -> None:
    result = dispatch("loadFasta", {"path": str(fasta)})

    assert result["sequenceCount"] == 4
    assert result["alignmentLength"] == 4
    assert result["headers"] == [
        "alpha_1|AU|Leptacis",
        "alpha_2|GB|Leptacis",
        "beta_1|FR|Synopeas",
        "beta_2|DE|Synopeas",
    ]
    assert result["duplicateHeaderCount"] == 0
    assert result["alphabet"] == ["A", "C", "G", "T"]
    assert json_roundtrip(result) == result


def test_load_fasta_rejects_an_unaligned_file(tmp_path: Path) -> None:
    path = tmp_path / "ragged.fasta"
    path.write_text(">a\nACGT\n>b\nACG\n", encoding="utf-8")

    with pytest.raises(ServiceError) as caught:
        dispatch("loadFasta", {"path": str(path)})

    assert caught.value.code == ErrorCode.FASTA_NOT_ALIGNED
    # The original wording is preserved for logs rather than discarded.
    assert "same length" in (caught.value.detail or "")


def test_load_fasta_rejects_an_empty_file(tmp_path: Path) -> None:
    path = tmp_path / "empty.fasta"
    path.write_text("", encoding="utf-8")

    with pytest.raises(ServiceError) as caught:
        dispatch("loadFasta", {"path": str(path)})

    assert caught.value.code == ErrorCode.FASTA_EMPTY


def test_load_fasta_rejects_a_missing_file(tmp_path: Path) -> None:
    with pytest.raises(ServiceError) as caught:
        dispatch("loadFasta", {"path": str(tmp_path / "nope.fasta")})

    assert caught.value.code == ErrorCode.FASTA_NOT_FOUND


def test_load_fasta_reports_duplicate_headers(tmp_path: Path) -> None:
    # parse_fasta builds a dict, so a repeated header collapses; say so.
    path = tmp_path / "dupes.fasta"
    path.write_text(">a\nACGT\n>a\nTGCA\n>b\nAAAA\n", encoding="utf-8")

    result = dispatch("loadFasta", {"path": str(path)})

    assert result["headerLineCount"] == 3
    assert result["sequenceCount"] == 2
    assert result["duplicateHeaderCount"] == 1


# ---------------------------------------------------------------------------
# validateFocalStrings
# ---------------------------------------------------------------------------


def test_validate_one_focal_selector(fasta: Path) -> None:
    result = dispatch(
        "validateFocalStrings", {"path": str(fasta), "focalStrings": ["Leptacis"]}
    )

    assert result["results"][0]["matched"] is True
    assert result["results"][0]["matchCount"] == 2
    assert result["unionMatchCount"] == 2
    assert json_roundtrip(result) == result


def test_validate_multiple_selectors_act_as_or(fasta: Path) -> None:
    result = dispatch(
        "validateFocalStrings",
        {"path": str(fasta), "focalStrings": ["Leptacis", "Synopeas"]},
    )

    assert [entry["matchCount"] for entry in result["results"]] == [2, 2]
    assert result["unionMatchCount"] == 4


def test_validate_overlapping_selectors_do_not_double_count(fasta: Path) -> None:
    result = dispatch(
        "validateFocalStrings",
        {"path": str(fasta), "focalStrings": ["Leptacis", "alpha_1"]},
    )

    assert [entry["matchCount"] for entry in result["results"]] == [2, 1]
    # alpha_1 is already inside the Leptacis set, so the union stays at 2.
    assert result["unionMatchCount"] == 2


def test_validate_selector_matching_nothing_is_red(fasta: Path) -> None:
    result = dispatch(
        "validateFocalStrings", {"path": str(fasta), "focalStrings": ["Nothosaurus"]}
    )

    assert result["results"][0]["matched"] is False
    assert result["results"][0]["matchCount"] == 0
    assert result["unionMatchCount"] == 0


def test_validate_is_case_sensitive(fasta: Path) -> None:
    result = dispatch(
        "validateFocalStrings", {"path": str(fasta), "focalStrings": ["leptacis"]}
    )

    assert result["results"][0]["matched"] is False


def test_validate_reports_duplicates_individually(fasta: Path) -> None:
    # The UI asks about the entries it holds, in order, including repeats.
    result = dispatch(
        "validateFocalStrings",
        {"path": str(fasta), "focalStrings": ["Leptacis", "Leptacis"]},
    )

    assert len(result["results"]) == 2
    assert all(entry["matched"] for entry in result["results"])
    assert result["unionMatchCount"] == 2


def test_validate_flags_an_empty_selector_without_matching_everything(fasta: Path) -> None:
    result = dispatch(
        "validateFocalStrings", {"path": str(fasta), "focalStrings": ["", "Leptacis"]}
    )

    assert result["results"][0]["invalid"] is True
    assert result["results"][0]["matchCount"] == 0
    # The empty entry must not drag every header into the focal set.
    assert result["unionMatchCount"] == 2


def test_validate_selector_containing_a_semicolon(tmp_path: Path) -> None:
    path = tmp_path / "semi.fasta"
    path.write_text(">wei;rd|AU|Leptacis\nAACC\n>plain|FR|Synopeas\nAGCT\n", encoding="utf-8")

    result = dispatch("validateFocalStrings", {"path": str(path), "focalStrings": ["wei;rd"]})

    assert result["results"][0]["matched"] is True
    assert result["results"][0]["matchCount"] == 1


# ---------------------------------------------------------------------------
# runMolecularDiagnosis
# ---------------------------------------------------------------------------


def test_run_molecular_diagnosis_end_to_end(fasta: Path, tmp_path: Path) -> None:
    out = tmp_path / "out"
    out.mkdir()

    result = dispatch(
        "runMolecularDiagnosis",
        {
            "fastaPath": str(fasta),
            "focalStrings": ["Leptacis"],
            "outputDirectory": str(out),
            "ignoreGaps": False,
            "giveBenefitOfDoubtToAmbiguousBases": False,
            "minCandidateSize": 1,
            "maxCandidateSize": 2,
        },
    )

    assert result["focalStrings"] == ["Leptacis"]
    assert Path(result["outputs"]["reportTxt"]).is_file()
    assert Path(result["outputs"]["workbookXlsx"]).is_file()
    assert Path(result["outputs"]["consensusTxt"]).is_file()
    assert result["dmc"]["singleSites"] == [1, 3]
    assert result["dmc"]["stopReason"] == "found_at_or_above_minimum_length"
    assert result["canContinue"] is False
    assert result["resume"] is None
    assert json_roundtrip(result) == result


def test_run_molecular_diagnosis_accepts_multiple_selectors(
    fasta: Path, tmp_path: Path
) -> None:
    out = tmp_path / "out"
    out.mkdir()

    result = dispatch(
        "runMolecularDiagnosis",
        {
            "fastaPath": str(fasta),
            "focalStrings": ["alpha_1", "alpha_2"],
            "outputDirectory": str(out),
        },
    )

    # Same focal group as the single "Leptacis" selector, reached via OR.
    assert result["dmc"]["singleSites"] == [1, 3]


def test_run_molecular_diagnosis_deduplicates_selectors(fasta: Path, tmp_path: Path) -> None:
    out = tmp_path / "out"
    out.mkdir()

    result = dispatch(
        "runMolecularDiagnosis",
        {
            "fastaPath": str(fasta),
            "focalStrings": ["Leptacis", "Leptacis"],
            "outputDirectory": str(out),
        },
    )

    assert result["focalStrings"] == ["Leptacis"]


def test_run_molecular_diagnosis_rejects_an_empty_selector(
    fasta: Path, tmp_path: Path
) -> None:
    out = tmp_path / "out"
    out.mkdir()

    with pytest.raises(ServiceError) as caught:
        dispatch(
            "runMolecularDiagnosis",
            {
                "fastaPath": str(fasta),
                "focalStrings": ["Leptacis", "   "],
                "outputDirectory": str(out),
            },
        )

    assert caught.value.code == ErrorCode.FOCAL_EMPTY


def test_run_molecular_diagnosis_surfaces_a_focal_miss_as_a_code(
    fasta: Path, tmp_path: Path
) -> None:
    out = tmp_path / "out"
    out.mkdir()

    with pytest.raises(ServiceError) as caught:
        dispatch(
            "runMolecularDiagnosis",
            {
                "fastaPath": str(fasta),
                "focalStrings": ["Nothosaurus"],
                "outputDirectory": str(out),
            },
        )

    assert caught.value.code == ErrorCode.FOCAL_NO_MATCH


def test_run_molecular_diagnosis_rejects_a_missing_output_directory(
    fasta: Path, tmp_path: Path
) -> None:
    with pytest.raises(ServiceError) as caught:
        dispatch(
            "runMolecularDiagnosis",
            {
                "fastaPath": str(fasta),
                "focalStrings": ["Leptacis"],
                "outputDirectory": str(tmp_path / "absent"),
            },
        )

    assert caught.value.code == ErrorCode.OUTPUT_DIR_MISSING


# ---------------------------------------------------------------------------
# Continuation / resume
# ---------------------------------------------------------------------------


def _continuable(tmp_path: Path) -> tuple[Path, Path]:
    """
    A dataset with no 1- or 2-site diagnosis, so the search exhausts the
    configured maximum and offers to continue.
    """
    path = tmp_path / "hard.fasta"
    path.write_text(
        ">focal_1\nAAAC\n"
        ">focal_2\nAAAC\n"
        ">other_1\nAAAT\n"
        ">other_2\nAATC\n"
        ">other_3\nATAC\n"
        ">other_4\nTAAC\n",
        encoding="utf-8",
    )
    out = tmp_path / "out"
    out.mkdir()
    return path, out


def test_search_offers_continuation_when_it_hits_the_maximum(tmp_path: Path) -> None:
    fasta_path, out = _continuable(tmp_path)

    result = dispatch(
        "runMolecularDiagnosis",
        {
            "fastaPath": str(fasta_path),
            "focalStrings": ["focal"],
            "outputDirectory": str(out),
            "minCandidateSize": 1,
            "maxCandidateSize": 1,
        },
    )

    assert result["dmc"]["stopReason"] == "reached_maximum_length"
    assert result["canContinue"] is True
    assert result["resume"]["startCombinationLength"] == 2
    assert json_roundtrip(result) == result


def test_continuation_resumes_instead_of_restarting(tmp_path: Path) -> None:
    fasta_path, out = _continuable(tmp_path)

    request = {
        "fastaPath": str(fasta_path),
        "focalStrings": ["focal"],
        "outputDirectory": str(out),
        "minCandidateSize": 1,
        "maxCandidateSize": 1,
    }
    first = dispatch("runMolecularDiagnosis", request)
    assert first["canContinue"] is True

    resumed = dispatch(
        "runMolecularDiagnosis",
        {**request, "maxCandidateSize": 3, "resume": first["resume"]},
    )

    # Picks up above the exhausted length rather than starting again at 1.
    assert resumed["dmc"]["startCombinationLength"] == 2
    tested = resumed["dmc"]["diagnostics"]["combinationsTestedByLength"]
    # The length-1 count is carried forward, not recomputed from zero.
    assert tested["1"] == first["dmc"]["diagnostics"]["combinationsTestedByLength"]["1"]
    assert "2" in tested


def test_continuation_rejects_a_maximum_below_the_resume_point(tmp_path: Path) -> None:
    fasta_path, out = _continuable(tmp_path)

    request = {
        "fastaPath": str(fasta_path),
        "focalStrings": ["focal"],
        "outputDirectory": str(out),
        "minCandidateSize": 1,
        "maxCandidateSize": 1,
    }
    first = dispatch("runMolecularDiagnosis", request)

    with pytest.raises(ServiceError) as caught:
        dispatch("runMolecularDiagnosis", {**request, "resume": first["resume"]})

    assert caught.value.code == ErrorCode.INVALID_PARAMETER


def test_resume_state_survives_a_json_round_trip() -> None:
    payload = resume_state_to_payload(
        stopped_at_length=2,
        diagnostic_combinations=[(1, 3), (0, 2, 5)],
        combinations_tested_by_length={1: 4, 2: 6},
    )
    revived = resume_state_from_payload(json_roundtrip(payload))

    assert revived == (3, [(1, 3), (0, 2, 5)], {1: 4, 2: 6})


def test_resume_state_of_none_is_no_continuation() -> None:
    assert resume_state_from_payload(None) is None


def test_malformed_resume_state_is_rejected_not_ignored() -> None:
    # Silently restarting from scratch would waste a long search.
    with pytest.raises(ServiceError):
        resume_state_from_payload({"startCombinationLength": 0})

    with pytest.raises(ServiceError):
        resume_state_from_payload({"startCombinationLength": 2, "diagnosticCombinations": [["x"]]})


# ---------------------------------------------------------------------------
# Protocol
# ---------------------------------------------------------------------------


def test_decode_request_reads_a_well_formed_line() -> None:
    request = decode_request('{"id":"7","method":"ping","params":{"a":1}}')

    assert request.id == "7"
    assert request.method == "ping"
    assert request.params == {"a": 1}


def test_decode_request_defaults_missing_params() -> None:
    assert decode_request('{"id":"7","method":"ping"}').params == {}


@pytest.mark.parametrize(
    "line",
    ['not json', '[]', '{"method":"ping"}', '{"id":"7"}', '{"id":"7","method":"ping","params":3}'],
)
def test_decode_request_rejects_malformed_lines(line: str) -> None:
    with pytest.raises(ServiceError) as caught:
        decode_request(line)

    assert caught.value.code == ErrorCode.INVALID_REQUEST


def test_encoded_success_and_failure_are_single_line_json() -> None:
    success = encode_success("7", {"value": 1})
    assert "\n" not in success
    assert json.loads(success) == {"id": "7", "ok": True, "result": {"value": 1}}

    failure = encode_failure("7", ServiceError(ErrorCode.UNKNOWN, "boom", detail="why"))
    assert "\n" not in failure
    decoded = json.loads(failure)
    assert decoded["ok"] is False
    assert decoded["error"]["code"] == ErrorCode.UNKNOWN
    assert decoded["error"]["detail"] == "why"


def test_unknown_method_is_a_coded_error() -> None:
    with pytest.raises(ServiceError) as caught:
        dispatch("noSuchMethod", {})

    assert caught.value.code == ErrorCode.UNKNOWN_METHOD


def test_to_service_error_preserves_the_original_exception() -> None:
    original = RuntimeError("something specific went wrong")
    wrapped = to_service_error(original)

    assert wrapped.code == ErrorCode.UNKNOWN
    assert wrapped.cause is original
    assert "something specific went wrong" in (wrapped.detail or "")
    assert "RuntimeError" in (wrapped.traceback_text or "")


# ---------------------------------------------------------------------------
# The process boundary itself
# ---------------------------------------------------------------------------


def test_stdio_process_serves_requests_and_keeps_stdout_clean(fasta: Path) -> None:
    """
    Drive the real child process the way Electron does, and confirm every
    stdout line is protocol JSON.
    """
    requests = (
        json.dumps({"id": "1", "method": "ping", "params": {}})
        + "\n"
        + json.dumps({"id": "2", "method": "loadFasta", "params": {"path": str(fasta)}})
        + "\n"
        + json.dumps({"id": "3", "method": "noSuchMethod", "params": {}})
        + "\n"
    )

    completed = subprocess.run(
        [sys.executable, "-m", "molecular_diagnosis.service"],
        input=requests,
        capture_output=True,
        text=True,
        timeout=120,
        cwd=str(Path(__file__).resolve().parent.parent),
    )

    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    messages = [json.loads(line) for line in lines]  # every line must be JSON

    assert messages[0]["id"] == "ready"
    by_id = {message["id"]: message for message in messages}

    assert by_id["1"]["ok"] is True
    assert by_id["2"]["ok"] is True
    assert by_id["2"]["result"]["sequenceCount"] == 4
    assert by_id["3"]["ok"] is False
    assert by_id["3"]["error"]["code"] == ErrorCode.UNKNOWN_METHOD
