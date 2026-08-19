"""
Old-vs-new Molecular Diagnosis parity.

Runs `parity_driver.py` twice — once against the current working tree and once
against a read-only git worktree of the pre-integration commit — and asserts the
two produce identical scientific results for the same single-selector inputs.

The claim under test: with exactly one focal selector, the refactor that
introduced `molecular_diagnosis.focal` and the service boundary changed nothing
about what the analysis computes or writes.

The baseline is created with `git worktree add --detach`, which never touches
the current working tree. It is removed again in the fixture teardown. If the
repository or the baseline commit is unavailable the whole module skips rather
than failing, so the suite still runs from a source export.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
DRIVER = REPO_ROOT / "tests" / "parity_driver.py"

#: Commit holding the pre-integration Python. Everything the backend-integration
#: pass changed in `molecular_diagnosis/` was uncommitted at the time this test
#: was written, so the tip of the branch is the correct baseline.
BASELINE_REF = "1b01284"


def _git(*args: str, cwd: Path = REPO_ROOT) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", *args], cwd=str(cwd), capture_output=True, text=True, timeout=180
    )


@pytest.fixture(scope="module")
def baseline_tree(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A read-only checkout of the pre-integration commit."""
    if shutil.which("git") is None:
        pytest.skip("git is not available")

    probe = _git("rev-parse", "--verify", f"{BASELINE_REF}^{{commit}}")
    if probe.returncode != 0:
        pytest.skip(f"baseline commit {BASELINE_REF} is not in this repository")

    target = tmp_path_factory.mktemp("baseline") / "tree"
    created = _git("worktree", "add", "--detach", str(target), BASELINE_REF)
    if created.returncode != 0:
        pytest.skip(f"could not create baseline worktree: {created.stderr.strip()}")

    # Sanity: the baseline must genuinely predate the integration work.
    assert not (target / "molecular_diagnosis" / "focal.py").exists()
    assert not (target / "molecular_diagnosis" / "service").exists()

    try:
        yield target
    finally:
        _git("worktree", "remove", "--force", str(target))


def _run_driver(tree: Path, spec: dict, work: Path, label: str) -> dict:
    """Execute the driver with `tree` on the import path, and read its dump."""
    spec_path = work / f"{label}-spec.json"
    out_path = work / f"{label}-out.json"
    output_dir = work / f"{label}-outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    spec = {**spec, "output_dir": str(output_dir)}
    spec_path.write_text(json.dumps(spec), encoding="utf-8")

    env = {**os.environ, "PYTHONPATH": str(tree), "PYTHONDONTWRITEBYTECODE": "1"}
    completed = subprocess.run(
        [sys.executable, str(DRIVER), str(spec_path), str(out_path)],
        capture_output=True,
        text=True,
        timeout=600,
        env=env,
        cwd=str(tree),
    )

    if completed.returncode != 0:
        raise AssertionError(
            f"parity driver failed for {label}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )

    return json.loads(out_path.read_text(encoding="utf-8"))


def run_both(baseline_tree: Path, spec: dict, work: Path) -> tuple[dict, dict]:
    old = _run_driver(baseline_tree, spec, work, "old")
    new = _run_driver(REPO_ROOT, spec, work, "new")
    return old, new


def write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    path.write_text(
        "".join(f">{header}\n{sequence}\n" for header, sequence in records),
        encoding="utf-8",
    )
    return path


def base_spec(fasta: Path, focal: str, **overrides) -> dict:
    spec = {
        "fasta_path": str(fasta),
        "focal": focal,
        "include_ambiguous_dmc_bd": False,
        "include_gappy_consensus_dmc_sites": False,
        "min_combination_length": 1,
        "max_combination_length": 2,
    }
    spec.update(overrides)
    return spec


# ---------------------------------------------------------------------------
# Datasets
# ---------------------------------------------------------------------------

#: A clean single-site diagnosis: focal fixed at sites 1 and 3, contrast varies.
SIMPLE = [
    ("focal_1|AU|Target", "AACC"),
    ("focal_2|GB|Target", "AACC"),
    ("other_1|FR|Contrast", "AGCT"),
    ("other_2|DE|Contrast", "CGTT"),
]

#: No single-site diagnosis, so a max of 1 exhausts and offers to continue.
HARD = [
    ("focal_1|AU|Target", "AAAC"),
    ("focal_2|GB|Target", "AAAC"),
    ("other_1|FR|Contrast", "AAAT"),
    ("other_2|DE|Contrast", "AATC"),
    ("other_3|ES|Contrast", "ATAC"),
    ("other_4|IT|Contrast", "TAAC"),
]

#: Ambiguity codes in focal and contrast, to exercise the BD option.
AMBIGUOUS = [
    ("focal_1|AU|Target", "AACCGT"),
    ("focal_2|GB|Target", "AACCRT"),
    ("focal_3|NL|Target", "AACCGT"),
    ("other_1|FR|Contrast", "AGCTGN"),
    ("other_2|DE|Contrast", "CGTTAT"),
    ("other_3|ES|Contrast", "AYCTGT"),
]

#: Gaps in focal and contrast columns, to exercise the gap option.
GAPPY = [
    ("focal_1|AU|Target", "AACC-GTA"),
    ("focal_2|GB|Target", "AA-CGGTA"),
    ("focal_3|NL|Target", "AACCGGTA"),
    ("other_1|FR|Contrast", "AGCTG-TA"),
    ("other_2|DE|Contrast", "CGTTA-TA"),
    ("other_3|ES|Contrast", "A-CTGGTA"),
]


# ---------------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------------

#: Fields whose difference would be a scientific change, checked one at a time
#: so a failure names the thing that drifted.
DMC_FIELDS = [
    "single",
    "pairs",
    "unique",
    "states",
    "diagnostic_combinations",
    "diagnostic_combinations_by_length",
    "min_combination_length",
    "max_combination_length",
    "start_combination_length",
    "stopped_at_length",
    "stop_reason",
    "combinations_tested_by_length",
    "total_combinations_tested",
    "fixed_count",
    "skipped_non_acgt",
    "globally_conserved_removed",
    "candidate_count",
    "pairs_tested",
    "ambiguous_bd_sites_included",
    "gappy_consensus_sites_included",
    "include_ambiguous_dmc_bd",
    "include_gappy_consensus_dmc_sites",
]


def assert_parity(old: dict, new: dict, *, case: str) -> None:
    """Every comparable aspect of the two runs must agree."""
    assert old["focal_headers"] == new["focal_headers"], f"{case}: focal selection"
    assert old["non_focal_headers"] == new["non_focal_headers"], f"{case}: non-focal selection"
    assert old["ref_id"] == new["ref_id"], f"{case}: reference sequence"
    assert old["alignment_length"] == new["alignment_length"], f"{case}: alignment length"

    for field in DMC_FIELDS:
        assert old["dmc"][field] == new["dmc"][field], f"{case}: dmc.{field}"

    assert old["pipeline_dmc"] == new["pipeline_dmc"], f"{case}: pipeline DMC result"
    assert old["five_site"] == new["five_site"], f"{case}: five-site optimisation"
    assert old["consensus"] == new["consensus"], f"{case}: consensus"

    assert old["outputs"] == new["outputs"], f"{case}: output file names"
    assert old["report_text"] == new["report_text"], f"{case}: DMC text report"
    assert old["consensus_text"] == new["consensus_text"], f"{case}: consensus text report"

    assert (
        old["workbook"]["sheet_names"] == new["workbook"]["sheet_names"]
    ), f"{case}: workbook sheets"
    for sheet in old["workbook"]["sheets"]:
        assert (
            old["workbook"]["sheets"][sheet] == new["workbook"]["sheets"][sheet]
        ), f"{case}: workbook sheet {sheet}"


# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------


def test_parity_default_settings(baseline_tree: Path, tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path / "simple.fasta", SIMPLE)
    old, new = run_both(baseline_tree, base_spec(fasta, "Target"), tmp_path)

    assert_parity(old, new, case="default settings")
    # Guard the fixture itself: this case must really exercise a 1-site find.
    assert new["dmc"]["single"], "fixture no longer produces a single-site diagnosis"
    assert new["dmc"]["stop_reason"] == "found_at_or_above_minimum_length"


def test_parity_single_site_diagnosis_details(baseline_tree: Path, tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path / "simple.fasta", SIMPLE)
    old, new = run_both(baseline_tree, base_spec(fasta, "Target"), tmp_path)

    assert old["dmc"]["single"] == new["dmc"]["single"]
    assert old["dmc"]["unique"] == new["dmc"]["unique"]
    # Early stop means pairs were never enumerated; that must still hold.
    assert new["dmc"]["pairs"] == []
    assert new["dmc"]["pairs_tested"] == 0
    assert old["dmc"]["pairs_tested"] == new["dmc"]["pairs_tested"]


def test_parity_reaches_configured_maximum(baseline_tree: Path, tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path / "hard.fasta", HARD)
    spec = base_spec(fasta, "Target", min_combination_length=1, max_combination_length=1)
    old, new = run_both(baseline_tree, spec, tmp_path)

    assert_parity(old, new, case="reaches maximum")
    assert new["dmc"]["stop_reason"] == "reached_maximum_length"
    assert new["dmc"]["stopped_at_length"] == 1


def test_parity_continuation_to_a_larger_maximum(baseline_tree: Path, tmp_path: Path) -> None:
    """
    The continuation the Tkinter loop performed: resume from the next length,
    carrying the previous combinations and tested-counts forward.
    """
    fasta = write_fasta(tmp_path / "hard.fasta", HARD)

    first_spec = base_spec(fasta, "Target", min_combination_length=1, max_combination_length=1)
    first_old, first_new = run_both(baseline_tree, first_spec, tmp_path / "pass1")
    assert_parity(first_old, first_new, case="continuation pass 1")
    assert first_new["dmc"]["stop_reason"] == "reached_maximum_length"

    resume = {
        "start_combination_length": first_new["dmc"]["stopped_at_length"] + 1,
        "initial_diagnostic_combinations": first_new["dmc"]["diagnostic_combinations"],
        "initial_combinations_tested_by_length": first_new["dmc"][
            "combinations_tested_by_length"
        ],
    }
    second_spec = base_spec(
        fasta,
        "Target",
        min_combination_length=1,
        max_combination_length=3,
        resume=resume,
    )
    second_old, second_new = run_both(baseline_tree, second_spec, tmp_path / "pass2")

    assert_parity(second_old, second_new, case="continuation pass 2")
    assert second_new["dmc"]["start_combination_length"] == 2
    # The carried-forward count survives rather than being recomputed.
    assert (
        second_new["dmc"]["combinations_tested_by_length"]["1"]
        == first_new["dmc"]["combinations_tested_by_length"]["1"]
    )


def test_parity_non_default_min_and_max(baseline_tree: Path, tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path / "hard.fasta", HARD)
    spec = base_spec(fasta, "Target", min_combination_length=2, max_combination_length=3)
    old, new = run_both(baseline_tree, spec, tmp_path)

    assert_parity(old, new, case="min=2 max=3")
    # The documented quirk: `min` does not raise the starting length.
    assert new["dmc"]["start_combination_length"] == 1
    assert "1" in new["dmc"]["combinations_tested_by_length"]


def test_parity_min_above_max_is_still_rejected(baseline_tree: Path, tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path / "simple.fasta", SIMPLE)
    spec = base_spec(fasta, "Target", min_combination_length=3, max_combination_length=2)

    for tree, label in ((baseline_tree, "old"), (REPO_ROOT, "new")):
        with pytest.raises(AssertionError) as caught:
            _run_driver(tree, spec, tmp_path, label)
        assert "Minimum combination length cannot exceed" in str(caught.value)


@pytest.mark.parametrize("ambiguous", [False, True])
def test_parity_ambiguous_base_option(
    baseline_tree: Path, tmp_path: Path, ambiguous: bool
) -> None:
    fasta = write_fasta(tmp_path / "ambiguous.fasta", AMBIGUOUS)
    spec = base_spec(fasta, "Target", include_ambiguous_dmc_bd=ambiguous)
    old, new = run_both(baseline_tree, spec, tmp_path)

    assert_parity(old, new, case=f"ambiguous BD={ambiguous}")


@pytest.mark.parametrize("gaps", [False, True])
def test_parity_gap_option(baseline_tree: Path, tmp_path: Path, gaps: bool) -> None:
    fasta = write_fasta(tmp_path / "gappy.fasta", GAPPY)
    spec = base_spec(fasta, "Target", include_gappy_consensus_dmc_sites=gaps)
    old, new = run_both(baseline_tree, spec, tmp_path)

    assert_parity(old, new, case=f"ignore gaps={gaps}")


def test_parity_both_options_together(baseline_tree: Path, tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path / "gappy.fasta", GAPPY)
    spec = base_spec(
        fasta,
        "Target",
        include_ambiguous_dmc_bd=True,
        include_gappy_consensus_dmc_sites=True,
        max_combination_length=3,
    )
    old, new = run_both(baseline_tree, spec, tmp_path)

    assert_parity(old, new, case="both options on")


def test_parity_five_site_optimisation(baseline_tree: Path, tmp_path: Path) -> None:
    """A dataset with enough unique sites that the 5-site search actually runs."""
    records = [
        ("focal_1|AU|Target", "AACCGTAC"),
        ("focal_2|GB|Target", "AACCGTAC"),
        ("other_1|FR|Contrast", "CGTTACGT"),
        ("other_2|DE|Contrast", "GTAGCATG"),
        ("other_3|ES|Contrast", "TCGATGCA"),
    ]
    fasta = write_fasta(tmp_path / "five.fasta", records)
    old, new = run_both(baseline_tree, base_spec(fasta, "Target"), tmp_path)

    assert_parity(old, new, case="five-site optimisation")
    assert new["five_site"]["total_combinations_tested"] > 0, "fixture did not reach the 5-site search"
    assert new["five_site"]["best_gap_sites"] is not None


def test_parity_on_the_real_trial_alignment(baseline_tree: Path, tmp_path: Path) -> None:
    """The repository's own aligned sample, if it is present on this machine."""
    sample = REPO_ROOT / "input" / "Leptacis_trial_alignment_output.fasta"
    if not sample.is_file():
        pytest.skip("input/Leptacis_trial_alignment_output.fasta is not present")

    spec = base_spec(sample, "phantasmatica", max_combination_length=3)
    old, new = run_both(baseline_tree, spec, tmp_path)

    assert_parity(old, new, case="real trial alignment")


# ---------------------------------------------------------------------------
# The integration path itself
# ---------------------------------------------------------------------------


def test_service_path_matches_the_baseline(baseline_tree: Path, tmp_path: Path) -> None:
    """
    Parity for the route Electron actually takes.

    The tests above compare the refactored library against the baseline. This
    one drives `service.dispatch('runMolecularDiagnosis', ...)` — the same entry
    point the desktop app calls — and checks its payload against the baseline
    run, so a regression introduced by the service layer itself would be caught
    rather than hidden behind the library comparison.
    """
    from molecular_diagnosis.service.handlers import dispatch

    fasta = write_fasta(tmp_path / "simple.fasta", SIMPLE)
    old = _run_driver(baseline_tree, base_spec(fasta, "Target"), tmp_path, "old")

    service_out = tmp_path / "service-outputs"
    service_out.mkdir()
    payload = dispatch(
        "runMolecularDiagnosis",
        {
            "fastaPath": str(fasta),
            "focalStrings": ["Target"],
            "outputDirectory": str(service_out),
            "ignoreGaps": False,
            "giveBenefitOfDoubtToAmbiguousBases": False,
            "minCandidateSize": 1,
            "maxCandidateSize": 2,
        },
    )

    dmc = payload["dmc"]
    diagnostics = dmc["diagnostics"]

    assert dmc["singleSites"] == old["dmc"]["single"]
    assert dmc["pairs"] == old["dmc"]["pairs"]
    assert dmc["uniqueSites"] == old["dmc"]["unique"]
    assert dmc["states"] == old["dmc"]["states"]
    assert dmc["combinationsByLength"] == old["dmc"]["diagnostic_combinations_by_length"]
    assert dmc["stopReason"] == old["dmc"]["stop_reason"]
    assert dmc["stoppedAtLength"] == old["dmc"]["stopped_at_length"]
    assert dmc["minCombinationLength"] == old["dmc"]["min_combination_length"]
    assert dmc["maxCombinationLength"] == old["dmc"]["max_combination_length"]
    assert dmc["startCombinationLength"] == old["dmc"]["start_combination_length"]

    assert diagnostics["fixedCount"] == old["dmc"]["fixed_count"]
    assert diagnostics["skippedSites"] == old["dmc"]["skipped_non_acgt"]
    assert diagnostics["globallyConservedRemoved"] == old["dmc"]["globally_conserved_removed"]
    assert diagnostics["candidateCount"] == old["dmc"]["candidate_count"]
    assert diagnostics["pairsTested"] == old["dmc"]["pairs_tested"]
    assert diagnostics["totalCombinationsTested"] == old["dmc"]["total_combinations_tested"]
    assert (
        diagnostics["combinationsTestedByLength"]
        == old["dmc"]["combinations_tested_by_length"]
    )
    assert diagnostics["ambiguousBdSitesIncluded"] == old["dmc"]["ambiguous_bd_sites_included"]
    assert (
        diagnostics["gappyConsensusSitesIncluded"]
        == old["dmc"]["gappy_consensus_sites_included"]
    )

    # The files the service produced must read the same as the baseline's.
    report = Path(payload["outputs"]["reportTxt"]).read_text(encoding="utf-8")
    report = report.replace(str(service_out), "<OUTPUT_DIR>").replace(
        str(service_out).replace("\\", "/"), "<OUTPUT_DIR>"
    )
    assert report == old["report_text"]

    consensus = Path(payload["outputs"]["consensusTxt"]).read_text(encoding="utf-8")
    assert consensus == old["consensus_text"]

    assert dump_workbook_via_driver(Path(payload["outputs"]["workbookXlsx"])) == old["workbook"]


def dump_workbook_via_driver(path: Path) -> dict:
    """Reuse the driver's workbook dump so both sides are read identically."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("parity_driver", DRIVER)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.dump_workbook(path)


def test_service_continuation_matches_the_baseline(
    baseline_tree: Path, tmp_path: Path
) -> None:
    """The service's resume blob must reproduce the baseline continuation."""
    from molecular_diagnosis.service.handlers import dispatch

    fasta = write_fasta(tmp_path / "hard.fasta", HARD)

    first_spec = base_spec(fasta, "Target", min_combination_length=1, max_combination_length=1)
    old_first = _run_driver(baseline_tree, first_spec, tmp_path / "p1", "old")

    resume_spec = base_spec(
        fasta,
        "Target",
        min_combination_length=1,
        max_combination_length=3,
        resume={
            "start_combination_length": old_first["dmc"]["stopped_at_length"] + 1,
            "initial_diagnostic_combinations": old_first["dmc"]["diagnostic_combinations"],
            "initial_combinations_tested_by_length": old_first["dmc"][
                "combinations_tested_by_length"
            ],
        },
    )
    old_second = _run_driver(baseline_tree, resume_spec, tmp_path / "p2", "old")

    service_out = tmp_path / "service-outputs"
    service_out.mkdir()
    request = {
        "fastaPath": str(fasta),
        "focalStrings": ["Target"],
        "outputDirectory": str(service_out),
        "minCandidateSize": 1,
        "maxCandidateSize": 1,
    }
    first = dispatch("runMolecularDiagnosis", request)
    assert first["canContinue"] is True

    second = dispatch(
        "runMolecularDiagnosis",
        {**request, "maxCandidateSize": 3, "resume": first["resume"]},
    )

    assert second["dmc"]["singleSites"] == old_second["dmc"]["single"]
    assert second["dmc"]["uniqueSites"] == old_second["dmc"]["unique"]
    assert second["dmc"]["combinationsByLength"] == old_second["dmc"][
        "diagnostic_combinations_by_length"
    ]
    assert second["dmc"]["stopReason"] == old_second["dmc"]["stop_reason"]
    assert second["dmc"]["stoppedAtLength"] == old_second["dmc"]["stopped_at_length"]
    assert second["dmc"]["startCombinationLength"] == old_second["dmc"][
        "start_combination_length"
    ]
    assert (
        second["dmc"]["diagnostics"]["combinationsTestedByLength"]
        == old_second["dmc"]["combinations_tested_by_length"]
    )
