from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DMCResult:
    single: list[int]
    pairs: list[tuple[int, int]]
    unique: list[int]
    fixed_count: int
    skipped_non_acgt: int
    globally_conserved_removed: int
    candidate_count: int
    pairs_tested: int


@dataclass(frozen=True)
class FiveSiteResult:
    total_combinations_tested: int
    best_gap_score: float | None
    best_gap_sites: tuple[int, ...] | None
    best_avg_score: float | None
    best_avg_sites: tuple[int, ...] | None


@dataclass(frozen=True)
class PipelineResult:
    txt_output_path: Path
    xlsx_output_path: Path