from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DMCResult:
    single: list[int]
    pairs: list[tuple[int, int]]
    unique: list[int]

    states: dict[int, str]

    diagnostic_combinations: list[tuple[int, ...]]
    diagnostic_combinations_by_length: dict[int, list[tuple[int, ...]]]

    min_combination_length: int
    max_combination_length: int
    start_combination_length: int
    stopped_at_length: int
    stop_reason: str

    combinations_tested_by_length: dict[int, int]
    total_combinations_tested: int

    fixed_count: int
    skipped_non_acgt: int
    globally_conserved_removed: int
    candidate_count: int
    pairs_tested: int

    ambiguous_bd_sites_included: list[int]
    gappy_consensus_sites_included: list[int]

    include_ambiguous_dmc_bd: bool
    include_gappy_consensus_dmc_sites: bool


@dataclass(frozen=True)
class FiveSiteResult:
    total_combinations_tested: int
    best_gap_score: float | None
    best_gap_sites: tuple[int, ...] | None
    best_avg_score: float | None
    best_avg_sites: tuple[int, ...] | None


@dataclass(frozen=True)
class PunishmentEvent:
    sequence_id: str
    site: int
    state: str
    category: str
    score: float
    empty_fraction: float
    direction: str | None = None
    note: str = ""

@dataclass(frozen=True)
class PunishmentResult:
    total_scores: dict[str, float]
    polymorphism_scores: dict[str, float]
    balancing_scores: dict[str, float]
    prolongation_scores: dict[str, float]
    insertion_scores: dict[str, float]
    empty_weight_scores: dict[str, float]
    bd_counts: dict[str, int]
    prl_counts: dict[str, int]
    ins_counts: dict[str, int]
    events: list[PunishmentEvent]

@dataclass(frozen=True)
class PunishmentPipelineResult:
    xlsx_output_path: Path
    

@dataclass(frozen=True)
class PipelineResult:
    txt_output_path: Path
    xlsx_output_path: Path
    dmc: DMCResult