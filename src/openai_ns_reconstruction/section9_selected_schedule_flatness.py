"""Fail-closed Section 9 admission for the closed actual residual-flatness schedule.

The pinned formal development contains
``LocalResidualFlatness.selected_schedule``.  That theorem instantiates the
repository's actual selected finite-stage estimate record, chooses one common
diagonal schedule, retains the three cut bounds, and proves all quantitative
residual jet rates for every derivative order and every nonnegative decay
power.

This module is only a formal-export handoff.  It does not choose the schedule
in Python, evaluate any stage field or residual, infer rates from samples, or
materialize the paper velocity.  Lean-only functions and existential witnesses
remain opaque stable identifiers and evidence other than an external
``lean-formal-export`` is rejected.
"""
from __future__ import annotations

from dataclasses import dataclass


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FORMAL_FILE = "NavierStokes/LocalResidualFlatness.lean"
PINNED_SELECTED_SCHEDULE_THEOREM = (
    "NavierStokes.LocalResidualFlatness.selected_schedule"
)
PINNED_EXISTS_SCHEDULE_THEOREM = (
    "NavierStokes.LocalResidualFlatness.exists_schedule_all_jetRates"
)
PINNED_ALL_RESIDUAL_JET_RATES = (
    "NavierStokes.LocalResidualFlatness.AllResidualJetRates"
)
PINNED_SELECTED_ESTIMATES = (
    "NavierStokes.LocalResidualFlatness.selectedEstimates"
)
PINNED_ACTUAL_ESTIMATES = "NavierStokes.ActualCandidateAssembly.estimates"
PINNED_SELECTED_POTENTIAL_STAGES = (
    "NavierStokes.ActualCandidateAssembly.selectedPotentialStages"
)
PINNED_SELECTED_DIRECT_STAGES = (
    "NavierStokes.ActualCandidateAssembly.selectedDirectStages"
)
PINNED_SELECTED_PRESSURE_STAGES = (
    "NavierStokes.ActualCandidateAssembly.selectedPressureStages"
)
PINNED_SELECTED_QBIG = "NavierStokes.ActualCandidateConstruction.selectedQbig"
PINNED_SELECTED_SCHEDULE_STRUCTURE = (
    "NavierStokes.MixedCandidateWitness.SelectedSchedule"
)
PINNED_THREE_CUT_BOUNDS = "NavierStokes.MixedDiagonalSchedule.ThreeCutBounds"
REQUIRED_FORMAL_SYMBOLS = (
    PINNED_SELECTED_SCHEDULE_THEOREM,
    PINNED_EXISTS_SCHEDULE_THEOREM,
    PINNED_ALL_RESIDUAL_JET_RATES,
    PINNED_SELECTED_ESTIMATES,
    PINNED_ACTUAL_ESTIMATES,
    PINNED_SELECTED_POTENTIAL_STAGES,
    PINNED_SELECTED_DIRECT_STAGES,
    PINNED_SELECTED_PRESSURE_STAGES,
    PINNED_SELECTED_QBIG,
    PINNED_SELECTED_SCHEDULE_STRUCTURE,
    PINNED_THREE_CUT_BOUNDS,
)
PRODUCER_KIND = "lean-formal-export"


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be formal-export certified true")
    return True


@dataclass(frozen=True)
class Section9SelectedScheduleFlatnessWitness:
    """Opaque metadata for one exported application of ``selected_schedule``."""

    selected_estimates_id: str
    selected_qbig_id: str
    potential_stages_id: str
    direct_stages_id: str
    pressure_stages_id: str
    selected_schedule_id: str
    three_cut_bounds_id: str
    all_residual_jet_rates_id: str
    application_id: str
    producer_kind: str
    provenance: str

    selected_estimates_identification_certified: bool
    selected_qbig_identification_certified: bool
    selected_stage_families_identification_certified: bool
    schedule_witness_export_certified: bool
    three_cut_bounds_certified: bool
    all_residual_jet_rates_certified: bool
    theorem_application_certified: bool

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_FORMAL_FILE
    theorem_symbol: str = PINNED_SELECTED_SCHEDULE_THEOREM
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    def __post_init__(self) -> None:
        for name in (
            "selected_estimates_id",
            "selected_qbig_id",
            "potential_stages_id",
            "direct_stages_id",
            "pressure_stages_id",
            "selected_schedule_id",
            "three_cut_bounds_id",
            "all_residual_jet_rates_id",
            "application_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        if self.producer_kind != PRODUCER_KIND:
            raise ValueError(
                "producer_kind must be lean-formal-export; "
                "sampled/fitted/numeric-scan evidence is rejected"
            )
        if self.formal_repository != PINNED_FORMAL_REPOSITORY:
            raise ValueError("formal_repository must match the pinned source exactly")
        if self.formal_commit != PINNED_FORMAL_COMMIT:
            raise ValueError("formal_commit must match the pinned source exactly")
        if self.formal_file != PINNED_FORMAL_FILE:
            raise ValueError("formal_file must match LocalResidualFlatness.lean exactly")
        if self.theorem_symbol != PINNED_SELECTED_SCHEDULE_THEOREM:
            raise ValueError("theorem_symbol must be LocalResidualFlatness.selected_schedule")
        if self.dependency_symbols != REQUIRED_FORMAL_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned selected-schedule chain exactly")

        for name in (
            "selected_estimates_identification_certified",
            "selected_qbig_identification_certified",
            "selected_stage_families_identification_certified",
            "schedule_witness_export_certified",
            "three_cut_bounds_certified",
            "all_residual_jet_rates_certified",
            "theorem_application_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section9SelectedScheduleFlatnessAdmission:
    """Validated formal handoff for the closed actual schedule-flatness result."""

    witness: Section9SelectedScheduleFlatnessWitness

    def __post_init__(self) -> None:
        if not isinstance(self.witness, Section9SelectedScheduleFlatnessWitness):
            raise TypeError("witness must be a Section9SelectedScheduleFlatnessWitness")

    @property
    def closed_selected_schedule_theorem_admitted(self) -> bool:
        return True

    @property
    def all_quantitative_residual_jet_rates_admitted(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def selected_schedule_theorem_machine_replayed(self) -> bool:
        return False

    @property
    def selected_schedule_values_materialized(self) -> bool:
        return False

    @property
    def finite_stage_fields_materialized(self) -> bool:
        return False

    @property
    def residual_jet_rate_values_materialized(self) -> bool:
        return False

    @property
    def section9_iteration_machine_materialized(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
