"""Fail-closed Section 9 admission for the literal finite-stage estimates.

The pinned formal development contains
``ActualStageEstimates.stageEstimates_of_representations``.  Starting from one
literal correction-cycle ``RunData`` together with its represented physical
stress/divergence/mean data and the checked analytic controls, that theorem
produces ``MixedCandidateAssembly.StageEstimates`` for the canonical finite
velocity, pressure, and forcing ledgers.  Its gain/loss schedules are fixed by
``q1 q = q^(1/2)``.

This module is deliberately only a formal-export handoff.  It does not build a
cycle, evaluate a stress, fit a rate, choose q numerically, or materialize any
paper field.  Lean-only objects stay opaque stable identifiers and evidence
other than an external ``lean-formal-export`` is rejected.
"""
from __future__ import annotations

from dataclasses import dataclass


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FORMAL_FILE = "NavierStokes/ActualStageEstimates.lean"
PINNED_STAGE_ESTIMATES_THEOREM = (
    "NavierStokes.ActualStageEstimates.stageEstimates_of_representations"
)
PINNED_LEDGER_COROLLARY = (
    "NavierStokes.ActualStageEstimates.stageEstimates_ledger"
)
PINNED_STAGE_ESTIMATES_STRUCTURE = (
    "NavierStokes.MixedCandidateAssembly.StageEstimates"
)
PINNED_FINITE_VELOCITY = "NavierStokes.ActualIterationLedger.finiteVelocity"
PINNED_FINITE_PRESSURE = "NavierStokes.ActualIterationLedger.finitePressure"
PINNED_FINITE_FORCING = "NavierStokes.ActualIterationLedger.finiteForcing"
REQUIRED_FORMAL_SYMBOLS = (
    PINNED_STAGE_ESTIMATES_THEOREM,
    PINNED_LEDGER_COROLLARY,
    PINNED_STAGE_ESTIMATES_STRUCTURE,
    PINNED_FINITE_VELOCITY,
    PINNED_FINITE_PRESSURE,
    PINNED_FINITE_FORCING,
)
PRODUCER_KIND = "lean-formal-export"
Q1_DEFINITION = "q ^ (1 / 2 : Real)"
GAIN_DEFINITION = "q1 q ^ (g : Real)"
LOSS_DEFINITION = "q1 q ^ (g : Real) * (1 + (g : Real))"


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be formal-export certified true")
    return True


@dataclass(frozen=True)
class Section9ActualStageEstimatesWitness:
    """Opaque metadata for one formal stage-estimates theorem application."""

    domain_id: str
    forcing_id: str
    q_id: str
    run_data_id: str
    physical_data_id: str
    positive_bump_id: str
    finite_velocity_family_id: str
    finite_pressure_family_id: str
    finite_forcing_family_id: str
    mixed_physical_data_id: str
    stage_estimates_id: str
    application_id: str
    producer_kind: str
    provenance: str

    q_open_unit_interval_certified: bool
    input_normalized_certified: bool
    output_strip_controlled_certified: bool
    output_mean_balanced_certified: bool
    stress_controlled_certified: bool
    div_osc_controlled_certified: bool
    velocity_increments_controlled_certified: bool
    pressure_increments_controlled_certified: bool
    forcing_increments_controlled_certified: bool
    forcing_increment_margin_certified: bool
    ledger_velocity_identification_certified: bool
    ledger_pressure_identification_certified: bool
    ledger_forcing_identification_certified: bool
    theorem_application_certified: bool

    q1_definition: str = Q1_DEFINITION
    gain_definition: str = GAIN_DEFINITION
    loss_definition: str = LOSS_DEFINITION
    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_FORMAL_FILE
    theorem_symbol: str = PINNED_STAGE_ESTIMATES_THEOREM
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    def __post_init__(self) -> None:
        for name in (
            "domain_id",
            "forcing_id",
            "q_id",
            "run_data_id",
            "physical_data_id",
            "positive_bump_id",
            "finite_velocity_family_id",
            "finite_pressure_family_id",
            "finite_forcing_family_id",
            "mixed_physical_data_id",
            "stage_estimates_id",
            "application_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        if self.producer_kind != PRODUCER_KIND:
            raise ValueError(
                "producer_kind must be lean-formal-export; "
                "sampled/fitted/numeric-scan evidence is rejected"
            )
        if self.q1_definition != Q1_DEFINITION:
            raise ValueError("q1_definition must match ActualStageEstimates.q1 exactly")
        if self.gain_definition != GAIN_DEFINITION:
            raise ValueError("gain_definition must match the pinned theorem exactly")
        if self.loss_definition != LOSS_DEFINITION:
            raise ValueError("loss_definition must match the pinned theorem exactly")
        if self.formal_repository != PINNED_FORMAL_REPOSITORY:
            raise ValueError("formal_repository must match the pinned source exactly")
        if self.formal_commit != PINNED_FORMAL_COMMIT:
            raise ValueError("formal_commit must match the pinned source exactly")
        if self.formal_file != PINNED_FORMAL_FILE:
            raise ValueError("formal_file must match ActualStageEstimates.lean exactly")
        if self.theorem_symbol != PINNED_STAGE_ESTIMATES_THEOREM:
            raise ValueError("theorem_symbol must be stageEstimates_of_representations")
        if self.dependency_symbols != REQUIRED_FORMAL_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned stage-estimate chain exactly")

        for name in (
            "q_open_unit_interval_certified",
            "input_normalized_certified",
            "output_strip_controlled_certified",
            "output_mean_balanced_certified",
            "stress_controlled_certified",
            "div_osc_controlled_certified",
            "velocity_increments_controlled_certified",
            "pressure_increments_controlled_certified",
            "forcing_increments_controlled_certified",
            "forcing_increment_margin_certified",
            "ledger_velocity_identification_certified",
            "ledger_pressure_identification_certified",
            "ledger_forcing_identification_certified",
            "theorem_application_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section9ActualStageEstimatesAdmission:
    """Validated formal handoff for literal finite-stage quantitative bounds."""

    witness: Section9ActualStageEstimatesWitness

    def __post_init__(self) -> None:
        if not isinstance(self.witness, Section9ActualStageEstimatesWitness):
            raise TypeError("witness must be a Section9ActualStageEstimatesWitness")

    @property
    def literal_finite_stage_estimates_admitted(self) -> bool:
        return True

    @property
    def canonical_ledger_families_admitted(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def stage_estimates_theorem_machine_replayed(self) -> bool:
        return False

    @property
    def actual_run_data_materialized(self) -> bool:
        return False

    @property
    def actual_physical_data_materialized(self) -> bool:
        return False

    @property
    def finite_stage_fields_materialized(self) -> bool:
        return False

    @property
    def stagewise_quantitative_bounds_materialized(self) -> bool:
        return False

    @property
    def section9_iteration_machine_materialized(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
