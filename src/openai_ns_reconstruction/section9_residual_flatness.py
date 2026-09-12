"""Fail-closed Section 9 residual-flatness closure contract.

The pinned formal development already contains the abstract convergence theorem
``DiagonalResidual.allJetsFlat_residual_of_stages``.  It says that a smooth
limit velocity/pressure has an all-jets-flat Navier--Stokes residual once a
family of finite stages satisfies a divergent gain schedule together with
uniformly formulated background, tail, pressure-tail, and stage-residual jet
rates.

This module does *not* construct those finite stages, infer any rate from
sampled data, or claim that the paper's Section 9 iteration has been replayed.
It only provides a strict machine handoff for a future external
``lean-formal-export`` application of that theorem.  Every mathematical object
that Python cannot faithfully encode -- filters, domains, fields, stage
families, and gain/loss functions -- remains an opaque stable identifier.
"""
from __future__ import annotations

from dataclasses import dataclass


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FORMAL_FILE = "NavierStokes/DiagonalResidual.lean"
PINNED_ALL_JETS_FLAT_THEOREM = (
    "NavierStokes.DiagonalResidual.allJetsFlat_residual_of_stages"
)
PINNED_STAGE_JET_RATE_THEOREM = (
    "NavierStokes.DiagonalResidual.residual_jetRate_of_stages"
)
PINNED_RESIDUAL_DIFFERENCE_THEOREM = (
    "NavierStokes.DiagonalResidual.residualDifference_jetRate"
)
REQUIRED_FORMAL_SYMBOLS = (
    PINNED_ALL_JETS_FLAT_THEOREM,
    PINNED_STAGE_JET_RATE_THEOREM,
    PINNED_RESIDUAL_DIFFERENCE_THEOREM,
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
class Section9ResidualFlatnessWitness:
    """Opaque metadata for one formal residual-flatness theorem application.

    The identifier fields name the exact objects used by the external Lean
    application.  They intentionally do not define a Python representation of
    the corresponding functions or filters.
    """

    domain_id: str
    filter_id: str
    scale_q_id: str
    limit_velocity_id: str
    limit_pressure_id: str
    velocity_stage_family_id: str
    pressure_stage_family_id: str
    gain_schedule_id: str
    background_loss_id: str
    tail_loss_id: str
    residual_loss_id: str
    application_id: str
    producer_kind: str
    provenance: str

    domain_open_certified: bool
    filter_eventually_in_domain_certified: bool
    scale_eventually_unit_interval_certified: bool
    limit_velocity_smooth_certified: bool
    limit_pressure_smooth_certified: bool
    velocity_stages_smooth_certified: bool
    pressure_stages_smooth_certified: bool
    gain_tendsto_at_top_certified: bool
    background_stage_jet_rates_certified: bool
    velocity_tail_jet_rates_certified: bool
    pressure_tail_jet_rates_certified: bool
    stage_residual_jet_rates_certified: bool
    theorem_application_certified: bool

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_FORMAL_FILE
    theorem_symbol: str = PINNED_ALL_JETS_FLAT_THEOREM
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    def __post_init__(self) -> None:
        for name in (
            "domain_id",
            "filter_id",
            "scale_q_id",
            "limit_velocity_id",
            "limit_pressure_id",
            "velocity_stage_family_id",
            "pressure_stage_family_id",
            "gain_schedule_id",
            "background_loss_id",
            "tail_loss_id",
            "residual_loss_id",
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
            raise ValueError("formal_file must match DiagonalResidual.lean exactly")
        if self.theorem_symbol != PINNED_ALL_JETS_FLAT_THEOREM:
            raise ValueError("theorem_symbol must be allJetsFlat_residual_of_stages")
        if self.dependency_symbols != REQUIRED_FORMAL_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned residual chain exactly")

        for name in (
            "domain_open_certified",
            "filter_eventually_in_domain_certified",
            "scale_eventually_unit_interval_certified",
            "limit_velocity_smooth_certified",
            "limit_pressure_smooth_certified",
            "velocity_stages_smooth_certified",
            "pressure_stages_smooth_certified",
            "gain_tendsto_at_top_certified",
            "background_stage_jet_rates_certified",
            "velocity_tail_jet_rates_certified",
            "pressure_tail_jet_rates_certified",
            "stage_residual_jet_rates_certified",
            "theorem_application_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section9ResidualFlatnessAdmission:
    """Validated formal handoff for the all-jets-flat residual conclusion."""

    witness: Section9ResidualFlatnessWitness

    def __post_init__(self) -> None:
        if not isinstance(self.witness, Section9ResidualFlatnessWitness):
            raise TypeError("witness must be a Section9ResidualFlatnessWitness")

    @property
    def all_jets_flat_residual_conclusion_admitted(self) -> bool:
        return True

    @property
    def stage_convergence_contract_admitted(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def all_jets_flat_theorem_machine_replayed(self) -> bool:
        return False

    @property
    def section9_iteration_machine_materialized(self) -> bool:
        return False

    @property
    def stage_fields_materialized(self) -> bool:
        return False

    @property
    def stage_rate_bounds_materialized(self) -> bool:
        return False

    @property
    def limit_velocity_materialized(self) -> bool:
        return False

    @property
    def navier_stokes_residual_values_materialized(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
