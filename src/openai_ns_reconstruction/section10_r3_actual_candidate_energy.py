"""Fail-closed handoff to the final whole-space candidate energy theorem.

The pinned OpenAI formalization contains
``NavierStokesR3.ActualCandidate.of_localized_fields``.  Starting from the
same periodic Section 10 candidate already used by ``R3CompactCandidate``, the
theorem keeps the localized whole-space velocity/pressure, makes the compact
force globally smooth with ``PositiveTimeForce.force``, proves the exact
forced equation, and invokes ``CompactEnergy.uniform_finite_energy`` to obtain
one kinetic-energy bound on the full pre-singular interval ``0 <= t < 1``.

This is stronger than applying compact-slab continuity separately on every
``[0,T]`` with ``T<1``: the theorem returns ``UniformFiniteEnergy (Ico 0 1)``.
Python does not replay Lean or materialize the fields.  This module therefore
admits only a pinned ``lean-formal-export`` whose identities agree with the
already-landed ``Section10R3CompactCandidateAdmission``.  Formal theorem
outputs are never promoted to runtime reconstruction truth.
"""
from __future__ import annotations

from dataclasses import dataclass

from .section10_mixed_periodic_candidate_force import PRODUCER_KIND
from .section10_r3_compact_candidate import Section10R3CompactCandidateAdmission


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FORMAL_FILE = "NavierStokes/R3/ActualCandidate.lean"
PINNED_ACTUAL_CANDIDATE_THEOREM = "NavierStokesR3.ActualCandidate.of_localized_fields"
PINNED_ENERGY_FILE = "NavierStokes/R3/CompactEnergy.lean"
PINNED_ENERGY_THEOREM = "NavierStokesR3.CompactEnergy.uniform_finite_energy"

REQUIRED_FORMAL_SYMBOLS = (
    "NavierStokesR3.ActualCandidate.localized_residual_zero_early",
    "NavierStokesR3.ActualCandidate.localized_velocity_tsupport",
    "NavierStokesR3.ActualCandidate.localized_pressure_tsupport",
    "NavierStokesR3.ActualCandidate.compactForce_smooth",
    "NavierStokesR3.ActualCandidate.compactForce_supported_all_times",
    "NavierStokes.R3CompactCandidate.compactForce",
    "NavierStokes.R3CompactCandidate.outerSupport",
    "NavierStokesR3.PositiveTimeForce.force",
    "NavierStokesR3.PositiveTimeForce.force_contDiff",
    "NavierStokesR3.PositiveTimeForce.force_compactPositiveTimeSupport",
    PINNED_ENERGY_THEOREM,
    PINNED_ACTUAL_CANDIDATE_THEOREM,
)


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be formal-export certified true")
    return True


@dataclass(frozen=True)
class Section10R3ActualCandidateEnergyWitness:
    """Opaque identities for one exact final R3 candidate theorem application."""

    compact_candidate: Section10R3CompactCandidateAdmission

    periodic_candidate_force_id: str
    compact_velocity_id: str
    compact_pressure_id: str
    compact_force_id: str
    positive_time_force_id: str
    support_cylinder_id: str
    outer_force_support_id: str
    application_id: str
    producer_kind: str
    provenance: str

    periodic_force_global_contdiff_certified: bool
    localized_velocity_tsupport_certified: bool
    localized_pressure_tsupport_certified: bool
    compact_force_global_contdiff_certified: bool
    compact_force_supported_all_times_certified: bool
    positive_time_force_global_contdiff_certified: bool
    positive_time_force_compact_positive_time_support_certified: bool
    early_residual_zero_certified: bool
    positive_time_force_equation_certified: bool
    compact_energy_uniform_finite_energy_certified: bool
    theorem_application_certified: bool

    output_candidate_properties_viscosity_one_certified: bool
    output_uniform_finite_energy_Ico_0_1_certified: bool
    output_force_smooth_certified: bool
    output_force_compact_positive_time_support_certified: bool
    output_divergence_free_certified: bool
    output_navier_stokes_certified: bool
    output_speed_unbounded_certified: bool

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_FORMAL_FILE
    theorem_symbol: str = PINNED_ACTUAL_CANDIDATE_THEOREM
    energy_file: str = PINNED_ENERGY_FILE
    energy_theorem_symbol: str = PINNED_ENERGY_THEOREM
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    def __post_init__(self) -> None:
        if not isinstance(self.compact_candidate, Section10R3CompactCandidateAdmission):
            raise TypeError(
                "compact_candidate must be a Section10R3CompactCandidateAdmission"
            )

        for name in (
            "periodic_candidate_force_id",
            "compact_velocity_id",
            "compact_pressure_id",
            "compact_force_id",
            "positive_time_force_id",
            "support_cylinder_id",
            "outer_force_support_id",
            "application_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        if self.producer_kind != PRODUCER_KIND:
            raise ValueError(
                "producer_kind must be lean-formal-export; sampled/fitted/numeric evidence is rejected"
            )
        if self.formal_repository != PINNED_FORMAL_REPOSITORY:
            raise ValueError("formal_repository must match the pinned source exactly")
        if self.formal_commit != PINNED_FORMAL_COMMIT:
            raise ValueError("formal_commit must match the pinned source exactly")
        if self.formal_file != PINNED_FORMAL_FILE:
            raise ValueError("formal_file must match R3/ActualCandidate.lean exactly")
        if self.theorem_symbol != PINNED_ACTUAL_CANDIDATE_THEOREM:
            raise ValueError("theorem_symbol must be ActualCandidate.of_localized_fields")
        if self.energy_file != PINNED_ENERGY_FILE:
            raise ValueError("energy_file must match R3/CompactEnergy.lean exactly")
        if self.energy_theorem_symbol != PINNED_ENERGY_THEOREM:
            raise ValueError("energy_theorem_symbol must be CompactEnergy.uniform_finite_energy")
        if self.dependency_symbols != REQUIRED_FORMAL_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned final-candidate chain")

        compact = self.compact_candidate.witness
        periodic = compact.periodic_candidate.witness
        exact_identity_pairs = (
            (
                "periodic_candidate_force_id",
                self.periodic_candidate_force_id,
                periodic.candidate_force_id,
            ),
            ("compact_velocity_id", self.compact_velocity_id, compact.compact_velocity_id),
            ("compact_pressure_id", self.compact_pressure_id, compact.compact_pressure_id),
            ("compact_force_id", self.compact_force_id, compact.compact_force_id),
            ("support_cylinder_id", self.support_cylinder_id, compact.support_cylinder_id),
            ("outer_force_support_id", self.outer_force_support_id, compact.outer_force_support_id),
        )
        for name, supplied, expected in exact_identity_pairs:
            if supplied != expected:
                raise ValueError(f"{name} must match the admitted compact candidate exactly")

        if periodic.force_global_contdiff_output_certified is not True:
            raise ValueError("the admitted periodic candidate must certify global smooth force")

        for name in (
            "periodic_force_global_contdiff_certified",
            "localized_velocity_tsupport_certified",
            "localized_pressure_tsupport_certified",
            "compact_force_global_contdiff_certified",
            "compact_force_supported_all_times_certified",
            "positive_time_force_global_contdiff_certified",
            "positive_time_force_compact_positive_time_support_certified",
            "early_residual_zero_certified",
            "positive_time_force_equation_certified",
            "compact_energy_uniform_finite_energy_certified",
            "theorem_application_certified",
            "output_candidate_properties_viscosity_one_certified",
            "output_uniform_finite_energy_Ico_0_1_certified",
            "output_force_smooth_certified",
            "output_force_compact_positive_time_support_certified",
            "output_divergence_free_certified",
            "output_navier_stokes_certified",
            "output_speed_unbounded_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section10R3ActualCandidateEnergyAdmission:
    """Validated formal handoff for the final R3 candidate energy theorem."""

    witness: Section10R3ActualCandidateEnergyWitness

    def __post_init__(self) -> None:
        if not isinstance(self.witness, Section10R3ActualCandidateEnergyWitness):
            raise TypeError("witness must be a Section10R3ActualCandidateEnergyWitness")

    @property
    def actual_candidate_theorem_admitted(self) -> bool:
        return True

    @property
    def uniform_finite_energy_theorem_output_admitted(self) -> bool:
        return True

    @property
    def force_smooth_theorem_output_admitted(self) -> bool:
        return True

    @property
    def compact_positive_time_force_support_theorem_output_admitted(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    # Formal theorem outputs stay separate from executable reconstruction truth.
    @property
    def theorem_machine_replayed(self) -> bool:
        return False

    @property
    def actual_fields_materialized(self) -> bool:
        return False

    @property
    def actual_section9_sequence_verified(self) -> bool:
        return False

    @property
    def section9_all_order_endpoint_limits_verified(self) -> bool:
        return False

    @property
    def residual_artifact_ready(self) -> bool:
        return False

    @property
    def forcing_artifact_ready(self) -> bool:
        return False

    @property
    def endpoint_residual_closure_verified(self) -> bool:
        return False

    @property
    def compact_support_runtime_verified(self) -> bool:
        return False

    @property
    def divergence_free_closure_verified(self) -> bool:
        return False

    @property
    def finite_energy_closure_verified(self) -> bool:
        return False

    @property
    def blow_up_closure_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False

    @property
    def full_reconstruction(self) -> bool:
        return False
