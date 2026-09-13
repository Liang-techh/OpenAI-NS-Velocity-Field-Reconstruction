"""Fail-closed handoff from the periodic Section 10 candidate to the compact R3 model.

This module binds an external Lean/formal export to the pinned theorem

    NavierStokes.R3CompactCandidate.of_localized_fields

at the repository-pinned OpenAI formalization commit.  The theorem reuses the
same raw potential/direct/pressure data as the periodic candidate, removes the
periodic copies by keeping the fixed Section 10 localized fields, applies the
paper's larger compact force cutoff, and returns the whole-space properties:
smoothness, compact spatial support, zero initial velocity, divergence-free
closure, the exact forced Navier--Stokes equation, and the speed-unbounded path.

The Python layer does not replay Lean or materialize any of those fields.  It
accepts only exact ``lean-formal-export`` metadata and requires the periodic
candidate identity to be the one already admitted by
``section10_mixed_periodic_candidate_force``.  Uniform finite energy is a
separate theorem-level consequence and intentionally remains unclaimed here.
"""
from __future__ import annotations

from dataclasses import dataclass

from .section10_mixed_periodic_candidate_force import (
    PRODUCER_KIND,
    Section10MixedPeriodicCandidateForceAdmission,
)


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FORMAL_FILE = "NavierStokes/R3CompactCandidate.lean"
PINNED_R3_COMPACT_THEOREM = "NavierStokes.R3CompactCandidate.of_localized_fields"

REQUIRED_FORMAL_SYMBOLS = (
    "NavierStokes.R3CompactCandidate.velocity",
    "NavierStokes.R3CompactCandidate.pressure",
    "NavierStokes.R3CompactCandidate.periodicVelocity",
    "NavierStokes.R3CompactCandidate.periodicPressure",
    "NavierStokes.R3CompactCandidate.velocity_supported",
    "NavierStokes.R3CompactCandidate.pressure_supported",
    "NavierStokes.R3CompactCandidate.velocity_locally_eq",
    "NavierStokes.R3CompactCandidate.pressure_locally_eq",
    "NavierStokes.R3CompactCandidate.compactForce",
    "NavierStokes.R3CompactCandidate.compactForce_smooth",
    "NavierStokes.R3CompactCandidate.compactForce_supported",
    "NavierStokes.R3CompactCandidate.compactForce_time_support",
    "NavierStokes.R3CompactCandidate.of_periodic_local_model",
    "NavierStokes.R3CompactCandidate.of_localized_fields",
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
class Section10R3CompactCandidateWitness:
    """Opaque identities for one exact ``of_localized_fields`` application."""

    periodic_candidate: Section10MixedPeriodicCandidateForceAdmission

    potential_field_id: str
    direct_field_id: str
    pressure_field_id: str
    periodic_candidate_velocity_id: str
    periodic_candidate_pressure_id: str
    periodic_candidate_force_id: str
    compact_velocity_id: str
    compact_pressure_id: str
    compact_force_id: str
    support_cylinder_id: str
    outer_force_support_id: str
    application_id: str
    producer_kind: str
    provenance: str

    candidate_properties_input_identity_certified: bool
    localized_velocity_identity_certified: bool
    localized_pressure_identity_certified: bool
    compact_force_identity_certified: bool
    support_cylinder_identity_certified: bool
    outer_force_support_identity_certified: bool
    velocity_supported_certified: bool
    pressure_supported_certified: bool
    velocity_locally_eq_periodic_certified: bool
    pressure_locally_eq_periodic_certified: bool
    theorem_application_certified: bool

    output_velocity_smooth_certified: bool
    output_pressure_smooth_certified: bool
    output_force_smooth_certified: bool
    output_velocity_compact_support_certified: bool
    output_pressure_compact_support_certified: bool
    output_force_compact_support_certified: bool
    output_zero_initial_velocity_certified: bool
    output_force_time_support_certified: bool
    output_divergence_free_certified: bool
    output_navier_stokes_certified: bool
    output_speed_unbounded_certified: bool

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_FORMAL_FILE
    theorem_symbol: str = PINNED_R3_COMPACT_THEOREM
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    def __post_init__(self) -> None:
        if not isinstance(
            self.periodic_candidate, Section10MixedPeriodicCandidateForceAdmission
        ):
            raise TypeError(
                "periodic_candidate must be a Section10MixedPeriodicCandidateForceAdmission"
            )

        for name in (
            "potential_field_id",
            "direct_field_id",
            "pressure_field_id",
            "periodic_candidate_velocity_id",
            "periodic_candidate_pressure_id",
            "periodic_candidate_force_id",
            "compact_velocity_id",
            "compact_pressure_id",
            "compact_force_id",
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
            raise ValueError("formal_file must match R3CompactCandidate.lean exactly")
        if self.theorem_symbol != PINNED_R3_COMPACT_THEOREM:
            raise ValueError("theorem_symbol must be R3CompactCandidate.of_localized_fields")
        if self.dependency_symbols != REQUIRED_FORMAL_SYMBOLS:
            raise ValueError("dependency_symbols must match the pinned compact-candidate chain")

        previous = self.periodic_candidate.witness
        exact_identity_pairs = (
            ("potential_field_id", self.potential_field_id, previous.potential_field_id),
            ("direct_field_id", self.direct_field_id, previous.direct_field_id),
            ("pressure_field_id", self.pressure_field_id, previous.pressure_field_id),
            (
                "periodic_candidate_velocity_id",
                self.periodic_candidate_velocity_id,
                previous.activated_velocity_id,
            ),
            (
                "periodic_candidate_pressure_id",
                self.periodic_candidate_pressure_id,
                previous.activated_pressure_id,
            ),
            (
                "periodic_candidate_force_id",
                self.periodic_candidate_force_id,
                previous.candidate_force_id,
            ),
        )
        for name, supplied, expected in exact_identity_pairs:
            if supplied != expected:
                raise ValueError(f"{name} must match the admitted periodic candidate exactly")

        for name in (
            "candidate_properties_input_identity_certified",
            "localized_velocity_identity_certified",
            "localized_pressure_identity_certified",
            "compact_force_identity_certified",
            "support_cylinder_identity_certified",
            "outer_force_support_identity_certified",
            "velocity_supported_certified",
            "pressure_supported_certified",
            "velocity_locally_eq_periodic_certified",
            "pressure_locally_eq_periodic_certified",
            "theorem_application_certified",
            "output_velocity_smooth_certified",
            "output_pressure_smooth_certified",
            "output_force_smooth_certified",
            "output_velocity_compact_support_certified",
            "output_pressure_compact_support_certified",
            "output_force_compact_support_certified",
            "output_zero_initial_velocity_certified",
            "output_force_time_support_certified",
            "output_divergence_free_certified",
            "output_navier_stokes_certified",
            "output_speed_unbounded_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section10R3CompactCandidateAdmission:
    """Validated formal handoff for the compact whole-space candidate theorem."""

    witness: Section10R3CompactCandidateWitness

    def __post_init__(self) -> None:
        if not isinstance(self.witness, Section10R3CompactCandidateWitness):
            raise TypeError("witness must be a Section10R3CompactCandidateWitness")

    @property
    def r3_compact_candidate_theorem_admitted(self) -> bool:
        return True

    @property
    def compact_support_theorem_output_admitted(self) -> bool:
        return True

    @property
    def divergence_free_theorem_output_admitted(self) -> bool:
        return True

    @property
    def navier_stokes_theorem_output_admitted(self) -> bool:
        return True

    @property
    def blow_up_theorem_output_admitted(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    # The formal theorem output is admitted, but no runtime field is materialized.
    @property
    def theorem_machine_replayed(self) -> bool:
        return False

    @property
    def r3_compact_fields_materialized(self) -> bool:
        return False

    @property
    def compact_force_field_materialized(self) -> bool:
        return False

    @property
    def residual_artifact_ready(self) -> bool:
        return False

    @property
    def forcing_artifact_ready(self) -> bool:
        return False

    @property
    def compact_support_runtime_verified(self) -> bool:
        return False

    @property
    def divergence_free_closure_verified(self) -> bool:
        return False

    @property
    def blow_up_closure_verified(self) -> bool:
        return False

    @property
    def finite_energy_closure_verified(self) -> bool:
        return False

    @property
    def endpoint_residual_closure_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
