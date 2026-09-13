"""Fail-closed admission for the paper's mixed periodic candidate-force theorem.

This module binds one external Lean/formal export to the pinned theorem

    NavierStokes.MixedPeriodicAssembly.exists_candidate_force

at the repository-pinned OpenAI formalization commit.  The theorem is the
paper-specific bridge from the actual mixed potential/direct/pressure inputs to
Section 10 periodization, time activation, residual-limit forcing, divergence
closure, and the blow-up path.  In particular it does materially more than the
repository's generic numerical ``reconstruct_forcing_numeric`` diagnostic.

The Python layer deliberately does not replay Lean, invent the incoming fields,
materialize the locally-uniform residual limits, or evaluate the constructed
force.  It accepts only exact ``lean-formal-export`` metadata for one theorem
application.  Therefore this remains a provenance gate and cannot promote the
runtime reconstruction to paper-exact velocity or forcing.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FORMAL_FILE = "NavierStokes/MixedPeriodicAssembly.lean"
PINNED_CANDIDATE_FORCE_THEOREM = (
    "NavierStokes.MixedPeriodicAssembly.exists_candidate_force"
)
PINNED_ENDPOINT_TIME = Fraction(1, 1)
PRODUCER_KIND = "lean-formal-export"

REQUIRED_FORMAL_SYMBOLS = (
    "NavierStokes.MixedPeriodicAssembly.velocity",
    "NavierStokes.MixedPeriodicAssembly.cutVelocity",
    "NavierStokes.MixedPeriodicAssembly.periodicVelocity",
    "NavierStokes.MixedPeriodicAssembly.originalResidual",
    "NavierStokes.MixedPeriodicAssembly.cutResidual",
    "NavierStokes.MixedPeriodicAssembly.periodicResidual",
    "NavierStokes.MixedPeriodicAssembly.boundaryLimits",
    "NavierStokes.MixedPeriodicAssembly.boundaryLimits_locallyUniform",
    "NavierStokes.MixedPeriodicAssembly.periodicVelocity_divergence_free",
    "NavierStokes.MixedPeriodicAssembly.periodicVelocity_speed_unbounded",
    "NavierStokes.CandidateFromLimits.force",
    "NavierStokes.CandidateFromLimits.force_smooth",
    "NavierStokes.CandidateFromLimits.force_boundary_jets",
    "NavierStokes.CandidateFromLimits.candidate_properties",
    "NavierStokes.TimeLocalization.activatedVelocity",
    "NavierStokes.TimeLocalization.activatedPressure",
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
class Section10MixedPeriodicCandidateForceWitness:
    """Opaque identities for one exact ``exists_candidate_force`` application.

    The identities name the actual incoming mixed fields and the exact objects
    constructed by the pinned formal theorem.  They are intentionally opaque:
    Python must not substitute a manufactured profile, fitted schedule, generic
    cutoff, or residual-defined force for any of these objects.
    """

    potential_field_id: str
    direct_field_id: str
    pressure_field_id: str
    periodic_velocity_id: str
    periodic_pressure_id: str
    activated_velocity_id: str
    activated_pressure_id: str
    boundary_limits_family_id: str
    candidate_force_id: str
    application_id: str
    producer_kind: str
    provenance: str

    endpoint_time: Fraction

    potential_open_past_smooth_certified: bool
    direct_open_past_smooth_certified: bool
    pressure_open_past_smooth_certified: bool
    localized_direct_divergence_zero_certified: bool
    original_residual_vanishing_joint_jets_certified: bool
    potential_away_extensions_certified: bool
    direct_away_extensions_certified: bool
    pressure_away_extensions_certified: bool
    axis_speed_tends_to_infinity_certified: bool

    periodic_velocity_identity_certified: bool
    periodic_pressure_identity_certified: bool
    activated_velocity_identity_certified: bool
    activated_pressure_identity_certified: bool
    boundary_limits_identity_certified: bool
    theorem_application_certified: bool
    candidate_properties_output_certified: bool
    force_global_contdiff_output_certified: bool
    force_endpoint_jets_equal_boundary_limits_certified: bool

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_FORMAL_FILE
    theorem_symbol: str = PINNED_CANDIDATE_FORCE_THEOREM
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    def __post_init__(self) -> None:
        for name in (
            "potential_field_id",
            "direct_field_id",
            "pressure_field_id",
            "periodic_velocity_id",
            "periodic_pressure_id",
            "activated_velocity_id",
            "activated_pressure_id",
            "boundary_limits_family_id",
            "candidate_force_id",
            "application_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        if self.producer_kind != PRODUCER_KIND:
            raise ValueError(
                "producer_kind must be lean-formal-export; "
                "sampled/fitted/numeric-scan evidence is rejected"
            )
        if not isinstance(self.endpoint_time, Fraction):
            raise TypeError("endpoint_time must be fractions.Fraction")
        if self.endpoint_time != PINNED_ENDPOINT_TIME:
            raise ValueError("endpoint_time must be the paper endpoint T=1")
        if self.formal_repository != PINNED_FORMAL_REPOSITORY:
            raise ValueError("formal_repository must match the pinned source exactly")
        if self.formal_commit != PINNED_FORMAL_COMMIT:
            raise ValueError("formal_commit must match the pinned source exactly")
        if self.formal_file != PINNED_FORMAL_FILE:
            raise ValueError("formal_file must match MixedPeriodicAssembly.lean exactly")
        if self.theorem_symbol != PINNED_CANDIDATE_FORCE_THEOREM:
            raise ValueError(
                "theorem_symbol must be MixedPeriodicAssembly.exists_candidate_force"
            )
        if self.dependency_symbols != REQUIRED_FORMAL_SYMBOLS:
            raise ValueError(
                "dependency_symbols must match the pinned mixed candidate-force chain exactly"
            )

        for name in (
            "potential_open_past_smooth_certified",
            "direct_open_past_smooth_certified",
            "pressure_open_past_smooth_certified",
            "localized_direct_divergence_zero_certified",
            "original_residual_vanishing_joint_jets_certified",
            "potential_away_extensions_certified",
            "direct_away_extensions_certified",
            "pressure_away_extensions_certified",
            "axis_speed_tends_to_infinity_certified",
            "periodic_velocity_identity_certified",
            "periodic_pressure_identity_certified",
            "activated_velocity_identity_certified",
            "activated_pressure_identity_certified",
            "boundary_limits_identity_certified",
            "theorem_application_certified",
            "candidate_properties_output_certified",
            "force_global_contdiff_output_certified",
            "force_endpoint_jets_equal_boundary_limits_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section10MixedPeriodicCandidateForceAdmission:
    """Validated formal handoff for the paper-specific candidate-force theorem."""

    witness: Section10MixedPeriodicCandidateForceWitness

    def __post_init__(self) -> None:
        if not isinstance(self.witness, Section10MixedPeriodicCandidateForceWitness):
            raise TypeError(
                "witness must be a Section10MixedPeriodicCandidateForceWitness"
            )

    @property
    def endpoint_time(self) -> Fraction:
        return self.witness.endpoint_time

    @property
    def mixed_periodic_candidate_force_theorem_admitted(self) -> bool:
        return True

    @property
    def candidate_properties_theorem_output_admitted(self) -> bool:
        return True

    @property
    def smooth_force_theorem_output_admitted(self) -> bool:
        return True

    @property
    def force_endpoint_jet_theorem_output_admitted(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    # Runtime truth remains fail-closed until the actual formal witnesses are
    # replayed/materialized end to end in this repository.
    @property
    def theorem_machine_replayed(self) -> bool:
        return False

    @property
    def actual_input_fields_materialized(self) -> bool:
        return False

    @property
    def residual_limit_values_materialized(self) -> bool:
        return False

    @property
    def candidate_force_field_materialized(self) -> bool:
        return False

    @property
    def actual_section9_sequence_verified(self) -> bool:
        return False

    @property
    def section9_field_smooth_extension_through_t1_constructed(self) -> bool:
        return False

    @property
    def residual_artifact_ready(self) -> bool:
        return False

    @property
    def forcing_artifact_ready(self) -> bool:
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
