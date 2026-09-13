"""Fail-closed admission for the paper's actual spacetime endpoint extension.

This module binds one external Lean/formal export to the pinned constructive
endpoint theorem

    NavierStokes.SpacetimeGluing.exists_smooth_periodic_extension_of_limits

at the paper endpoint ``T = 1``.  The theorem does materially more than an
arbitrary time window: from an actual joint derivative family on ``t < T`` and
locally-uniform limits of every derivative, it constructs a jointly smooth
periodic spacetime extension, preserves the open-past field, preserves every
mixed endpoint jet, and vanishes after ``T + 1``.

The Python layer deliberately does not replay Lean, invent endpoint limits, or
materialize the extension field.  It accepts only exact ``lean-formal-export``
metadata for the theorem application and therefore remains a provenance gate,
not a paper-exact velocity/forcing constructor.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FORMAL_FILE = "NavierStokes/SpacetimeGluing.lean"
PINNED_EXTENSION_THEOREM = (
    "NavierStokes.SpacetimeGluing.exists_smooth_periodic_extension_of_limits"
)
PINNED_SMOOTH_EXTENSION = "NavierStokes.SpacetimeGluing.smoothExtension"
PINNED_SMOOTH_EXTENSION_CONTDIFF = (
    "NavierStokes.SpacetimeGluing.smoothExtension_contDiff"
)
PINNED_SMOOTH_EXTENSION_EQ_PAST = (
    "NavierStokes.SpacetimeGluing.smoothExtension_eqOn_past"
)
PINNED_SMOOTH_EXTENSION_ZERO_FROM = (
    "NavierStokes.SpacetimeGluing.smoothExtension_zero_from"
)
PINNED_SMOOTH_EXTENSION_JETS = (
    "NavierStokes.SpacetimeGluing.smoothExtension_iteratedFDeriv"
)
PINNED_ENDPOINT_JETS = "NavierStokes.SpacetimeEndpoint.boundary_jets_eq_limits"
PINNED_JOINT_LEFT_EXTENSION = (
    "NavierStokes.SpacetimeEndpoint.contDiffOn_joint_extension"
)
PINNED_BOREL_RIGHT_EXTENSION = "NavierStokes.SpatialBorelExtension.rightExtension"
PINNED_BOREL_RIGHT_JETS = (
    "NavierStokes.SpatialBorelExtension.rightExtension_right_jets"
)
PINNED_ENDPOINT_TIME = Fraction(1, 1)
PRODUCER_KIND = "lean-formal-export"

REQUIRED_FORMAL_SYMBOLS = (
    PINNED_EXTENSION_THEOREM,
    PINNED_SMOOTH_EXTENSION,
    PINNED_SMOOTH_EXTENSION_CONTDIFF,
    PINNED_SMOOTH_EXTENSION_EQ_PAST,
    PINNED_SMOOTH_EXTENSION_ZERO_FROM,
    PINNED_SMOOTH_EXTENSION_JETS,
    PINNED_ENDPOINT_JETS,
    PINNED_JOINT_LEFT_EXTENSION,
    PINNED_BOREL_RIGHT_EXTENSION,
    PINNED_BOREL_RIGHT_JETS,
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
class Section10SpacetimeEndpointExtensionWitness:
    """Opaque identities for one exact theorem application at ``T=1``.

    ``section9_field_id`` must identify the actual left/open-past field supplied
    to the formal theorem application.  ``joint_jet_family_id`` identifies the
    theorem's ``J`` family and ``endpoint_limit_jet_id`` its locally-uniform
    limit family ``L``.  The remaining identities name the constructed
    closed-past and global extension objects exported by the formal runner.
    """

    section9_field_id: str
    joint_jet_family_id: str
    endpoint_limit_jet_id: str
    closed_past_extension_id: str
    global_extension_id: str
    application_id: str
    producer_kind: str
    provenance: str

    endpoint_time: Fraction
    actual_section9_field_identification_certified: bool
    open_past_value_identity_certified: bool
    derivative_recurrence_certified: bool
    locally_uniform_left_limits_certified: bool
    unit_spatial_periods_on_open_past_certified: bool
    theorem_application_certified: bool
    global_contdiff_output_certified: bool
    agrees_with_open_past_output_certified: bool
    global_unit_spatial_periods_output_certified: bool
    right_tail_zero_from_t_plus_one_certified: bool
    endpoint_mixed_jets_equal_limits_certified: bool

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_FORMAL_FILE
    theorem_symbol: str = PINNED_EXTENSION_THEOREM
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    def __post_init__(self) -> None:
        for name in (
            "section9_field_id",
            "joint_jet_family_id",
            "endpoint_limit_jet_id",
            "closed_past_extension_id",
            "global_extension_id",
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
            raise ValueError("formal_file must match SpacetimeGluing.lean exactly")
        if self.theorem_symbol != PINNED_EXTENSION_THEOREM:
            raise ValueError(
                "theorem_symbol must be SpacetimeGluing.exists_smooth_periodic_extension_of_limits"
            )
        if self.dependency_symbols != REQUIRED_FORMAL_SYMBOLS:
            raise ValueError(
                "dependency_symbols must match the pinned constructive endpoint chain exactly"
            )

        for name in (
            "actual_section9_field_identification_certified",
            "open_past_value_identity_certified",
            "derivative_recurrence_certified",
            "locally_uniform_left_limits_certified",
            "unit_spatial_periods_on_open_past_certified",
            "theorem_application_certified",
            "global_contdiff_output_certified",
            "agrees_with_open_past_output_certified",
            "global_unit_spatial_periods_output_certified",
            "right_tail_zero_from_t_plus_one_certified",
            "endpoint_mixed_jets_equal_limits_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section10SpacetimeEndpointExtensionAdmission:
    """Validated formal handoff for the constructive spacetime extension theorem."""

    witness: Section10SpacetimeEndpointExtensionWitness

    def __post_init__(self) -> None:
        if not isinstance(self.witness, Section10SpacetimeEndpointExtensionWitness):
            raise TypeError(
                "witness must be a Section10SpacetimeEndpointExtensionWitness"
            )

    @property
    def endpoint_time(self) -> Fraction:
        return self.witness.endpoint_time

    @property
    def certified_zero_from_time(self) -> Fraction:
        """The pinned theorem gives zero for all ``t >= T+1``; here ``T=1``."""
        return self.endpoint_time + 1

    @property
    def constructive_spacetime_extension_theorem_admitted(self) -> bool:
        return True

    @property
    def actual_left_endpoint_limit_chain_admitted(self) -> bool:
        return True

    @property
    def global_smooth_periodic_extension_output_admitted(self) -> bool:
        return True

    @property
    def endpoint_mixed_jet_preservation_admitted(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def endpoint_extension_theorem_machine_replayed(self) -> bool:
        return False

    @property
    def section9_endpoint_limit_values_materialized(self) -> bool:
        return False

    @property
    def global_extension_field_materialized(self) -> bool:
        return False

    @property
    def actual_section9_sequence_verified(self) -> bool:
        return False

    @property
    def section9_field_smooth_extension_through_t1_constructed(self) -> bool:
        """Remain fail-closed until the actual field/export is materialized end to end."""
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
    def paper_exact_velocity_available(self) -> bool:
        return False
