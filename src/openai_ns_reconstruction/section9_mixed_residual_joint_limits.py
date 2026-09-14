"""Bind the actual mixed Section 9 residual to all-order endpoint limits.

This module is deliberately narrower than the generic Section 10 endpoint/Borel
infrastructure.  It records the exact Lean theorem chain that can turn the
*paper mixed residual* into the all-order one-sided endpoint limits consumed by
Section 10:

1. ``MixedDiagonalResidual.residual_eq_originalResidual`` identifies the mixed
   physical residual with ``MixedPeriodicAssembly.originalResidual``.
2. ``MixedDiagonalResidual.physical_vanishingJointJets`` proves all-order joint
   jet vanishing for that actual residual from the paper schedule/rate
   hypotheses.
3. ``MixedPeriodicAssembly.cutResidual_vanishingJointJets`` transports the
   all-order vanishing through the fixed Section 10 spatial localization.
4. ``MixedPeriodicAssembly.cutResidual_awayExtensions`` supplies away-from-origin
   one-sided extensions from the actual A/B/P field extensions.
5. ``MixedPeriodicAssembly.boundaryLimits_locallyUniform`` then yields the
   locally-uniform endpoint limits for every natural derivative order.

Admission is intentionally fail-closed.  It accepts only one pinned
``lean-formal-export`` with universal natural-number quantifiers.  Finite stage
prefixes, finite derivative ladders, sampled/numerical evidence, or residual
identities from another source cannot satisfy this interface.

Passing the bridge is still only ``formal-structure`` in this Python repository:
it does not materialize/replay the Lean proof, the infinite Section 9 sequence,
the residual values, the analytic majorant/Borel family required by Section 10,
or the final force.  Consequently all runtime reconstruction truth flags remain
false.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


FORMAL_REPOSITORY = "openai/NavierStokesAndEuler"
FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
MIXED_DIAGONAL_SOURCE = "NavierStokes/MixedDiagonalResidual.lean"
MIXED_PERIODIC_SOURCE = "NavierStokes/MixedPeriodicAssembly.lean"

RESIDUAL_IDENTITY_THEOREM = (
    "NavierStokes.MixedDiagonalResidual.residual_eq_originalResidual"
)
PHYSICAL_VANISHING_JOINT_JETS_THEOREM = (
    "NavierStokes.MixedDiagonalResidual.physical_vanishingJointJets"
)
CUT_RESIDUAL_VANISHING_JOINT_JETS_THEOREM = (
    "NavierStokes.MixedPeriodicAssembly.cutResidual_vanishingJointJets"
)
CUT_RESIDUAL_AWAY_EXTENSIONS_THEOREM = (
    "NavierStokes.MixedPeriodicAssembly.cutResidual_awayExtensions"
)
BOUNDARY_LIMITS_LOCALLY_UNIFORM_THEOREM = (
    "NavierStokes.MixedPeriodicAssembly.boundaryLimits_locallyUniform"
)

SECTION10_ENDPOINT_EXACT = Fraction(1, 1)
_UNIVERSAL_STAGE_RATE_QUANTIFIER = "forall-natural-J-m"
_UNIVERSAL_DERIVATIVE_ORDER_QUANTIFIER = "forall-natural-derivative-orders"
_ACCEPTED_EVIDENCE = "lean-formal-export"


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


@dataclass(frozen=True)
class Section9MixedResidualJointLimitsWitness:
    """Pinned theorem export for one actual mixed Section 9 residual."""

    source_id: str
    source_revision: str
    scale_schedule_id: str
    physical_q_id: str
    potential_family_id: str
    direct_family_id: str
    pressure_family_id: str
    mixed_residual_id: str
    original_residual_id: str
    cut_residual_id: str
    boundary_limits_family_id: str
    provenance: str
    evidence_kind: str = _ACCEPTED_EVIDENCE
    endpoint: Fraction = SECTION10_ENDPOINT_EXACT

    formal_repository: str = FORMAL_REPOSITORY
    formal_commit: str = FORMAL_COMMIT
    mixed_diagonal_source: str = MIXED_DIAGONAL_SOURCE
    mixed_periodic_source: str = MIXED_PERIODIC_SOURCE
    residual_identity_theorem: str = RESIDUAL_IDENTITY_THEOREM
    physical_vanishing_joint_jets_theorem: str = (
        PHYSICAL_VANISHING_JOINT_JETS_THEOREM
    )
    cut_residual_vanishing_joint_jets_theorem: str = (
        CUT_RESIDUAL_VANISHING_JOINT_JETS_THEOREM
    )
    cut_residual_away_extensions_theorem: str = (
        CUT_RESIDUAL_AWAY_EXTENSIONS_THEOREM
    )
    boundary_limits_locally_uniform_theorem: str = (
        BOUNDARY_LIMITS_LOCALLY_UNIFORM_THEOREM
    )

    stage_rate_quantifier: str = _UNIVERSAL_STAGE_RATE_QUANTIFIER
    derivative_order_quantifier: str = _UNIVERSAL_DERIVATIVE_ORDER_QUANTIFIER
    max_stage_index: int | None = None
    max_derivative_order: int | None = None
    finite: bool = False

    physical_h_range_certified: bool = False
    scale_schedule_tends_to_infinity_certified: bool = False
    preterminal_open_neighborhood_certified: bool = False
    all_stage_fields_smooth_certified: bool = False
    gauge_monotone_unbounded_certified: bool = False
    cut_stage_bounds_potential_certified: bool = False
    cut_stage_bounds_direct_certified: bool = False
    cut_stage_bounds_pressure_certified: bool = False
    uncut_velocity_jet_rate_all_J_m_certified: bool = False
    actual_residual_jet_rate_all_J_m_certified: bool = False

    residual_identity_theorem_applied: bool = False
    physical_vanishing_joint_jets_theorem_applied: bool = False
    potential_away_extensions_certified: bool = False
    direct_away_extensions_certified: bool = False
    pressure_away_extensions_certified: bool = False
    cut_residual_vanishing_joint_jets_theorem_applied: bool = False
    cut_residual_away_extensions_theorem_applied: bool = False
    boundary_limits_locally_uniform_theorem_applied: bool = False
    same_actual_field_source_certified: bool = False

    def __post_init__(self) -> None:
        for name in (
            "source_id",
            "source_revision",
            "scale_schedule_id",
            "physical_q_id",
            "potential_family_id",
            "direct_family_id",
            "pressure_family_id",
            "mixed_residual_id",
            "original_residual_id",
            "cut_residual_id",
            "boundary_limits_family_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        if self.evidence_kind != _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be lean-formal-export; "
                "sampled/fitted/numerical evidence is rejected"
            )
        if type(self.endpoint) is not Fraction or self.endpoint != SECTION10_ENDPOINT_EXACT:
            raise ValueError("endpoint must be the exact rational T=1")

        pinned_text = {
            "formal_repository": FORMAL_REPOSITORY,
            "formal_commit": FORMAL_COMMIT,
            "mixed_diagonal_source": MIXED_DIAGONAL_SOURCE,
            "mixed_periodic_source": MIXED_PERIODIC_SOURCE,
            "residual_identity_theorem": RESIDUAL_IDENTITY_THEOREM,
            "physical_vanishing_joint_jets_theorem": (
                PHYSICAL_VANISHING_JOINT_JETS_THEOREM
            ),
            "cut_residual_vanishing_joint_jets_theorem": (
                CUT_RESIDUAL_VANISHING_JOINT_JETS_THEOREM
            ),
            "cut_residual_away_extensions_theorem": (
                CUT_RESIDUAL_AWAY_EXTENSIONS_THEOREM
            ),
            "boundary_limits_locally_uniform_theorem": (
                BOUNDARY_LIMITS_LOCALLY_UNIFORM_THEOREM
            ),
        }
        for name, expected in pinned_text.items():
            if getattr(self, name) != expected:
                raise ValueError(f"{name} must equal pinned formal source value {expected!r}")

        if self.stage_rate_quantifier != _UNIVERSAL_STAGE_RATE_QUANTIFIER:
            raise ValueError("stage_rate_quantifier must universally quantify all natural J,m")
        if self.derivative_order_quantifier != _UNIVERSAL_DERIVATIVE_ORDER_QUANTIFIER:
            raise ValueError(
                "derivative_order_quantifier must universally quantify all natural orders"
            )
        if self.max_stage_index is not None:
            raise ValueError("max_stage_index must be None; finite stage prefixes are rejected")
        if self.max_derivative_order is not None:
            raise ValueError(
                "max_derivative_order must be None; finite endpoint ladders are rejected"
            )
        if self.finite is not False:
            raise ValueError("finite must be False for the all-order theorem bridge")

        for name in (
            "physical_h_range_certified",
            "scale_schedule_tends_to_infinity_certified",
            "preterminal_open_neighborhood_certified",
            "all_stage_fields_smooth_certified",
            "gauge_monotone_unbounded_certified",
            "cut_stage_bounds_potential_certified",
            "cut_stage_bounds_direct_certified",
            "cut_stage_bounds_pressure_certified",
            "uncut_velocity_jet_rate_all_J_m_certified",
            "actual_residual_jet_rate_all_J_m_certified",
            "residual_identity_theorem_applied",
            "physical_vanishing_joint_jets_theorem_applied",
            "potential_away_extensions_certified",
            "direct_away_extensions_certified",
            "pressure_away_extensions_certified",
            "cut_residual_vanishing_joint_jets_theorem_applied",
            "cut_residual_away_extensions_theorem_applied",
            "boundary_limits_locally_uniform_theorem_applied",
            "same_actual_field_source_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def field_ids(self) -> tuple[str, str, str]:
        return self.potential_family_id, self.direct_family_id, self.pressure_family_id


@dataclass(frozen=True)
class Section9MixedResidualJointLimitsCertificate:
    """Formal theorem-chain binding; runtime reconstruction remains unavailable."""

    witness: Section9MixedResidualJointLimitsWitness
    expected_source_key: tuple[str, str]
    expected_field_ids: tuple[str, str, str]
    expected_mixed_residual_id: str
    expected_original_residual_id: str
    expected_cut_residual_id: str
    status: str = "formal-structure"

    actual_mixed_residual_identity_formally_bound: bool = True
    original_residual_vanishing_joint_jets_formally_bound: bool = True
    cut_residual_boundary_limits_locally_uniform_formally_bound: bool = True

    actual_section9_sequence_verified: bool = False
    section9_all_order_endpoint_limits_verified: bool = False
    residual_artifact_ready: bool = False
    forcing_artifact_ready: bool = False
    endpoint_residual_closure_verified: bool = False
    paper_exact_velocity_available: bool = False
    full_reconstruction: bool = False

    @property
    def formal_joint_limits_bridge_ready(self) -> bool:
        return (
            self.witness.source_key == self.expected_source_key
            and self.witness.field_ids == self.expected_field_ids
            and self.witness.mixed_residual_id == self.expected_mixed_residual_id
            and self.witness.original_residual_id == self.expected_original_residual_id
            and self.witness.cut_residual_id == self.expected_cut_residual_id
            and self.actual_mixed_residual_identity_formally_bound
            and self.original_residual_vanishing_joint_jets_formally_bound
            and self.cut_residual_boundary_limits_locally_uniform_formally_bound
            and not self.actual_section9_sequence_verified
            and not self.section9_all_order_endpoint_limits_verified
            and not self.residual_artifact_ready
            and not self.forcing_artifact_ready
            and not self.endpoint_residual_closure_verified
            and not self.paper_exact_velocity_available
            and not self.full_reconstruction
        )


def admit_section9_mixed_residual_joint_limits(
    witness: Section9MixedResidualJointLimitsWitness,
    *,
    expected_source_id: str,
    expected_source_revision: str,
    expected_potential_family_id: str,
    expected_direct_family_id: str,
    expected_pressure_family_id: str,
    expected_mixed_residual_id: str,
    expected_original_residual_id: str,
    expected_cut_residual_id: str,
) -> Section9MixedResidualJointLimitsCertificate:
    """Bind one selected actual mixed residual to the pinned endpoint theorem chain."""

    if not isinstance(witness, Section9MixedResidualJointLimitsWitness):
        raise TypeError("witness must be a Section9MixedResidualJointLimitsWitness")

    expected_source_key = (
        _nonempty_text(expected_source_id, "expected_source_id"),
        _nonempty_text(expected_source_revision, "expected_source_revision"),
    )
    expected_field_ids = (
        _nonempty_text(expected_potential_family_id, "expected_potential_family_id"),
        _nonempty_text(expected_direct_family_id, "expected_direct_family_id"),
        _nonempty_text(expected_pressure_family_id, "expected_pressure_family_id"),
    )
    expected_mixed = _nonempty_text(
        expected_mixed_residual_id, "expected_mixed_residual_id"
    )
    expected_original = _nonempty_text(
        expected_original_residual_id, "expected_original_residual_id"
    )
    expected_cut = _nonempty_text(expected_cut_residual_id, "expected_cut_residual_id")

    if witness.source_key != expected_source_key:
        raise ValueError("mixed residual theorem export source identity mismatch")
    if witness.field_ids != expected_field_ids:
        raise ValueError("mixed residual theorem export field identity mismatch")
    if witness.mixed_residual_id != expected_mixed:
        raise ValueError("mixed residual identity mismatch")
    if witness.original_residual_id != expected_original:
        raise ValueError("original residual identity mismatch")
    if witness.cut_residual_id != expected_cut:
        raise ValueError("cut residual identity mismatch")

    certificate = Section9MixedResidualJointLimitsCertificate(
        witness=witness,
        expected_source_key=expected_source_key,
        expected_field_ids=expected_field_ids,
        expected_mixed_residual_id=expected_mixed,
        expected_original_residual_id=expected_original,
        expected_cut_residual_id=expected_cut,
    )
    if not certificate.formal_joint_limits_bridge_ready:
        raise ArithmeticError("mixed residual joint-limits bridge invariant failed")
    return certificate
