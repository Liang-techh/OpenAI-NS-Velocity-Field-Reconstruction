"""Bind the pinned Section 9 physical schedule theorem to joint residual flatness.

This module records one upstream theorem-shaped contract, deliberately before
Section 10 endpoint/Borel plumbing.  The pinned Lean theorem
``MixedDiagonalResidual.exists_physical_schedule_residual_zero`` selects one
common rapidly growing schedule from the three genuine raw A/B/P stage-bound
families and proves that the corresponding mixed physical residual has
``VanishingJointJets`` at ``(t,x)=(1,0)``.

The contract is fail-closed: it accepts only one source/revision, one A/B/P
family, universal ``J,m : Nat`` rate hypotheses, exact rational scalar metadata,
and a ``lean-formal-export`` of the pinned theorem application.  Finite stage
prefixes, finite derivative ladders, sampled/numerical residuals, or schedules
selected independently for the three components are rejected.

Passing this bridge is still only ``formal-structure`` here.  It does not replay
Lean, materialize the infinite Section 9 fields, construct AwayExtensions,
prove the official late-plateau analytic majorants, materialize residual/force
artifacts, or close energy/blow-up.  Runtime reconstruction truth therefore
remains false.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


FORMAL_REPOSITORY = "openai/NavierStokesAndEuler"
FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
FORMAL_SOURCE = "NavierStokes/MixedDiagonalResidual.lean"
THEOREM = "NavierStokes.MixedDiagonalResidual.exists_physical_schedule_residual_zero"
ENDPOINT = Fraction(1, 1)
HALF = Fraction(1, 2)
_ACCEPTED_EVIDENCE = "lean-formal-export"
_STAGE_RATE_QUANTIFIER = "forall-natural-J-m"
_STAGE_FAMILY_QUANTIFIER = "forall-natural-stage"


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a natural number")
    return value


def _exact_positive_fraction(value: object, name: str) -> Fraction:
    if type(value) is not Fraction:
        raise ValueError(f"{name} must be an exact Fraction; floats are rejected")
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


@dataclass(frozen=True)
class Section9PhysicalScheduleResidualZeroWitness:
    """One pinned formal export of the physical schedule/residual-zero theorem."""

    source_id: str
    source_revision: str
    potential_family_id: str
    direct_family_id: str
    pressure_family_id: str
    raw_potential_bounds_id: str
    raw_direct_bounds_id: str
    raw_pressure_bounds_id: str
    background_rate_family_id: str
    uncut_residual_rate_family_id: str
    selected_schedule_id: str
    physical_q_id: str
    support_region_id: str
    mixed_residual_id: str
    provenance: str

    h: Fraction
    qbig: Fraction
    lower_stage: int
    endpoint: Fraction = ENDPOINT
    evidence_kind: str = _ACCEPTED_EVIDENCE

    formal_repository: str = FORMAL_REPOSITORY
    formal_commit: str = FORMAL_COMMIT
    formal_source: str = FORMAL_SOURCE
    theorem: str = THEOREM

    stage_family_quantifier: str = _STAGE_FAMILY_QUANTIFIER
    stage_rate_quantifier: str = _STAGE_RATE_QUANTIFIER
    max_stage_index: int | None = None
    max_derivative_order: int | None = None
    finite: bool = False

    support_open_certified: bool = False
    support_preterminal_certified: bool = False
    support_eventually_near_endpoint_origin_certified: bool = False
    all_stage_fields_smooth_on_physical_sublevel_certified: bool = False
    raw_potential_stage_bounds_certified: bool = False
    raw_direct_stage_bounds_certified: bool = False
    raw_pressure_stage_bounds_certified: bool = False
    gauge_nonnegative_at_zero_certified: bool = False
    gauge_positive_after_zero_certified: bool = False
    gauge_monotone_certified: bool = False
    gauge_tends_to_infinity_certified: bool = False
    background_jet_rate_all_J_m_certified: bool = False
    uncut_residual_jet_rate_all_J_m_certified: bool = False
    theorem_applied: bool = False

    selected_schedule_positive_certified: bool = False
    selected_schedule_doubling_certified: bool = False
    selected_schedule_strict_mono_certified: bool = False
    selected_schedule_tends_to_infinity_certified: bool = False
    selected_schedule_reciprocal_sublevel_certified: bool = False
    three_cut_bounds_certified: bool = False
    three_smooth_sums_certified: bool = False
    vanishing_joint_jets_certified: bool = False
    same_actual_field_source_certified: bool = False

    def __post_init__(self) -> None:
        for name in (
            "source_id",
            "source_revision",
            "potential_family_id",
            "direct_family_id",
            "pressure_family_id",
            "raw_potential_bounds_id",
            "raw_direct_bounds_id",
            "raw_pressure_bounds_id",
            "background_rate_family_id",
            "uncut_residual_rate_family_id",
            "selected_schedule_id",
            "physical_q_id",
            "support_region_id",
            "mixed_residual_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        h = _exact_positive_fraction(self.h, "h")
        if not h < HALF:
            raise ValueError("h must satisfy the pinned theorem range 0 < h < 1/2")
        _exact_positive_fraction(self.qbig, "qbig")
        _natural(self.lower_stage, "lower_stage")
        if type(self.endpoint) is not Fraction or self.endpoint != ENDPOINT:
            raise ValueError("endpoint must be the exact rational T=1")
        if self.evidence_kind != _ACCEPTED_EVIDENCE:
            raise ValueError("evidence_kind must be lean-formal-export")

        pinned_text = {
            "formal_repository": FORMAL_REPOSITORY,
            "formal_commit": FORMAL_COMMIT,
            "formal_source": FORMAL_SOURCE,
            "theorem": THEOREM,
        }
        for name, expected in pinned_text.items():
            if getattr(self, name) != expected:
                raise ValueError(f"{name} must equal pinned formal source value {expected!r}")

        if self.stage_family_quantifier != _STAGE_FAMILY_QUANTIFIER:
            raise ValueError("stage_family_quantifier must cover every natural stage")
        if self.stage_rate_quantifier != _STAGE_RATE_QUANTIFIER:
            raise ValueError("stage_rate_quantifier must universally quantify all natural J,m")
        if self.max_stage_index is not None:
            raise ValueError("max_stage_index must be None; finite prefixes are rejected")
        if self.max_derivative_order is not None:
            raise ValueError("max_derivative_order must be None; finite ladders are rejected")
        if self.finite is not False:
            raise ValueError("finite must be False for the physical all-order theorem bridge")

        for name in (
            "support_open_certified",
            "support_preterminal_certified",
            "support_eventually_near_endpoint_origin_certified",
            "all_stage_fields_smooth_on_physical_sublevel_certified",
            "raw_potential_stage_bounds_certified",
            "raw_direct_stage_bounds_certified",
            "raw_pressure_stage_bounds_certified",
            "gauge_nonnegative_at_zero_certified",
            "gauge_positive_after_zero_certified",
            "gauge_monotone_certified",
            "gauge_tends_to_infinity_certified",
            "background_jet_rate_all_J_m_certified",
            "uncut_residual_jet_rate_all_J_m_certified",
            "theorem_applied",
            "selected_schedule_positive_certified",
            "selected_schedule_doubling_certified",
            "selected_schedule_strict_mono_certified",
            "selected_schedule_tends_to_infinity_certified",
            "selected_schedule_reciprocal_sublevel_certified",
            "three_cut_bounds_certified",
            "three_smooth_sums_certified",
            "vanishing_joint_jets_certified",
            "same_actual_field_source_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def field_ids(self) -> tuple[str, str, str]:
        return self.potential_family_id, self.direct_family_id, self.pressure_family_id

    @property
    def raw_bound_ids(self) -> tuple[str, str, str]:
        return self.raw_potential_bounds_id, self.raw_direct_bounds_id, self.raw_pressure_bounds_id


@dataclass(frozen=True)
class Section9PhysicalScheduleResidualZeroCertificate:
    witness: Section9PhysicalScheduleResidualZeroWitness
    expected_source_key: tuple[str, str]
    expected_field_ids: tuple[str, str, str]
    expected_raw_bound_ids: tuple[str, str, str]
    expected_mixed_residual_id: str
    status: str = "formal-structure"

    common_physical_schedule_formally_selected: bool = True
    actual_mixed_residual_vanishing_joint_jets_formally_bound: bool = True

    actual_section9_sequence_verified: bool = False
    section9_all_order_endpoint_limits_verified: bool = False
    residual_artifact_ready: bool = False
    forcing_artifact_ready: bool = False
    endpoint_residual_closure_verified: bool = False
    paper_exact_velocity_available: bool = False
    full_reconstruction: bool = False

    @property
    def formal_physical_schedule_residual_zero_ready(self) -> bool:
        return (
            self.witness.source_key == self.expected_source_key
            and self.witness.field_ids == self.expected_field_ids
            and self.witness.raw_bound_ids == self.expected_raw_bound_ids
            and self.witness.mixed_residual_id == self.expected_mixed_residual_id
            and self.common_physical_schedule_formally_selected
            and self.actual_mixed_residual_vanishing_joint_jets_formally_bound
            and not self.actual_section9_sequence_verified
            and not self.section9_all_order_endpoint_limits_verified
            and not self.residual_artifact_ready
            and not self.forcing_artifact_ready
            and not self.endpoint_residual_closure_verified
            and not self.paper_exact_velocity_available
            and not self.full_reconstruction
        )


def admit_section9_physical_schedule_residual_zero(
    witness: Section9PhysicalScheduleResidualZeroWitness,
    *,
    expected_source_id: str,
    expected_source_revision: str,
    expected_potential_family_id: str,
    expected_direct_family_id: str,
    expected_pressure_family_id: str,
    expected_raw_potential_bounds_id: str,
    expected_raw_direct_bounds_id: str,
    expected_raw_pressure_bounds_id: str,
    expected_mixed_residual_id: str,
) -> Section9PhysicalScheduleResidualZeroCertificate:
    """Bind one actual A/B/P source to the theorem-selected physical schedule."""
    if not isinstance(witness, Section9PhysicalScheduleResidualZeroWitness):
        raise TypeError("witness must be a Section9PhysicalScheduleResidualZeroWitness")

    expected_source_key = (
        _nonempty_text(expected_source_id, "expected_source_id"),
        _nonempty_text(expected_source_revision, "expected_source_revision"),
    )
    expected_field_ids = (
        _nonempty_text(expected_potential_family_id, "expected_potential_family_id"),
        _nonempty_text(expected_direct_family_id, "expected_direct_family_id"),
        _nonempty_text(expected_pressure_family_id, "expected_pressure_family_id"),
    )
    expected_raw_bound_ids = (
        _nonempty_text(expected_raw_potential_bounds_id, "expected_raw_potential_bounds_id"),
        _nonempty_text(expected_raw_direct_bounds_id, "expected_raw_direct_bounds_id"),
        _nonempty_text(expected_raw_pressure_bounds_id, "expected_raw_pressure_bounds_id"),
    )
    expected_mixed = _nonempty_text(expected_mixed_residual_id, "expected_mixed_residual_id")

    if witness.source_key != expected_source_key:
        raise ValueError("physical schedule theorem export source identity mismatch")
    if witness.field_ids != expected_field_ids:
        raise ValueError("physical schedule theorem export field identity mismatch")
    if witness.raw_bound_ids != expected_raw_bound_ids:
        raise ValueError("physical schedule theorem export raw-bound identity mismatch")
    if witness.mixed_residual_id != expected_mixed:
        raise ValueError("physical schedule theorem export mixed-residual identity mismatch")

    certificate = Section9PhysicalScheduleResidualZeroCertificate(
        witness=witness,
        expected_source_key=expected_source_key,
        expected_field_ids=expected_field_ids,
        expected_raw_bound_ids=expected_raw_bound_ids,
        expected_mixed_residual_id=expected_mixed,
    )
    if not certificate.formal_physical_schedule_residual_zero_ready:
        raise ArithmeticError("physical schedule residual-zero bridge invariant failed")
    return certificate
