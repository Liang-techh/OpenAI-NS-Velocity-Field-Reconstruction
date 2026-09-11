"""Exact physical-slab support test for one admitted Section 9 correction stage.

This module composes two fail-closed pieces already present on ``main``:

* :mod:`section9_correction_extension_admission`, which admits one positive
  Eq. (9.21) stage only after theorem/proof witnesses certify coherent
  ``A_j/B_j/p_j`` data, its actual schedule scale ``a_j``, the paper cutoff
  support rule, and smooth zero extension; and
* :mod:`section9_physical_window`, which maps the fixed Section 10 support on a
  closed pre-endpoint slab into an exact rational q-strip using Eq. (4.1).

For an admitted stage, Eq. (9.21) multiplies every correction component by
``chi(a_j q)``.  The paper cutoff satisfies ``chi(s)=0`` for ``s>=1``.  Hence
on a physical slab with certified lower q-bound ``q_min``, the entire stage is
identically zero throughout the fixed Section 10 support whenever

    a_j * q_min >= 1.

This is a rigorous support consequence that needs no correction samples and can
be sharper than the generic doubling-schedule bound because it uses the
witnessed *actual* ``a_j``.  Failure of the inequality is deliberately not
interpreted as proof that the stage is active.

The bridge remains pre-endpoint only.  It does not construct correction values,
the Eq. (9.21) sum, residual majorants, or a smooth extension through ``t=1``.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .section9_correction_extension_admission import (
    Section9CorrectionStageAdmissionCertificate,
)
from .section9_local_sum_admission import Section9LocalSumHypothesis
from .section9_physical_window import (
    Section9PhysicalSlabLocalSumCertificate,
    certify_section10_physical_slab_local_finiteness,
)


@dataclass(frozen=True)
class Section9AdmittedStagePhysicalSupportCertificate:
    """Support-only consequence for one admitted stage on one physical slab."""

    stage: int
    scale: int
    t_lower: Fraction
    t_upper: Fraction
    q_lower_bound: Fraction
    q_upper_bound: Fraction
    cutoff_argument_lower: Fraction
    cutoff_argument_upper: Fraction
    physical_slab: Section9PhysicalSlabLocalSumCertificate
    uniformly_zero_on_slab: bool
    generic_doubling_bound_excludes_stage: bool
    actual_scale_sharpens_generic_bound: bool
    paper_cutoff_support_rule_verified: bool = True
    smooth_zero_extension_verified: bool = True
    actual_correction_field_values_verified: bool = False
    stage_activity_verified: bool = False
    endpoint_covered: bool = False
    eq_9_21_sum_constructed: bool = False
    endpoint_uniform_residual_majorants_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def zero_short_circuit_certified(self) -> bool:
        """Whether the whole admitted stage can be omitted on this slab."""
        return (
            self.uniformly_zero_on_slab
            and self.paper_cutoff_support_rule_verified
            and self.smooth_zero_extension_verified
        )


def certify_admitted_stage_physical_support(
    hypothesis: Section9LocalSumHypothesis,
    stage_certificate: Section9CorrectionStageAdmissionCertificate,
    *,
    t_lower: Fraction | float | int,
    t_upper: Fraction | float | int,
) -> Section9AdmittedStagePhysicalSupportCertificate:
    """Certify whether an admitted Eq. (9.21) stage vanishes on a whole slab.

    The supplied stage certificate must be coherent with the same global
    ``Section9LocalSumHypothesis`` used for the physical slab.  In particular,
    this checks the common ``q_big`` and recomputes the exact schedule lower
    bound from ``first_scale`` and the admitted stage number.  This prevents a
    certificate produced from a different schedule from being silently mixed
    into the Section 10 bridge.

    ``uniformly_zero_on_slab`` is proved exactly from ``a_j*q_min >= 1`` and
    the witnessed paper cutoff support rule.  When the inequality fails, the
    result says only that this support argument cannot exclude the stage; it
    never upgrades that failure to a positive activity statement.
    """
    if not isinstance(hypothesis, Section9LocalSumHypothesis):
        raise TypeError("hypothesis must be a Section9LocalSumHypothesis")
    if not isinstance(stage_certificate, Section9CorrectionStageAdmissionCertificate):
        raise TypeError(
            "stage_certificate must be a Section9CorrectionStageAdmissionCertificate"
        )
    if not stage_certificate.formal_admission_ready:
        raise ValueError("stage certificate is not a complete fail-closed formal admission")
    if stage_certificate.q_big != hypothesis.q_big:
        raise ValueError("stage certificate and physical slab must use one common q_big")

    expected_lower_bound = hypothesis.first_scale * (1 << (stage_certificate.stage - 1))
    if stage_certificate.schedule_lower_bound != expected_lower_bound:
        raise ValueError(
            "stage certificate was not derived from the supplied schedule first_scale"
        )
    if stage_certificate.scale < expected_lower_bound:
        raise ArithmeticError("admitted stage scale violates the supplied schedule bound")

    physical = certify_section10_physical_slab_local_finiteness(
        hypothesis,
        t_lower=t_lower,
        t_upper=t_upper,
    )
    q_lower = physical.q_lower_bound
    q_upper = physical.q_upper_bound
    scale = stage_certificate.scale
    arg_lower = scale * q_lower
    arg_upper = scale * q_upper

    uniformly_zero = arg_lower >= 1
    generic_excludes = (
        stage_certificate.stage >= physical.first_guaranteed_inactive_stage
    )

    # The generic local-finiteness certificate uses only
    # a_j >= a_1*2^(j-1).  If it excludes a stage, the actual witnessed scale
    # must also make the support argument vanish.  Keep this as a cross-module
    # invariant so schedule/strict-inequality drift fails loudly.
    if generic_excludes and not uniformly_zero:
        raise ArithmeticError(
            "generic doubling bound excludes the stage but the actual-scale support test does not"
        )

    sharpens = uniformly_zero and not generic_excludes
    return Section9AdmittedStagePhysicalSupportCertificate(
        stage=stage_certificate.stage,
        scale=scale,
        t_lower=physical.t_lower,
        t_upper=physical.t_upper,
        q_lower_bound=q_lower,
        q_upper_bound=q_upper,
        cutoff_argument_lower=arg_lower,
        cutoff_argument_upper=arg_upper,
        physical_slab=physical,
        uniformly_zero_on_slab=uniformly_zero,
        generic_doubling_bound_excludes_stage=generic_excludes,
        actual_scale_sharpens_generic_bound=sharpens,
        paper_cutoff_support_rule_verified=stage_certificate.cutoff_support_verified,
        smooth_zero_extension_verified=stage_certificate.smooth_zero_extension_verified,
    )
