from dataclasses import replace
from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.coordinates import solve_q
from openai_ns_reconstruction.section9_correction_extension_admission import (
    Section9CorrectionExtensionWitness,
    admit_section9_correction_stage,
)
from openai_ns_reconstruction.section9_local_sum_admission import Section9LocalSumHypothesis
from openai_ns_reconstruction.section9_stage_physical_support import (
    certify_admitted_stage_physical_support,
)


def _hypothesis(q_big=Fraction(1, 4)):
    return Section9LocalSumHypothesis(
        q_big=q_big,
        first_scale=8,
        kind="paper-derived",
        provenance="Lemma 9.7 + Lemma 5.4 infinite schedule",
    )


def _stage_certificate(*, stage=1, scale=12, q_big=Fraction(1, 4)):
    hypothesis = _hypothesis(q_big)
    witnesses = [
        Section9CorrectionExtensionWitness(
            stage=stage,
            component=component,
            q_big=q_big,
            scale=scale,
            kind="paper-derived",
            provenance=f"Eq. (9.21) stage {stage} {component} extension theorem",
            common_domain_verified=True,
            schedule_membership_verified=True,
            cutoff_support_verified=True,
            smooth_zero_extension_verified=True,
        )
        for component in ("A", "B", "p")
    ]
    return admit_section9_correction_stage(hypothesis, witnesses)


def test_actual_admitted_scale_can_sharpen_generic_doubling_support_bound():
    hypothesis = _hypothesis()
    stage = _stage_certificate(stage=1, scale=12)
    cert = certify_admitted_stage_physical_support(
        hypothesis,
        stage,
        t_lower=Fraction(7, 8),
        t_upper=Fraction(11, 12),
    )

    # Eq. (4.1) physical-slab bridge gives q >= 1-t_upper = 1/12.
    # The witnessed actual a_1=12 therefore gives a_1*q >= 1 exactly, so the
    # paper support rule chi(s)=0 for s>=1 kills this whole stage.  The generic
    # schedule knows only a_1>=8, for which 8*(1/12)<1, so it cannot exclude it.
    assert cert.q_lower_bound == Fraction(1, 12)
    assert cert.q_upper_bound == Fraction(3, 16)
    assert cert.cutoff_argument_lower == 1
    assert cert.cutoff_argument_upper == Fraction(9, 4)
    assert cert.uniformly_zero_on_slab
    assert cert.zero_short_circuit_certified
    assert not cert.generic_doubling_bound_excludes_stage
    assert cert.actual_scale_sharpens_generic_bound

    assert cert.paper_cutoff_support_rule_verified
    assert cert.smooth_zero_extension_verified
    assert not cert.actual_correction_field_values_verified
    assert not cert.stage_activity_verified
    assert not cert.endpoint_covered
    assert not cert.eq_9_21_sum_constructed
    assert not cert.endpoint_uniform_residual_majorants_verified
    assert not cert.paper_exact_velocity_available


def test_failure_to_prove_uniform_zero_is_not_promoted_to_stage_activity():
    hypothesis = _hypothesis()
    stage = _stage_certificate(stage=1, scale=8)
    cert = certify_admitted_stage_physical_support(
        hypothesis,
        stage,
        t_lower=Fraction(7, 8),
        t_upper=Fraction(11, 12),
    )

    assert cert.cutoff_argument_lower == Fraction(2, 3)
    assert not cert.uniformly_zero_on_slab
    assert not cert.zero_short_circuit_certified
    assert not cert.generic_doubling_bound_excludes_stage
    assert not cert.actual_scale_sharpens_generic_bound
    assert not cert.stage_activity_verified


def test_generic_local_finiteness_exclusion_agrees_with_actual_scale_support_test():
    hypothesis = _hypothesis()
    stage = _stage_certificate(stage=3, scale=40)
    cert = certify_admitted_stage_physical_support(
        hypothesis,
        stage,
        t_lower=Fraction(7, 8),
        t_upper=Fraction(11, 12),
    )

    assert cert.physical_slab.first_guaranteed_inactive_stage == 2
    assert cert.stage == 3
    assert cert.generic_doubling_bound_excludes_stage
    assert cert.uniformly_zero_on_slab
    assert not cert.actual_scale_sharpens_generic_bound


def test_independent_eq_4_1_solver_confirms_certified_uniform_zero():
    hypothesis = _hypothesis()
    stage = _stage_certificate(stage=1, scale=12)
    cert = certify_admitted_stage_physical_support(
        hypothesis,
        stage,
        t_lower=Fraction(7, 8),
        t_upper=Fraction(11, 12),
    )
    assert cert.uniformly_zero_on_slab

    # Independent regression: solve Eq. (4.1) at sampled physical points rather
    # than reusing the rational support algebra from the production certificate.
    for h in (0.005, 0.1, 0.2):
        for t in np.linspace(7.0 / 8.0, 11.0 / 12.0, 7):
            for z in np.linspace(-0.25, 0.25, 9):
                q = solve_q(float(z), float(t), h)
                assert stage.scale * q >= 1.0 - 2e-12


def test_rejects_stage_certificate_from_a_different_common_domain_or_schedule():
    hypothesis = _hypothesis()

    foreign_domain = _stage_certificate(stage=1, scale=12, q_big=Fraction(1, 3))
    with pytest.raises(ValueError, match="one common q_big"):
        certify_admitted_stage_physical_support(
            hypothesis,
            foreign_domain,
            t_lower=Fraction(7, 8),
            t_upper=Fraction(11, 12),
        )

    stage = _stage_certificate(stage=1, scale=12)
    corrupted_schedule = replace(stage, schedule_lower_bound=7)
    assert corrupted_schedule.formal_admission_ready
    with pytest.raises(ValueError, match="supplied schedule first_scale"):
        certify_admitted_stage_physical_support(
            hypothesis,
            corrupted_schedule,
            t_lower=Fraction(7, 8),
            t_upper=Fraction(11, 12),
        )


def test_endpoint_remains_fail_closed():
    with pytest.raises(ValueError, match="strictly before t=1"):
        certify_admitted_stage_physical_support(
            _hypothesis(),
            _stage_certificate(stage=1, scale=12),
            t_lower=Fraction(7, 8),
            t_upper=1,
        )
