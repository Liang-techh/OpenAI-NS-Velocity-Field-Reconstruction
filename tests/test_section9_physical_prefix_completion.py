from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.section9_correction_extension_admission import (
    Section9CorrectionExtensionWitness,
    admit_section9_correction_stage,
)
from openai_ns_reconstruction.section9_finite_prefix import (
    Section9CorrectionStageValue,
    evaluate_admitted_eq_9_21_prefix,
)
from openai_ns_reconstruction.section9_local_sum_admission import Section9LocalSumHypothesis
from openai_ns_reconstruction.section9_physical_prefix_completion import (
    certify_section9_physical_prefix_completion,
)
from openai_ns_reconstruction.section9_physical_window import (
    certify_section10_physical_slab_local_finiteness,
)


def _hypothesis():
    return Section9LocalSumHypothesis(
        q_big=Fraction(1, 2),
        first_scale=8,
        kind="paper-derived",
        provenance="Lemma 9.7 common domain + Lemma 5.4 infinite schedule",
    )


def _stage_certificate(hypothesis, stage, scale):
    return admit_section9_correction_stage(
        hypothesis,
        [
            Section9CorrectionExtensionWitness(
                stage=stage,
                component=component,
                q_big=hypothesis.q_big,
                scale=scale,
                kind="paper-derived",
                provenance=f"stage {stage} {component} extension theorem",
                common_domain_verified=True,
                schedule_membership_verified=True,
                cutoff_support_verified=True,
                smooth_zero_extension_verified=True,
            )
            for component in ("A", "B", "p")
        ],
    )


def _value(stage):
    return Section9CorrectionStageValue(
        stage=stage,
        A=np.array([stage, -stage, stage + 0.5], dtype=float),
        B=float(stage),
        p=-float(stage),
        provenance="manufactured regression payload; not paper data",
    )


def _test_cutoff(s):
    if s <= 0.5:
        return 1.0
    if s >= 1.0:
        raise AssertionError("support-zero contribution should short-circuit")
    return 2.0 * (1.0 - s)


def test_prefix_reaching_physical_active_bound_has_uniform_zero_omitted_tail():
    hypothesis = _hypothesis()
    slab = certify_section10_physical_slab_local_finiteness(
        hypothesis,
        t_lower=Fraction(3, 4),
        t_upper=Fraction(15, 16),
    )
    assert slab.q_lower_bound == Fraction(1, 16)
    assert slab.max_potentially_active_stage == 1

    c1 = _stage_certificate(hypothesis, 1, 8)
    evaluation = evaluate_admitted_eq_9_21_prefix(
        hypothesis,
        [c1],
        [_value(1)],
        q=Fraction(1, 10),
        cutoff_evaluator=_test_cutoff,
        cutoff_evaluator_provenance="test-only cutoff; not manuscript pointwise chi",
    )
    cert = certify_section9_physical_prefix_completion(
        hypothesis,
        slab,
        evaluation,
        [c1],
    )

    assert cert.prefix_order == 1
    assert cert.max_potentially_active_stage == 1
    assert cert.first_guaranteed_inactive_stage == 2
    assert cert.first_omitted_stage == 2
    assert cert.omitted_tail_support_zero_on_whole_slab
    assert cert.formal_slab_prefix_ready

    # Independent arithmetic oracle: every omitted j>=2 has
    # a_j*q >= 8*2^(j-1)*(1/16) >= 1 on the whole slab.
    for stage in range(cert.first_omitted_stage, cert.first_omitted_stage + 6):
        schedule_lower_bound = hypothesis.first_scale * (1 << (stage - 1))
        assert schedule_lower_bound * slab.q_lower_bound >= 1

    assert not cert.endpoint_covered
    assert not cert.actual_correction_field_values_verified
    assert not cert.paper_fixed_cutoff_pointwise_evaluator_verified
    assert not cert.eq_9_21_infinite_sum_constructed
    assert not cert.endpoint_uniform_residual_majorants_verified
    assert not cert.proposition_9_9_verified
    assert not cert.paper_exact_velocity_available


def test_rejects_prefix_too_short_for_a_slab_closer_to_endpoint():
    hypothesis = _hypothesis()
    slab = certify_section10_physical_slab_local_finiteness(
        hypothesis,
        t_lower=Fraction(7, 8),
        t_upper=Fraction(127, 128),
    )
    assert slab.q_lower_bound == Fraction(1, 128)
    assert slab.max_potentially_active_stage == 4

    c1 = _stage_certificate(hypothesis, 1, 8)
    evaluation = evaluate_admitted_eq_9_21_prefix(
        hypothesis,
        [c1],
        [_value(1)],
        q=Fraction(1, 64),
        cutoff_evaluator=_test_cutoff,
        cutoff_evaluator_provenance="test only",
    )
    with pytest.raises(ValueError, match="too short for whole-slab support completeness"):
        certify_section9_physical_prefix_completion(hypothesis, slab, evaluation, [c1])


def test_rejects_pointwise_evaluation_outside_physical_q_enclosure():
    hypothesis = _hypothesis()
    slab = certify_section10_physical_slab_local_finiteness(
        hypothesis,
        t_lower=Fraction(3, 4),
        t_upper=Fraction(15, 16),
    )
    c1 = _stage_certificate(hypothesis, 1, 8)
    evaluation = evaluate_admitted_eq_9_21_prefix(
        hypothesis,
        [c1],
        [_value(1)],
        q=Fraction(1, 100),
        cutoff_evaluator=_test_cutoff,
        cutoff_evaluator_provenance="test only",
    )
    with pytest.raises(ValueError, match="outside the physical slab q enclosure"):
        certify_section9_physical_prefix_completion(hypothesis, slab, evaluation, [c1])


def test_rejects_admission_payload_that_does_not_match_evaluated_scale():
    hypothesis = _hypothesis()
    slab = certify_section10_physical_slab_local_finiteness(
        hypothesis,
        t_lower=Fraction(3, 4),
        t_upper=Fraction(15, 16),
    )
    c1_eval = _stage_certificate(hypothesis, 1, 8)
    c1_other = _stage_certificate(hypothesis, 1, 12)
    evaluation = evaluate_admitted_eq_9_21_prefix(
        hypothesis,
        [c1_eval],
        [_value(1)],
        q=Fraction(1, 10),
        cutoff_evaluator=_test_cutoff,
        cutoff_evaluator_provenance="test only",
    )
    with pytest.raises(ValueError, match="scale does not match"):
        certify_section9_physical_prefix_completion(
            hypothesis, slab, evaluation, [c1_other]
        )


def test_rejects_physical_slab_from_different_hypothesis_provenance():
    hypothesis = _hypothesis()
    other = Section9LocalSumHypothesis(
        q_big=Fraction(1, 2),
        first_scale=8,
        kind="paper-derived",
        provenance="different independent theorem ledger",
    )
    slab = certify_section10_physical_slab_local_finiteness(
        other,
        t_lower=Fraction(3, 4),
        t_upper=Fraction(15, 16),
    )
    c1 = _stage_certificate(hypothesis, 1, 8)
    evaluation = evaluate_admitted_eq_9_21_prefix(
        hypothesis,
        [c1],
        [_value(1)],
        q=Fraction(1, 10),
        cutoff_evaluator=_test_cutoff,
        cutoff_evaluator_provenance="test only",
    )
    with pytest.raises(ValueError, match="supplied local-sum hypothesis"):
        certify_section9_physical_prefix_completion(hypothesis, slab, evaluation, [c1])
