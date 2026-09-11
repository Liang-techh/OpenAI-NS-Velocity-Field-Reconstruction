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


def _hypothesis(*, first_scale=8):
    return Section9LocalSumHypothesis(
        q_big=Fraction(1, 4),
        first_scale=first_scale,
        kind="paper-derived",
        provenance="Proposition 9.9 + Lemma 5.4 infinite schedule",
    )


def _stage_certificate(stage, scale, *, hypothesis=None):
    hypothesis = _hypothesis() if hypothesis is None else hypothesis
    rows = [
        Section9CorrectionExtensionWitness(
            stage=stage,
            component=component,
            q_big=hypothesis.q_big,
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
    return admit_section9_correction_stage(hypothesis, rows)


def _value(stage):
    return Section9CorrectionStageValue(
        stage=stage,
        A=np.array([stage, stage + 1, -stage], dtype=float),
        B=10.0 * stage,
        p=-float(stage),
        provenance=f"manufactured regression payload stage {stage}; not paper data",
    )


def _independent_linear_cutoff(s):
    # Test-only pointwise evaluator. It has the same plateau/support endpoints
    # used by the paper but is deliberately not claimed to be the manuscript chi.
    if s <= 0.5:
        return 1.0
    if s >= 1.0:
        raise AssertionError("support-zero stages must be short-circuited before evaluation")
    return 2.0 * (1.0 - s)


def test_evaluates_contiguous_admitted_prefix_and_short_circuits_support_zero():
    hypothesis = _hypothesis()
    certificates = [
        _stage_certificate(1, 8),
        _stage_certificate(2, 16),
        _stage_certificate(3, 40),
        _stage_certificate(4, 80),
    ]
    result = evaluate_admitted_eq_9_21_prefix(
        hypothesis,
        certificates,
        [_value(1), _value(2), _value(3), _value(4)],
        q=Fraction(1, 50),
        cutoff_evaluator=_independent_linear_cutoff,
        cutoff_evaluator_provenance="test-only linear cutoff; not paper-exact",
    )

    # Independent arithmetic oracle at q=1/50:
    # a_j q = 4/25, 8/25, 4/5, 8/5.  Hence weights are
    # 1, 1, 2*(1-4/5)=2/5, 0 (the last from the support theorem).
    assert result.stages == (1, 2, 3, 4)
    assert [c.cutoff_argument for c in result.contributions] == [
        Fraction(4, 25),
        Fraction(8, 25),
        Fraction(4, 5),
        Fraction(8, 5),
    ]
    assert np.allclose([c.cutoff_weight for c in result.contributions], [1.0, 1.0, 0.4, 0.0])
    assert np.allclose(result.A_correction_prefix, [4.2, 6.6, -4.2])
    assert result.B_correction_prefix == pytest.approx(42.0)
    assert result.p_correction_prefix == pytest.approx(-4.2)
    assert [c.support_zero_short_circuit for c in result.contributions] == [
        False,
        False,
        False,
        True,
    ]
    assert result.formal_prefix_ready

    # A finite evaluation is intentionally not the paper's infinite sum.
    assert not result.actual_correction_field_values_verified
    assert not result.paper_fixed_cutoff_pointwise_evaluator_verified
    assert not result.infinite_schedule_constructed_here
    assert not result.eq_9_21_infinite_sum_constructed
    assert not result.proposition_9_9_verified
    assert not result.paper_exact_velocity_available


def test_rejects_sparse_or_duplicate_stage_certificates_and_inexact_value_coverage():
    hypothesis = _hypothesis()
    c1 = _stage_certificate(1, 8)
    c2 = _stage_certificate(2, 16)
    c3 = _stage_certificate(3, 40)

    with pytest.raises(ValueError, match="contiguous positive-stage prefix"):
        evaluate_admitted_eq_9_21_prefix(
            hypothesis,
            [c1, c3],
            [_value(1), _value(3)],
            q=Fraction(1, 50),
            cutoff_evaluator=_independent_linear_cutoff,
            cutoff_evaluator_provenance="test only",
        )

    with pytest.raises(ValueError, match="duplicate admitted Section 9 stage"):
        evaluate_admitted_eq_9_21_prefix(
            hypothesis,
            [c1, c1],
            [_value(1)],
            q=Fraction(1, 50),
            cutoff_evaluator=_independent_linear_cutoff,
            cutoff_evaluator_provenance="test only",
        )

    with pytest.raises(ValueError, match="cover exactly the admitted prefix"):
        evaluate_admitted_eq_9_21_prefix(
            hypothesis,
            [c1, c2, c3],
            [_value(1), _value(2)],
            q=Fraction(1, 50),
            cutoff_evaluator=_independent_linear_cutoff,
            cutoff_evaluator_provenance="test only",
        )

    with pytest.raises(ValueError, match="duplicate Section 9 stage value"):
        evaluate_admitted_eq_9_21_prefix(
            hypothesis,
            [c1, c2],
            [_value(1), _value(1), _value(2)],
            q=Fraction(1, 50),
            cutoff_evaluator=_independent_linear_cutoff,
            cutoff_evaluator_provenance="test only",
        )


def test_rejects_mixed_schedule_or_q_domain_and_bad_cutoff_values():
    original_hypothesis = _hypothesis()
    c1 = _stage_certificate(1, 8, hypothesis=original_hypothesis)

    # Same q_big but a different first scale: do not mix a certificate from
    # another admitted infinite schedule into this prefix.
    different_schedule = _hypothesis(first_scale=6)
    with pytest.raises(ValueError, match="different schedule first_scale"):
        evaluate_admitted_eq_9_21_prefix(
            different_schedule,
            [c1],
            [_value(1)],
            q=Fraction(1, 50),
            cutoff_evaluator=_independent_linear_cutoff,
            cutoff_evaluator_provenance="test only",
        )

    with pytest.raises(ValueError, match="0 < q < q_big"):
        evaluate_admitted_eq_9_21_prefix(
            original_hypothesis,
            [c1],
            [_value(1)],
            q=Fraction(1, 4),
            cutoff_evaluator=_independent_linear_cutoff,
            cutoff_evaluator_provenance="test only",
        )

    with pytest.raises(ValueError, match=r"value in \[0,1\]"):
        evaluate_admitted_eq_9_21_prefix(
            original_hypothesis,
            [c1],
            [_value(1)],
            q=Fraction(1, 100),
            cutoff_evaluator=lambda _: 1.5,
            cutoff_evaluator_provenance="invalid test evaluator",
        )

    with pytest.raises(ValueError, match="must be finite"):
        evaluate_admitted_eq_9_21_prefix(
            original_hypothesis,
            [c1],
            [_value(1)],
            q=Fraction(1, 100),
            cutoff_evaluator=lambda _: float("nan"),
            cutoff_evaluator_provenance="invalid test evaluator",
        )


def test_stage_value_payload_fails_closed_on_shape_finiteness_and_provenance():
    with pytest.raises(ValueError, match="finite three-vector"):
        Section9CorrectionStageValue(
            stage=1, A=[1.0, 2.0], B=0.0, p=0.0, provenance="test"
        )
    with pytest.raises(ValueError, match="B must be finite"):
        Section9CorrectionStageValue(
            stage=1, A=[1.0, 2.0, 3.0], B=float("inf"), p=0.0, provenance="test"
        )
    with pytest.raises(ValueError, match="value provenance must be nonempty"):
        Section9CorrectionStageValue(
            stage=1, A=[1.0, 2.0, 3.0], B=0.0, p=0.0, provenance="   "
        )
    with pytest.raises(ValueError, match="stage must be a positive integer"):
        Section9CorrectionStageValue(
            stage=0, A=[1.0, 2.0, 3.0], B=0.0, p=0.0, provenance="test"
        )
