from dataclasses import replace
from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.section9_correction_extension_admission import (
    Section9CorrectionStageAdmissionCertificate,
)
from openai_ns_reconstruction.section9_eq921_prefix_jet import (
    assemble_eq_9_21_prefix_jet,
    required_spacetime_multiindices,
)
from openai_ns_reconstruction.section9_finite_prefix import (
    Section9FinitePrefixEvaluationCertificate,
    Section9FinitePrefixStageContribution,
)
from openai_ns_reconstruction.section9_weighted_correction_jet import (
    Section9CutoffWeightJet,
    Section9UnweightedCorrectionStageJet,
    derive_weighted_correction_stage_jet,
)


ZERO = (0, 0, 0, 0)
TX = (1, 1, 0, 0)
T = (1, 0, 0, 0)
X = (0, 1, 0, 0)


def _admission() -> Section9CorrectionStageAdmissionCertificate:
    result = Section9CorrectionStageAdmissionCertificate(
        stage=1,
        q_big=Fraction(1, 2),
        scale=4,
        schedule_lower_bound=4,
        components=("A", "B", "p"),
        evidence_kinds=("paper-derived", "paper-derived", "paper-derived"),
        evidence_provenance=("paper:A:1", "paper:B:1", "paper:p:1"),
    )
    assert result.formal_admission_ready
    return result


def _evaluation() -> Section9FinitePrefixEvaluationCertificate:
    contribution = Section9FinitePrefixStageContribution(
        stage=1,
        scale=4,
        cutoff_argument=Fraction(1, 4),
        cutoff_weight=0.5,
        A_term=np.array([1.0, 2.0, 3.0]),
        B_term=4.0,
        p_term=5.0,
        value_provenance="fixture:stage-1-values",
        support_zero_short_circuit=False,
    )
    result = Section9FinitePrefixEvaluationCertificate(
        q=Fraction(1, 16),
        q_big=Fraction(1, 2),
        prefix_order=1,
        stages=(1,),
        contributions=(contribution,),
        A_correction_prefix=np.array([1.0, 2.0, 3.0]),
        B_correction_prefix=4.0,
        p_correction_prefix=5.0,
        cutoff_evaluator_provenance="fixture:cutoff-jet",
    )
    assert result.formal_prefix_ready
    return result


def _blank_vector_jet(order: int):
    return {
        alpha: np.zeros(3, dtype=float)
        for alpha in required_spacetime_multiindices(order)
    }


def _blank_scalar_jet(order: int):
    return {alpha: 0.0 for alpha in required_spacetime_multiindices(order)}


def _inputs():
    order = 2
    admission = _admission()
    evaluation = _evaluation()

    A = _blank_vector_jet(order)
    B = _blank_scalar_jet(order)
    p = _blank_scalar_jet(order)
    A[ZERO] = np.array([2.0, 4.0, 6.0])
    A[T] = np.array([1.0, -2.0, 3.0])
    A[X] = np.array([4.0, 1.0, -1.0])
    A[TX] = np.array([2.0, 5.0, 7.0])
    B[ZERO], B[T], B[X], B[TX] = 8.0, 2.0, -1.0, 6.0
    p[ZERO], p[T], p[X], p[TX] = 10.0, -4.0, 3.0, 2.0

    correction = Section9UnweightedCorrectionStageJet(
        stage=1,
        scale=4,
        q=Fraction(1, 16),
        derivative_order=order,
        A_derivatives=A,
        B_derivatives=B,
        p_derivatives=p,
        evidence_kinds=admission.evidence_kinds,
        evidence_provenance=admission.evidence_provenance,
        value_provenance="fixture:stage-1-values",
        provider_id="fixture-correction-provider",
        provider_revision="r1",
        provider_provenance="fixture:unweighted-correction-jet",
    )

    w = _blank_scalar_jet(order)
    w[ZERO] = 0.5
    w[T] = 2.0
    w[X] = -3.0
    w[TX] = 5.0
    cutoff = Section9CutoffWeightJet(
        stage=1,
        scale=4,
        q=Fraction(1, 16),
        derivative_order=order,
        derivatives=w,
        provider_id="fixture-cutoff-provider",
        provider_revision="chi-r1",
        provider_provenance="fixture:cutoff-jet",
    )
    return evaluation, admission, correction, cutoff


def test_product_rule_derives_weighted_jet_and_feeds_prefix_aggregator():
    evaluation, admission, correction, cutoff = _inputs()
    weighted = derive_weighted_correction_stage_jet(
        evaluation,
        admission,
        correction,
        cutoff,
    )

    np.testing.assert_array_equal(weighted.A_derivatives[ZERO], evaluation.A_correction_prefix)
    assert weighted.B_derivatives[ZERO] == evaluation.B_correction_prefix
    assert weighted.p_derivatives[ZERO] == evaluation.p_correction_prefix

    # Independent literal mixed-derivative formula:
    # D_tx(w F) = w_tx F + w_t F_x + w_x F_t + w F_tx.
    expected_A_tx = (
        5.0 * correction.A_derivatives[ZERO]
        + 2.0 * correction.A_derivatives[X]
        - 3.0 * correction.A_derivatives[T]
        + 0.5 * correction.A_derivatives[TX]
    )
    np.testing.assert_array_equal(weighted.A_derivatives[TX], expected_A_tx)
    expected_B_tx = 5.0 * 8.0 + 2.0 * (-1.0) - 3.0 * 2.0 + 0.5 * 6.0
    expected_p_tx = 5.0 * 10.0 + 2.0 * 3.0 - 3.0 * (-4.0) + 0.5 * 2.0
    assert weighted.B_derivatives[TX] == expected_B_tx
    assert weighted.p_derivatives[TX] == expected_p_tx

    prefix = assemble_eq_9_21_prefix_jet(evaluation, (admission,), (weighted,))
    assert prefix.formal_prefix_jet_ready
    assert not prefix.analytic_derivatives_machine_derived_from_actual_corrections
    assert not prefix.paper_fixed_cutoff_derivatives_machine_verified
    assert not prefix.residual_artifact_ready
    assert not prefix.paper_exact_velocity_available


def test_product_rule_rejects_cutoff_zero_order_drift():
    evaluation, admission, correction, cutoff = _inputs()
    derivatives = dict(cutoff.derivatives)
    derivatives[ZERO] = 0.25
    bad = replace(cutoff, derivatives=derivatives)

    with pytest.raises(ValueError, match="zero-order value"):
        derive_weighted_correction_stage_jet(evaluation, admission, correction, bad)


def test_product_rule_rejects_value_provenance_drift():
    evaluation, admission, correction, cutoff = _inputs()
    bad = replace(correction, value_provenance="other-stage-values")

    with pytest.raises(ValueError, match="value provenance"):
        derive_weighted_correction_stage_jet(evaluation, admission, bad, cutoff)


def test_product_rule_rejects_cutoff_evaluator_provenance_drift():
    evaluation, admission, correction, cutoff = _inputs()
    bad = replace(cutoff, provider_provenance="different-cutoff-provider")

    with pytest.raises(ValueError, match="cutoff jet provenance"):
        derive_weighted_correction_stage_jet(evaluation, admission, correction, bad)


def test_unweighted_jet_requires_complete_multiindex_coverage():
    evaluation, admission, correction, _ = _inputs()
    A = dict(correction.A_derivatives)
    del A[TX]

    with pytest.raises(ValueError, match="cover exactly all derivatives"):
        Section9UnweightedCorrectionStageJet(
            stage=1,
            scale=4,
            q=evaluation.q,
            derivative_order=2,
            A_derivatives=A,
            B_derivatives=correction.B_derivatives,
            p_derivatives=correction.p_derivatives,
            evidence_kinds=admission.evidence_kinds,
            evidence_provenance=admission.evidence_provenance,
            value_provenance="fixture:stage-1-values",
            provider_id="fixture-correction-provider",
            provider_revision="r1",
            provider_provenance="fixture:unweighted-correction-jet",
        )
