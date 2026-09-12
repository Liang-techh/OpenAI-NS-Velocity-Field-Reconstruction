from dataclasses import replace
from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.section9_correction_extension_admission import (
    Section9CorrectionStageAdmissionCertificate,
)
from openai_ns_reconstruction.section9_eq921_prefix_jet import (
    Section9WeightedCorrectionStageJet,
    assemble_eq_9_21_prefix_jet,
    required_spacetime_multiindices,
)
from openai_ns_reconstruction.section9_finite_prefix import (
    Section9FinitePrefixEvaluationCertificate,
    Section9FinitePrefixStageContribution,
)


ZERO = (0, 0, 0, 0)


def _stage_cert(stage: int, scale: int) -> Section9CorrectionStageAdmissionCertificate:
    return Section9CorrectionStageAdmissionCertificate(
        stage=stage,
        q_big=Fraction(1, 2),
        scale=scale,
        schedule_lower_bound=4 * (1 << (stage - 1)),
        components=("A", "B", "p"),
        evidence_kinds=("paper-derived", "paper-derived", "paper-derived"),
        evidence_provenance=(
            f"paper:A:{stage}",
            f"paper:B:{stage}",
            f"paper:p:{stage}",
        ),
    )


def _evaluation() -> Section9FinitePrefixEvaluationCertificate:
    contributions = (
        Section9FinitePrefixStageContribution(
            stage=1,
            scale=4,
            cutoff_argument=Fraction(1, 4),
            cutoff_weight=0.5,
            A_term=np.array([1.0, 2.0, 3.0]),
            B_term=4.0,
            p_term=5.0,
            value_provenance="fixture:stage-1",
            support_zero_short_circuit=False,
        ),
        Section9FinitePrefixStageContribution(
            stage=2,
            scale=8,
            cutoff_argument=Fraction(1, 2),
            cutoff_weight=0.25,
            A_term=np.array([-1.0, 0.5, 2.0]),
            B_term=-2.0,
            p_term=1.5,
            value_provenance="fixture:stage-2",
            support_zero_short_circuit=False,
        ),
    )
    result = Section9FinitePrefixEvaluationCertificate(
        q=Fraction(1, 16),
        q_big=Fraction(1, 2),
        prefix_order=2,
        stages=(1, 2),
        contributions=contributions,
        A_correction_prefix=np.array([0.0, 2.5, 5.0]),
        B_correction_prefix=2.0,
        p_correction_prefix=6.5,
        cutoff_evaluator_provenance="fixture:cutoff",
    )
    assert result.formal_prefix_ready
    return result


def _jet(
    stage: int,
    scale: int,
    A0: np.ndarray,
    B0: float,
    p0: float,
    *,
    revision: str = "r1",
) -> Section9WeightedCorrectionStageJet:
    order = 1
    alphas = required_spacetime_multiindices(order)
    A = {}
    B = {}
    p = {}
    for i, alpha in enumerate(alphas):
        if alpha == ZERO:
            A[alpha] = A0
            B[alpha] = B0
            p[alpha] = p0
        else:
            A[alpha] = np.array([stage + i, stage + 2 * i, stage - i], dtype=float)
            B[alpha] = float(10 * stage + i)
            p[alpha] = float(-10 * stage - i)
    cert = _stage_cert(stage, scale)
    return Section9WeightedCorrectionStageJet(
        stage=stage,
        scale=scale,
        q=Fraction(1, 16),
        derivative_order=order,
        A_derivatives=A,
        B_derivatives=B,
        p_derivatives=p,
        evidence_kinds=cert.evidence_kinds,
        evidence_provenance=cert.evidence_provenance,
        provider_id="fixture-provider",
        provider_revision=revision,
        provider_provenance="analytic fixture only",
    )


def _inputs():
    evaluation = _evaluation()
    certs = (_stage_cert(1, 4), _stage_cert(2, 8))
    jets = (
        _jet(1, 4, np.array([1.0, 2.0, 3.0]), 4.0, 5.0),
        _jet(2, 8, np.array([-1.0, 0.5, 2.0]), -2.0, 1.5),
    )
    return evaluation, certs, jets


def test_prefix_jet_sums_complete_analytic_derivatives_and_preserves_truth_boundary():
    evaluation, certs, jets = _inputs()
    result = assemble_eq_9_21_prefix_jet(evaluation, certs, jets)

    assert result.formal_prefix_jet_ready
    assert result.derivative_order == 1
    assert set(result.A_derivatives) == set(required_spacetime_multiindices(1))
    np.testing.assert_array_equal(
        result.A_derivatives[ZERO],
        evaluation.A_correction_prefix,
    )
    assert result.B_derivatives[ZERO] == evaluation.B_correction_prefix
    assert result.p_derivatives[ZERO] == evaluation.p_correction_prefix

    alpha = (1, 0, 0, 0)
    np.testing.assert_array_equal(
        result.A_derivatives[alpha],
        jets[0].A_derivatives[alpha] + jets[1].A_derivatives[alpha],
    )
    assert result.B_derivatives[alpha] == (
        jets[0].B_derivatives[alpha] + jets[1].B_derivatives[alpha]
    )
    assert result.p_derivatives[alpha] == (
        jets[0].p_derivatives[alpha] + jets[1].p_derivatives[alpha]
    )

    assert not result.analytic_derivatives_machine_derived_from_actual_corrections
    assert not result.paper_fixed_cutoff_derivatives_machine_verified
    assert not result.actual_section9_sequence_verified
    assert not result.residual_artifact_ready
    assert not result.paper_exact_velocity_available


def test_stage_jet_requires_complete_multiindex_coverage():
    cert = _stage_cert(1, 4)
    alphas = required_spacetime_multiindices(1)
    A = {alpha: np.ones(3) for alpha in alphas}
    B = {alpha: 1.0 for alpha in alphas}
    p = {alpha: 1.0 for alpha in alphas}
    del A[(1, 0, 0, 0)]

    with pytest.raises(ValueError, match="cover exactly all derivatives"):
        Section9WeightedCorrectionStageJet(
            stage=1,
            scale=4,
            q=Fraction(1, 16),
            derivative_order=1,
            A_derivatives=A,
            B_derivatives=B,
            p_derivatives=p,
            evidence_kinds=cert.evidence_kinds,
            evidence_provenance=cert.evidence_provenance,
            provider_id="fixture-provider",
            provider_revision="r1",
            provider_provenance="fixture",
        )


def test_prefix_jet_rejects_zero_order_switch_from_evaluated_stage_value():
    evaluation, certs, jets = _inputs()
    bad_A = dict(jets[1].A_derivatives)
    bad_A[ZERO] = np.array([999.0, 0.5, 2.0])
    bad = replace(jets[1], A_derivatives=bad_A)

    with pytest.raises(ValueError, match="zero-order A term"):
        assemble_eq_9_21_prefix_jet(evaluation, certs, (jets[0], bad))


def test_prefix_jet_rejects_stage_provenance_swap():
    evaluation, certs, jets = _inputs()
    bad = replace(
        jets[1],
        evidence_provenance=("paper:A:1", "paper:B:2", "paper:p:2"),
    )

    with pytest.raises(ValueError, match="evidence provenance"):
        assemble_eq_9_21_prefix_jet(evaluation, certs, (jets[0], bad))


def test_prefix_jet_rejects_mixed_provider_revisions():
    evaluation, certs, jets = _inputs()
    bad = replace(jets[1], provider_revision="r2")

    with pytest.raises(ValueError, match="one provider"):
        assemble_eq_9_21_prefix_jet(evaluation, certs, (jets[0], bad))
