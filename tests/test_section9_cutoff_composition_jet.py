from dataclasses import replace
from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.section9_correction_extension_admission import (
    Section9CorrectionStageAdmissionCertificate,
)
from openai_ns_reconstruction.section9_cutoff_composition_jet import (
    Section9ScalarCutoffDerivativeJet,
    Section9SimilarityCoordinateJet,
    derive_section9_cutoff_weight_jet,
)
from openai_ns_reconstruction.section9_eq921_prefix_jet import (
    required_spacetime_multiindices,
)
from openai_ns_reconstruction.section9_finite_prefix import (
    Section9FinitePrefixEvaluationCertificate,
    Section9FinitePrefixStageContribution,
)


ZERO = (0, 0, 0, 0)
T = (1, 0, 0, 0)
X = (0, 1, 0, 0)
TX = (1, 1, 0, 0)


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
        cutoff_evaluator_provenance="fixture:scalar-cutoff",
    )
    assert result.formal_prefix_ready
    return result


def _q_jet(order: int = 2) -> Section9SimilarityCoordinateJet:
    derivatives = {
        alpha: 0.0 for alpha in required_spacetime_multiindices(order)
    }
    derivatives[ZERO] = 1.0 / 16.0
    derivatives[T] = 0.1
    derivatives[X] = -0.2
    derivatives[TX] = 0.3
    return Section9SimilarityCoordinateJet(
        q=Fraction(1, 16),
        derivative_order=order,
        derivatives=derivatives,
        provider_id="fixture-q-provider",
        provider_revision="q-r1",
        provider_provenance="fixture:eq4.1-q-jet-not-yet-machine-derived",
    )


def _cutoff_jet(order: int = 2) -> Section9ScalarCutoffDerivativeJet:
    # At s=a_1 q=1/4, use manufactured chi, chi', chi'' data.
    return Section9ScalarCutoffDerivativeJet(
        argument=Fraction(1, 4),
        derivative_order=order,
        derivatives={0: 0.5, 1: 2.0, 2: 3.0},
        provider_id="fixture-chi-provider",
        provider_revision="chi-r1",
        provider_provenance="fixture:scalar-cutoff",
    )


def test_normalized_taylor_composition_derives_mixed_spacetime_derivative():
    certificate = derive_section9_cutoff_weight_jet(
        _evaluation(), _admission(), _q_jet(), _cutoff_jet()
    )
    jet = certificate.weight_jet

    assert certificate.formal_composition_ready
    assert not certificate.similarity_coordinate_jet_machine_derived_from_eq_4_1
    assert not certificate.paper_fixed_cutoff_derivatives_machine_verified
    assert not certificate.paper_exact_velocity_available
    assert jet.derivatives[ZERO] == 0.5

    # Independent literal chain-rule checks with a=4:
    # w_t = chi'(aq) a q_t
    # w_tx = chi''(aq) (a q_t)(a q_x) + chi'(aq) a q_tx.
    assert jet.derivatives[T] == pytest.approx(2.0 * 4.0 * 0.1)
    assert jet.derivatives[X] == pytest.approx(2.0 * 4.0 * -0.2)
    expected_tx = 3.0 * (4.0 * 0.1) * (4.0 * -0.2) + 2.0 * 4.0 * 0.3
    assert jet.derivatives[TX] == pytest.approx(expected_tx)


def test_composition_rejects_q_identity_drift():
    original = _q_jet()
    derivatives = dict(original.derivatives)
    derivatives[ZERO] = 1.0 / 15.0
    q_jet = Section9SimilarityCoordinateJet(
        q=Fraction(1, 15),
        derivative_order=original.derivative_order,
        derivatives=derivatives,
        provider_id=original.provider_id,
        provider_revision=original.provider_revision,
        provider_provenance=original.provider_provenance,
    )
    with pytest.raises(ValueError, match="finite-prefix evaluation q"):
        derive_section9_cutoff_weight_jet(
            _evaluation(), _admission(), q_jet, _cutoff_jet()
        )


def test_composition_rejects_scalar_argument_drift():
    cutoff = replace(_cutoff_jet(), argument=Fraction(1, 5))
    with pytest.raises(ValueError, match="argument must equal"):
        derive_section9_cutoff_weight_jet(
            _evaluation(), _admission(), _q_jet(), cutoff
        )


def test_composition_rejects_scalar_cutoff_provenance_drift():
    cutoff = replace(_cutoff_jet(), provider_provenance="other-cutoff")
    with pytest.raises(ValueError, match="provenance"):
        derive_section9_cutoff_weight_jet(
            _evaluation(), _admission(), _q_jet(), cutoff
        )


def test_q_jet_requires_complete_multiindex_coverage():
    jet = _q_jet()
    derivatives = dict(jet.derivatives)
    del derivatives[TX]
    with pytest.raises(ValueError, match="cover exactly every spacetime derivative"):
        Section9SimilarityCoordinateJet(
            q=jet.q,
            derivative_order=jet.derivative_order,
            derivatives=derivatives,
            provider_id=jet.provider_id,
            provider_revision=jet.provider_revision,
            provider_provenance=jet.provider_provenance,
        )
