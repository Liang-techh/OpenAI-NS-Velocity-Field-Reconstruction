from fractions import Fraction
import math

import numpy as np
import pytest

from openai_ns_reconstruction.section9_correction_extension_admission import (
    Section9CorrectionStageAdmissionCertificate,
)
from openai_ns_reconstruction.section9_cutoff_composition_jet import (
    Section9ScalarCutoffDerivativeJet,
)
from openai_ns_reconstruction.section9_eq41_similarity_jet import (
    derive_section9_cutoff_weight_jet_from_eq41,
    derive_section9_similarity_coordinate_jet_from_eq41,
)
from openai_ns_reconstruction.section9_finite_prefix import (
    Section9FinitePrefixEvaluationCertificate,
    Section9FinitePrefixStageContribution,
)


ZERO = (0, 0, 0, 0)
T = (1, 0, 0, 0)
X = (0, 1, 0, 0)
Y = (0, 0, 1, 0)
Z = (0, 0, 0, 1)
TT = (2, 0, 0, 0)
TX = (1, 1, 0, 0)


def test_eq41_normalized_taylor_recurrence_matches_independent_derivative_identities():
    tau = 0.1
    z = 0.2
    h = 0.25
    certificate = derive_section9_similarity_coordinate_jet_from_eq41(
        tau=tau,
        z=z,
        h=h,
        derivative_order=3,
        provider_id="paper-eq4.1",
        provider_revision="fixture-r1",
        provider_provenance="paper:Eq. (4.1) implicit relation",
    )
    jet = certificate.q_jet
    q = float(jet.q)
    exponent = 2.0 * h
    jacobian = 1.0 - exponent * z * z * q ** (exponent - 1.0)

    assert certificate.formal_eq41_jet_ready
    assert certificate.similarity_coordinate_jet_machine_derived_from_eq_4_1
    assert not certificate.used_finite_differences
    assert not certificate.paper_fixed_cutoff_derivatives_machine_verified
    assert not certificate.paper_exact_velocity_available
    assert q - z * z * q**exponent == pytest.approx(tau, rel=5e-12, abs=5e-12)
    assert certificate.implicit_jacobian == pytest.approx(jacobian)

    # Independent implicit-differentiation identities for
    # F=q-z^2 q^(2h)-(1-t)=0.
    q_t = -1.0 / jacobian
    q_z = 2.0 * z * q**exponent / jacobian
    f_qq = -exponent * (exponent - 1.0) * z * z * q ** (exponent - 2.0)
    q_tt = -f_qq * q_t * q_t / jacobian
    assert jet.derivatives[T] == pytest.approx(q_t)
    assert jet.derivatives[Z] == pytest.approx(q_z)
    assert jet.derivatives[TT] == pytest.approx(q_tt)

    # Eq. (4.1) has no x/y dependence, including mixed derivatives.
    assert jet.derivatives[X] == 0.0
    assert jet.derivatives[Y] == 0.0
    assert jet.derivatives[TX] == 0.0


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


def _cutoff_jet() -> Section9ScalarCutoffDerivativeJet:
    return Section9ScalarCutoffDerivativeJet(
        argument=Fraction(1, 4),
        derivative_order=2,
        derivatives={0: 0.5, 1: 2.0, 2: 3.0},
        provider_id="fixture-chi-provider",
        provider_revision="chi-r1",
        provider_provenance="fixture:scalar-cutoff",
    )


def test_eq41_jet_feeds_existing_cutoff_composition_without_upgrading_fixed_cutoff_truth():
    # z=0 makes Eq. (4.1) exactly q=tau=1/16, matching the independently
    # admitted finite-prefix evaluation identity without numeric tolerance glue.
    certificate = derive_section9_cutoff_weight_jet_from_eq41(
        _evaluation(),
        _admission(),
        _cutoff_jet(),
        tau=1.0 / 16.0,
        z=0.0,
        h=0.25,
        provider_id="paper-eq4.1",
        provider_revision="fixture-r1",
        provider_provenance="paper:Eq. (4.1) implicit relation",
    )

    assert certificate.formal_composition_ready
    assert certificate.similarity_coordinate_jet_machine_derived_from_eq_4_1
    assert not certificate.paper_fixed_cutoff_derivatives_machine_verified
    assert not certificate.paper_exact_velocity_available
    assert certificate.composition.formal_composition_ready

    # At z=0, q_t=-1.  With a_1=4 and chi'=2, the independent chain rule gives
    # d_t chi(a_1 q)=chi' * a_1 * q_t=-8.
    assert certificate.composition.weight_jet.derivatives[T] == pytest.approx(-8.0)


def test_eq41_adapter_rejects_nonphysical_endpoint_and_invalid_order():
    kwargs = dict(
        z=0.0,
        h=0.25,
        provider_id="paper-eq4.1",
        provider_revision="fixture-r1",
        provider_provenance="paper:Eq. (4.1) implicit relation",
    )
    with pytest.raises(ValueError, match="tau must be positive"):
        derive_section9_similarity_coordinate_jet_from_eq41(
            tau=0.0, derivative_order=2, **kwargs
        )
    with pytest.raises(ValueError, match="derivative_order"):
        derive_section9_similarity_coordinate_jet_from_eq41(
            tau=0.1, derivative_order=-1, **kwargs
        )
