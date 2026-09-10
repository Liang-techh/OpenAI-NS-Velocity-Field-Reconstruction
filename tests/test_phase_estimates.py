import math
import pytest

from openai_ns_reconstruction.phase_estimates import (
    PhaseScaleCertificate,
    phase_constant,
    phase_error,
)


def test_phase_error_and_constant_match_pinned_lean_definitions():
    S, eps, k, M = 10.0, 0.001, 1000.0, 2.0
    assert math.isclose(
        phase_error(S, eps, k),
        1 / S + S * eps**2 + S / k + eps * S,
        rel_tol=0.0,
        abs_tol=0.0,
    )
    assert phase_constant(M) == 8 * M**3 + 2 * M**4


def test_phase_error_four_over_S_certificate():
    # Exact admissible boundary pattern: epsilon*S^2 = 1 and S^2/k = 1.
    S = 20.0
    cert = PhaseScaleCertificate(M=3.0, S=S, epsilon=1 / S**2, k=S**2)
    assert cert.scale_hypotheses() == {
        "S2_epsilon2_le_1": True,
        "S2_over_k_le_1": True,
        "epsilon_S2_le_1": True,
    }
    assert cert.verify_scalar_inequality()
    assert cert.error <= cert.four_over_S
    assert cert.rounded_normal_bound <= cert.simplified_normal_bound


def test_certificate_fails_closed_when_any_scale_hypothesis_is_missing():
    S = 10.0
    with pytest.raises(ValueError, match="epsilon_S2_le_1"):
        PhaseScaleCertificate(M=2.0, S=S, epsilon=0.02, k=1000.0)
    with pytest.raises(ValueError, match="S2_over_k_le_1"):
        PhaseScaleCertificate(M=2.0, S=S, epsilon=1e-4, k=50.0)


def test_certificate_does_not_accept_invalid_domains():
    for kwargs in (
        {"M": 0.5, "S": 10.0, "epsilon": 1e-4, "k": 1000.0},
        {"M": 2.0, "S": 0.5, "epsilon": 1e-4, "k": 1000.0},
        {"M": 2.0, "S": 10.0, "epsilon": -1e-4, "k": 1000.0},
        {"M": 2.0, "S": 10.0, "epsilon": 1e-4, "k": 0.5},
        {"M": 2.0, "S": float("nan"), "epsilon": 1e-4, "k": 1000.0},
    ):
        with pytest.raises(ValueError):
            PhaseScaleCertificate(**kwargs)


def test_negative_slack_is_rejected():
    cert = PhaseScaleCertificate(M=2.0, S=10.0, epsilon=1e-4, k=1000.0)
    with pytest.raises(ValueError):
        cert.verify_scalar_inequality(slack=-1.0)
