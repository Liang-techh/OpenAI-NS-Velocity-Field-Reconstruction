import math
from fractions import Fraction

import pytest

from openai_ns_reconstruction.phase_estimates import phase_constant
from openai_ns_reconstruction.phase_large_band_scale import (
    LargeBandPhaseScaleCertificate,
    exact_large_band_hypotheses,
)


def test_default_h_large_band_crosses_carrier_sufficient_gate_at_23204():
    """Independent log check of the nontrivial carrier sufficient boundary."""

    h = Fraction(1, 200)
    below = exact_large_band_hypotheses(23203, h)
    at = exact_large_band_hypotheses(23204, h)

    assert below["S2_epsilon2_le_1"]
    assert not below["S2_over_carrier_le_1_via_carrier_lower_bound"]
    assert below["epsilon_S2_le_1"]
    assert all(at.values())

    # Independently rewrite the sufficient condition as
    #   ell^4 <= 2^(ell/400)
    # and compare its two sides in log2 scale, without reusing integer powers.
    assert 4.0 * math.log2(23203) > 23203.0 / 400.0
    assert 4.0 * math.log2(23204) < 23204.0 / 400.0


def test_certificate_works_beyond_binary64_Q_range_without_forming_Q():
    cert = LargeBandPhaseScaleCertificate(ell=23204, h=0.005, M=4.0)

    assert cert.h_fraction == Fraction(1, 200)
    assert cert.S_star == 23204**2
    assert cert.epsilon_log2 == Fraction(-23204, 200)
    assert cert.inverse_sqrt_epsilon_log2 == Fraction(23204, 400)
    assert all(cert.scale_hypotheses().values())

    # The existing chart path intentionally cannot represent Q this deep.
    assert math.ldexp(1.0, -23204) == 0.0
    assert cert.simplified_phase_error_bound == pytest.approx(4.0 / 23204**2)
    assert cert.simplified_rounded_normal_bound == pytest.approx(
        phase_constant(4.0) * 4.0 / 23204**2
    )


def test_certificate_preserves_fail_closed_paper_exact_boundaries():
    cert = LargeBandPhaseScaleCertificate(ell=23204, h=Fraction(1, 200), M=2.0)

    assert cert.status == "formal-structure"
    assert cert.actual_base_fields_verified is False
    assert cert.uniform_eq_7_9_to_7_11_verified is False
    assert cert.paper_exact_velocity_available is False


def test_certificate_rejects_band_below_sufficient_threshold():
    with pytest.raises(ValueError, match="S2_over_carrier_le_1"):
        LargeBandPhaseScaleCertificate(ell=23203, h=Fraction(1, 200), M=2.0)


def test_large_band_inputs_fail_closed():
    for ell in (0, -1, True, 1.5):
        with pytest.raises((TypeError, ValueError)):
            exact_large_band_hypotheses(ell, Fraction(1, 200))

    for h in (0.0, 0.5, -0.01, float("nan"), True):
        with pytest.raises((TypeError, ValueError)):
            exact_large_band_hypotheses(23204, h)

    for M in (0.5, float("inf"), float("nan"), True):
        with pytest.raises((TypeError, ValueError)):
            LargeBandPhaseScaleCertificate(ell=23204, h=Fraction(1, 200), M=M)
