from decimal import Decimal, localcontext
from fractions import Fraction

import pytest

from openai_ns_reconstruction.axis_coefficient_rational_pressure import (
    actual_schedule_reference_rational_pressure,
)
from openai_ns_reconstruction.axis_coefficient_wide_natural_pressure import (
    actual_schedule_reference_wide_pressure_log_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


@pytest.fixture(scope="module")
def amplitude():
    from openai_ns_reconstruction.axis_coefficient_amplitude import (
        actual_schedule_amplitude_log_state,
    )

    return actual_schedule_amplitude_log_state(
        TailData(
            OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
            h=0.01,
        ),
        0.05,
    )


@pytest.fixture(scope="module")
def rational_pressure(amplitude):
    return actual_schedule_reference_rational_pressure(amplitude)


def test_exact_low_radial_source_and_pressure_rows(rational_pressure) -> None:
    rd = rational_pressure.amplitude.rational_data
    chi = rd.jet_fraction("chi", 0, 0, 0.0)
    assert rational_pressure.source_factor(0, 0, 0.0) == Fraction(1)
    assert rational_pressure.source_factor(1, 0, 0.0) == -chi / 2
    assert rational_pressure.source_factor(2, 0, 0.0) == 5 * chi * chi / 48
    assert rational_pressure.pressure_factor(0, 0, 0.0) == 0
    assert rational_pressure.pressure_factor(1, 0, 0.0) == 0

    source = rational_pressure.source_log_enclosure(0, 0, 0.0)
    pressure = rational_pressure.pressure_log_enclosure(0, 0, 0.0)
    assert source.exact_factor == 1
    assert source.sign == 1
    assert source.total_log_lower == 2 * source.source.enclosure.lower
    assert source.total_log_upper == 2 * source.source.enclosure.upper
    assert pressure.sign == 0
    assert pressure.log_factor is None
    assert pressure.total_log_lower is None
    assert pressure.total_log_upper is None


def test_independent_low_order_pressure_identity_and_diagnostic(amplitude, rational_pressure) -> None:
    rd = amplitude.rational_data
    h = rd.h
    j = rd.j
    sigma = rd.sigma
    for eta in (0.0, 0.01):
        eta_fraction = Fraction.from_float(eta)
        # H and H' are calculated directly from the defining polynomial,
        # independently of the rational source/pressure implementation.
        D = Fraction(1, 2) - h
        H = D * eta_fraction + (1 - eta_fraction * eta_fraction) * (4 * eta_fraction + j)
        H_prime = D + 4 - 12 * eta_fraction * eta_fraction - 2 * j * eta_fraction
        denominator = H * H + sigma * sigma
        L = Fraction(1) - 2 * h * eta_fraction * eta_fraction
        L_prime = -4 * h * eta_fraction
        g = -L * H / denominator
        g_prime = (
            -(L_prime * H + L * H_prime) / denominator
            + L * H * (2 * H * H_prime) / (denominator * denominator)
        )
        B1 = 2 * Fraction(amplitude.Lambda) * g
        B2 = 2 * Fraction(amplitude.Lambda) * g_prime + B1 * B1
        d = rd.jet_fraction("d", 0, 0, eta)
        d_prime = rd.jet_fraction("d", 0, 1, eta)
        expected_pressure_zero = (d * B1 - (4 * rd.A + 2) * eta_fraction) / 4
        # At eta = 0 this is the compact form requested; at nonzero eta the
        # eta*B1 product-rule term is retained explicitly.
        expected_pressure_first = (
            d_prime * B1
            + d * B2
            - (4 * rd.A + 2) * (1 + eta_fraction * B1)
        ) / 4
        assert rational_pressure.source_factor(0, 2, eta) == B2
        assert rational_pressure.pressure_factor(2, 0, eta) == expected_pressure_zero
        assert rational_pressure.pressure_factor(2, 1, eta) == expected_pressure_first
        assert B2 != 0

    # The n = 3, m = 1 identity checks the radial/eta operator indices using
    # already-computed exact source factors, rather than reproving the source.
    eta = 0.01
    eta_fraction = Fraction.from_float(eta)
    S10 = rational_pressure.source_factor(1, 0, eta)
    S11 = rational_pressure.source_factor(1, 1, eta)
    S12 = rational_pressure.source_factor(1, 2, eta)
    n3_m1 = (
        rd.jet_fraction("d", 0, 0, eta) * S12
        - (4 * rd.A + 4) * rd.jet_fraction("eta", 0, 0, eta) * S11
        + rd.jet_fraction("d", 0, 1, eta) * S11
        - (4 * rd.A + 4) * rd.jet_fraction("eta", 0, 1, eta) * S10
    ) / 18
    assert rational_pressure.pressure_factor(3, 1, eta) == n3_m1

    # Existing wide state is only a rounded diagnostic comparison.
    wide = actual_schedule_reference_wide_pressure_log_state(
        TailData(
            OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
            h=0.01,
        ),
        0.05,
    )
    with localcontext() as context:
        context.prec = 96
        exact_decimal = Decimal(
            rational_pressure.pressure_factor(2, 0, eta).numerator
        ) / Decimal(rational_pressure.pressure_factor(2, 0, eta).denominator)
        rounded_decimal = wide.normalized_factor(2, 0, eta)
        assert abs(exact_decimal - rounded_decimal) <= (
            max(abs(exact_decimal), abs(rounded_decimal)) * Decimal("1e-12")
            + Decimal("1e-12")
        )


def test_log_bounds_sign_and_exact_total_width(rational_pressure) -> None:
    result = rational_pressure.source_log_enclosure(0, 1, 0.01)
    assert result.exact_factor != 0
    assert result.log_factor is not None
    assert result.sign == (1 if result.exact_factor > 0 else -1)
    assert result.total_log_lower == 2 * result.source.enclosure.lower + result.log_factor.lower
    assert result.total_log_upper == 2 * result.source.enclosure.upper + result.log_factor.upper
    assert result.total_log_width == (
        2 * result.source.enclosure.width + result.log_factor.width
    )


def test_invalid_inputs_and_caps_fail_before_zero(rational_pressure) -> None:
    with pytest.raises(ValueError, match="absolute_log_factor_tolerance must be positive"):
        rational_pressure.pressure_log_enclosure(
            0, 0, 0.0, absolute_log_factor_tolerance=Fraction(0)
        )
    with pytest.raises(ValueError, match="max_terms must be a positive integer"):
        rational_pressure.pressure_log_enclosure(0, 0, 0.0, max_terms=True)
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        rational_pressure.source_log_enclosure(-1, 0, 0.0)
    with pytest.raises(ValueError, match="kind must be"):
        rational_pressure.log_enclosure("other", 0, 0, 0.0)
    with pytest.raises(ArithmeticError, match="max_terms"):
        rational_pressure.source_log_enclosure(
            0, 1, 0.01, absolute_log_factor_tolerance=Fraction(1, 10**96), max_terms=1
        )
