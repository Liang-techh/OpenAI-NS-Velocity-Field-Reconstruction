from decimal import localcontext, Decimal
from fractions import Fraction

import pytest

from openai_ns_reconstruction.axis_amplitude_derivative_enclosure import (
    AmplitudeDerivativeLogEnclosure,
    rational_log_enclosure,
)
from openai_ns_reconstruction.axis_coefficient_amplitude import (
    actual_schedule_amplitude_log_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _amplitude():
    return actual_schedule_amplitude_log_state(
        TailData(
            OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
            h=0.01,
        ),
        0.05,
    )


def test_positive_logs_are_bounded_against_independent_decimal_reference() -> None:
    tolerance = Fraction(1, 10**96)
    with localcontext() as context:
        context.prec = 180
        for value in (Fraction(1), Fraction(2), Fraction(1, 2), Fraction(3, 2)):
            enclosure = rational_log_enclosure(
                value,
                absolute_tolerance=tolerance,
            )
            reference = (
                Decimal(value.numerator) / Decimal(value.denominator)
            ).ln()
            assert enclosure.width <= tolerance
            assert enclosure.lower <= Fraction(reference) <= enclosure.upper


def test_width_cap_and_reciprocal_intervals() -> None:
    tolerance = Fraction(1, 10**30)
    positive = rational_log_enclosure(
        Fraction(3, 2), absolute_tolerance=tolerance
    )
    reciprocal = rational_log_enclosure(
        Fraction(2, 3), absolute_tolerance=tolerance
    )
    assert positive.width <= tolerance
    assert reciprocal.width <= tolerance
    negated = (-positive.upper, -positive.lower)
    assert max(reciprocal.lower, negated[0]) <= min(reciprocal.upper, negated[1])


def test_invalid_inputs_and_series_cap_fail_closed() -> None:
    with pytest.raises(ValueError, match="value must be positive"):
        rational_log_enclosure(Fraction(0))
    with pytest.raises(TypeError, match="value must be a Fraction"):
        rational_log_enclosure(2)
    with pytest.raises(ValueError, match="absolute_tolerance must be positive"):
        rational_log_enclosure(Fraction(1), absolute_tolerance=Fraction(0))
    with pytest.raises(TypeError, match="absolute_tolerance must be a Fraction"):
        rational_log_enclosure(Fraction(1), absolute_tolerance=0.1)
    with pytest.raises(ValueError, match="max_terms must be a positive integer"):
        rational_log_enclosure(Fraction(1), max_terms=True)
    with pytest.raises(ArithmeticError, match="max_terms"):
        rational_log_enclosure(
            Fraction(3, 2),
            absolute_tolerance=Fraction(1, 10**96),
            max_terms=1,
        )


def test_actual_amplitude_bell_factors_and_full_bounds() -> None:
    amplitude = _amplitude()
    source = amplitude.log_amplitude_source(0.0)
    h = amplitude.rational_data.h
    j = amplitude.rational_data.j
    sigma = amplitude.rational_data.sigma

    # At eta = 0, H = j, H' = (1/2 - h) + 4, L = 1 and L' = 0.
    H0 = j
    H1 = Fraction(1, 2) - h + 4
    L0 = Fraction(1)
    L1 = Fraction(0)
    denominator = H0 * H0 + sigma * sigma
    g0 = -L0 * H0 / denominator
    g1 = -(L1 * H0 + L0 * H1) / denominator + (
        L0 * H0 * (2 * H0 * H1) / (denominator * denominator)
    )

    for power in (1, 2):
        lambda_power = power * Fraction(amplitude.Lambda)
        expected = (
            Fraction(1),
            lambda_power * g0,
            lambda_power * g1 + (lambda_power * g0) ** 2,
        )
        for order, expected_factor in enumerate(expected):
            enclosure = amplitude.derivative_log_enclosure(power, order, 0.0)
            assert enclosure.exact_bell_factor == expected_factor
            assert enclosure.sign == (1 if expected_factor > 0 else -1)
            assert enclosure.log_factor is not None
            assert enclosure.total_log_lower == (
                power * source.enclosure.lower + enclosure.log_factor.lower
            )
            assert enclosure.total_log_upper == (
                power * source.enclosure.upper + enclosure.log_factor.upper
            )
            assert enclosure.total_log_width == (
                power * source.enclosure.width + enclosure.log_factor.width
            )


def test_zero_factor_has_no_log_interval_and_flags_are_explicit() -> None:
    source = _amplitude().log_amplitude_source(0.0)
    enclosure = AmplitudeDerivativeLogEnclosure(
        source=source,
        power=1,
        order=3,
        exact_bell_factor=Fraction(0),
        log_factor=None,
    )
    assert enclosure.sign == 0
    assert enclosure.total_log_lower is None
    assert enclosure.total_log_upper is None
    assert enclosure.total_log_width is None
    assert enclosure.paper_exact is False
    assert enclosure.downstream_arithmetic_certified is False
