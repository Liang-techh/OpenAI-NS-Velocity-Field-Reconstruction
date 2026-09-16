from decimal import (
    Decimal,
    ROUND_CEILING,
    ROUND_FLOOR,
    ROUND_HALF_EVEN,
    localcontext,
)
from fractions import Fraction
import math

import pytest

from openai_ns_reconstruction.axis_coefficient_amplitude import (
    actual_schedule_amplitude_log_state,
)
from openai_ns_reconstruction.axis_coefficient_formal_solver import (
    formal_axis_coefficient_solver,
)
from openai_ns_reconstruction.axis_coefficient_rational_data import (
    RationalAxisCoefficientData,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_natural_source import (
    actual_schedule_reference_natural_source_log_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def selected():
    data = _schedule_data()
    amplitude = actual_schedule_amplitude_log_state(data, 0.05, phase_samples=101)
    rational = RationalAxisCoefficientData(amplitude.data)
    source = actual_schedule_reference_natural_source_log_state(
        data,
        0.05,
        phase_samples=101,
    )
    x1 = actual_schedule_wide_first_picard_state(data, 0.05, phase_samples=101)
    return rational, amplitude, source, x1


def _bell_formula(
    rational: RationalAxisCoefficientData,
    power: int,
    order: int,
    eta: float,
    Lambda: Decimal,
) -> Fraction:
    if order == 0:
        return Fraction(1)
    t = power * Fraction(Lambda)
    g0 = rational.jet_fraction("normalizedGradient", 0, 0, eta)
    if order == 1:
        return t * g0
    g1 = rational.jet_fraction("normalizedGradient", 0, 1, eta)
    if order == 2:
        return t * g1 + t * t * g0 * g0
    assert order == 3
    g2 = rational.jet_fraction("normalizedGradient", 0, 2, eta)
    return t * g2 + 3 * t * t * g0 * g1 + t**3 * g0**3


def _nearest96(value: Fraction) -> Decimal:
    with localcontext() as context:
        context.prec = 96
        context.rounding = ROUND_HALF_EVEN
        return +(Decimal(value.numerator) / Decimal(value.denominator))


@pytest.mark.parametrize("power", (1, 2, 3))
def test_selected_bell_orders_match_independent_closed_forms(selected, power) -> None:
    rational, amplitude, _, _ = selected
    eta = 0.23
    for order in range(4):
        expected = _bell_formula(rational, power, order, eta, amplitude.Lambda)
        assert rational.amplitude_power_bell_fraction(
            power, order, eta, amplitude.Lambda
        ) == expected


def test_zero_power_and_square_bell_leibniz_convolution(selected) -> None:
    rational, amplitude, _, _ = selected
    eta = -0.17
    for order in range(1, 5):
        assert rational.amplitude_power_bell_fraction(
            0, order, eta, amplitude.Lambda
        ) == 0
    assert rational.amplitude_power_bell_fraction(0, 0, eta, amplitude.Lambda) == 1

    first = [
        rational.amplitude_power_bell_fraction(1, order, eta, amplitude.Lambda)
        for order in range(5)
    ]
    for order in range(5):
        expected = sum(
            (math.comb(order, k) * first[k] * first[order - k] for k in range(order + 1)),
            Fraction(0),
        )
        assert rational.amplitude_power_bell_fraction(
            2, order, eta, amplitude.Lambda
        ) == expected


def test_decimal_enclosure_and_actual_consumers_are_context_independent(selected) -> None:
    rational, amplitude, source, x1 = selected
    eta = 0.23
    order = 4
    exact = rational.amplitude_power_bell_fraction(2, order, eta, amplitude.Lambda)
    expected_decimal = _nearest96(exact)

    with localcontext() as context:
        context.prec = 8
        context.rounding = ROUND_FLOOR
        low_decimal = rational.amplitude_power_bell_decimal(
            2, order, eta, amplitude.Lambda
        )
        low_enclosure = rational.amplitude_power_bell_enclosure(
            2, order, eta, amplitude.Lambda
        )
    with localcontext() as context:
        context.prec = 120
        context.rounding = ROUND_CEILING
        high_decimal = rational.amplitude_power_bell_decimal(
            2, order, eta, amplitude.Lambda
        )
        high_enclosure = rational.amplitude_power_bell_enclosure(
            2, order, eta, amplitude.Lambda
        )

    assert low_decimal == high_decimal == expected_decimal
    assert low_enclosure.lower == high_enclosure.lower
    assert low_enclosure.upper == high_enclosure.upper
    assert Fraction(low_enclosure.lower) <= exact <= Fraction(low_enclosure.upper)

    for order in range(5):
        assert amplitude._bell_factor(order, eta) == rational.amplitude_power_bell_decimal(
            1, order, eta, amplitude.Lambda
        )
        assert source._amplitude_square_bell_factor(
            order, eta
        ) == rational.amplitude_power_bell_decimal(
            2, order, eta, amplitude.Lambda
        )

    # The formal pressure graph consumes the same actual q=2 source channel.
    formal = formal_axis_coefficient_solver(x1)
    pressure = formal.profile_jet_prefix(1, 0, eta)[1][2]
    assert pressure.channels.get((2, 0)) == Decimal(1)

