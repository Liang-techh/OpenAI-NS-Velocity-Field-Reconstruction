from fractions import Fraction
from math import factorial

import pytest

from openai_ns_reconstruction.schedule_axis_pressure_kernel_enclosure import (
    RationalInterval,
    pressure_kernel_derivative_enclosure,
    pressure_kernel_normalized_taylor_enclosure,
    pressure_kernel_value_enclosure,
)


def test_a_zero_is_exact_for_interval_eta():
    eta = RationalInterval(Fraction(-3), Fraction(5))
    jet = pressure_kernel_normalized_taylor_enclosure(0, eta, 6)
    assert jet.normalized_coefficients[0] == RationalInterval.point(1)
    assert all(x == RationalInterval.point(0) for x in jet.normalized_coefficients[1:])
    assert jet.full_derivatives[0] == RationalInterval.point(1)
    assert all(x == RationalInterval.point(0) for x in jet.full_derivatives[1:])


def test_a_one_point_value_first_and_second_derivatives_are_exact():
    eta = Fraction(2, 3)
    jet = pressure_kernel_normalized_taylor_enclosure(1, eta, 2)
    assert jet.normalized_coefficients[0] == RationalInterval.point(Fraction(81, 169))
    assert jet.full_derivatives[1] == RationalInterval.point(Fraction(-1944, 2197))
    assert jet.full_derivatives[2] == RationalInterval.point(Fraction(32076, 28561))
    assert jet.normalized_coefficients[2] == RationalInterval.point(Fraction(16038, 28561))


def test_normalized_coefficients_are_distinct_from_full_derivatives():
    jet = pressure_kernel_normalized_taylor_enclosure(1, Fraction(2, 3), 5)
    for m, normalized in enumerate(jet.normalized_coefficients):
        assert jet.full_derivatives[m] == normalized.scale(factorial(m))


def test_interval_exponent_is_not_replaced_by_a_midpoint_and_refinement_is_nested():
    a = RationalInterval(Fraction(1, 4), Fraction(3, 4))
    eta = RationalInterval(Fraction(1, 3), Fraction(2, 3))
    coarse = pressure_kernel_normalized_taylor_enclosure(
        a, eta, 3, exponent_cells=1, eta_cells=1
    )
    refined = pressure_kernel_normalized_taylor_enclosure(
        a, eta, 3, exponent_cells=4, eta_cells=4
    )
    assert refined.exponent == a
    assert refined.eta == eta
    for c, r in zip(coarse.normalized_coefficients, refined.normalized_coefficients):
        assert c.lower <= r.lower <= r.upper <= c.upper

    for aa in (a.lower, a.upper):
        for ee in (eta.lower, eta.upper):
            point = pressure_kernel_normalized_taylor_enclosure(aa, ee, 3)
            for enclosing, exact_point in zip(
                refined.normalized_coefficients, point.normalized_coefficients
            ):
                assert enclosing.lower <= exact_point.lower
                assert exact_point.upper <= enclosing.upper


def test_eta_interval_crossing_zero_preserves_even_value_and_derivative_sign_range():
    eta = RationalInterval(Fraction(-1, 2), Fraction(1, 2))
    jet = pressure_kernel_normalized_taylor_enclosure(1, eta, 1, eta_cells=4)
    assert jet.normalized_coefficients[0].contains(1)
    assert jet.full_derivatives[1].contains(0)
    endpoint_value = Fraction(16, 25)
    endpoint_slope = Fraction(128, 125)
    assert jet.normalized_coefficients[0].contains(endpoint_value)
    assert jet.full_derivatives[1].contains(-endpoint_slope)
    assert jet.full_derivatives[1].contains(endpoint_slope)


def test_value_and_derivative_entrypoints_match_jet():
    a = Fraction(2, 5)
    eta = Fraction(3, 7)
    jet = pressure_kernel_normalized_taylor_enclosure(a, eta, 3)
    assert pressure_kernel_value_enclosure(a, eta) == jet.normalized_coefficients[0]
    assert pressure_kernel_derivative_enclosure(a, eta, 3) == jet.full_derivatives[3]


@pytest.mark.parametrize(
    "call",
    [
        lambda: pressure_kernel_value_enclosure(0.5, Fraction(1, 3)),
        lambda: pressure_kernel_value_enclosure(Fraction(1, 2), 0.25),
        lambda: pressure_kernel_value_enclosure(Fraction(-1, 10), 0),
        lambda: pressure_kernel_value_enclosure(Fraction(11, 10), 0),
        lambda: pressure_kernel_normalized_taylor_enclosure(Fraction(1, 2), 0, -1),
        lambda: pressure_kernel_normalized_taylor_enclosure(
            Fraction(1, 2), 0, 5, max_order=4
        ),
        lambda: pressure_kernel_normalized_taylor_enclosure(
            RationalInterval(Fraction(0), Fraction(1)),
            RationalInterval(Fraction(-1), Fraction(1)),
            1,
            exponent_cells=100,
            eta_cells=100,
            max_cells=16,
        ),
        lambda: pressure_kernel_normalized_taylor_enclosure(
            Fraction(1, 2), Fraction(10**20), 0, max_range_reductions=1
        ),
        lambda: pressure_kernel_value_enclosure(complex(0.5, 0.0), 0),
    ],
)
def test_fail_closed_types_domains_orders_and_caps(call):
    with pytest.raises((TypeError, ValueError, ArithmeticError)):
        call()
