from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_base_heat_initial_jet import (
    heat_initial_derivative,
    heat_initial_jet,
    ode_taylor_recurrence_residual,
    rising_factorial,
    supplied_recurrence_residual,
)


def test_rising_factorial_is_exact() -> None:
    assert rising_factorial(Fraction(1, 8), 4) == Fraction(3825, 4096)
    assert rising_factorial(Fraction(9, 8), 4) == Fraction(126225, 4096)


def test_public_gamma_moments_give_exact_initial_jet() -> None:
    h = Fraction(1, 8)
    assert heat_initial_jet(h, 4) == (
        Fraction(1),
        Fraction(-9, 64),
        Fraction(1377, 4096),
        Fraction(-585225, 262144),
        Fraction(482810625, 16777216),
    )


def test_normalization_and_first_derivative_are_the_ode_zero_order_relation() -> None:
    h = Fraction(1, 8)
    h0 = heat_initial_derivative(h, 0)
    h1 = heat_initial_derivative(h, 1)
    assert h0 == 1
    assert h1 == -h * (1 + h)
    assert h1 + h * (1 + h) * h0 == 0


def test_integral_moment_jet_satisfies_all_checked_ode_taylor_recurrences() -> None:
    h = Fraction(1, 113)
    assert [ode_taylor_recurrence_residual(h, k) for k in range(8)] == [0] * 8


def test_one_bit_scale_mutation_is_detected_exactly() -> None:
    h = Fraction(1, 8)
    h1 = heat_initial_derivative(h, 1)
    h2 = heat_initial_derivative(h, 2)
    mutation = Fraction(1, 2**40)
    assert supplied_recurrence_residual(h, 1, h1, h2 + mutation) == mutation


def test_positive_h_jet_has_alternating_signs() -> None:
    jet = heat_initial_jet(Fraction(1, 100), 7)
    assert all((value > 0) == (k % 2 == 0) for k, value in enumerate(jet))


def test_invalid_or_inexact_inputs_fail_closed() -> None:
    with pytest.raises(TypeError):
        heat_initial_derivative(0.01, 1)
    with pytest.raises(TypeError):
        heat_initial_derivative(Fraction(1, 100), True)
    with pytest.raises(ValueError):
        heat_initial_derivative(0, 1)
    with pytest.raises(ValueError):
        heat_initial_derivative(Fraction(1, 100), -1)
