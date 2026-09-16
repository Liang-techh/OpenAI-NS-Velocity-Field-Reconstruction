from fractions import Fraction as F

import pytest

from openai_ns_reconstruction.kokuno_mean_pressure_moment_replay import (
    full_reconstruction,
    mc24_fifth_row_linear,
    mc29_full_jz,
    mc29_linear_jz,
    mc29_quadratic_defect,
    multiply,
    paper_exact_velocity_available,
    polynomial,
    pressure_source_linear_moment,
    radial_shift,
    integrate_unit_interval,
    verify_pressure_moment_sign,
)


G = polynomial([F(2, 3), F(-1, 5), F(3, 7)])
d = polynomial([F(-3, 8), F(4, 9)])
V_over_R = polynomial([F(5, 6), F(-2, 7), F(1, 4)])
u = polynomial([F(7, 10), F(1, 3)])
RG_NONLINEAR = polynomial([F(2, 11), F(-5, 13), F(1, 6)])


def test_pressure_source_linear_moment_is_exact_negative_rvu_term():
    source_value = pressure_source_linear_moment(V_over_R, u)
    direct = -integrate_unit_interval(radial_shift(multiply(V_over_R, u), 2))
    assert source_value == direct
    assert source_value == F(-3071, 12600)


def test_mc29_linear_jz_matches_mc24_fifth_row_exactly():
    lhs = mc29_linear_jz(G, d, V_over_R, u)
    rhs = mc24_fifth_row_linear(G, d, V_over_R, u)
    assert lhs == rhs == F(-66863, 226800)
    assert verify_pressure_moment_sign(G, d, V_over_R, u)


def test_wrong_pressure_source_sign_fails_closed():
    correct = mc29_linear_jz(G, d, V_over_R, u)
    pressure = pressure_source_linear_moment(V_over_R, u)
    wrong = correct - pressure - pressure
    expected = mc24_fifth_row_linear(G, d, V_over_R, u)
    assert wrong != expected
    assert wrong - expected == F(3071, 6300)


def test_quadratic_mc29_defect_is_retained_not_zeroed():
    defect = mc29_quadratic_defect(d, RG_NONLINEAR)
    assert defect == F(71843, 7413120)
    assert defect != 0
    assert mc29_full_jz(G, d, V_over_R, u, RG_NONLINEAR) - mc29_linear_jz(G, d, V_over_R, u) == defect


def test_tiny_sign_sensitive_mutation_is_exactly_detected():
    mutated = list(V_over_R)
    mutated[0] += F(1, 2**40)
    original = pressure_source_linear_moment(V_over_R, u)
    changed = pressure_source_linear_moment(mutated, u)
    assert changed != original
    assert changed - original == F(-19, 60 * 2**40)


def test_approximate_coefficients_are_rejected():
    with pytest.raises(TypeError):
        polynomial([0.5])
    with pytest.raises(TypeError):
        polynomial([True])


def test_truth_boundary_remains_fail_closed():
    assert paper_exact_velocity_available is False
    assert full_reconstruction is False
