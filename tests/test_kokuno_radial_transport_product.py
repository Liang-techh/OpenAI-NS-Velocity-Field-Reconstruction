from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_radial_transport_product import (
    FULL_RECONSTRUCTION,
    PAPER_EXACT_VELOCITY_AVAILABLE,
    corrected_polynomial_bound,
    exact_product_derivative,
    multi_binomial,
    transport_defect,
)


def test_multiindex_leibniz_rule_exact():
    I = (1, 1, 0, 0, 0)
    zero = (0, 0, 0, 0, 0)
    r = (1, 0, 0, 0, 0)
    z = (0, 1, 0, 0, 0)
    rz = (1, 1, 0, 0, 0)
    b = {zero: Fraction(2, 7), r: Fraction(-3, 11), z: Fraction(5, 13), rz: Fraction(7, 17)}
    n = {zero: Fraction(11, 19), r: Fraction(13, 23), z: Fraction(-17, 29), rz: Fraction(19, 31)}
    expected = b[zero] * n[rz] + b[r] * n[z] + b[z] * n[r] + b[rz] * n[zero]
    assert multi_binomial(I, r) == 1
    assert exact_product_derivative(I, b, n) == expected


def test_corrected_bound_retains_polynomial_s_star_factor():
    I = (1, 0, 0, 0, 0)
    zero = (0, 0, 0, 0, 0)
    r = I
    epsilon = Fraction(1, 1024)
    s_star = Fraction(16)
    cb = {zero: 1, r: 0}
    cn = {zero: 1, r: 1}
    degrees = {zero: 0, r: 1}

    exact, envelope, max_degree = corrected_polynomial_bound(
        I, epsilon, s_star, cb, cn, degrees
    )
    assert exact == Fraction(1, 64)
    assert envelope == Fraction(1, 64)
    assert max_degree == 1
    assert exact > epsilon


def test_max_degree_envelope_dominates_full_finite_sum():
    I = (2, 0, 0, 0, 0)
    zero = (0, 0, 0, 0, 0)
    one = (1, 0, 0, 0, 0)
    two = I
    exact, envelope, max_degree = corrected_polynomial_bound(
        I,
        Fraction(1, 128),
        9,
        {zero: 2, one: 3, two: 5},
        {zero: 7, one: 11, two: 13},
        {zero: 0, one: 1, two: 2},
    )
    assert exact == Fraction(1, 128) * (2 * 13 * 81 + 2 * 3 * 11 * 9 + 5 * 7)
    assert max_degree == 2
    assert envelope >= exact


def test_transport_defect_keeps_radial_product_term():
    eps = Fraction(1, 32)
    first = transport_defect(eps, 3, 5, 2, 1, Fraction(7, 64), Fraction(11, 13))
    without_radial = eps * 3 * (5 - 2)
    assert first - without_radial == Fraction(7, 64) * Fraction(11, 13)
    assert first != without_radial


def test_fail_closed_on_approximate_or_invalid_bound_inputs():
    I = (1, 0, 0, 0, 0)
    zero = (0, 0, 0, 0, 0)
    with pytest.raises(TypeError):
        corrected_polynomial_bound(I, 0.01, 4, {zero: 1, I: 1}, {zero: 1, I: 1}, {zero: 0, I: 1})
    with pytest.raises(ValueError):
        corrected_polynomial_bound(I, Fraction(1, 100), Fraction(1, 2), {zero: 1, I: 1}, {zero: 1, I: 1}, {zero: 0, I: 1})


def test_truth_boundary_stays_false():
    assert PAPER_EXACT_VELOCITY_AVAILABLE is False
    assert FULL_RECONSTRUCTION is False
