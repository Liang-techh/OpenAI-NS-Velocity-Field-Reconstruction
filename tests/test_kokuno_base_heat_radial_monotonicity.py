from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_base_heat_radial_monotonicity import (
    combined_factor_direct,
    combined_factor_identity_residual,
    combined_factor_positive_form,
    heat_similarity_exponent,
    local_log_slope,
    radial_sign_certificate,
)


def test_similarity_exponent_is_exact() -> None:
    assert heat_similarity_exponent(Fraction(1, 113)) == Fraction(115, 226)


def test_public_combined_factor_reduces_to_manifestly_positive_form() -> None:
    h = Fraction(1, 113)
    z = Fraction(7, 5)
    v = Fraction(11, 3)
    assert combined_factor_identity_residual(h, z, v) == 0
    assert combined_factor_direct(h, z, v) == combined_factor_positive_form(h, z, v)


def test_pointwise_factor_has_exact_half_margin() -> None:
    h = Fraction(1, 113)
    z = Fraction(7, 5)
    v = Fraction(11, 3)
    cert = radial_sign_certificate(h, z, v)
    # Zv = 77/15, hence 1+Zv = 92/15.
    assert cert.half_margin == Fraction(15, 10396)
    assert cert.factor == Fraction(5213, 10396)
    assert cert.outward_bracket == Fraction(-5213, 10396)


def test_zero_similarity_or_kernel_node_keeps_strict_margin() -> None:
    h = Fraction(1, 200)
    for z, v in ((Fraction(0), Fraction(9, 4)), (Fraction(5, 7), Fraction(0))):
        cert = radial_sign_certificate(h, z, v)
        assert cert.factor == Fraction(101, 200)
        assert cert.half_margin == h
        assert cert.outward_bracket < 0


def test_local_log_slope_is_strictly_below_h_for_finite_exact_inputs() -> None:
    h = Fraction(1, 113)
    slope = local_log_slope(h, Fraction(10**6), Fraction(10**6))
    assert 0 <= slope < h
    assert h - slope == Fraction(1, 113 * (10**12 + 1))


def test_one_bit_mutation_of_A_is_detected_exactly() -> None:
    h = Fraction(1, 113)
    z = Fraction(7, 5)
    v = Fraction(11, 3)
    mutation = Fraction(1, 2**40)
    mutated_direct = combined_factor_direct(h, z, v) + mutation
    assert mutated_direct - combined_factor_positive_form(h, z, v) == mutation


def test_invalid_or_inexact_inputs_fail_closed() -> None:
    with pytest.raises(TypeError):
        radial_sign_certificate(0.005, Fraction(1), Fraction(1))
    with pytest.raises(TypeError):
        radial_sign_certificate(True, Fraction(1), Fraction(1))
    with pytest.raises(ValueError):
        radial_sign_certificate(0, Fraction(1), Fraction(1))
    with pytest.raises(ValueError):
        radial_sign_certificate(Fraction(1, 113), Fraction(-1), Fraction(1))
    with pytest.raises(ValueError):
        radial_sign_certificate(Fraction(1, 113), Fraction(1), Fraction(-1))
