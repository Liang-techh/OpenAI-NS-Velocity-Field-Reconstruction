from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_base_heat_initial_radial import (
    base_swirl_radial_core,
    heat_r_radial_core,
    heat_s_radial_core,
    initial_radial_seam_residual,
)


def test_public_chain_rule_fixture_is_exact() -> None:
    A = Fraction(5, 8)
    Z = Fraction(1, 3)
    H = Fraction(5, 7)
    H_prime = Fraction(-3, 11)

    assert heat_s_radial_core(A, Z, H, H_prime) == Fraction(-219, 616)
    assert heat_r_radial_core(A, Z, H, H_prime) == Fraction(-219, 308)


def test_r_chain_rule_factor_is_exactly_two() -> None:
    args = (Fraction(5, 8), Fraction(1, 3), Fraction(5, 7), Fraction(-3, 11))
    assert heat_r_radial_core(*args) == 2 * heat_s_radial_core(*args)


def test_omitting_z_hprime_has_detectable_exact_defect() -> None:
    A = Fraction(5, 8)
    Z = Fraction(1, 3)
    H = Fraction(5, 7)
    H_prime = Fraction(-3, 11)

    correct = heat_r_radial_core(A, Z, H, H_prime)
    incorrectly_drop_z_hprime = -2 * A * H
    assert correct - incorrectly_drop_z_hprime == Fraction(2, 11)


def test_initial_heat_core_matches_base_swirl_homogeneity() -> None:
    h = Fraction(1, 8)
    A = Fraction(1, 2) + h

    assert heat_r_radial_core(A, 0, 1, 0) == Fraction(-5, 4)
    assert base_swirl_radial_core(h) == Fraction(-5, 4)
    assert initial_radial_seam_residual(h) == 0


def test_initial_normalization_mutation_fails_closed() -> None:
    h = Fraction(1, 8)
    mutated_h0 = 1 + Fraction(1, 2**40)
    assert initial_radial_seam_residual(h, mutated_h0) == Fraction(
        -5, 4398046511104
    )


def test_inexact_inputs_are_rejected() -> None:
    with pytest.raises(TypeError):
        heat_r_radial_core(Fraction(5, 8), 0.0, 1, 0)
    with pytest.raises(TypeError):
        initial_radial_seam_residual(True)
