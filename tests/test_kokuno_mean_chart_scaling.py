from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_mean_chart_scaling import (
    mean_chart_scaling,
    source_mean_moment_exponents,
    verify_mc27_ratio_exponents,
    verify_negative_dz_potential,
)


def test_mc15_generic_weighted_moment_and_primitive_exponents():
    s = mean_chart_scaling(h=Fraction(1, 200), moment_a=Fraction(7, 5), moment_e=2)
    assert s.weighted_moment_exponent == Fraction(-1, 10)
    assert s.radial_primitive_exponent == Fraction(9, 10)


def test_source_mean_moment_exponents_are_exact():
    m_theta, m_z = source_mean_moment_exponents(Fraction(1, 200))
    assert m_theta == Fraction(-199, 200)
    assert m_z == Fraction(-99, 200)


def test_field_normalization_exponents_and_negative_dz_sign():
    s = mean_chart_scaling(h=Fraction(1, 200), moment_a=Fraction(1), moment_e=0)
    assert s.A == Fraction(101, 200)
    assert s.D == Fraction(99, 200)
    assert s.velocity_exponent == Fraction(101, 200)
    assert s.residual_exponent == Fraction(151, 100)
    assert s.pressure_exponent == Fraction(101, 100)
    assert s.potential_exponent == Fraction(1, 200)
    assert s.normalized_negative_dz_potential_exponent == Fraction(1, 200)
    assert s.normalized_negative_dz_potential_sign == -1


def test_wrong_sign_and_2pow40_exponent_mutation_fail_closed():
    h = Fraction(1, 200)
    assert verify_negative_dz_potential(
        h=h, reported_sign=-1, reported_exponent=Fraction(1, 200)
    )
    assert not verify_negative_dz_potential(
        h=h, reported_sign=1, reported_exponent=Fraction(1, 200)
    )
    assert not verify_negative_dz_potential(
        h=h,
        reported_sign=-1,
        reported_exponent=Fraction(1, 200) + Fraction(1, 2**40),
    )


def test_mc27_fixed_q_ratio_exponents():
    h = Fraction(1, 200)
    assert verify_mc27_ratio_exponents(
        h=h, amplitude_exponent=Fraction(101, 200), argument_exponent=Fraction(1, 2)
    )
    assert not verify_mc27_ratio_exponents(
        h=h, amplitude_exponent=Fraction(101, 200), argument_exponent=-Fraction(1, 2)
    )


def test_rejects_approximate_and_out_of_domain_inputs():
    with pytest.raises(TypeError):
        mean_chart_scaling(h=0.005, moment_a=Fraction(1), moment_e=0)
    with pytest.raises(TypeError):
        mean_chart_scaling(h=Fraction(1, 200), moment_a=1.0, moment_e=0)
    with pytest.raises(TypeError):
        mean_chart_scaling(h=Fraction(1, 200), moment_a=Fraction(1), moment_e=True)
    with pytest.raises(ValueError):
        mean_chart_scaling(h=Fraction(1, 100), moment_a=Fraction(1), moment_e=0)
    with pytest.raises(ValueError):
        mean_chart_scaling(h=Fraction(1, 200), moment_a=Fraction(1), moment_e=-1)


def test_exact_identities_have_zero_symbolic_residual():
    h = Fraction(1, 200)
    s = mean_chart_scaling(h=h, moment_a=Fraction(7, 5), moment_e=2)
    assert s.normalized_negative_dz_potential_exponent - h == 0
    assert s.potential_exponent - h == 0
    assert s.mc27_amplitude_ratio_exponent - s.A == 0
    assert s.mc27_argument_ratio_exponent - Fraction(1, 2) == 0
