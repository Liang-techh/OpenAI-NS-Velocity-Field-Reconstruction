from decimal import Decimal
from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_base_heat_radial_replay import (
    FULL_RECONSTRUCTION,
    PAPER_EXACT_VELOCITY_AVAILABLE,
    HeatJet,
    cylindrical_swirl_core,
    heat_equation_residual,
    heat_ode_residual,
    replay_identity,
    time_core,
)


def fixture() -> HeatJet:
    return HeatJet.exact(
        Fraction(1, 8),
        Fraction(2, 3),
        Fraction(5, 7),
        Fraction(-3, 11),
        Fraction(25785, 19712),
    )


def test_exact_public_heat_operator_bridge():
    jet = fixture()
    result = replay_identity(jet)
    assert heat_ode_residual(jet) == 0
    assert time_core(jet) == Fraction(6, 11)
    assert cylindrical_swirl_core(jet) == Fraction(6, 11)
    assert heat_equation_residual(jet) == 0
    assert result["bridge_residual"] == 0
    assert result["passes"] is True


def test_bridge_is_twice_ode_residual_for_arbitrary_exact_jet():
    jet = HeatJet.exact(
        Fraction(3, 20),
        Fraction(7, 5),
        Fraction(4, 9),
        Fraction(-2, 13),
        Fraction(11, 17),
    )
    assert heat_equation_residual(jet) == 2 * heat_ode_residual(jet)


def test_h2_mutation_fails_closed_at_zero_tolerance():
    jet = fixture()
    eps = Fraction(1, 2**40)
    mutated = HeatJet.exact(jet.h, jet.Z, jet.H, jet.H1, jet.H2 + eps)
    expected = Fraction(1, 9 * 2**37)
    assert heat_equation_residual(mutated) == expected
    assert replay_identity(mutated)["passes"] is False


def test_approximate_inputs_are_rejected():
    with pytest.raises(TypeError):
        HeatJet.exact(0.125, Fraction(2, 3), 1, 1, 1)
    with pytest.raises(TypeError):
        HeatJet.exact(Fraction(1, 8), Decimal("0.666"), 1, 1, 1)


def test_truth_boundary_stays_false():
    assert PAPER_EXACT_VELOCITY_AVAILABLE is False
    assert FULL_RECONSTRUCTION is False
