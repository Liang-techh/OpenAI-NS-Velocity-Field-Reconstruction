from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_osc_production_margin import (
    replay_oscillatory_production_budget,
    require_half_leading_production,
)


def fixture(ell: int = 500):
    return dict(
        h=Fraction(1, 200),
        ell=ell,
        g_plus=Fraction(2),
        g_minus=Fraction(1),
        c_minus=Fraction(1, 2),
        c_t=Fraction(3),
        c_b=Fraction(1),
        l_g=Fraction(1),
        c_box=Fraction(1),
        s_max=Fraction(2),
        gamma_max=Fraction(2),
        c_x=Fraction(3, 4),
    )


def test_exact_corrected_budget_at_passing_dyadic_stage():
    out = require_half_leading_production(**fixture(500))
    assert out.dyadic_exponent == 5
    assert out.borel_actual_shear_error == Fraction(1, 32)
    assert out.frozen_frame_error == Fraction(3, 125000)
    assert out.box_actual_shear_error == Fraction(1, 15625000000000000)
    assert out.stress_factor == Fraction(1000003, 250000)
    assert out.total_error == Fraction(488376464843751000003, 3906250000000000000000)
    assert out.half_leading_budget == Fraction(1, 4)
    assert out.half_leading_margin == Fraction(488186035156248999997, 3906250000000000000000)
    assert out.lower_production_coefficient == Fraction(9, 64)
    assert out.certifies_half_leading_bound


def test_corrected_actual_shear_term_prevents_false_positive():
    out = replay_oscillatory_production_budget(**fixture(400))
    assert out.dyadic_exponent == 4
    assert out.borel_actual_shear_error == Fraction(1, 16)
    assert out.half_leading_margin == Fraction(
        -25344000000640003, 655360000000000000000
    )
    assert not out.certifies_half_leading_bound

    # If the Delta-g contribution were incorrectly dropped, the frozen-frame
    # term alone would fit easily inside the same half-leading budget.
    leading_only_margin = out.half_leading_budget - out.frozen_frame_error
    assert leading_only_margin > 0
    with pytest.raises(ValueError, match="actual-shear error"):
        require_half_leading_production(**fixture(400))


def test_single_bit_increase_in_borel_constant_changes_margin_exactly():
    base = replay_oscillatory_production_budget(**fixture(500))
    mutated_args = fixture(500)
    mutated_args["c_b"] += Fraction(1, 2**40)
    mutated = replay_oscillatory_production_budget(**mutated_args)
    expected_drift = base.stress_factor * Fraction(1, 2**45)
    assert base.half_leading_margin - mutated.half_leading_margin == expected_drift
    assert expected_drift > 0


def test_non_integral_dyadic_exponent_fails_closed():
    args = fixture(501)
    with pytest.raises(ValueError, match=r"integral 2\*h\*ell"):
        replay_oscillatory_production_budget(**args)


def test_source_h_boundary_fails_closed():
    args = fixture()
    args["h"] = Fraction(1, 100)
    with pytest.raises(ValueError, match="0 < h < 1/100"):
        replay_oscillatory_production_budget(**args)


def test_float_theorem_input_fails_closed():
    args = fixture()
    args["g_plus"] = 2.0
    with pytest.raises(TypeError, match="fractions.Fraction"):
        replay_oscillatory_production_budget(**args)


def test_nonpositive_lower_bound_inputs_fail_closed():
    args = fixture()
    args["g_minus"] = Fraction(0)
    with pytest.raises(ValueError, match="g_minus must be positive"):
        replay_oscillatory_production_budget(**args)
