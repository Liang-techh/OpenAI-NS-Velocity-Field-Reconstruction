from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_profile_stress_primitive_replay import (
    replay_profile_stress_primitives,
)


H = (Fraction(0), Fraction(1), Fraction(1, 3))
S_Q = (Fraction(2, 5), Fraction(-3, 7), Fraction(5, 11))
S_N = (Fraction(-1, 4), Fraction(2, 9), Fraction(7, 13))
X = Fraction(3, 5)


def test_regular_primitives_satisfy_both_odes_exactly():
    replay = replay_profile_stress_primitives(X, H, S_Q, S_N)
    assert replay.q_ode_residual == 0
    assert replay.n_ode_residual == 0
    assert replay.tolerance == 0
    assert replay.passed


def test_axis_limits_match_public_regular_extension():
    replay = replay_profile_stress_primitives(X, H, S_Q, S_N)
    assert replay.q_axis_limit == Fraction(1, 5)
    assert replay.q_axis_expected == Fraction(1, 5)
    assert replay.n_axis_limit == Fraction(-1, 4)
    assert replay.n_axis_expected == Fraction(-1, 4)
    assert replay.axis_limits_exact


def test_nonzero_q_integration_constant_is_rejected_even_though_ode_still_holds():
    replay = replay_profile_stress_primitives(
        X, H, S_Q, S_N, q_integration_constant=Fraction(1, 2**40)
    )
    assert replay.q_ode_residual == 0
    assert not replay.regular_axis_constants
    assert not replay.passed


def test_nonzero_n_integration_constant_is_rejected_even_though_ode_still_holds():
    replay = replay_profile_stress_primitives(
        X, H, S_Q, S_N, n_integration_constant=Fraction(1, 2**40)
    )
    assert replay.n_ode_residual == 0
    assert not replay.regular_axis_constants
    assert not replay.passed


def test_float_and_non_axis_profile_inputs_fail_closed():
    with pytest.raises(TypeError):
        replay_profile_stress_primitives(0.6, H, S_Q, S_N)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        replay_profile_stress_primitives(X, (Fraction(0), 1.0), S_Q, S_N)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        replay_profile_stress_primitives(X, (Fraction(1), Fraction(1)), S_Q, S_N)


def test_fixture_values_are_pinned_exactly():
    replay = replay_profile_stress_primitives(X, H, S_Q, S_N)
    assert replay.h_value == Fraction(18, 25)
    assert replay.h_prime == Fraction(7, 5)
    assert replay.l_value == Fraction(7, 6)
    assert replay.s_q == Fraction(118, 385)
    assert replay.s_n == Fraction(301, 3900)
