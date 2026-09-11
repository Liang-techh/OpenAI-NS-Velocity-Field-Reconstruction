from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_residual_decay import (
    first_stage_for_increment_power,
    first_stage_for_residual_power,
    section9_increment_derivative_loss,
    section9_increment_gain,
    section9_residual_gain,
    section9_residual_power_record,
    section9_sigma,
)


def test_section9_schedule_is_exact_affine_schedule():
    assert [section9_sigma(j) for j in range(4)] == [
        Fraction(1, 5),
        Fraction(3, 10),
        Fraction(2, 5),
        Fraction(1, 2),
    ]
    assert all(section9_sigma(j + 1) - section9_sigma(j) == Fraction(1, 10)
               for j in range(20))


def test_eq_917_loss_and_first_power_threshold_are_exact():
    # Independent substitution into ell_m = 2A + (m+1)(1+3h/2).
    h = Fraction(1, 200)
    m = 2
    expected_loss = Fraction(1613, 400)
    assert section9_increment_derivative_loss(m, h) == expected_loss

    # g_j = h*j/10.  Requiring g_j-ell_m >= 1 gives j >= 10065.
    j = first_stage_for_increment_power(h, m, 1)
    assert j == 10065
    assert section9_increment_gain(j, h) - expected_loss == 1
    assert section9_increment_gain(j - 1, h) - expected_loss < 1


def test_eq_918_first_residual_power_threshold_is_minimal():
    # Independent closed-form check of h*(1/5+j/10)-K >= N.
    h = Fraction(1, 200)
    K = Fraction(3, 2)
    N = 2
    j = first_stage_for_residual_power(h, K, N)
    assert j == 6998

    exponent = h * (Fraction(1, 5) + Fraction(j, 10)) - K
    previous = h * (Fraction(1, 5) + Fraction(j - 1, 10)) - K
    assert exponent == N
    assert previous < N
    assert section9_residual_gain(j, h) - K == N


def test_residual_record_stays_explicitly_fail_closed():
    record = section9_residual_power_record(
        derivative_order=4,
        h=Fraction(1, 200),
        derivative_loss=Fraction(3, 2),
        target_power=2,
    )
    assert record.first_stage == 6998
    assert record.leading_exponent == 2
    assert record.status == "formal-structure"
    assert record.actual_residual_bound_verified is False
    assert record.flat_remainder_verified is False
    assert record.paper_exact_velocity_available is False


def test_section9_ledger_rejects_invalid_or_nonpaper_inputs():
    with pytest.raises(ValueError):
        section9_sigma(-1)
    with pytest.raises(ValueError):
        section9_sigma(True)
    with pytest.raises(ValueError):
        section9_increment_gain(0, 0)
    with pytest.raises(ValueError):
        section9_increment_gain(0, Fraction(1, 2))
    with pytest.raises(ValueError):
        section9_increment_derivative_loss(-1, Fraction(1, 200))
    with pytest.raises(ValueError):
        first_stage_for_residual_power(Fraction(1, 200), -1, 1)
    with pytest.raises(ValueError):
        first_stage_for_residual_power(Fraction(1, 200), 1, -1)
    with pytest.raises(ValueError):
        first_stage_for_increment_power(float("nan"), 0, 1)
