from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_stage9_physical_loss_replay import (
    Stage9ReplayMismatch,
    replay_stage9_physical_loss,
    underestimated_time_loss_defect,
)


def test_exact_source_fixture() -> None:
    r = replay_stage9_physical_loss(Fraction(1, 200), 3, 17)
    assert r.A == Fraction(101, 200)
    assert r.physical_residual_loss == Fraction(151, 100)
    assert r.spatial_derivative_loss == Fraction(201, 400)
    assert r.temporal_derivative_loss == Fraction(403, 400)
    assert r.fixed_loss_Km == Fraction(1813, 400)
    assert r.sigma_j == Fraction(19, 10)
    assert r.stage_gain == Fraction(19, 2000)
    assert r.next_stage_gain == Fraction(1, 100)
    assert r.gain_increment == Fraction(1, 2000)
    assert r.net_power == Fraction(-4523, 1000)


def test_fixed_loss_is_independent_of_stage() -> None:
    r0 = replay_stage9_physical_loss(Fraction(1, 200), 4, 0)
    r1 = replay_stage9_physical_loss(Fraction(1, 200), 4, 12345)
    assert r0.fixed_loss_Km == r1.fixed_loss_Km
    assert r1.stage_gain - r0.stage_gain == Fraction(12345, 2000)


def test_gain_progression_is_exact_h_over_ten() -> None:
    r = replay_stage9_physical_loss(Fraction(3, 1000), 2, 91)
    assert r.gain_increment == Fraction(3, 10000)


def test_perturbed_A_fails_closed() -> None:
    h = Fraction(1, 200)
    expected_A = Fraction(1, 2) + h
    delta = Fraction(1, 2**40)
    with pytest.raises(Stage9ReplayMismatch, match="A mismatch"):
        replay_stage9_physical_loss(h, 3, 17, claimed_A=expected_A + delta)
    # If the bad A were used in Q^(-2A-1/2), its loss would drift by 2*delta.
    assert 2 * delta == Fraction(1, 2**39)


def test_underestimating_time_derivative_loss_is_detected_exactly() -> None:
    h = Fraction(1, 200)
    claimed = Fraction(1) + h
    expected = Fraction(1) + 3 * h / 2
    assert expected - claimed == Fraction(1, 400)
    assert underestimated_time_loss_defect(h, 3) == Fraction(3, 400)
    with pytest.raises(Stage9ReplayMismatch, match="temporal derivative loss mismatch"):
        replay_stage9_physical_loss(h, 3, 17, claimed_temporal_loss=claimed)


def test_exact_claims_are_accepted() -> None:
    h = Fraction(1, 250)
    replay_stage9_physical_loss(
        h,
        5,
        7,
        claimed_A=Fraction(1, 2) + h,
        claimed_spatial_loss=Fraction(1, 2) + h / 2,
        claimed_temporal_loss=Fraction(1) + 3 * h / 2,
    )


def test_approximate_and_out_of_domain_inputs_fail_closed() -> None:
    with pytest.raises(TypeError):
        replay_stage9_physical_loss(0.005, 3, 17)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        replay_stage9_physical_loss(Fraction(1, 100), 3, 17)
    with pytest.raises(TypeError):
        replay_stage9_physical_loss(Fraction(1, 200), True, 17)  # type: ignore[arg-type]
