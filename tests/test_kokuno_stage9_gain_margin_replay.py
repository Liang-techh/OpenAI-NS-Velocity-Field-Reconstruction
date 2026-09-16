from decimal import Decimal
from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_stage9_gain_margin_replay import (
    RADIAL_GAIN_FLOOR,
    SOURCE_KAPPA_S,
    replay_stage9_gain_margins,
)


def test_stage_zero_exact_published_margins() -> None:
    replay = replay_stage9_gain_margins(0)
    assert replay.sigma_j == Fraction(1, 5)
    assert replay.b_j == Fraction(7, 10)
    assert replay.c_j_star == Fraction(6, 5)
    assert replay.h_j == Fraction(59_999, 50_000)
    assert replay.transport_radial_gain == Fraction(8_999, 50_000)
    assert replay.curl_radial_gain == Fraction(49_997, 100_000)
    assert replay.self_product_radial_gain == Fraction(19_997, 100_000)
    assert replay.mean_to_wave_gain == Fraction(24_999, 50_000)
    assert replay.auxiliary_gain == Fraction(39_999, 100_000)
    assert replay.pressure_gain == Fraction(24_999, 25_000)
    assert replay.high_gain == Fraction(22_499, 25_000)
    assert replay.tolerance == 0
    assert replay.passed


def test_stage_relations_hold_for_later_finite_stage() -> None:
    replay = replay_stage9_gain_margins(7)
    assert replay.sigma_j == Fraction(9, 10)
    assert replay.b_j == Fraction(7, 5)
    assert replay.c_j_star == Fraction(19, 10)
    assert replay.stage_step == Fraction(1, 10)
    assert replay.self_product_radial_gain == Fraction(89_997, 100_000)
    assert replay.stage_relations_exact
    assert replay.radial_margin_floor_holds
    assert replay.common_gain_holds
    assert replay.passed


def test_stage_zero_is_worst_self_product_margin() -> None:
    stage_zero = replay_stage9_gain_margins(0)
    stage_one = replay_stage9_gain_margins(1)
    assert stage_zero.self_product_radial_gain == Fraction(19_997, 100_000)
    assert stage_zero.self_product_radial_gain > RADIAL_GAIN_FLOOR
    assert stage_one.self_product_radial_gain > stage_zero.self_product_radial_gain


def test_tiny_kappa_mutation_fails_exact_replay() -> None:
    mutated = replay_stage9_gain_margins(
        0, kappa_s=SOURCE_KAPPA_S + Fraction(1, 2**40)
    )
    assert mutated.transport_radial_gain != Fraction(8_999, 50_000)
    assert mutated.published_constants_exact is False
    assert mutated.passed is False


def test_radial_floor_is_not_inferred_when_kappa_is_too_large() -> None:
    mutated = replay_stage9_gain_margins(0, kappa_s=Fraction(1, 100))
    assert mutated.transport_radial_gain == Fraction(4, 25)
    assert mutated.transport_radial_gain < RADIAL_GAIN_FLOOR
    assert mutated.radial_margin_floor_holds is False
    assert mutated.passed is False


@pytest.mark.parametrize("bad", [1e-5, Decimal("0.00001"), 1])
def test_approximate_or_untyped_kappa_rejected(bad: object) -> None:
    with pytest.raises(TypeError):
        replay_stage9_gain_margins(0, kappa_s=bad)  # type: ignore[arg-type]


def test_invalid_stage_rejected() -> None:
    with pytest.raises(TypeError):
        replay_stage9_gain_margins(-1)
    with pytest.raises(TypeError):
        replay_stage9_gain_margins(True)
