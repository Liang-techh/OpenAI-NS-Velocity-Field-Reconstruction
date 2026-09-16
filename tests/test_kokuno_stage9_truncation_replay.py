from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_stage9_truncation_replay import (
    Stage9TruncationData,
    gain,
    minimal_truncation_stage,
    sigma,
    stage_is_admissible,
    truncation_margins,
)


def fixture() -> Stage9TruncationData:
    return Stage9TruncationData(
        h=Fraction(1, 200),
        derivative_order=1,
        target_power=2,
        phase_loss_m_plus_2=Fraction(1, 3),
        residual_loss_m_plus_2=Fraction(1, 4),
        residual_loss_m=Fraction(1, 5),
    )


def test_public_stage_sequences_are_exact() -> None:
    data = fixture()
    assert sigma(0) == Fraction(1, 5)
    assert sigma(37) == Fraction(39, 10)
    assert gain(data, 37) == Fraction(37, 2000)


def test_minimal_stage_satisfies_both_strict_inequalities() -> None:
    data = fixture()
    stage = minimal_truncation_stage(data, max_stage=20000)
    assert stage == 18333
    first, second = truncation_margins(data, stage)
    assert first == Fraction(1, 6000)
    assert second == Fraction(2387, 400)
    assert first > 0 and second > 0


def test_predecessor_is_not_admissible() -> None:
    data = fixture()
    stage = 18333
    first, second = truncation_margins(data, stage - 1)
    assert first == Fraction(-1, 12000)
    assert second > 0
    assert not stage_is_admissible(data, stage - 1)
    assert stage_is_admissible(data, stage)


def test_exact_equality_does_not_pass_a_strict_bound() -> None:
    h = Fraction(1, 200)
    j = 10000
    n = 0
    k2 = Fraction(0)
    phase_loss = h * Fraction(j + 1, 20) - (n + k2 + 2)
    data = Stage9TruncationData(
        h=h,
        derivative_order=0,
        target_power=n,
        phase_loss_m_plus_2=phase_loss,
        residual_loss_m_plus_2=k2,
        residual_loss_m=Fraction(0),
    )
    first, second = truncation_margins(data, j)
    assert first == 0
    assert second > 0
    assert not stage_is_admissible(data, j)


def test_search_cap_exhaustion_fails_closed() -> None:
    with pytest.raises(RuntimeError, match="search cap exhausted"):
        minimal_truncation_stage(fixture(), max_stage=18332)


def test_construction_range_and_exact_types_fail_closed() -> None:
    with pytest.raises(ValueError, match="0 < h < 1/100"):
        Stage9TruncationData(
            h=Fraction(1, 100),
            derivative_order=0,
            target_power=0,
            phase_loss_m_plus_2=Fraction(0),
            residual_loss_m_plus_2=Fraction(0),
            residual_loss_m=Fraction(0),
        )
    with pytest.raises(TypeError, match="Fraction"):
        Stage9TruncationData(  # type: ignore[arg-type]
            h=0.005,
            derivative_order=0,
            target_power=0,
            phase_loss_m_plus_2=Fraction(0),
            residual_loss_m_plus_2=Fraction(0),
            residual_loss_m=Fraction(0),
        )


def test_larger_target_power_cannot_reduce_selected_stage() -> None:
    low = fixture()
    high = Stage9TruncationData(
        h=low.h,
        derivative_order=low.derivative_order,
        target_power=low.target_power + 1,
        phase_loss_m_plus_2=low.phase_loss_m_plus_2,
        residual_loss_m_plus_2=low.residual_loss_m_plus_2,
        residual_loss_m=low.residual_loss_m,
    )
    j_low = minimal_truncation_stage(low, max_stage=30000)
    j_high = minimal_truncation_stage(high, max_stage=30000)
    assert j_high > j_low
