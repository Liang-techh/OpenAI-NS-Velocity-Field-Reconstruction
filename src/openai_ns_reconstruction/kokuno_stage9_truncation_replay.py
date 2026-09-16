"""Clean-room exact replay of the public finite Stage-9 truncation choice.

This module independently implements only the displayed exponent arithmetic.
It does not encode or certify the analytic estimates that supply the losses.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


_ONE_HUNDREDTH = Fraction(1, 100)
_ONE_FIFTH = Fraction(1, 5)
_ONE_TENTH = Fraction(1, 10)


def _require_fraction(name: str, value: object) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be fractions.Fraction")
    return value


def _require_nonnegative_int(name: str, value: object) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise TypeError(f"{name} must be a nonnegative int")
    return value


@dataclass(frozen=True)
class Stage9TruncationData:
    """Exact data used by the finite truncation inequalities."""

    h: Fraction
    derivative_order: int
    target_power: int
    phase_loss_m_plus_2: Fraction
    residual_loss_m_plus_2: Fraction
    residual_loss_m: Fraction

    def __post_init__(self) -> None:
        h = _require_fraction("h", self.h)
        _require_nonnegative_int("derivative_order", self.derivative_order)
        _require_nonnegative_int("target_power", self.target_power)
        losses = (
            ("phase_loss_m_plus_2", self.phase_loss_m_plus_2),
            ("residual_loss_m_plus_2", self.residual_loss_m_plus_2),
            ("residual_loss_m", self.residual_loss_m),
        )
        if not (Fraction(0) < h < _ONE_HUNDREDTH):
            raise ValueError("h must satisfy the construction range 0 < h < 1/100")
        for name, value in losses:
            value = _require_fraction(name, value)
            if value < 0:
                raise ValueError(f"{name} must be nonnegative")


def sigma(stage: int) -> Fraction:
    stage = _require_nonnegative_int("stage", stage)
    return _ONE_FIFTH + stage * _ONE_TENTH


def gain(data: Stage9TruncationData, stage: int) -> Fraction:
    stage = _require_nonnegative_int("stage", stage)
    return data.h * stage * _ONE_TENTH


def truncation_margins(data: Stage9TruncationData, stage: int) -> tuple[Fraction, Fraction]:
    """Return the two strict source margins; both must be positive."""

    stage = _require_nonnegative_int("stage", stage)
    first = (
        gain(data, stage + 1) / 2
        - data.phase_loss_m_plus_2
        - (
            data.target_power
            + data.residual_loss_m_plus_2
            + 2
        )
    )
    second = (
        data.h * sigma(stage)
        - data.residual_loss_m
        - (data.target_power + 1)
    )
    return first, second


def stage_is_admissible(data: Stage9TruncationData, stage: int) -> bool:
    if stage < data.derivative_order + 2:
        return False
    first, second = truncation_margins(data, stage)
    return first > 0 and second > 0


def minimal_truncation_stage(
    data: Stage9TruncationData,
    *,
    max_stage: int,
) -> int:
    """Find the first finite stage satisfying both strict inequalities.

    ``max_stage`` is an explicit computational cap; exhaustion is an error rather
    than evidence that no mathematical stage exists.
    """

    max_stage = _require_nonnegative_int("max_stage", max_stage)
    start = data.derivative_order + 2
    if max_stage < start:
        raise ValueError("max_stage lies below the mandatory J >= m+2 boundary")
    for stage in range(start, max_stage + 1):
        if stage_is_admissible(data, stage):
            return stage
    raise RuntimeError("finite search cap exhausted before both strict margins became positive")
