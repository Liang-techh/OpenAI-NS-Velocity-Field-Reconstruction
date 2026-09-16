"""Shared exact-rational interval primitives for certified enclosures.

This module is intentionally small. It provides the arithmetic needed by the
outgoing-schedule enclosure code without attaching theorem status to the
result. All endpoints are :class:`fractions.Fraction` values and all rounding
is outward.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction


def positive_fraction(value: object, name: str) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be a Fraction")
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def positive_cap(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def dyadic_step(bound: Fraction) -> Fraction:
    """Return the largest dyadic step no larger than positive ``bound``."""

    if not isinstance(bound, Fraction) or bound <= 0:
        raise ValueError("dyadic step bound must be positive")
    exponent = max(0, bound.denominator.bit_length() - bound.numerator.bit_length())
    step = Fraction(1, 1 << exponent)
    while step > bound:
        exponent += 1
        step = Fraction(1, 1 << exponent)
    while exponent > 0 and Fraction(1, 1 << (exponent - 1)) <= bound:
        exponent -= 1
        step = Fraction(1, 1 << exponent)
    return step


def floor_grid(value: Fraction, step: Fraction) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError("value must be a Fraction")
    step = positive_fraction(step, "step")
    return (value // step) * step


def ceil_grid(value: Fraction, step: Fraction) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError("value must be a Fraction")
    step = positive_fraction(step, "step")
    quotient = value / step
    return (-((-quotient.numerator) // quotient.denominator)) * step


@dataclass(frozen=True)
class RationalInterval:
    """An ordered interval whose endpoints are exact rational numbers."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if not isinstance(self.lower, Fraction):
            raise TypeError("lower must be a Fraction")
        if not isinstance(self.upper, Fraction):
            raise TypeError("upper must be a Fraction")
        if self.lower > self.upper:
            raise ValueError("interval bounds must be ordered")

    @property
    def width(self) -> Fraction:
        return self.upper - self.lower


def interval_add(left: RationalInterval, right: RationalInterval) -> RationalInterval:
    return RationalInterval(left.lower + right.lower, left.upper + right.upper)


def interval_subtract(left: RationalInterval, right: RationalInterval) -> RationalInterval:
    return RationalInterval(left.lower - right.upper, left.upper - right.lower)


def interval_multiply(left: RationalInterval, right: RationalInterval) -> RationalInterval:
    products = (
        left.lower * right.lower,
        left.lower * right.upper,
        left.upper * right.lower,
        left.upper * right.upper,
    )
    return RationalInterval(min(products), max(products))


def interval_divide_positive(
    numerator: RationalInterval,
    denominator: RationalInterval,
) -> RationalInterval:
    """Divide by an interval proved strictly positive."""

    if denominator.lower <= 0:
        raise ValueError("denominator interval must be strictly positive")
    reciprocal = RationalInterval(
        Fraction(1, 1) / denominator.upper,
        Fraction(1, 1) / denominator.lower,
    )
    return interval_multiply(numerator, reciprocal)


def map_monotone(
    interval: RationalInterval,
    function: Callable[[Fraction], Fraction],
    *,
    increasing: bool = True,
) -> RationalInterval:
    """Propagate exact endpoints through a monotone rational-valued function."""

    if not callable(function):
        raise TypeError("function must be callable")
    first = function(interval.lower)
    second = function(interval.upper)
    if not isinstance(first, Fraction) or not isinstance(second, Fraction):
        raise TypeError("monotone function must return Fractions")
    lower, upper = (first, second) if increasing else (second, first)
    if lower > upper:
        raise ValueError("monotonicity direction is inconsistent with endpoint values")
    return RationalInterval(lower, upper)


def outward_round_dyadic(
    interval: RationalInterval,
    step: Fraction,
) -> RationalInterval:
    """Round both endpoints outward to an exact positive dyadic grid step."""

    step = positive_fraction(step, "step")
    return RationalInterval(
        floor_grid(interval.lower, step),
        ceil_grid(interval.upper, step),
    )


__all__ = [
    "RationalInterval",
    "ceil_grid",
    "dyadic_step",
    "floor_grid",
    "interval_add",
    "interval_divide_positive",
    "interval_multiply",
    "interval_subtract",
    "map_monotone",
    "outward_round_dyadic",
    "positive_cap",
    "positive_fraction",
]
