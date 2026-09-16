"""Exact rational enclosures for the pinned outgoing smooth step.

The schedule uses

    ``sigma(x) = edge(x) / (edge(x) + edge(1 - x))``

with ``edge(x) = exp(-1/x**2)`` for positive ``x``.  This module evaluates
that particular cutoff and its primitive with rational interval arithmetic.
It is deliberately separate from :mod:`cutoffs`: ``cutoffs.smooth_step``
uses the different exponent ``exp(-1/x)`` and is not the schedule function.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


_DEFAULT_SIGMA_TOLERANCE = Fraction(1, 1024)
_DEFAULT_MAX_TERMS = 1024
_DEFAULT_MAX_SQUARINGS = 4096
_DEFAULT_MAX_CELLS = 4096


def _positive_fraction(value: object, name: str) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be a Fraction")
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _positive_cap(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _dyadic_step(bound: Fraction) -> Fraction:
    """Return a dyadic step no larger than the positive rational ``bound``."""

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


def _floor_grid(value: Fraction, step: Fraction) -> Fraction:
    return (value // step) * step


def _ceil_grid(value: Fraction, step: Fraction) -> Fraction:
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


def _binary_scale_to_unit(value: Fraction) -> tuple[int, Fraction]:
    """Return minimal ``s >= 0`` and ``y = value / 2**s <= 1``."""

    if value <= 0:
        return 0, value
    if value <= 1:
        return 0, value

    numerator = value.numerator
    denominator = value.denominator
    floor_log = numerator.bit_length() - denominator.bit_length()
    if floor_log < 0:
        floor_log = 0
    elif value < Fraction(1 << floor_log):
        floor_log -= 1
    exact_power = value == Fraction(1 << floor_log)
    squarings = floor_log if exact_power else floor_log + 1
    return squarings, value / (1 << squarings)


def _exp_positive_sum(y: Fraction, degree: int) -> Fraction:
    """Return ``sum(y**n/n!, n=0..degree)`` exactly."""

    total = Fraction(1)
    term = Fraction(1)
    for n in range(1, degree + 1):
        term = term * y / n
        total += term
    return total


def _exp_positive_tail(y: Fraction, degree: int, next_term: Fraction) -> Fraction:
    """Bound the positive exponential tail after ``degree`` exactly."""

    ratio = y / (degree + 2)
    if ratio >= 1:
        raise ArithmeticError("exponential Taylor tail ratio is not below one")
    return next_term / (1 - ratio)


def validated_exp_negative(
    x: Fraction,
    *,
    absolute_tolerance: Fraction,
    max_terms: int = _DEFAULT_MAX_TERMS,
    max_squarings: int = _DEFAULT_MAX_SQUARINGS,
) -> RationalInterval:
    """Enclose ``exp(-x)`` for a nonnegative rational ``x``.

    The positive Taylor series is evaluated at ``y = x / 2**s <= 1`` and the
    resulting reciprocal interval is squared ``s`` times.  Every endpoint is
    rounded outward to one fixed dyadic grid, keeping the denominators bounded
    during repeated squaring.  Invalid controls are checked before the exact
    ``x == 0`` branch.
    """

    if not isinstance(x, Fraction):
        raise TypeError("x must be a Fraction")
    if x < 0:
        raise ValueError("x must be nonnegative")
    tolerance = _positive_fraction(absolute_tolerance, "absolute_tolerance")
    max_terms = _positive_cap(max_terms, "max_terms")
    max_squarings = _positive_cap(max_squarings, "max_squarings")

    epsilon = min(tolerance, Fraction(1))
    if x == 0:
        return RationalInterval(Fraction(1), Fraction(1))

    squarings, y = _binary_scale_to_unit(x)
    if squarings > max_squarings:
        raise ArithmeticError("exponential enclosure exceeded max_squarings")
    if not Fraction(0) < y <= Fraction(1):
        raise ArithmeticError("binary exponential reduction did not produce 0 < y <= 1")

    tail_target = epsilon / (4 * (1 << squarings))
    degree = 0
    term = Fraction(1)
    while degree < max_terms:
        next_term = term * y / (degree + 1)
        tail = _exp_positive_tail(y, degree, next_term)
        if tail <= tail_target:
            break
        term = next_term
        degree += 1
    else:
        raise ArithmeticError("exponential enclosure exceeded max_terms")

    finite_sum = _exp_positive_sum(y, degree)
    # exp(y) lies in [finite_sum, finite_sum + tail], so reciprocals reverse
    # the endpoints and enclose exp(-y).
    lower = Fraction(1, 1) / (finite_sum + tail)
    upper = Fraction(1, 1) / finite_sum

    delta = _dyadic_step(epsilon / (8 * (1 << squarings)))

    def rounded_interval(left: Fraction, right: Fraction) -> RationalInterval:
        left = max(Fraction(0), _floor_grid(left, delta))
        right = min(Fraction(1), _ceil_grid(right, delta))
        if left > right:
            raise ArithmeticError("rounded exponential interval became unordered")
        return RationalInterval(left, right)

    interval = rounded_interval(lower, upper)
    for _ in range(squarings):
        interval = rounded_interval(interval.lower * interval.lower, interval.upper * interval.upper)

    if interval.width > epsilon:
        raise ArithmeticError("exponential enclosure exceeded requested tolerance")
    return interval


def _validate_sigma_controls(
    absolute_tolerance: object,
    max_terms: object,
    max_squarings: object,
) -> tuple[Fraction, int, int]:
    tolerance = _positive_fraction(absolute_tolerance, "absolute_tolerance")
    terms = _positive_cap(max_terms, "max_terms")
    squarings = _positive_cap(max_squarings, "max_squarings")
    return tolerance, terms, squarings


def validated_outgoing_sigma(
    x: Fraction,
    *,
    absolute_tolerance: Fraction = _DEFAULT_SIGMA_TOLERANCE,
    max_terms: int = _DEFAULT_MAX_TERMS,
    max_squarings: int = _DEFAULT_MAX_SQUARINGS,
) -> RationalInterval:
    """Enclose the pinned outgoing ``sigma(x)`` using exact rationals."""

    if not isinstance(x, Fraction):
        raise TypeError("x must be a Fraction")
    tolerance, max_terms, max_squarings = _validate_sigma_controls(
        absolute_tolerance,
        max_terms,
        max_squarings,
    )

    if x <= 0:
        return RationalInterval(Fraction(0), Fraction(0))
    if x >= 1:
        return RationalInterval(Fraction(1), Fraction(1))
    if x == Fraction(1, 2):
        return RationalInterval(Fraction(1, 2), Fraction(1, 2))

    if x > Fraction(1, 2):
        reflected = validated_outgoing_sigma(
            1 - x,
            absolute_tolerance=tolerance,
            max_terms=max_terms,
            max_squarings=max_squarings,
        )
        return RationalInterval(
            1 - reflected.upper,
            1 - reflected.lower,
        )

    t = Fraction(1, x * x) - Fraction(1, (1 - x) * (1 - x))
    ratio = validated_exp_negative(
        t,
        absolute_tolerance=tolerance,
        max_terms=max_terms,
        max_squarings=max_squarings,
    )
    lower = ratio.lower / (1 + ratio.lower)
    upper = ratio.upper / (1 + ratio.upper)
    interval = RationalInterval(lower, upper)
    if interval.width > tolerance:
        raise ArithmeticError("outgoing sigma enclosure exceeded requested tolerance")
    return interval


def _round_sigma_interval(interval: RationalInterval, step: Fraction) -> RationalInterval:
    lower = max(Fraction(0), _floor_grid(interval.lower, step))
    upper = min(Fraction(1), _ceil_grid(interval.upper, step))
    if lower > upper:
        raise ArithmeticError("rounded sigma interval became unordered")
    return RationalInterval(lower, upper)


def validated_sigma_primitive(
    x: Fraction,
    *,
    absolute_tolerance: Fraction = _DEFAULT_SIGMA_TOLERANCE,
    max_cells: int = _DEFAULT_MAX_CELLS,
    max_terms: int = _DEFAULT_MAX_TERMS,
    max_squarings: int = _DEFAULT_MAX_SQUARINGS,
) -> RationalInterval:
    """Enclose ``integral_0^x outgoing_sigma(t) dt`` with Darboux sums."""

    if not isinstance(x, Fraction):
        raise TypeError("x must be a Fraction")
    tolerance, max_terms, max_squarings = _validate_sigma_controls(
        absolute_tolerance,
        max_terms,
        max_squarings,
    )
    max_cells = _positive_cap(max_cells, "max_cells")

    # Validate every control before exact outside-domain branches.
    epsilon = min(tolerance, Fraction(1))
    if x <= 0:
        return RationalInterval(Fraction(0), Fraction(0))
    if x >= 1:
        return RationalInterval(x - Fraction(1, 2), x - Fraction(1, 2))
    if x > Fraction(1, 2):
        reflected = validated_sigma_primitive(
            1 - x,
            absolute_tolerance=tolerance,
            max_cells=max_cells,
            max_terms=max_terms,
            max_squarings=max_squarings,
        )
        return RationalInterval(
            x - Fraction(1, 2) + reflected.lower,
            x - Fraction(1, 2) + reflected.upper,
        )

    endpoint_tolerance = epsilon / (16 * x)
    common_step = _dyadic_step(min(Fraction(1), endpoint_tolerance))
    endpoint_cache: dict[Fraction, RationalInterval] = {}

    def endpoint(value: Fraction) -> RationalInterval:
        if value not in endpoint_cache:
            raw = validated_outgoing_sigma(
                value,
                absolute_tolerance=endpoint_tolerance,
                max_terms=max_terms,
                max_squarings=max_squarings,
            )
            endpoint_cache[value] = _round_sigma_interval(raw, common_step)
        return endpoint_cache[value]

    cells = 1
    while cells <= max_cells:
        cell_width = x / cells
        lower = Fraction(0)
        upper = Fraction(0)
        for index in range(cells):
            left = cell_width * index
            right = left + cell_width
            left_value = endpoint(left)
            right_value = endpoint(right)
            lower += cell_width * left_value.lower
            upper += cell_width * right_value.upper
        interval = RationalInterval(lower, upper)
        if interval.width <= tolerance:
            return interval
        cells *= 2

    raise ArithmeticError("sigma primitive did not reach requested tolerance before max_cells")


__all__ = [
    "RationalInterval",
    "validated_exp_negative",
    "validated_outgoing_sigma",
    "validated_sigma_primitive",
]
