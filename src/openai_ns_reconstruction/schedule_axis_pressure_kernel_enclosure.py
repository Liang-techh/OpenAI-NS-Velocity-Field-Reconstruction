"""Strict real interval enclosures for the SchedulePressure kernel and finite eta jets.

This is the rigorous counterpart of
``schedule_axis_pressure_jets.pressure_kernel_normalized_taylor`` for the
paper/Lean kernel

    k_a(eta) = (1 + eta**2)**(-2*a),   0 <= a <= 1.

Only real exact-rational parameter intervals are accepted.  The implementation
never substitutes interval midpoints.  For a box in ``(a, eta)`` it encloses
the value by exact rational log/exp bounds and propagates every eta derivative
through the differential identity

    (1+eta^2) k' + 4 a eta k = 0.

Subdivisions refine the same shared ``a`` and ``eta`` boxes.  They do not turn
the result into an independence assumption.  Normalized Taylor coefficients
``k^(m)/m!`` and full derivatives ``k^(m)`` are exposed separately.

This module makes no complex-domain claim and does not certify the outer
SchedulePressure y-integral.  Those are separate obligations.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import factorial
from typing import Iterable


ExactScalar = int | Fraction


def _fraction(value: ExactScalar, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an int or Fraction")
    return Fraction(value)


def _nonnegative_index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _positive_index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


@dataclass(frozen=True)
class RationalInterval:
    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        lower = _fraction(self.lower, "lower")
        upper = _fraction(self.upper, "upper")
        if lower > upper:
            raise ValueError("interval lower endpoint exceeds upper endpoint")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)

    @classmethod
    def point(cls, value: ExactScalar) -> "RationalInterval":
        q = _fraction(value, "value")
        return cls(q, q)

    @property
    def width(self) -> Fraction:
        return self.upper - self.lower

    @property
    def is_point(self) -> bool:
        return self.lower == self.upper

    def contains(self, value: ExactScalar) -> bool:
        q = _fraction(value, "value")
        return self.lower <= q <= self.upper

    def hull(self, other: "RationalInterval") -> "RationalInterval":
        if not isinstance(other, RationalInterval):
            raise TypeError("other must be RationalInterval")
        return RationalInterval(min(self.lower, other.lower), max(self.upper, other.upper))

    def __neg__(self) -> "RationalInterval":
        return RationalInterval(-self.upper, -self.lower)

    def __add__(self, other: "RationalInterval") -> "RationalInterval":
        if not isinstance(other, RationalInterval):
            return NotImplemented
        return RationalInterval(self.lower + other.lower, self.upper + other.upper)

    def __sub__(self, other: "RationalInterval") -> "RationalInterval":
        if not isinstance(other, RationalInterval):
            return NotImplemented
        return self + (-other)

    def __mul__(self, other: "RationalInterval") -> "RationalInterval":
        if not isinstance(other, RationalInterval):
            return NotImplemented
        products = (
            self.lower * other.lower,
            self.lower * other.upper,
            self.upper * other.lower,
            self.upper * other.upper,
        )
        return RationalInterval(min(products), max(products))

    def reciprocal(self) -> "RationalInterval":
        if self.lower <= 0 <= self.upper:
            raise ZeroDivisionError("interval reciprocal crosses zero")
        return RationalInterval(1 / self.upper, 1 / self.lower)

    def __truediv__(self, other: "RationalInterval") -> "RationalInterval":
        if not isinstance(other, RationalInterval):
            return NotImplemented
        return self * other.reciprocal()

    def scale(self, value: ExactScalar) -> "RationalInterval":
        return self * RationalInterval.point(value)


def _coerce_interval(value: ExactScalar | RationalInterval, name: str) -> RationalInterval:
    if isinstance(value, RationalInterval):
        return value
    return RationalInterval.point(_fraction(value, name))


def _square_interval(x: RationalInterval) -> RationalInterval:
    if x.lower <= 0 <= x.upper:
        lower = Fraction(0)
    else:
        lower = min(x.lower * x.lower, x.upper * x.upper)
    upper = max(x.lower * x.lower, x.upper * x.upper)
    return RationalInterval(lower, upper)


def _split_interval(x: RationalInterval, cells: int) -> tuple[RationalInterval, ...]:
    cells = _positive_index(cells, "cells")
    if x.is_point:
        return (x,)
    step = x.width / cells
    return tuple(
        RationalInterval(x.lower + i * step, x.lower + (i + 1) * step)
        for i in range(cells)
    )


def _log_atanh_enclosure(
    x: Fraction,
    target_width: Fraction,
    max_terms: int,
) -> RationalInterval:
    if x < 1:
        raise ValueError("atanh log helper requires x >= 1")
    if target_width <= 0:
        raise ValueError("target_width must be positive")
    max_terms = _positive_index(max_terms, "max_terms")
    if x == 1:
        return RationalInterval.point(0)

    z = (x - 1) / (x + 1)
    if not 0 <= z < 1:
        raise ArithmeticError("invalid atanh reduction")
    z2 = z * z
    power = z
    partial = Fraction(0)
    for j in range(max_terms):
        denominator = 2 * j + 1
        partial += 2 * power / denominator
        next_power = power * z2
        next_denominator = denominator + 2
        tail = 2 * next_power / (next_denominator * (1 - z2))
        if tail <= target_width:
            return RationalInterval(partial, partial + tail)
        power = next_power
    raise ArithmeticError("log series term cap reached before requested precision")


def _log_fraction_enclosure(
    x: Fraction,
    *,
    precision_bits: int,
    max_terms: int,
    max_range_reductions: int,
) -> RationalInterval:
    if x < 1:
        raise ValueError("kernel log input must be at least one")
    if x == 1:
        return RationalInterval.point(0)
    precision_bits = _positive_index(precision_bits, "precision_bits")
    max_range_reductions = _nonnegative_index(max_range_reductions, "max_range_reductions")
    target = Fraction(1, 1 << precision_bits)

    y = x
    k = 0
    while y >= 2:
        if k >= max_range_reductions:
            raise ArithmeticError("log range-reduction cap reached")
        y /= 2
        k += 1

    if k == 0:
        return _log_atanh_enclosure(y, target, max_terms)

    log2_budget = target / (2 * k)
    reduced_budget = target / 2
    log2 = _log_atanh_enclosure(Fraction(2), log2_budget, max_terms).scale(k)
    reduced = _log_atanh_enclosure(y, reduced_budget, max_terms)
    out = log2 + reduced
    if out.width > target:
        raise ArithmeticError("log enclosure exceeded requested approximation width")
    return out


def _exp_neg_fraction_enclosure(
    y: Fraction,
    *,
    precision_bits: int,
    max_terms: int,
    max_scale_reductions: int,
) -> RationalInterval:
    if y < 0:
        raise ValueError("exp(-y) helper requires y >= 0")
    if y == 0:
        return RationalInterval.point(1)
    precision_bits = _positive_index(precision_bits, "precision_bits")
    max_terms = _positive_index(max_terms, "max_terms")
    max_scale_reductions = _nonnegative_index(max_scale_reductions, "max_scale_reductions")

    q = y
    scale = 0
    while q > Fraction(1, 2):
        if scale >= max_scale_reductions:
            raise ArithmeticError("exp range-reduction cap reached")
        q /= 2
        scale += 1

    target = Fraction(1, 1 << precision_bits)
    base_target = target / (1 << scale)

    partial = Fraction(1)
    upper = Fraction(1)
    lower = Fraction(0)
    term = Fraction(1)
    found = False
    for n in range(1, max_terms + 1):
        term = term * q / n
        if n % 2:
            partial -= term
            lower = partial
        else:
            partial += term
            upper = partial
        if upper < lower:
            raise ArithmeticError("alternating exponential bracket inverted")
        if upper - lower <= base_target:
            found = True
            break
    if not found:
        raise ArithmeticError("exp series term cap reached before requested precision")

    out = RationalInterval(lower, upper)
    for _ in range(scale):
        out = RationalInterval(out.lower * out.lower, out.upper * out.upper)
    if out.width > target:
        raise ArithmeticError("exp enclosure exceeded requested approximation width")
    return out


def _log_interval(
    x: RationalInterval,
    *,
    precision_bits: int,
    max_terms: int,
    max_range_reductions: int,
) -> RationalInterval:
    if x.lower < 1:
        raise ValueError("kernel logarithm interval must be at least one")
    lo = _log_fraction_enclosure(
        x.lower,
        precision_bits=precision_bits,
        max_terms=max_terms,
        max_range_reductions=max_range_reductions,
    )
    hi = _log_fraction_enclosure(
        x.upper,
        precision_bits=precision_bits,
        max_terms=max_terms,
        max_range_reductions=max_range_reductions,
    )
    return RationalInterval(lo.lower, hi.upper)


def _exp_neg_interval(
    y: RationalInterval,
    *,
    precision_bits: int,
    max_terms: int,
    max_scale_reductions: int,
) -> RationalInterval:
    if y.lower < 0:
        raise ValueError("exp(-y) interval must be nonnegative")
    at_upper = _exp_neg_fraction_enclosure(
        y.upper,
        precision_bits=precision_bits,
        max_terms=max_terms,
        max_scale_reductions=max_scale_reductions,
    )
    at_lower = _exp_neg_fraction_enclosure(
        y.lower,
        precision_bits=precision_bits,
        max_terms=max_terms,
        max_scale_reductions=max_scale_reductions,
    )
    return RationalInterval(at_upper.lower, at_lower.upper)


def _kernel_value_box(
    a: RationalInterval,
    eta: RationalInterval,
    *,
    precision_bits: int,
    max_terms: int,
    max_range_reductions: int,
) -> RationalInterval:
    if a.lower < 0 or a.upper > 1:
        raise ValueError("schedule exponent interval must lie in [0, 1]")
    eta2 = _square_interval(eta)
    base = RationalInterval.point(1) + eta2

    if a.upper == 0 or base.upper == 1:
        return RationalInterval.point(1)

    if a.lower == a.upper == 1:
        return RationalInterval(
            Fraction(1, 1) / (base.upper * base.upper),
            Fraction(1, 1) / (base.lower * base.lower),
        )

    log_base = _log_interval(
        base,
        precision_bits=precision_bits,
        max_terms=max_terms,
        max_range_reductions=max_range_reductions,
    )
    exponent = a.scale(2) * log_base
    return _exp_neg_interval(
        exponent,
        precision_bits=precision_bits,
        max_terms=max_terms,
        max_scale_reductions=max_range_reductions,
    )


def _normalized_jet_box(
    a: RationalInterval,
    eta: RationalInterval,
    order: int,
    *,
    precision_bits: int,
    max_terms: int,
    max_range_reductions: int,
) -> tuple[RationalInterval, ...]:
    base = RationalInterval.point(1) + _square_interval(eta)
    c0 = _kernel_value_box(
        a,
        eta,
        precision_bits=precision_bits,
        max_terms=max_terms,
        max_range_reductions=max_range_reductions,
    )
    out = [c0]
    if order == 0:
        return tuple(out)

    four_a = a.scale(4)
    c1 = -(four_a * eta * c0) / base
    out.append(c1)

    for n in range(1, order):
        first_factor = RationalInterval.point(2 * n) + four_a
        second_factor = RationalInterval.point(n - 1) + four_a
        numerator = first_factor * eta * out[n] + second_factor * out[n - 1]
        denominator = base.scale(n + 1)
        out.append((-numerator) / denominator)
    return tuple(out)


def _hull_rows(rows: Iterable[tuple[RationalInterval, ...]]) -> tuple[RationalInterval, ...]:
    iterator = iter(rows)
    try:
        first = list(next(iterator))
    except StopIteration as exc:
        raise ValueError("at least one enclosure row is required") from exc
    for row in iterator:
        if len(row) != len(first):
            raise ArithmeticError("jet rows have inconsistent lengths")
        first = [left.hull(right) for left, right in zip(first, row)]
    return tuple(first)


@dataclass(frozen=True)
class PressureKernelJetEnclosure:
    """Rigorous finite real eta-jet for one exact-rational parameter box."""

    exponent: RationalInterval
    eta: RationalInterval
    normalized_coefficients: tuple[RationalInterval, ...]
    exponent_cells: int
    eta_cells: int
    precision_bits: int

    @property
    def order(self) -> int:
        return len(self.normalized_coefficients) - 1

    @property
    def full_derivatives(self) -> tuple[RationalInterval, ...]:
        return tuple(
            coefficient.scale(factorial(m))
            for m, coefficient in enumerate(self.normalized_coefficients)
        )


def pressure_kernel_normalized_taylor_enclosure(
    exponent: ExactScalar | RationalInterval,
    eta: ExactScalar | RationalInterval,
    order: int,
    *,
    precision_bits: int = 96,
    exponent_cells: int = 1,
    eta_cells: int = 1,
    max_order: int = 64,
    max_cells: int = 4096,
    max_series_terms: int = 256,
    max_range_reductions: int = 4096,
) -> PressureKernelJetEnclosure:
    """Enclose ``k_a^(m)(eta)/m!`` for every ``0 <= m <= order``.

    ``exponent`` and ``eta`` must be exact rational points or
    :class:`RationalInterval` objects.  No float/Decimal/complex coercion is
    accepted.  Optional subdivisions preserve the same shared parameter boxes
    and hull their rigorous results; they never sample a midpoint.
    """

    order = _nonnegative_index(order, "order")
    max_order = _nonnegative_index(max_order, "max_order")
    if order > max_order:
        raise ValueError("requested order exceeds max_order")
    precision_bits = _positive_index(precision_bits, "precision_bits")
    if precision_bits < 8 or precision_bits > 4096:
        raise ValueError("precision_bits must lie in [8, 4096]")
    exponent_cells = _positive_index(exponent_cells, "exponent_cells")
    eta_cells = _positive_index(eta_cells, "eta_cells")
    max_cells = _positive_index(max_cells, "max_cells")
    max_series_terms = _positive_index(max_series_terms, "max_series_terms")
    max_range_reductions = _nonnegative_index(max_range_reductions, "max_range_reductions")

    a_box = _coerce_interval(exponent, "exponent")
    eta_box = _coerce_interval(eta, "eta")
    if a_box.lower < 0 or a_box.upper > 1:
        raise ValueError("schedule exponent interval must lie in [0, 1]")

    a_parts = _split_interval(a_box, exponent_cells)
    eta_parts = _split_interval(eta_box, eta_cells)
    if len(a_parts) * len(eta_parts) > max_cells:
        raise ValueError("requested subdivision exceeds max_cells")

    rows = []
    for a_part in a_parts:
        for eta_part in eta_parts:
            rows.append(
                _normalized_jet_box(
                    a_part,
                    eta_part,
                    order,
                    precision_bits=precision_bits,
                    max_terms=max_series_terms,
                    max_range_reductions=max_range_reductions,
                )
            )

    normalized = _hull_rows(rows)
    return PressureKernelJetEnclosure(
        exponent=a_box,
        eta=eta_box,
        normalized_coefficients=normalized,
        exponent_cells=len(a_parts),
        eta_cells=len(eta_parts),
        precision_bits=precision_bits,
    )


def pressure_kernel_value_enclosure(
    exponent: ExactScalar | RationalInterval,
    eta: ExactScalar | RationalInterval,
    **kwargs: object,
) -> RationalInterval:
    """Return a strict real enclosure of ``(1+eta^2)^(-2*exponent)``."""

    return pressure_kernel_normalized_taylor_enclosure(
        exponent, eta, 0, **kwargs
    ).normalized_coefficients[0]


def pressure_kernel_derivative_enclosure(
    exponent: ExactScalar | RationalInterval,
    eta: ExactScalar | RationalInterval,
    order: int,
    **kwargs: object,
) -> RationalInterval:
    """Return a strict enclosure of the full eta derivative ``k_a^(order)``."""

    jet = pressure_kernel_normalized_taylor_enclosure(
        exponent, eta, order, **kwargs
    )
    return jet.full_derivatives[order]
