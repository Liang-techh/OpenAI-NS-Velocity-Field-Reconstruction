"""Theorem-shaped fixed-prefix tail exponents for the Section 5 slow Borel sum.

The pinned ``SlowBorelBase`` formalization proves, for an admissible coefficient
hierarchy and cutoff schedule, that the order-``J`` uncut prefix has chart-tail
jet bound

    2**(-J) * q**(h*(J+1) - m),              m <= J+3,

and that after the physical pullback/power factor the common derivative-order
``M`` tail has exponent

    h*(J+1) + b - 2*M.

This module exposes only that exponent arithmetic and a finite-schedule gate.
It does *not* certify the missing coefficient identities, admissibility of a
paper-exact hierarchy, the infinite cutoff schedule, or the Navier--Stokes
residual itself.  Those remain upstream requirements before Proposition 5.3
can be claimed.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
import math

from .background_cutoff_schedule import SlowBorelCutoffSchedule


def _nonnegative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _finite_fraction(value: float, name: str) -> Fraction:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return Fraction.from_float(value)


def _height_fraction(value: float) -> Fraction:
    value = float(value)
    if not math.isfinite(value) or not 0.0 < value < 0.5:
        raise ValueError("h must be finite with 0 < h < 1/2")
    return Fraction.from_float(value)


def _ceil_fraction(value: Fraction) -> int:
    """Exact integer ceiling for a rational with positive denominator."""
    return -((-value.numerator) // value.denominator)


def chart_tail_exponent_exact(h: float, prefix_order: int, derivative_order: int) -> Fraction:
    """Return ``h*(J+1)-m`` using the exact supplied binary64 value of ``h``."""
    hq = _height_fraction(h)
    J = _nonnegative_int(prefix_order, "prefix_order")
    m = _nonnegative_int(derivative_order, "derivative_order")
    if m > J + 3:
        raise ValueError("uncut_fixed_prefix_bound requires m <= J+3")
    return hq * (J + 1) - m


def physical_tail_exponent_exact(
    h: float,
    prefix_order: int,
    max_derivative_order: int,
    leading_power: float,
) -> Fraction:
    """Return ``h*(J+1)+b-2*M`` from ``powered_fixed_prefix_bound`` exactly."""
    hq = _height_fraction(h)
    J = _nonnegative_int(prefix_order, "prefix_order")
    M = _nonnegative_int(max_derivative_order, "max_derivative_order")
    bq = _finite_fraction(leading_power, "leading_power")
    if M > J + 3:
        raise ValueError("powered_fixed_prefix_bound requires M <= J+3")
    return hq * (J + 1) + bq - 2 * M


@dataclass(frozen=True)
class FixedPrefixTailOrderCertificate:
    """Finite theorem-side exponent certificate for one chosen uncut prefix.

    The certificate records only the quantitative conclusion available *if*
    the pinned theorem hypotheses (smooth actual coefficients and admissible
    scales on the relevant compact set) are independently certified.
    """

    h: float
    prefix_order: int
    max_derivative_order: int
    leading_power: float = 0.0

    def __post_init__(self) -> None:
        _height_fraction(self.h)
        J = _nonnegative_int(self.prefix_order, "prefix_order")
        M = _nonnegative_int(self.max_derivative_order, "max_derivative_order")
        b = float(self.leading_power)
        _finite_fraction(b, "leading_power")
        if M > J + 3:
            raise ValueError("fixed-prefix derivative budget requires M <= J+3")
        object.__setattr__(self, "h", float(self.h))
        object.__setattr__(self, "prefix_order", J)
        object.__setattr__(self, "max_derivative_order", M)
        object.__setattr__(self, "leading_power", b)

    @property
    def dyadic_prefactor_exact(self) -> Fraction:
        return Fraction(1, 2**self.prefix_order)

    def chart_exponent_exact(self, derivative_order: int) -> Fraction:
        derivative_order = _nonnegative_int(derivative_order, "derivative_order")
        if derivative_order > self.max_derivative_order:
            raise ValueError("derivative_order exceeds the certified derivative budget")
        return chart_tail_exponent_exact(self.h, self.prefix_order, derivative_order)

    @property
    def worst_chart_exponent_exact(self) -> Fraction:
        return self.chart_exponent_exact(self.max_derivative_order)

    @property
    def physical_exponent_exact(self) -> Fraction:
        return physical_tail_exponent_exact(
            self.h,
            self.prefix_order,
            self.max_derivative_order,
            self.leading_power,
        )

    def certifies_chart_target(self, target_exponent: float) -> bool:
        target = _finite_fraction(target_exponent, "target_exponent")
        return self.worst_chart_exponent_exact >= target

    def certifies_physical_target(self, target_exponent: float) -> bool:
        target = _finite_fraction(target_exponent, "target_exponent")
        return self.physical_exponent_exact >= target


def minimal_prefix_for_chart_target(
    h: float,
    max_derivative_order: int,
    target_exponent: float,
    *,
    minimum_prefix_order: int = 0,
) -> FixedPrefixTailOrderCertificate:
    """Choose the smallest Lean-existence-style ``J`` for a chart-tail target.

    ``exists_ordinary_uncut_tail`` chooses ``J >= M`` and then requires
    ``h*(J+1)-M >= P``.  This function solves that inequality with exact
    rational arithmetic for the supplied binary64 inputs.
    """
    hq = _height_fraction(h)
    M = _nonnegative_int(max_derivative_order, "max_derivative_order")
    Jmin = _nonnegative_int(minimum_prefix_order, "minimum_prefix_order")
    target = _finite_fraction(target_exponent, "target_exponent")
    required = (target + M) / hq - 1
    J = max(M, Jmin, _ceil_fraction(required))
    cert = FixedPrefixTailOrderCertificate(float(h), J, M, 0.0)
    if not cert.certifies_chart_target(float(target)):
        raise RuntimeError("exact prefix selection failed to meet the requested chart exponent")
    return cert


def minimal_prefix_for_physical_target(
    h: float,
    max_derivative_order: int,
    leading_power: float,
    target_exponent: float,
    *,
    minimum_prefix_order: int = 0,
) -> FixedPrefixTailOrderCertificate:
    """Choose the smallest ``J>=M`` meeting the powered physical-tail target.

    This mirrors ``exists_powered_physical_tail_finite`` through the exact
    inequality ``h*(J+1)+b-2*M >= P``.  It selects an order only; it does not
    manufacture the theorem's compact-set constants or admissibility proof.
    """
    hq = _height_fraction(h)
    M = _nonnegative_int(max_derivative_order, "max_derivative_order")
    Jmin = _nonnegative_int(minimum_prefix_order, "minimum_prefix_order")
    bq = _finite_fraction(leading_power, "leading_power")
    target = _finite_fraction(target_exponent, "target_exponent")
    required = (target - bq + 2 * M) / hq - 1
    J = max(M, Jmin, _ceil_fraction(required))
    cert = FixedPrefixTailOrderCertificate(float(h), J, M, float(leading_power))
    if not cert.certifies_physical_target(float(target)):
        raise RuntimeError("exact prefix selection failed to meet the requested physical exponent")
    return cert


def certify_schedule_physical_tail(
    schedule: SlowBorelCutoffSchedule,
    max_derivative_order: int,
    leading_power: float,
    target_exponent: float,
    *,
    minimum_prefix_order: int = 0,
) -> FixedPrefixTailOrderCertificate:
    """Require a constructed finite schedule to reach the selected tail order.

    This is intentionally fail-closed: if the supplied finite schedule stops
    before the theorem-side prefix order, no statement is extrapolated to the
    unconstructed infinite schedule.
    """
    if not isinstance(schedule, SlowBorelCutoffSchedule):
        raise TypeError("schedule must be a SlowBorelCutoffSchedule")
    cert = minimal_prefix_for_physical_target(
        schedule.h,
        max_derivative_order,
        leading_power,
        target_exponent,
        minimum_prefix_order=minimum_prefix_order,
    )
    if cert.prefix_order > schedule.max_order:
        raise ValueError(
            "constructed cutoff schedule does not reach the required fixed-prefix order"
        )
    return cert
