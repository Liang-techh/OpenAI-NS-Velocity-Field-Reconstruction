"""Exact rational log enclosures for finite amplitude derivatives.

The amplitude coefficient is represented as ``a(eta) = exp(f(eta))`` with
``f`` supplied by :class:`AmplitudeLogSource`.  The rational coefficient
provider supplies the exact Bell factor ``B_m`` in

    d_eta^m a(eta)^q = a(eta)^q B_m.

This module encloses ``log(abs(B_m))`` without converting to a floating point
logarithm.  Every endpoint used by the enclosure is a :class:`Fraction`;
the resulting interval is conditional on the selected finite kernel and
scalar data and does not certify downstream coefficient arithmetic.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from functools import lru_cache

from .axis_amplitude_log_scale import AmplitudeLogSource


_DEFAULT_LOG_TOLERANCE = Fraction(1, 10**96)


def _require_positive_fraction(value: object, name: str) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be a Fraction")
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _require_positive_cap(value: object, name: str = "max_terms") -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _dyadic_step(bound: Fraction) -> Fraction:
    """Return a dyadic step no larger than the positive rational ``bound``."""

    # The initial estimate only chooses a nearby exponent.  The exact loop
    # below is what establishes the required inequality, including for very
    # large or very small rationals.
    exponent = max(4, bound.denominator.bit_length() - bound.numerator.bit_length())
    step = Fraction(1, 1 << exponent)
    while step > bound:
        exponent += 1
        step = Fraction(1, 1 << exponent)
    while exponent > 4 and Fraction(1, 1 << (exponent - 1)) <= bound:
        exponent -= 1
        step = Fraction(1, 1 << exponent)
    return step


def _floor_dyadic(value: Fraction, step: Fraction) -> Fraction:
    return Fraction(value.numerator // value.denominator) * step if step == 1 else (
        value // step
    ) * step


def _ceil_dyadic(value: Fraction, step: Fraction) -> Fraction:
    quotient = value / step
    return (-((-quotient.numerator) // quotient.denominator)) * step


def _binary_reduce(value: Fraction) -> tuple[int, Fraction]:
    """Return ``k, m`` with ``value = 2**k * m`` and ``1 <= m < 2``."""

    numerator = value.numerator
    denominator = value.denominator
    k = numerator.bit_length() - denominator.bit_length()
    if k >= 0:
        if numerator < (denominator << k):
            k -= 1
    else:
        if (numerator << (-k)) < denominator:
            k -= 1
    if k < 0:
        shift = -k
        mantissa = Fraction(numerator << shift, denominator)
    else:
        mantissa = Fraction(numerator, denominator << k)
    if not Fraction(1) <= mantissa < Fraction(2):
        raise ArithmeticError("binary reduction failed to produce a normalized mantissa")
    return k, mantissa


def _finite_atanh_sum(value: Fraction, terms: int) -> Fraction:
    """Return ``2 * sum(value**(2n+1)/(2n+1), n < terms)`` exactly."""

    total = Fraction(0)
    for n in range(terms):
        exponent = 2 * n + 1
        total += Fraction(2) * value**exponent / exponent
    return total


def _tail_bound(value: Fraction, terms: int) -> Fraction:
    """Bound the omitted positive tail after ``terms`` terms."""

    exponent = 2 * terms + 1
    return Fraction(2) * value**exponent / (
        exponent * (Fraction(1) - value * value)
    )


def _atanh_log_interval(
    t: Fraction,
    budget: Fraction,
    max_terms: int,
) -> "RationalLogInterval":
    """Enclose ``2 * atanh(t)`` for ``0 <= t < 1/2``."""

    if not Fraction(0) <= t < Fraction(1, 2):
        raise ValueError("atanh argument must lie in [0, 1/2)")
    if t == 0:
        return RationalLogInterval(Fraction(0), Fraction(0), 0)

    # Adjacent grid points have width d.  Since the derivative is at most
    # 8/3 on [0, 1/2], the input-rounding contribution is at most b/6.
    step = _dyadic_step(min(Fraction(1, 16), budget / 16))
    lower_t = _floor_dyadic(t, step)
    upper_t = _ceil_dyadic(t, step)
    if upper_t > Fraction(1, 2):
        # t < 1/2 and d <= 1/16 make this unreachable, but retain the
        # invariant explicitly if the grid policy is changed later.
        raise ArithmeticError("dyadic atanh grid escaped its convergence interval")

    terms = 0
    tail = None
    while terms < max_terms:
        terms += 1
        tail = _tail_bound(upper_t, terms)
        if tail <= budget / 2:
            break
    else:
        raise ArithmeticError(
            "rational logarithm enclosure did not reach the requested tolerance "
            "before max_terms"
        )

    assert tail is not None
    lower = _finite_atanh_sum(lower_t, terms)
    upper = _finite_atanh_sum(upper_t, terms) + tail

    # Round both endpoints outward to the same dyadic grid.  Each endpoint
    # moves by less than d, so the total rounding contribution is at most
    # b/8.  The exact interval remains wider than the unrounded one only by
    # this controlled amount.
    lower = _floor_dyadic(lower, step)
    upper = _ceil_dyadic(upper, step)
    interval = RationalLogInterval(lower, upper, terms)
    if interval.width > budget:
        raise ArithmeticError("rational logarithm component exceeded its exact budget")
    return interval


@lru_cache(maxsize=64)
def _cached_ln2(budget: Fraction, max_terms: int) -> "RationalLogInterval":
    # ln(2) = 2 atanh(1/3), evaluated using exactly the same rational bounds.
    return _atanh_log_interval(Fraction(1, 3), budget, max_terms)


@dataclass(frozen=True)
class RationalLogInterval:
    """An ordered exact rational enclosure of one logarithm."""

    lower: Fraction
    upper: Fraction
    terms_used: int

    def __post_init__(self) -> None:
        if not isinstance(self.lower, Fraction):
            raise TypeError("lower must be a Fraction")
        if not isinstance(self.upper, Fraction):
            raise TypeError("upper must be a Fraction")
        if self.lower > self.upper:
            raise ValueError("log interval bounds must be ordered")
        if isinstance(self.terms_used, bool) or not isinstance(self.terms_used, int):
            raise ValueError("terms_used must be a nonnegative integer")
        if self.terms_used < 0:
            raise ValueError("terms_used must be a nonnegative integer")

    @property
    def width(self) -> Fraction:
        return self.upper - self.lower


def rational_log_enclosure(
    value: Fraction,
    *,
    absolute_tolerance: Fraction = _DEFAULT_LOG_TOLERANCE,
    max_terms: int = 1024,
) -> RationalLogInterval:
    """Return an exact interval for ``ln(value)``.

    Binary reduction bounds the mantissa in ``[1, 2)``.  The mantissa log and
    ``ln(2)`` are evaluated with the positive atanh series, all in Fraction
    arithmetic.  The per-component budget is ``tol / (abs(k) + 1)`` so the
    final weighted interval is guaranteed to fit the requested tolerance.
    """

    if not isinstance(value, Fraction):
        raise TypeError("value must be a Fraction")
    if value <= 0:
        raise ValueError("value must be positive")
    tolerance = _require_positive_fraction(absolute_tolerance, "absolute_tolerance")
    max_terms = _require_positive_cap(max_terms)

    # Keep validation ahead of this fast path: invalid tolerance/cap values
    # must not be hidden by the exact logarithm of one.
    if value == 1:
        return RationalLogInterval(Fraction(0), Fraction(0), 0)

    k, mantissa = _binary_reduce(value)
    budget = tolerance / (abs(k) + 1)
    mantissa_t = (mantissa - 1) / (mantissa + 1)
    mantissa_interval = _atanh_log_interval(mantissa_t, budget, max_terms)

    if k == 0:
        lower = mantissa_interval.lower
        upper = mantissa_interval.upper
        terms_used = mantissa_interval.terms_used
    else:
        ln2 = _cached_ln2(budget, max_terms)
        if k > 0:
            lower = mantissa_interval.lower + k * ln2.lower
            upper = mantissa_interval.upper + k * ln2.upper
        else:
            lower = mantissa_interval.lower + k * ln2.upper
            upper = mantissa_interval.upper + k * ln2.lower
        terms_used = max(mantissa_interval.terms_used, ln2.terms_used)

    result = RationalLogInterval(lower, upper, terms_used)
    if result.width > tolerance:
        raise ArithmeticError("rational logarithm enclosure exceeded the requested tolerance")
    return result


def _require_power(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value not in (1, 2):
        raise ValueError("power must be 1 or 2")
    return value


def _require_order(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("order must be a nonnegative integer")
    return value


@dataclass(frozen=True)
class AmplitudeDerivativeLogEnclosure:
    """Full log enclosure for ``d_eta^order (a(eta)^power)``.

    ``exact_bell_factor`` retains the exact signed rational factor.  A zero
    factor is represented by ``sign == 0`` and no logarithmic interval.  For
    a nonzero factor, ``total_log_lower`` and ``total_log_upper`` enclose the
    log magnitude of the complete derivative, including the common amplitude
    power from ``source.enclosure``.
    """

    source: AmplitudeLogSource
    power: int
    order: int
    exact_bell_factor: Fraction
    log_factor: RationalLogInterval | None
    sign: int = field(init=False)
    total_log_lower: Fraction | None = field(init=False)
    total_log_upper: Fraction | None = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.source, AmplitudeLogSource):
            raise TypeError("source must be an AmplitudeLogSource")
        power = _require_power(self.power)
        order = _require_order(self.order)
        if not isinstance(self.exact_bell_factor, Fraction):
            raise TypeError("exact_bell_factor must be a Fraction")
        if self.exact_bell_factor == 0:
            if self.log_factor is not None:
                raise ValueError("zero Bell factor must not carry a logarithm interval")
            sign = 0
            lower = None
            upper = None
        else:
            if not isinstance(self.log_factor, RationalLogInterval):
                raise TypeError("nonzero Bell factor requires a RationalLogInterval")
            sign = 1 if self.exact_bell_factor > 0 else -1
            lower = power * self.source.enclosure.lower + self.log_factor.lower
            upper = power * self.source.enclosure.upper + self.log_factor.upper
            if lower > upper:
                raise ArithmeticError("full amplitude derivative log bounds are unordered")

        object.__setattr__(self, "power", power)
        object.__setattr__(self, "order", order)
        object.__setattr__(self, "sign", sign)
        object.__setattr__(self, "total_log_lower", lower)
        object.__setattr__(self, "total_log_upper", upper)

    @property
    def total_log_width(self) -> Fraction | None:
        if self.total_log_lower is None or self.total_log_upper is None:
            return None
        return self.total_log_upper - self.total_log_lower

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def downstream_arithmetic_certified(self) -> bool:
        return False


__all__ = [
    "AmplitudeDerivativeLogEnclosure",
    "RationalLogInterval",
    "rational_log_enclosure",
]
