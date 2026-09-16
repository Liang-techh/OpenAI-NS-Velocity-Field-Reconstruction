"""Exact affine phase-to-log-amplitude interval transport.

Given a rational interval for a phase value and explicitly selected scalar
values ``Lambda`` and ``log_C``, the natural-axis relation

    log(a) = Lambda * phase - log_C

is transported coefficientwise into an exact rational interval.  The phase
interval and scalar choices are caller-supplied inputs: this constructor does
not evaluate a phase quadrature, certify the parameter selection, or establish
the global reconstruction.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import (
    Context,
    Decimal,
    MAX_EMAX,
    MIN_EMIN,
    ROUND_CEILING,
    ROUND_FLOOR,
    ROUND_HALF_EVEN,
    localcontext,
)
from fractions import Fraction


def _positive_precision(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError("precision must be a positive integer")
    return value


def _decimal_from_fraction(
    value: Fraction,
    precision: int,
    rounding: str,
) -> Decimal:
    """Round one rational with a fresh, wide context and no inherited traps."""

    context = Context(
        prec=precision,
        rounding=rounding,
        Emax=MAX_EMAX,
        Emin=MIN_EMIN,
    )
    for signal in context.traps:
        context.traps[signal] = False
    with localcontext(context):
        result = +(Decimal(value.numerator) / Decimal(value.denominator))
    if not result.is_finite():
        raise ArithmeticError("Decimal bound is outside the implementation range")
    return result


@dataclass(frozen=True)
class RationalLogAmplitudeEnclosure:
    """Exact interval for ``Lambda * phase - log_C``.

    The enclosure is conditional on the supplied rational phase interval and
    explicitly supplied scalar values.  It does not certify the phase
    interval, parameter selection, global norm bounds, or the full velocity
    reconstruction.
    """

    phase_lower: Fraction
    phase_upper: Fraction
    Lambda: Decimal
    log_C: Decimal
    lower: Fraction = field(init=False)
    upper: Fraction = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.phase_lower, Fraction):
            raise TypeError("phase_lower must be a Fraction")
        if not isinstance(self.phase_upper, Fraction):
            raise TypeError("phase_upper must be a Fraction")
        if self.phase_lower > self.phase_upper:
            raise ValueError("phase bounds must be ordered")
        if not isinstance(self.Lambda, Decimal):
            raise TypeError("Lambda must be a Decimal")
        if not self.Lambda.is_finite():
            raise ValueError("Lambda must be finite")
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")
        if not isinstance(self.log_C, Decimal):
            raise TypeError("log_C must be a Decimal")
        if not self.log_C.is_finite():
            raise ValueError("log_C must be finite")

        lambda_fraction = Fraction(self.Lambda)
        log_C_fraction = Fraction(self.log_C)
        lower = lambda_fraction * self.phase_lower - log_C_fraction
        upper = lambda_fraction * self.phase_upper - log_C_fraction
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)

    @property
    def width(self) -> Fraction:
        """Return the exact scaled phase width as a rational number."""

        return self.upper - self.lower

    @property
    def paper_exact(self) -> bool:
        """The supplied interval and scalars are not a paper-exact proof."""

        return False

    def decimal_bounds(self, precision: int = 96) -> tuple[Decimal, Decimal]:
        """Return directed Decimal endpoints at the requested precision.

        A fresh context uses the requested precision and direction while
        widening exponent limits and disabling inherited traps.  The returned
        endpoints are checked against the exact rational bounds before they
        are exposed.
        """

        precision = _positive_precision(precision)
        lower = _decimal_from_fraction(self.lower, precision, ROUND_FLOOR)
        upper = _decimal_from_fraction(self.upper, precision, ROUND_CEILING)
        lower_fraction = Fraction(lower)
        upper_fraction = Fraction(upper)
        if lower_fraction > self.lower or upper_fraction < self.upper:
            raise ArithmeticError("directed Decimal bounds failed to enclose the exact interval")
        return lower, upper

    def midpoint_decimal(self, precision: int = 96) -> Decimal:
        """Return the exact interval midpoint rounded to nearest Decimal."""

        precision = _positive_precision(precision)
        midpoint = (self.lower + self.upper) / 2
        return _decimal_from_fraction(midpoint, precision, ROUND_HALF_EVEN)


__all__ = ["RationalLogAmplitudeEnclosure"]
