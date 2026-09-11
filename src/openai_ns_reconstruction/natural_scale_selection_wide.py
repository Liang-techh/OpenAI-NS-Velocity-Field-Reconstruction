"""Wide/log-domain natural-profile scale selection.

The pinned Stage-1 theorem chain chooses

``Lambda = max(1 + B + L, 1 + (14000/9) B)``

and then ``C = exp(Lambda * realPartSup(axisPhase, compactSet))``.  The narrow
adapter in :mod:`natural_scale_selection` stores every quantity in binary64,
which is intentionally fail-closed but cannot consume the wide Decimal
remainder majorants produced by :mod:`axis_remainder_wide_bounds`.

This module changes only the scalar representation.  Positive algebra is
performed in high-precision Decimal arithmetic rounded toward ``+infinity``.
The normalization constant is kept as the exact symbolic exponential
``exp(exponent_upper)`` rather than evaluating an astronomically large real
number.  That symbolic representation is enough to carry the theorem's
Lambda/C *selection* past a finite-float representation barrier; it is not a
materialized coefficient-space fixed point and is not a paper-exact profile.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, ROUND_CEILING, localcontext
import math
from typing import Union

DecimalLike = Union[Decimal, float, int]
_DECIMAL_PRECISION = 96


def _finite_decimal(value: DecimalLike, name: str) -> Decimal:
    if isinstance(value, Decimal):
        result = value
    elif isinstance(value, int):
        result = Decimal(value)
    else:
        numeric = float(value)
        if not math.isfinite(numeric):
            raise ValueError(f"{name} must be finite")
        result = Decimal.from_float(numeric)
    if not result.is_finite():
        raise ValueError(f"{name} must be finite")
    return result


def _nonnegative_decimal(value: DecimalLike, name: str) -> Decimal:
    result = _finite_decimal(value, name)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _positive_decimal(value: DecimalLike, name: str) -> Decimal:
    result = _finite_decimal(value, name)
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _add_up(*values: DecimalLike) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.rounding = ROUND_CEILING
        total = Decimal(0)
        for value in values:
            total += _nonnegative_decimal(value, "summand")
        return +total


def _mul_up(*values: DecimalLike) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.rounding = ROUND_CEILING
        product = Decimal(1)
        for value in values:
            product *= _nonnegative_decimal(value, "factor")
        return +product


def _mul_ratio_up(value: DecimalLike, numerator: int, denominator: int) -> Decimal:
    if numerator < 0 or denominator <= 0:
        raise ValueError("ratio must be nonnegative with positive denominator")
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.rounding = ROUND_CEILING
        result = _nonnegative_decimal(value, "ratio input") * Decimal(numerator)
        result /= Decimal(denominator)
        return +result


@dataclass(frozen=True)
class SymbolicExponentialThreshold:
    """Representation of a positive threshold ``exp(exponent_upper)``.

    We deliberately store the exponent, not the exponential value.  For the
    actual landed schedule the latter can have a decimal exponent far beyond
    ordinary floating-point range even though the theorem-side real number is
    perfectly well-defined.
    """

    exponent_upper: Decimal

    def __post_init__(self) -> None:
        _nonnegative_decimal(self.exponent_upper, "exponent_upper")

    def admits_exponential_exponent(self, candidate_exponent: DecimalLike) -> bool:
        """Use monotonicity of exp to compare another symbolic positive C."""

        try:
            candidate = _nonnegative_decimal(candidate_exponent, "candidate_exponent")
        except ValueError:
            return False
        return candidate >= self.exponent_upper


@dataclass(frozen=True)
class WideNaturalScaleSelection:
    """Wide representation of the pinned deterministic Lambda/C choice."""

    remainder_bound: Decimal
    remainder_lipschitz: Decimal
    phase_real_part_sup: Decimal
    contraction: Decimal
    stability: Decimal
    Lambda: Decimal
    C: SymbolicExponentialThreshold

    def __post_init__(self) -> None:
        _nonnegative_decimal(self.remainder_bound, "remainder_bound")
        _nonnegative_decimal(self.remainder_lipschitz, "remainder_lipschitz")
        _nonnegative_decimal(self.phase_real_part_sup, "phase_real_part_sup")
        _positive_decimal(self.contraction, "contraction")
        _positive_decimal(self.stability, "stability")
        _positive_decimal(self.Lambda, "Lambda")

    def admits_symbolic(self, Lambda: DecimalLike, C_exponent: DecimalLike) -> bool:
        """Check theorem inequalities when C is represented by its exp exponent."""

        try:
            candidate_lambda = _positive_decimal(Lambda, "Lambda")
            candidate_c_exponent = _nonnegative_decimal(C_exponent, "C_exponent")
        except ValueError:
            return False
        if candidate_lambda < self.contraction or candidate_lambda < self.stability:
            return False
        required_exponent = _mul_up(candidate_lambda, self.phase_real_part_sup)
        return candidate_c_exponent >= required_exponent


def contraction_threshold_wide(
    remainder_bound: DecimalLike,
    remainder_lipschitz: DecimalLike,
) -> Decimal:
    """Wide image of pinned ``1 + remainderBound + remainderLip``."""

    bound = _nonnegative_decimal(remainder_bound, "remainder_bound")
    lip = _nonnegative_decimal(remainder_lipschitz, "remainder_lipschitz")
    return _add_up(1, bound, lip)


def stability_scale_wide(remainder_bound: DecimalLike) -> Decimal:
    """Wide image of pinned ``1 + (14000/9) * remainderBound``."""

    bound = _nonnegative_decimal(remainder_bound, "remainder_bound")
    return _add_up(1, _mul_ratio_up(bound, 14000, 9))


def select_natural_scale_wide(
    *,
    remainder_bound: DecimalLike,
    remainder_lipschitz: DecimalLike,
    phase_real_part_sup: DecimalLike,
) -> WideNaturalScaleSelection:
    """Carry the exact theorem-side Lambda/C algebra without float conversion.

    ``C`` is returned as the symbolic threshold ``exp(exponent_upper)`` with
    ``exponent_upper = Lambda * phase_real_part_sup``.  This avoids mistaking an
    enormous but finite theorem quantity for numerical infinity.
    """

    bound = _nonnegative_decimal(remainder_bound, "remainder_bound")
    lip = _nonnegative_decimal(remainder_lipschitz, "remainder_lipschitz")
    phase_sup = _nonnegative_decimal(phase_real_part_sup, "phase_real_part_sup")

    contraction = contraction_threshold_wide(bound, lip)
    stability = stability_scale_wide(bound)
    Lambda = max(contraction, stability)
    exponent = _mul_up(Lambda, phase_sup)

    return WideNaturalScaleSelection(
        remainder_bound=bound,
        remainder_lipschitz=lip,
        phase_real_part_sup=phase_sup,
        contraction=contraction,
        stability=stability,
        Lambda=Lambda,
        C=SymbolicExponentialThreshold(exponent),
    )
