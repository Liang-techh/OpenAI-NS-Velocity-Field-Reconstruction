"""Conditional radial-tail bounds for the natural-axis coefficient weights.

The bound in this module is an analytic consequence of a caller-supplied
global ``AxisCoefficientSpace`` norm bound.  It does not inspect finite
coefficient rows, infer a norm from them, or certify the formal solver's
uncomputed tail.  Decimal inputs are converted to exact ``Fraction`` values
before the positive majorant is evaluated; the final result is rounded upward
in a local 96-digit Decimal context.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, ROUND_CEILING, localcontext
from fractions import Fraction
import math


_DECIMAL_PRECISION = 96
_TWENTY = Decimal(20)


def _decimal(value: Decimal, name: str) -> Decimal:
    if isinstance(value, bool) or not isinstance(value, Decimal):
        raise TypeError(f"{name} must be a Decimal")
    if not value.is_finite():
        raise ValueError(f"{name} must be finite")
    return value


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _fraction_to_decimal_ceiling(value: Fraction) -> Decimal:
    if value < 0:
        raise ValueError("tail majorant must be nonnegative")
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.rounding = ROUND_CEILING
        return +(Decimal(value.numerator) / Decimal(value.denominator))


def _negative_binomial_tail(q: Fraction, J: int, K: int) -> Fraction:
    """Return ``sum(k >= J) binom(k+K,K) q^k`` for ``0 <= q < 1``."""

    one_minus_q = Fraction(1, 1) - q
    if J == 0:
        return one_minus_q ** (-(K + 1))
    if q == 0:
        return Fraction(0, 1)
    polynomial = sum(
        (Fraction(math.comb(J - 1 + k, k), 1) * one_minus_q**k)
        for k in range(K + 1)
    )
    return q**J * one_minus_q ** (-(K + 1)) * polynomial


@dataclass(frozen=True)
class AxisCoefficientRadialTailBound:
    """One conditional majorant for omitted radial coefficient terms."""

    upper_bound: Decimal
    norm_upper: Decimal
    epsilon: Decimal
    Y: Decimal
    max_n: int
    eta_order: int = 0
    radial_order: int = 0
    average: bool = False
    conditional_on_global_norm: bool = True
    includes_coefficient_roundoff: bool = False
    paper_exact: bool = False

    def __post_init__(self) -> None:
        for name in ("upper_bound", "norm_upper", "epsilon", "Y"):
            _decimal(getattr(self, name), name)
        if self.upper_bound < 0:
            raise ValueError("upper_bound must be nonnegative")
        if self.norm_upper < 0:
            raise ValueError("norm_upper must be nonnegative")
        if self.epsilon <= 0:
            raise ValueError("epsilon must be positive")
        if self.Y.copy_abs() >= _TWENTY:
            raise ValueError("Y must satisfy abs(Y) < 20")
        _index(self.max_n, "max_n")
        _index(self.eta_order, "eta_order")
        _index(self.radial_order, "radial_order")
        if not isinstance(self.average, bool):
            raise TypeError("average must be a boolean")
        if self.conditional_on_global_norm is not True:
            raise ValueError("the bound is conditional on a global norm upper bound")
        if self.includes_coefficient_roundoff is not False:
            raise ValueError("coefficient roundoff is not included in this bound")
        if self.paper_exact is not False:
            raise ValueError("this analytic majorant is not paper-exact")


def axis_coefficient_radial_tail_bound(
    norm_upper: Decimal,
    epsilon: Decimal,
    Y: Decimal,
    max_n: int,
    eta_order: int = 0,
    radial_order: int = 0,
    average: bool = False,
) -> AxisCoefficientRadialTailBound:
    """Bound the omitted radial tail under a supplied global norm bound.

    The omitted rows begin at ``L = max(max_n + 1, radial_order)``.  The
    estimate uses the pinned weight

    ``20^-n epsilon^-m m! binom(n+m,m) / ((n+1)^2 (m+1)^2)``

    and the absolute radial factor ``Y^(n-r) n!/(n-r)!``.  When ``average`` is
    true, the additional radial-average divisor ``1/(n+1)`` is included.
    """

    norm_upper = _decimal(norm_upper, "norm_upper")
    epsilon = _decimal(epsilon, "epsilon")
    Y = _decimal(Y, "Y")
    if norm_upper < 0:
        raise ValueError("norm_upper must be nonnegative")
    if epsilon <= 0:
        raise ValueError("epsilon must be positive")
    if Y.copy_abs() >= _TWENTY:
        raise ValueError("Y must satisfy abs(Y) < 20")
    max_n = _index(max_n, "max_n")
    eta_order = _index(eta_order, "eta_order")
    radial_order = _index(radial_order, "radial_order")
    if not isinstance(average, bool):
        raise TypeError("average must be a boolean")

    # Every operation below is exact rational arithmetic.  In particular,
    # no Decimal division is used before the final directed conversion.
    M = Fraction(norm_upper)
    epsilon_exact = Fraction(epsilon)
    q = Fraction(Y.copy_abs()) / Fraction(20, 1)
    L = max(max_n + 1, radial_order)
    J = L - radial_order
    K = eta_order + radial_order
    denominator_power = 2 + int(average)
    H = _negative_binomial_tail(q, J, K)
    rational_bound = (
        M
        * Fraction(math.factorial(K), 20**radial_order)
        / (epsilon_exact**eta_order)
        / ((eta_order + 1) ** 2)
        / ((L + 1) ** denominator_power)
        * H
    )
    upper_bound = _fraction_to_decimal_ceiling(rational_bound)
    return AxisCoefficientRadialTailBound(
        upper_bound=upper_bound,
        norm_upper=norm_upper,
        epsilon=epsilon,
        Y=Y,
        max_n=max_n,
        eta_order=eta_order,
        radial_order=radial_order,
        average=average,
    )


__all__ = [
    "AxisCoefficientRadialTailBound",
    "axis_coefficient_radial_tail_bound",
]
