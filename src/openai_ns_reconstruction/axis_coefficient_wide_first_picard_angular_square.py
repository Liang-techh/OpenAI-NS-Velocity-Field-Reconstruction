"""Pinned coefficient product ``phi1 * phi1`` at the genuine first Picard state.

The landed first Picard angular field is intentionally represented as

    phi1 = r + b / Lambda + c / Lambda^2,

where ``r`` is the actual SchedulePressure reference coefficient and ``b,c``
are the theorem-derived correction numerators.  A direct binary64 projection of
``phi1`` would erase the nonzero theorem-scale corrections before the next
nonlinear remainder evaluation.

This module lifts exactly one pinned ``AxisOperators.product`` seam to that
mixed-scale representation: the angular square ``phi1 * phi1``.  For each
radial/parameter jet it applies the official radial convolution and eta-Leibniz
sum and keeps the resulting powers ``Lambda^0`` through ``Lambda^-4`` separate.
No caller-supplied Lambda, coefficient table, surrogate angular field, or
product cutoff is accepted.

This is only one nonlinear primitive needed for ``naturalRemainder(x1)``.  It
does not materialize the complete remainder, ``x2``, convergence to the fixed
point, global weighted ``AxisSpace`` membership, or a paper-exact velocity.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math

from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState


_DECIMAL_PRECISION = 96
_MAX_POWER = 4


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _eta_in_window(value: float) -> float:
    eta = float(value)
    if not math.isfinite(eta) or not WINDOW_LEFT <= eta <= WINDOW_RIGHT:
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]")
    return eta


def _finite_decimal(value: Decimal, name: str) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError(f"{name} must be a finite Decimal")
    return value


@dataclass(frozen=True)
class MixedScaleFirstPicardAngularSquareCoefficientJet:
    """One coefficient jet of ``phi1 * phi1`` with all inverse-Lambda powers split.

    The represented value is

    ``p0 + p1/Lambda + p2/Lambda^2 + p3/Lambda^3 + p4/Lambda^4``.
    """

    reference: Decimal
    inverse_lambda_numerator: Decimal
    inverse_lambda_squared_numerator: Decimal
    inverse_lambda_cubed_numerator: Decimal
    inverse_lambda_fourth_numerator: Decimal
    Lambda: Decimal

    def __post_init__(self) -> None:
        _finite_decimal(self.reference, "reference")
        _finite_decimal(self.inverse_lambda_numerator, "inverse_lambda_numerator")
        _finite_decimal(
            self.inverse_lambda_squared_numerator,
            "inverse_lambda_squared_numerator",
        )
        _finite_decimal(
            self.inverse_lambda_cubed_numerator,
            "inverse_lambda_cubed_numerator",
        )
        _finite_decimal(
            self.inverse_lambda_fourth_numerator,
            "inverse_lambda_fourth_numerator",
        )
        _finite_decimal(self.Lambda, "Lambda")
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    def correction_terms_decimal(self) -> tuple[Decimal, Decimal, Decimal, Decimal]:
        """Return the four non-reference terms without summing them into ``p0``."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            inverse = +(Decimal(1) / self.Lambda)
            inverse2 = +(inverse * inverse)
            return (
                +(self.inverse_lambda_numerator * inverse),
                +(self.inverse_lambda_squared_numerator * inverse2),
                +(self.inverse_lambda_cubed_numerator * inverse2 * inverse),
                +(self.inverse_lambda_fourth_numerator * inverse2 * inverse2),
            )


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardAngularSquareState:
    """The exact mixed-scale pinned product ``phi1 * phi1`` on actual schedule data."""

    x1: ActualScheduleWideFirstPicardState

    def __post_init__(self) -> None:
        if not isinstance(self.x1, ActualScheduleWideFirstPicardState):
            raise TypeError("x1 must be ActualScheduleWideFirstPicardState")
        if not self.x1.picard_x1_materialized:
            raise ValueError("the genuine first Picard iterate must be materialized")
        if self.x1.fixed_point_materialized:
            raise ValueError("x1 must not be mislabeled as the final fixed point")

    @property
    def Lambda(self) -> Decimal:
        return self.x1.Lambda

    @property
    def epsilon(self) -> float:
        return self.x1.epsilon

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def mixed_scale_angular_square_materialized(self) -> bool:
        return True

    @property
    def natural_remainder_x1_materialized(self) -> bool:
        return False

    @property
    def fixed_point_materialized(self) -> bool:
        return False

    def jet(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> MixedScaleFirstPicardAngularSquareCoefficientJet:
        """Evaluate the pinned radial-convolution/eta-Leibniz square jet.

        If ``phi1 = sum_{p=0}^2 phi_p Lambda^-p`` coefficientwise, then the
        product law is applied before collecting equal inverse-Lambda powers:

        ``J[n,m] = sum_{i+j=n} sum_{k+l=m} binom(m,k) J_i,k J_j,l``.
        """

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)

        # Force the landed first-Picard guards at the requested jet before the
        # finite convolution.  Every summand below remains tied to this same x1.
        self.x1.jet_pair(n, m, eta)

        totals = [Decimal(0) for _ in range(_MAX_POWER + 1)]
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for i in range(n + 1):
                j = n - i
                for k in range(m + 1):
                    l = m - k
                    left, _ = self.x1.jet_pair(i, k, eta)
                    right, _ = self.x1.jet_pair(j, l, eta)
                    weight = Decimal(math.comb(m, k))
                    left_terms = (
                        left.reference,
                        left.inverse_lambda_numerator,
                        left.inverse_lambda_squared_numerator,
                    )
                    right_terms = (
                        right.reference,
                        right.inverse_lambda_numerator,
                        right.inverse_lambda_squared_numerator,
                    )
                    for left_power, left_value in enumerate(left_terms):
                        for right_power, right_value in enumerate(right_terms):
                            totals[left_power + right_power] += (
                                weight * left_value * right_value
                            )
            totals = [+value for value in totals]

        return MixedScaleFirstPicardAngularSquareCoefficientJet(
            reference=totals[0],
            inverse_lambda_numerator=totals[1],
            inverse_lambda_squared_numerator=totals[2],
            inverse_lambda_cubed_numerator=totals[3],
            inverse_lambda_fourth_numerator=totals[4],
            Lambda=self.Lambda,
        )


def wide_first_picard_angular_square_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardAngularSquareState:
    """Lift the pinned product to ``phi1 * phi1`` without scale collapse."""

    return ActualScheduleWideFirstPicardAngularSquareState(x1=x1)
