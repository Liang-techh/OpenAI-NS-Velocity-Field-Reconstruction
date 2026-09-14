"""Pinned angular-linear remainder ``lin1(x1)`` at the genuine first Picard state.

The pinned natural-axis contraction defines

    angularLinearCoefficient(d) = wStar + h*one - 2*h*product(eta,uStar)
    lin1(x) = j2(product(angularLinearCoefficient(d), phi))
              + dot2(wStar, phi)
              + param2(phi, hStar).

For the landed first Picard state, ``phi1`` is represented without scale collapse
as three ordinary coefficient families: ``Lambda^0``, ``Lambda^-1`` and
``Lambda^-2``.  This module propagates each family independently through the
three pinned linear branches.  The fixed AxisData fields and coefficient
operators come only from the same actual SchedulePressure first-Picard state.

No pressure scale is introduced here: the pinned ``lin1`` branch depends only on
``phi1`` and fixed AxisData.  No caller-supplied coefficient table, Lambda,
cutoff, derivative table, fitted state, or absent-row default is accepted.
This closes only ``lin1(x1)``; ``slow1(x1)``, ``lin2(x1)``, the complete
``naturalRemainder(x1)``, x2, a fixed point, global AxisCoefficientSpace
membership, NaturalProfileAssembly, and paper-exact velocity remain open.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math

from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState


_DECIMAL_PRECISION = 96
_MAX_POWER = 2


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


def _decimal_from_float(value: float, name: str) -> Decimal:
    value = float(value)
    if not math.isfinite(value):
        raise ArithmeticError(f"{name} must be finite")
    return Decimal.from_float(value)


@dataclass(frozen=True)
class MixedScaleFirstPicardLin1CoefficientJet:
    """One exact inverse-Lambda decomposition of the genuine ``lin1(x1)`` jet."""

    reference: Decimal
    inverse_lambda_numerator: Decimal
    inverse_lambda_squared_numerator: Decimal
    Lambda: Decimal

    def __post_init__(self) -> None:
        _finite_decimal(self.reference, "reference")
        _finite_decimal(self.inverse_lambda_numerator, "inverse_lambda_numerator")
        _finite_decimal(
            self.inverse_lambda_squared_numerator,
            "inverse_lambda_squared_numerator",
        )
        _finite_decimal(self.Lambda, "Lambda")
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    def factors_decimal(self) -> tuple[Decimal, Decimal, Decimal]:
        """Return factors multiplying ``1, Lambda^-1, Lambda^-2``."""

        return (
            self.reference,
            self.inverse_lambda_numerator,
            self.inverse_lambda_squared_numerator,
        )

    def correction_terms_decimal(self) -> tuple[Decimal, Decimal]:
        """Return the two non-reference terms without adding them to O(1)."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            inverse = +(Decimal(1) / self.Lambda)
            return (
                +(self.inverse_lambda_numerator * inverse),
                +(self.inverse_lambda_squared_numerator * inverse * inverse),
            )


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardLin1State:
    """Typed genuine ``lin1(x1)`` bound to one actual-schedule first Picard state."""

    x1: ActualScheduleWideFirstPicardState

    def __post_init__(self) -> None:
        if not isinstance(self.x1, ActualScheduleWideFirstPicardState):
            raise TypeError("x1 must be ActualScheduleWideFirstPicardState")
        if not self.x1.picard_x1_materialized:
            raise ValueError("the genuine first Picard iterate must be materialized")
        if self.x1.fixed_point_materialized:
            raise ValueError("x1 must not be mislabeled as the final fixed point")
        if self.axis_data.epsilon != self.x1.epsilon:
            raise ValueError("AxisData epsilon must match x1")
        if self.operators.epsilon != self.x1.epsilon:
            raise ValueError("coefficientOperators epsilon must match x1")

    @property
    def Lambda(self) -> Decimal:
        return self.x1.Lambda

    @property
    def epsilon(self) -> float:
        return self.x1.epsilon

    @property
    def axis_data(self):
        return self.x1.remainder.angular.axis_data

    @property
    def operators(self):
        return self.x1.remainder.angular.operators

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def lin1_x1_materialized(self) -> bool:
        return True

    @property
    def natural_remainder_x1_materialized(self) -> bool:
        return False

    @property
    def picard_x2_materialized(self) -> bool:
        return False

    @property
    def fixed_point_materialized(self) -> bool:
        return False

    def _phi_terms(self, n: int, m: int, eta: float) -> tuple[Decimal, Decimal, Decimal]:
        phi, _ = self.x1.jet_pair(n, m, eta)
        if phi.Lambda != self.Lambda:
            raise ValueError("first Picard angular jet Lambda mismatch")
        return (
            phi.reference,
            phi.inverse_lambda_numerator,
            phi.inverse_lambda_squared_numerator,
        )

    def _angular_linear_jet(self, n: int, m: int, eta: float) -> Decimal:
        """Replay the pinned fixed coefficient ``wStar+h-2h eta*uStar``."""

        d = self.axis_data
        eta_u_star = self.operators.product(d.eta, d.uStar)
        w_value = _decimal_from_float(d.wStar.jet(n, m, eta), "wStar jet")
        one_value = _decimal_from_float(d.one.jet(n, m, eta), "one jet")
        eta_u_value = _decimal_from_float(
            eta_u_star.jet(n, m, eta),
            "product(eta,uStar) jet",
        )
        h = _decimal_from_float(d.h, "h")
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(w_value + h * one_value - Decimal(2) * h * eta_u_value)

    def jet(self, n: int, m: int, eta: float) -> MixedScaleFirstPicardLin1CoefficientJet:
        """Evaluate all three pinned ``lin1(x1)`` branches without scale collapse."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        factors = [Decimal(0) for _ in range(_MAX_POWER + 1)]

        # j2, dot2 and param2 are all zero-datum regular inverses; their row 0
        # is exactly zero.  The input eta/order guards above remain active.
        if n == 0:
            return MixedScaleFirstPicardLin1CoefficientJet(
                reference=factors[0],
                inverse_lambda_numerator=factors[1],
                inverse_lambda_squared_numerator=factors[2],
                Lambda=self.Lambda,
            )

        d = self.axis_data
        radial_degree = n - 1
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            divisor = Decimal(n) * Decimal(n + 1)

            # j2(product(angularLinearCoefficient(d), phi1))
            for i in range(radial_degree + 1):
                j = radial_degree - i
                for k in range(m + 1):
                    l = m - k
                    weight = Decimal(math.comb(m, k))
                    fixed = self._angular_linear_jet(i, k, eta)
                    if fixed == 0:
                        continue
                    for power, phi_value in enumerate(self._phi_terms(j, l, eta)):
                        factors[power] += weight * fixed * phi_value

            # dot2(wStar, phi1): radial Euler derivative acts on phi1.
            for i in range(radial_degree + 1):
                j = radial_degree - i
                if j == 0:
                    continue
                for k in range(m + 1):
                    l = m - k
                    weight = Decimal(math.comb(m, k) * j)
                    fixed = _decimal_from_float(
                        d.wStar.jet(i, k, eta),
                        "wStar jet",
                    )
                    if fixed == 0:
                        continue
                    for power, phi_value in enumerate(self._phi_terms(j, l, eta)):
                        factors[power] += weight * fixed * phi_value

            # param2(phi1, hStar): eta derivative acts on the first argument.
            for i in range(radial_degree + 1):
                j = radial_degree - i
                for k in range(m + 1):
                    l = m - k
                    weight = Decimal(math.comb(m, k))
                    fixed = _decimal_from_float(
                        d.hStar.jet(j, l, eta),
                        "hStar jet",
                    )
                    if fixed == 0:
                        continue
                    for power, phi_value in enumerate(self._phi_terms(i, k + 1, eta)):
                        factors[power] += weight * phi_value * fixed

            factors = [+(value / divisor) for value in factors]

        return MixedScaleFirstPicardLin1CoefficientJet(
            reference=factors[0],
            inverse_lambda_numerator=factors[1],
            inverse_lambda_squared_numerator=factors[2],
            Lambda=self.Lambda,
        )


def wide_first_picard_lin1_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardLin1State:
    """Bind the pinned angular-linear branch to one genuine first Picard state."""

    return ActualScheduleWideFirstPicardLin1State(x1=x1)
