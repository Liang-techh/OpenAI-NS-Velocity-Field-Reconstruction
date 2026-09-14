"""Pinned angular quadratic remainder ``quad1(x1)`` at genuine first Picard state.

The official natural-axis remainder contains

    quad1 = - j2 (product (primitive phi) phi).

The landed ``ActualScheduleWideFirstPicardState`` keeps

    phi1 = r + b / Lambda + c / Lambda^2

as three coefficient families so the theorem-scale corrections are never
collapsed into the O(1) reference.  This module applies the pinned primitive,
product, and j2 row rules directly to that mixed-scale representation and keeps
the resulting Lambda^0 through Lambda^-4 families separate.

No caller-supplied coefficient table, cutoff, fitted field, or default-zero
extension is accepted.  The only zero rows below are the literal row-zero
semantics of ``primitive`` and ``j2`` in ``AxisOperators``.

This closes ``quad1(x1)`` only.  It does not materialize the remaining lin1,
slow1, lin2 recombination, complete naturalRemainder(x1), x2, a fixed point,
global AxisCoefficientSpace membership, NaturalProfileAssembly, or paper-exact
velocity.
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
class MixedScaleFirstPicardQuad1CoefficientJet:
    """One coefficient jet of ``quad1(x1)`` with inverse-Lambda powers split."""

    reference: Decimal
    inverse_lambda_numerator: Decimal
    inverse_lambda_squared_numerator: Decimal
    inverse_lambda_cubed_numerator: Decimal
    inverse_lambda_fourth_numerator: Decimal
    Lambda: Decimal

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            _finite_decimal(value, name)
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    def correction_terms_decimal(self) -> tuple[Decimal, Decimal, Decimal, Decimal]:
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
class ActualScheduleWideFirstPicardQuad1State:
    """Typed genuine ``quad1(x1)`` bound to one actual-schedule first Picard state."""

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
    def quad1_x1_materialized(self) -> bool:
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

    def _primitive_phi_terms(
        self, n: int, m: int, eta: float
    ) -> tuple[Decimal, Decimal, Decimal]:
        """Apply the pinned ``primitiveScale``/``Nat.pred`` row rule."""

        if n == 0:
            return (Decimal(0), Decimal(0), Decimal(0))
        predecessor = self._phi_terms(n - 1, m, eta)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            divisor = Decimal(n)
            return tuple(+(value / divisor) for value in predecessor)

    def jet(self, n: int, m: int, eta: float) -> MixedScaleFirstPicardQuad1CoefficientJet:
        """Evaluate ``-j2(product(primitive(phi1), phi1))`` exactly by pinned rows."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)

        # j2 is regularInverse with r=2 and therefore has a literal zero row.
        if n == 0:
            return MixedScaleFirstPicardQuad1CoefficientJet(
                reference=Decimal(0),
                inverse_lambda_numerator=Decimal(0),
                inverse_lambda_squared_numerator=Decimal(0),
                inverse_lambda_cubed_numerator=Decimal(0),
                inverse_lambda_fourth_numerator=Decimal(0),
                Lambda=self.Lambda,
            )

        product_row = n - 1
        totals = [Decimal(0) for _ in range(_MAX_POWER + 1)]
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for i in range(product_row + 1):
                j = product_row - i
                for k in range(m + 1):
                    l = m - k
                    left_terms = self._primitive_phi_terms(i, k, eta)
                    right_terms = self._phi_terms(j, l, eta)
                    weight = Decimal(math.comb(m, k))
                    for left_power, left_value in enumerate(left_terms):
                        for right_power, right_value in enumerate(right_terms):
                            totals[left_power + right_power] += (
                                weight * left_value * right_value
                            )

            # For output radial row n, j2 divides the predecessor by
            # radialDivisor(2,n-1) = n*(n+1).  quad1 contributes the minus sign.
            divisor = Decimal(n) * Decimal(n + 1)
            totals = [+(value / -divisor) for value in totals]

        return MixedScaleFirstPicardQuad1CoefficientJet(
            reference=totals[0],
            inverse_lambda_numerator=totals[1],
            inverse_lambda_squared_numerator=totals[2],
            inverse_lambda_cubed_numerator=totals[3],
            inverse_lambda_fourth_numerator=totals[4],
            Lambda=self.Lambda,
        )


def wide_first_picard_quad1_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardQuad1State:
    """Bind the pinned quadratic angular branch to one genuine first Picard state."""

    return ActualScheduleWideFirstPicardQuad1State(x1=x1)
