"""Mixed-scale angular ``lin1(x1)`` for the actual first Picard state.

The pinned natural-axis remainder contains the linear angular branch

    lin1(phi) = j2(angularLinear * phi)
                 + dot2(wStar, phi)
                 + param2(phi, hStar),

where the radial-zero coefficient family is

    angularLinear = wStar + h * one - 2 * h * (eta * uStar).

The first Picard angular coefficient is represented without collapsing its
theorem-scale pieces:

    phi1 = Phi0 + Phi1 / Lambda + Phi2 / Lambda^2.

Since ``lin1`` is linear in its angular input, this module applies the three
pinned operator terms to those three channels and returns the existing
``MixedScaleFirstPicardAngularCoefficientJet`` split type.  No outer
``inverseL`` or angular resolvent is applied here.

This remains a formal-structure adapter.  It does not materialize the other
``naturalRemainder(x1)`` branches, a second Picard iterate, a fixed point, or a
paper-exact velocity profile.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, localcontext
import math

from .axis_coefficient_data import (
    ActualScheduleAxisCoefficientData,
    actual_schedule_axis_coefficient_data,
)
from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_wide_first_picard import (
    ActualScheduleWideFirstPicardState,
    MixedScaleFirstPicardAngularCoefficientJet,
)


_DECIMAL_PRECISION = 96


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
        raise ValueError(f"{name} must be finite Decimal")
    return value


def _decimal_from_float(value: float, name: str) -> Decimal:
    value = float(value)
    if not math.isfinite(value):
        raise ArithmeticError(f"{name} must be finite")
    return Decimal.from_float(value)


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardLin1State:
    """Raw ``lin1(x1)`` on one genuine actual-schedule first Picard state."""

    x1: ActualScheduleWideFirstPicardState
    axis_data: ActualScheduleAxisCoefficientData = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.x1, ActualScheduleWideFirstPicardState):
            raise TypeError("x1 must be ActualScheduleWideFirstPicardState")
        if not self.x1.picard_x1_materialized:
            raise ValueError("the genuine first Picard iterate must be materialized")
        if self.x1.fixed_point_materialized:
            raise ValueError("x1 must not be mislabeled as the final fixed point")

        axis_data = actual_schedule_axis_coefficient_data(self.x1.reference)
        if axis_data.epsilon != self.x1.epsilon:
            raise ValueError("AxisData epsilon must match the first Picard state")
        _finite_decimal(self.x1.Lambda, "x1 Lambda")
        if self.x1.Lambda <= 0:
            raise ValueError("x1 Lambda must be positive")
        object.__setattr__(self, "axis_data", axis_data)

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
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def lin1_x1_materialized(self) -> bool:
        return True

    @property
    def angular_lin1_materialized(self) -> bool:
        return True

    @property
    def lin1_complete(self) -> bool:
        return True

    @property
    def natural_remainder_x1_materialized(self) -> bool:
        return False

    @property
    def fixed_point_materialized(self) -> bool:
        return False

    @property
    def fixed_point_convergence_certified(self) -> bool:
        return False

    def _row_zero_product(
        self,
        left,
        right,
        m: int,
        eta: float,
    ) -> Decimal:
        """Evaluate one radial-zero product jet in Decimal arithmetic."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            total = Decimal(0)
            for k in range(m + 1):
                total += (
                    Decimal(math.comb(m, k))
                    * _decimal_from_float(left.jet(0, k, eta), "AxisData left jet")
                    * _decimal_from_float(
                        right.jet(0, m - k, eta),
                        "AxisData right jet",
                    )
                )
            return +total

    def _angular_linear_jet(self, m: int, eta: float) -> Decimal:
        """Return ``wStar + h*one - 2*h*(eta*uStar)`` at eta order ``m``."""

        data = self.axis_data
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            w_star = _decimal_from_float(
                data.wStar.jet(0, m, eta),
                "AxisData.wStar jet",
            )
            one = _decimal_from_float(data.one.jet(0, m, eta), "AxisData.one jet")
            eta_u_star = self._row_zero_product(data.eta, data.uStar, m, eta)
            h = _decimal_from_float(data.h, "AxisData.h")
            return +(w_star + h * one - Decimal(2) * h * eta_u_star)

    def _source_component(
        self,
        component: str,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        """Read one split x1 angular component at a coefficient coordinate."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        angular, _ = self.x1.jet_pair(n, m, eta)
        if component == "reference":
            return _finite_decimal(angular.reference, "x1 angular reference jet")
        if component == "inverse_lambda":
            return _finite_decimal(
                angular.inverse_lambda_numerator,
                "x1 angular inverse-Lambda numerator",
            )
        if component == "inverse_lambda_squared":
            return _finite_decimal(
                angular.inverse_lambda_squared_numerator,
                "x1 angular inverse-Lambda-squared numerator",
            )
        raise KeyError(f"unknown x1 lin1 component {component!r}")

    def _linear_action(
        self,
        component: str,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        """Apply the three pinned linear operators to one split component."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)

        # Preserve the source validation performed by the pinned zero-datum
        # inverses even though every output component at row zero is zero.
        if n == 0:
            self._angular_linear_jet(m, eta)
            self.axis_data.wStar.jet(0, m, eta)
            self.axis_data.hStar.jet(0, m + 1, eta)
            self._source_component(component, 0, m + 1, eta)
            return Decimal(0)

        source_row = n - 1
        divisor = Decimal(n * (n + 1))
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            total = Decimal(0)

            # j2((angularLinear) * source).  AxisData is radial degree zero,
            # so the product source row is exactly n-1.
            for k in range(m + 1):
                total += (
                    Decimal(math.comb(m, k))
                    * self._angular_linear_jet(k, eta)
                    * self._source_component(component, source_row, m - k, eta)
                )

            # dot2(wStar, source): the second operand carries Euler factor
            # j = n-1 in the only nonzero radial convolution slot.
            radial_factor = Decimal(source_row)
            for k in range(m + 1):
                total += (
                    radial_factor
                    * Decimal(math.comb(m, k))
                    * _decimal_from_float(
                        self.axis_data.wStar.jet(0, k, eta),
                        "AxisData.wStar jet",
                    )
                    * self._source_component(component, source_row, m - k, eta)
                )

            # param2(source, hStar): parameter differentiation acts on the
            # source (the left operand), hence the k+1 source jet.
            for k in range(m + 1):
                total += (
                    Decimal(math.comb(m, k))
                    * self._source_component(component, source_row, k + 1, eta)
                    * _decimal_from_float(
                        self.axis_data.hStar.jet(0, m - k, eta),
                        "AxisData.hStar jet",
                    )
                )

            return +(total / divisor)

    def jet(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> MixedScaleFirstPicardAngularCoefficientJet:
        """Return one split coefficient jet of raw ``lin1(x1)``."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        return MixedScaleFirstPicardAngularCoefficientJet(
            reference=self._linear_action("reference", n, m, eta),
            inverse_lambda_numerator=self._linear_action(
                "inverse_lambda",
                n,
                m,
                eta,
            ),
            inverse_lambda_squared_numerator=self._linear_action(
                "inverse_lambda_squared",
                n,
                m,
                eta,
            ),
            Lambda=self.Lambda,
        )


def wide_first_picard_lin1_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardLin1State:
    """Lift raw ``lin1`` to the genuine theorem-scale first Picard state."""

    return ActualScheduleWideFirstPicardLin1State(x1=x1)


__all__ = [
    "ActualScheduleWideFirstPicardLin1State",
    "MixedScaleFirstPicardAngularCoefficientJet",
    "wide_first_picard_lin1_state",
]
