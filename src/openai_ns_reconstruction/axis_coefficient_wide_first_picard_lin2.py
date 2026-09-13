"""Mixed-scale axial ``lin2(x1)`` for the actual first Picard state.

The pinned natural-axis remainder contains the linear axial branch

    lin2(u) = j1(axialLinear * u) + dot1(wStar, u) + param1(u, hStar),

where ``axialLinear`` is the radial-zero coefficient family

    A * one - 4 * A * (eta * uStar) + d * uStarEta.

The landed first Picard axial field is represented without collapsing its
theorem-scale pieces:

    u1 = U0 + U1 / Lambda + U2 / Lambda^2 + a^2 * P1 / Lambda.

This module applies the linear ``lin2`` branch to those four pieces and keeps
the result split as ``Decimal`` ordinary numerators plus one normalized
signed-log pressure numerator.  ``P1`` is taken from the actual wide pressure
state's *full pressure derivative jet*; it is not obtained by differentiating a
zeroth-order normalized pressure value.  No outer ``inverseL`` is applied here;
that multiplication belongs to the later axial natural-remainder assembly.

This remains a formal-structure adapter.  It does not materialize the other
``naturalRemainder(x1)`` branches, a second Picard iterate, a fixed point, or a
paper-exact velocity profile.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, localcontext
import math

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_data import (
    ActualScheduleAxisCoefficientData,
    actual_schedule_axis_coefficient_data,
)
from .axis_coefficient_reference_state import (
    WINDOW_LEFT,
    WINDOW_RIGHT,
)
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState


_DECIMAL_PRECISION = 96
_TWO = Decimal(2)


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


def _signed_log_term(
    numerator: Decimal,
    *,
    amplitude_log: Decimal,
    Lambda: Decimal,
) -> SignedLogCoefficientJet:
    """Encode ``a^2 * numerator / Lambda`` without narrowing its magnitude."""

    _finite_decimal(numerator, "pressure numerator")
    _finite_decimal(amplitude_log, "amplitude_log")
    _finite_decimal(Lambda, "Lambda")
    if Lambda <= 0:
        raise ValueError("Lambda must be positive")
    if numerator == 0:
        return SignedLogCoefficientJet.zero()

    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        log_scale = +(Decimal(2) * amplitude_log)
        log_factor = +(abs(numerator).ln() - Lambda.ln())
    return SignedLogCoefficientJet(
        sign=1 if numerator > 0 else -1,
        log_scale=log_scale,
        log_factor=log_factor,
    )


@dataclass(frozen=True)
class MixedScaleFirstPicardLin2CoefficientJet:
    """One theorem-scale coefficient jet of the raw ``lin2(x1)`` branch.

    The represented value is

    ``ordinary_reference + ordinary_inverse_lambda_numerator / Lambda``
    ``+ ordinary_inverse_lambda_squared_numerator / Lambda^2``
    ``+ a^2 * pressure_linear_inverse_lambda_numerator / Lambda``.

    The pressure term is exposed only through ``pressure_linear_term_log`` so a
    nonzero term cannot be erased by a binary64 projection.
    """

    ordinary_reference: Decimal
    ordinary_inverse_lambda_numerator: Decimal
    ordinary_inverse_lambda_squared_numerator: Decimal
    pressure_linear_inverse_lambda_numerator: Decimal
    Lambda: Decimal
    amplitude_log: Decimal

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            _finite_decimal(value, name)
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    def ordinary_correction_terms_decimal(self) -> tuple[Decimal, Decimal]:
        """Return the separate ``Lambda^-1`` and ``Lambda^-2`` terms."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            inverse = +(Decimal(1) / self.Lambda)
            return (
                +(self.ordinary_inverse_lambda_numerator * inverse),
                +(
                    self.ordinary_inverse_lambda_squared_numerator
                    * inverse
                    * inverse
                ),
            )

    def pressure_linear_term_log(self) -> SignedLogCoefficientJet:
        """Return the normalized pressure term ``a^2 P1 / Lambda``."""

        return _signed_log_term(
            self.pressure_linear_inverse_lambda_numerator,
            amplitude_log=self.amplitude_log,
            Lambda=self.Lambda,
        )

    def pressure_linear_terms_log(self) -> tuple[SignedLogCoefficientJet]:
        """Tuple alias matching the multi-term slow-branch API."""

        return (self.pressure_linear_term_log(),)


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardLin2State:
    """Raw ``lin2(x1)`` on one genuine actual-schedule first Picard state."""

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
        if self.x1.Lambda != self.x1.remainder.axial.wide_pressure.amplitude.Lambda:
            raise ValueError("x1 Lambda must match the actual pressure amplitude scale")
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
    def lin2_x1_materialized(self) -> bool:
        return True

    @property
    def axial_lin2_materialized(self) -> bool:
        return True

    @property
    def lin2_complete(self) -> bool:
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

    @property
    def _A(self) -> Decimal:
        return _decimal_from_float(self.axis_data.A, "AxisData.A")

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

    def _axial_linear_jet(self, m: int, eta: float) -> Decimal:
        """Return the radial-zero ``A*one - 4*A*eta*uStar + d*uStarEta`` jet."""

        data = self.axis_data
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            one = _decimal_from_float(data.one.jet(0, m, eta), "AxisData.one jet")
            eta_u_star = self._row_zero_product(data.eta, data.uStar, m, eta)
            d_u_star_eta = self._row_zero_product(data.d, data.uStarEta, m, eta)
            return +(
                self._A * one
                - Decimal(4) * self._A * eta_u_star
                + d_u_star_eta
            )

    def _source_component(
        self,
        component: str,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        """Read one split x1 source component at a coefficient coordinate."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        if component == "pressure":
            # This is the full pressure derivative jet normalized by a^2.  The
            # wide-pressure API already includes the eta derivatives of the
            # full pressure, including amplitude derivatives; do not
            # differentiate a zeroth-order normalized value here.
            factor = self.x1.remainder.axial.pressure_normalized_factor(n, m, eta)
            with localcontext() as ctx:
                ctx.prec = _DECIMAL_PRECISION
                return +(factor / _TWO)

        _, axial = self.x1.jet_pair(n, m, eta)
        if component == "reference":
            return _finite_decimal(axial.reference, "x1 reference jet")
        if component == "inverse_lambda":
            return _finite_decimal(
                axial.inverse_lambda_numerator,
                "x1 inverse-Lambda numerator",
            )
        if component == "inverse_lambda_squared":
            return _finite_decimal(
                axial.inverse_lambda_squared_numerator,
                "x1 inverse-Lambda-squared numerator",
            )
        raise KeyError(f"unknown x1 lin2 component {component!r}")

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
            self._axial_linear_jet(m, eta)
            self.axis_data.wStar.jet(0, m, eta)
            self.axis_data.hStar.jet(0, m + 1, eta)
            self._source_component(component, 0, m + 1, eta)
            return Decimal(0)

        source_row = n - 1
        divisor = Decimal(n * n)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            total = Decimal(0)

            # j1((axialLinear) * source).  AxisData is radial degree zero, so
            # the product source row is exactly n-1.
            for k in range(m + 1):
                total += (
                    Decimal(math.comb(m, k))
                    * self._axial_linear_jet(k, eta)
                    * self._source_component(component, source_row, m - k, eta)
                )

            # dot1(wStar, source): the second operand carries Euler factor
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

            # param1(source, hStar): parameter differentiation acts on the
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
    ) -> MixedScaleFirstPicardLin2CoefficientJet:
        """Return one split coefficient jet of raw ``lin2(x1)``."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        amplitude_log = _finite_decimal(
            self.x1.remainder.axial.wide_pressure.amplitude.log_amplitude(eta),
            "actual amplitude log",
        )
        return MixedScaleFirstPicardLin2CoefficientJet(
            ordinary_reference=self._linear_action("reference", n, m, eta),
            ordinary_inverse_lambda_numerator=self._linear_action(
                "inverse_lambda",
                n,
                m,
                eta,
            ),
            ordinary_inverse_lambda_squared_numerator=self._linear_action(
                "inverse_lambda_squared",
                n,
                m,
                eta,
            ),
            pressure_linear_inverse_lambda_numerator=self._linear_action(
                "pressure",
                n,
                m,
                eta,
            ),
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )


def wide_first_picard_lin2_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardLin2State:
    """Lift raw ``lin2`` to the genuine theorem-scale first Picard state."""

    return ActualScheduleWideFirstPicardLin2State(x1=x1)


__all__ = [
    "ActualScheduleWideFirstPicardLin2State",
    "MixedScaleFirstPicardLin2CoefficientJet",
    "wide_first_picard_lin2_state",
]
