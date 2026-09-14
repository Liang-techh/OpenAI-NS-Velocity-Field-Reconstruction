"""Pinned axial-linear remainder ``lin2(x1)`` at the genuine first Picard state.

The pinned natural-axis contraction defines

    axialLinearCoefficient(d) =
        A*one - 4*A*product(eta,uStar) + product(d,uStarEta)
    lin2(x) = j1(product(axialLinearCoefficient(d), u))
              + dot1(wStar, u)
              + param1(u, hStar).

For the landed first Picard state, ``u1`` contains ordinary ``Lambda^0``,
``Lambda^-1`` and ``Lambda^-2`` coefficient families plus the genuine
pressure-derived ``a^2/Lambda`` family.  This module propagates all four
families independently through the three pinned linear branches.  The pressure
family remains normalized by ``a^2`` and is exposed with the actual wide
pressure amplitude log rather than collapsed to binary64.

No caller-supplied coefficient table, amplitude, Lambda, cutoff, fitted field,
or absent-row default is accepted.  This closes only ``lin2(x1)``; ``slow1(x1)``,
the complete ``naturalRemainder(x1)``, x2, fixed point, global
AxisCoefficientSpace membership, NaturalProfileAssembly, and paper-exact
velocity remain open.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState


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
        raise ValueError(f"{name} must be a finite Decimal")
    return value


def _decimal_from_float(value: float, name: str) -> Decimal:
    value = float(value)
    if not math.isfinite(value):
        raise ArithmeticError(f"{name} must be finite")
    return Decimal.from_float(value)


def _signed_log_pressure_term(
    numerator: Decimal,
    *,
    amplitude_log: Decimal,
    Lambda: Decimal,
) -> SignedLogCoefficientJet:
    """Restore ``a^2 * numerator / Lambda`` without narrowing ``a^2``."""

    _finite_decimal(numerator, "numerator")
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
    """One scale-separated coefficient jet of the genuine ``lin2(x1)``."""

    reference: Decimal
    inverse_lambda_numerator: Decimal
    inverse_lambda_squared_numerator: Decimal
    pressure_inverse_lambda_numerator: Decimal
    Lambda: Decimal
    amplitude_log: Decimal

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            _finite_decimal(value, name)
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    def ordinary_factors_decimal(self) -> tuple[Decimal, Decimal, Decimal]:
        """Return ordinary factors multiplying ``1, Lambda^-1, Lambda^-2``."""

        return (
            self.reference,
            self.inverse_lambda_numerator,
            self.inverse_lambda_squared_numerator,
        )

    def pressure_factor_decimal(self) -> Decimal:
        """Return the normalized factor multiplying ``a^2/Lambda``."""

        return self.pressure_inverse_lambda_numerator

    def ordinary_correction_terms_decimal(self) -> tuple[Decimal, Decimal]:
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            inverse = +(Decimal(1) / self.Lambda)
            return (
                +(self.inverse_lambda_numerator * inverse),
                +(self.inverse_lambda_squared_numerator * inverse * inverse),
            )

    def pressure_term_log(self) -> SignedLogCoefficientJet:
        return _signed_log_pressure_term(
            self.pressure_inverse_lambda_numerator,
            amplitude_log=self.amplitude_log,
            Lambda=self.Lambda,
        )


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardLin2State:
    """Typed genuine ``lin2(x1)`` bound to one actual-schedule first Picard state."""

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
        if self.x1.remainder.axial.wide_pressure.epsilon != self.x1.epsilon:
            raise ValueError("wide pressure epsilon must match x1")

    @property
    def Lambda(self) -> Decimal:
        return self.x1.Lambda

    @property
    def epsilon(self) -> float:
        return self.x1.epsilon

    @property
    def axis_data(self):
        return self.x1.remainder.axial.axis_data

    @property
    def operators(self):
        return self.x1.remainder.axial.operators

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def lin2_x1_materialized(self) -> bool:
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

    def _u_ordinary_terms(self, n: int, m: int, eta: float) -> tuple[Decimal, Decimal, Decimal]:
        _, u = self.x1.jet_pair(n, m, eta)
        if u.Lambda != self.Lambda:
            raise ValueError("first Picard axial jet Lambda mismatch")
        return (
            u.reference,
            u.inverse_lambda_numerator,
            u.inverse_lambda_squared_numerator,
        )

    def _u_pressure_inverse_lambda_numerator(self, n: int, m: int, eta: float) -> Decimal:
        """Return ``p`` in the genuine x1 pressure term ``a^2 * p / Lambda``."""

        factor = self.x1.remainder.axial.pressure_normalized_factor(n, m, eta)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(factor / Decimal(2))

    def _axial_linear_jet(self, n: int, m: int, eta: float) -> Decimal:
        """Replay the pinned fixed ``axialLinearCoefficient`` coefficient jet."""

        d = self.axis_data
        eta_u_star = self.operators.product(d.eta, d.uStar)
        d_u_star_eta = self.operators.product(d.d, d.uStarEta)
        one_value = _decimal_from_float(d.one.jet(n, m, eta), "one jet")
        eta_u_value = _decimal_from_float(
            eta_u_star.jet(n, m, eta),
            "product(eta,uStar) jet",
        )
        d_u_eta_value = _decimal_from_float(
            d_u_star_eta.jet(n, m, eta),
            "product(d,uStarEta) jet",
        )
        A = _decimal_from_float(d.A, "A")
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(A * one_value - Decimal(4) * A * eta_u_value + d_u_eta_value)

    def jet(self, n: int, m: int, eta: float) -> MixedScaleFirstPicardLin2CoefficientJet:
        """Evaluate all three pinned ``lin2(x1)`` branches without scale collapse."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        amplitude = self.x1.remainder.axial.wide_pressure.amplitude
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            amplitude_log = +amplitude.log_amplitude(eta)

        ordinary = [Decimal(0), Decimal(0), Decimal(0)]
        pressure = Decimal(0)

        # j1, dot1 and param1 are zero-datum r=1 radial inverses.
        if n == 0:
            return MixedScaleFirstPicardLin2CoefficientJet(
                reference=ordinary[0],
                inverse_lambda_numerator=ordinary[1],
                inverse_lambda_squared_numerator=ordinary[2],
                pressure_inverse_lambda_numerator=pressure,
                Lambda=self.Lambda,
                amplitude_log=amplitude_log,
            )

        d = self.axis_data
        radial_degree = n - 1
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION

            # j1(product(axialLinearCoefficient(d), u1))
            for i in range(radial_degree + 1):
                j = radial_degree - i
                for k in range(m + 1):
                    l = m - k
                    weight = Decimal(math.comb(m, k))
                    fixed = self._axial_linear_jet(i, k, eta)
                    if fixed == 0:
                        continue
                    for power, u_value in enumerate(self._u_ordinary_terms(j, l, eta)):
                        ordinary[power] += weight * fixed * u_value
                    pressure += (
                        weight
                        * fixed
                        * self._u_pressure_inverse_lambda_numerator(j, l, eta)
                    )

            # dot1(wStar, u1): radial Euler derivative acts on the second argument.
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
                    for power, u_value in enumerate(self._u_ordinary_terms(j, l, eta)):
                        ordinary[power] += weight * fixed * u_value
                    pressure += (
                        weight
                        * fixed
                        * self._u_pressure_inverse_lambda_numerator(j, l, eta)
                    )

            # param1(u1, hStar): eta derivative acts on the first argument.
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
                    for power, u_value in enumerate(self._u_ordinary_terms(i, k + 1, eta)):
                        ordinary[power] += weight * u_value * fixed
                    pressure += (
                        weight
                        * self._u_pressure_inverse_lambda_numerator(i, k + 1, eta)
                        * fixed
                    )

            # radialDivisor(1,n-1) = n*n.
            divisor = Decimal(n) * Decimal(n)
            ordinary = [+(value / divisor) for value in ordinary]
            pressure = +(pressure / divisor)

        return MixedScaleFirstPicardLin2CoefficientJet(
            reference=ordinary[0],
            inverse_lambda_numerator=ordinary[1],
            inverse_lambda_squared_numerator=ordinary[2],
            pressure_inverse_lambda_numerator=pressure,
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )


def wide_first_picard_lin2_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardLin2State:
    """Bind the pinned axial-linear branch to one genuine first Picard state."""

    return ActualScheduleWideFirstPicardLin2State(x1=x1)
