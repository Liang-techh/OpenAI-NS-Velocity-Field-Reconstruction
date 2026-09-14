"""Pinned angular quadratic remainder ``quad1(x1)`` at genuine first Picard state.

The pinned natural-axis contraction defines

    angularQuadraticCoefficient(d) = product(d.d, d.normalizedGradient)
    quad1(x) = j2(product(product(angularQuadraticCoefficient(d), u), phi)).

For the landed first Picard state, ``phi1`` has ordinary Lambda powers 0,-1,-2,
while ``u1`` has those same ordinary powers plus the genuine pressure-derived
``a^2/(2*Lambda)`` contribution.  This module propagates *both* pieces through
the exact finite radial convolution, eta-Leibniz rule, and pinned ``j2`` row
map.  The ordinary result is kept as Lambda^0 through Lambda^-4; the pressure
part is kept as ``a^2`` times Lambda^-1 through Lambda^-3 normalized factors.

No caller-supplied coefficient table, amplitude, Lambda, cutoff, fitted field,
or default-zero extension is accepted.  This closes only ``quad1(x1)``; the
complete ``naturalRemainder(x1)``, x2, fixed point, global AxisCoefficientSpace
membership, NaturalProfileAssembly, and paper-exact velocity remain open.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState


_DECIMAL_PRECISION = 96
_ORDINARY_MAX_POWER = 4
_PRESSURE_MAX_POWER = 3


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
    inverse_lambda_power: int,
) -> SignedLogCoefficientJet:
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
        log_factor = +(
            abs(numerator).ln()
            - Decimal(inverse_lambda_power) * Lambda.ln()
        )
    return SignedLogCoefficientJet(
        sign=1 if numerator > 0 else -1,
        log_scale=log_scale,
        log_factor=log_factor,
    )


@dataclass(frozen=True)
class MixedScaleFirstPicardQuad1CoefficientJet:
    """One exact scale decomposition of the genuine ``quad1(x1)`` coefficient jet."""

    reference: Decimal
    inverse_lambda_numerator: Decimal
    inverse_lambda_squared_numerator: Decimal
    inverse_lambda_cubed_numerator: Decimal
    inverse_lambda_fourth_numerator: Decimal
    pressure_inverse_lambda_numerator: Decimal
    pressure_inverse_lambda_squared_numerator: Decimal
    pressure_inverse_lambda_cubed_numerator: Decimal
    Lambda: Decimal
    amplitude_log: Decimal

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            _finite_decimal(value, name)
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    def ordinary_factors_decimal(self) -> tuple[Decimal, Decimal, Decimal, Decimal, Decimal]:
        return (
            self.reference,
            self.inverse_lambda_numerator,
            self.inverse_lambda_squared_numerator,
            self.inverse_lambda_cubed_numerator,
            self.inverse_lambda_fourth_numerator,
        )

    def pressure_factors_decimal(self) -> tuple[Decimal, Decimal, Decimal]:
        """Return factors of ``a^2/Lambda^p`` for p=1,2,3."""

        return (
            self.pressure_inverse_lambda_numerator,
            self.pressure_inverse_lambda_squared_numerator,
            self.pressure_inverse_lambda_cubed_numerator,
        )

    def ordinary_correction_terms_decimal(self) -> tuple[Decimal, Decimal, Decimal, Decimal]:
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

    def pressure_terms_log(self) -> tuple[SignedLogCoefficientJet, SignedLogCoefficientJet, SignedLogCoefficientJet]:
        """Restore the three pressure-derived scales without narrowing ``a^2``."""

        return tuple(
            _signed_log_term(
                numerator,
                amplitude_log=self.amplitude_log,
                Lambda=self.Lambda,
                inverse_lambda_power=power,
            )
            for power, numerator in enumerate(self.pressure_factors_decimal(), start=1)
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
        if self.x1.remainder.angular.axis_data.epsilon != self.x1.epsilon:
            raise ValueError("AxisData epsilon must match x1")
        if self.x1.remainder.angular.operators.epsilon != self.x1.epsilon:
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

    @property
    def angular_quadratic(self):
        """Pinned ``product(d.d, d.normalizedGradient)`` from actual AxisData."""

        O = self.x1.remainder.angular.operators
        d = self.x1.remainder.angular.axis_data
        return O.product(d.d, d.normalizedGradient)

    def _phi_terms(self, n: int, m: int, eta: float) -> tuple[Decimal, Decimal, Decimal]:
        phi, _ = self.x1.jet_pair(n, m, eta)
        if phi.Lambda != self.Lambda:
            raise ValueError("first Picard angular jet Lambda mismatch")
        return (
            phi.reference,
            phi.inverse_lambda_numerator,
            phi.inverse_lambda_squared_numerator,
        )

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
        """Return the factor ``p`` in the x1 term ``a^2 * p / Lambda``.

        ``remainder.axial.pressure_normalized_factor`` is the actual derivative
        of the wide remainder pressure divided by ``a^2``.  The outer Picard
        map contributes exactly ``1/(2*Lambda)``.
        """

        factor = self.x1.remainder.axial.pressure_normalized_factor(n, m, eta)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(factor / Decimal(2))

    def jet(self, n: int, m: int, eta: float) -> MixedScaleFirstPicardQuad1CoefficientJet:
        """Evaluate pinned ``j2(q * u1 * phi1)`` with every genuine x1 scale retained."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        amplitude = self.x1.remainder.axial.wide_pressure.amplitude
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            amplitude_log = +amplitude.log_amplitude(eta)

        ordinary = [Decimal(0) for _ in range(_ORDINARY_MAX_POWER + 1)]
        pressure = [Decimal(0) for _ in range(_PRESSURE_MAX_POWER + 1)]
        if n == 0:
            return MixedScaleFirstPicardQuad1CoefficientJet(
                reference=ordinary[0],
                inverse_lambda_numerator=ordinary[1],
                inverse_lambda_squared_numerator=ordinary[2],
                inverse_lambda_cubed_numerator=ordinary[3],
                inverse_lambda_fourth_numerator=ordinary[4],
                pressure_inverse_lambda_numerator=pressure[1],
                pressure_inverse_lambda_squared_numerator=pressure[2],
                pressure_inverse_lambda_cubed_numerator=pressure[3],
                Lambda=self.Lambda,
                amplitude_log=amplitude_log,
            )

        input_row = n - 1
        q = self.angular_quadratic
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for q_row in range(input_row + 1):
                for u_row in range(input_row - q_row + 1):
                    phi_row = input_row - q_row - u_row
                    for q_eta in range(m + 1):
                        remaining = m - q_eta
                        for u_eta in range(remaining + 1):
                            phi_eta = remaining - u_eta
                            weight = Decimal(
                                math.comb(m, q_eta) * math.comb(remaining, u_eta)
                            )
                            q_value = _decimal_from_float(
                                q.jet(q_row, q_eta, eta),
                                "angularQuadraticCoefficient jet",
                            )
                            if q_value == 0:
                                continue
                            u_terms = self._u_ordinary_terms(u_row, u_eta, eta)
                            phi_terms = self._phi_terms(phi_row, phi_eta, eta)
                            for u_power, u_value in enumerate(u_terms):
                                for phi_power, phi_value in enumerate(phi_terms):
                                    ordinary[u_power + phi_power] += (
                                        weight * q_value * u_value * phi_value
                                    )

                            wide_u = self._u_pressure_inverse_lambda_numerator(
                                u_row, u_eta, eta
                            )
                            for phi_power, phi_value in enumerate(phi_terms):
                                pressure[1 + phi_power] += (
                                    weight * q_value * wide_u * phi_value
                                )

            divisor = Decimal(n) * Decimal(n + 1)
            ordinary = [+(value / divisor) for value in ordinary]
            pressure = [+(value / divisor) for value in pressure]

        return MixedScaleFirstPicardQuad1CoefficientJet(
            reference=ordinary[0],
            inverse_lambda_numerator=ordinary[1],
            inverse_lambda_squared_numerator=ordinary[2],
            inverse_lambda_cubed_numerator=ordinary[3],
            inverse_lambda_fourth_numerator=ordinary[4],
            pressure_inverse_lambda_numerator=pressure[1],
            pressure_inverse_lambda_squared_numerator=pressure[2],
            pressure_inverse_lambda_cubed_numerator=pressure[3],
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )


def wide_first_picard_quad1_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardQuad1State:
    """Bind the pinned quadratic angular branch to one genuine first Picard state."""

    return ActualScheduleWideFirstPicardQuad1State(x1=x1)
