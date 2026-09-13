"""Pinned ``dot1`` average/u branch of ``slow2`` at the genuine first Picard state.

One of the remaining terms of ``AxisContraction.naturalRemainder.slow2`` is

    dot1 (product averageCoefficient (average u)) u,

where ``averageCoefficient = 2 * D * eta`` and ``dot1`` is the pinned
``inverseDotProduct 1 0``.  The landed first Picard state keeps

    u1 = U0 + U1/Lambda + U2/Lambda^2 + a^2 P/Lambda

without collapsing the theorem scales.  This module propagates exactly that
split state through ``average``, multiplication by ``2*D*eta``, and ``dot1``.
The result remains separated into

* ordinary powers Lambda^0 through Lambda^-4;
* pressure-linear ``a^2`` powers Lambda^-1 through Lambda^-3; and
* the pressure-square ``a^4`` power Lambda^-2.

For output row ``n>=1``, pinned ``dot1`` reads pairs ``i+j=n-1`` with factor
``j / radialDivisor(1,n-1) = j/n^2`` and the usual eta-Leibniz binomial weight.
Row zero is exactly zero.  No caller-supplied D, Lambda, amplitude, pressure or
coefficient table, derivative table, or radial cutoff is accepted.

This is one mixed-scale ``slow2(x1)`` branch only.  It is not complete
``slow2(x1)``, ``naturalRemainder(x1)``, ``x2``, a fixed-point convergence
certificate, global weighted ``AxisSpace`` membership, or paper-exact velocity.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_data import actual_schedule_axis_coefficient_data
from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState


_DECIMAL_PRECISION = 96
_TWO = Decimal(2)
_MAX_ORDINARY_POWER = 4


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


def _signed_log_term(
    numerator: Decimal,
    *,
    amplitude_log: Decimal,
    amplitude_power: int,
    Lambda: Decimal,
    inverse_lambda_power: int,
) -> SignedLogCoefficientJet:
    _finite_decimal(numerator, "numerator")
    _finite_decimal(amplitude_log, "amplitude_log")
    _finite_decimal(Lambda, "Lambda")
    if Lambda <= 0:
        raise ValueError("Lambda must be positive")
    if amplitude_power <= 0 or inverse_lambda_power < 0:
        raise ValueError("invalid signed-log scale powers")
    if numerator == 0:
        return SignedLogCoefficientJet.zero()

    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        log_scale = +(Decimal(amplitude_power) * amplitude_log)
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
class MixedScaleFirstPicardSlow2AverageDotCoefficientJet:
    """One coefficient jet of ``dot1((2*D*eta)*average(u1), u1)``."""

    ordinary_reference: Decimal
    ordinary_inverse_lambda_numerator: Decimal
    ordinary_inverse_lambda_squared_numerator: Decimal
    ordinary_inverse_lambda_cubed_numerator: Decimal
    ordinary_inverse_lambda_fourth_numerator: Decimal
    pressure_linear_inverse_lambda_numerator: Decimal
    pressure_linear_inverse_lambda_squared_numerator: Decimal
    pressure_linear_inverse_lambda_cubed_numerator: Decimal
    pressure_square_inverse_lambda_squared_numerator: Decimal
    Lambda: Decimal
    amplitude_log: Decimal

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            _finite_decimal(value, name)
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    def ordinary_correction_terms_decimal(
        self,
    ) -> tuple[Decimal, Decimal, Decimal, Decimal]:
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            inverse = +(Decimal(1) / self.Lambda)
            inverse2 = +(inverse * inverse)
            return (
                +(self.ordinary_inverse_lambda_numerator * inverse),
                +(self.ordinary_inverse_lambda_squared_numerator * inverse2),
                +(self.ordinary_inverse_lambda_cubed_numerator * inverse2 * inverse),
                +(self.ordinary_inverse_lambda_fourth_numerator * inverse2 * inverse2),
            )

    def pressure_linear_terms_log(
        self,
    ) -> tuple[SignedLogCoefficientJet, SignedLogCoefficientJet, SignedLogCoefficientJet]:
        return tuple(
            _signed_log_term(
                numerator,
                amplitude_log=self.amplitude_log,
                amplitude_power=2,
                Lambda=self.Lambda,
                inverse_lambda_power=power,
            )
            for power, numerator in enumerate(
                (
                    self.pressure_linear_inverse_lambda_numerator,
                    self.pressure_linear_inverse_lambda_squared_numerator,
                    self.pressure_linear_inverse_lambda_cubed_numerator,
                ),
                start=1,
            )
        )

    def pressure_square_term_log(self) -> SignedLogCoefficientJet:
        return _signed_log_term(
            self.pressure_square_inverse_lambda_squared_numerator,
            amplitude_log=self.amplitude_log,
            amplitude_power=4,
            Lambda=self.Lambda,
            inverse_lambda_power=2,
        )


_FIELDS = (
    "ordinary_reference",
    "ordinary_inverse_lambda_numerator",
    "ordinary_inverse_lambda_squared_numerator",
    "ordinary_inverse_lambda_cubed_numerator",
    "ordinary_inverse_lambda_fourth_numerator",
    "pressure_linear_inverse_lambda_numerator",
    "pressure_linear_inverse_lambda_squared_numerator",
    "pressure_linear_inverse_lambda_cubed_numerator",
    "pressure_square_inverse_lambda_squared_numerator",
)


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardSlow2AverageDotState:
    """Pinned average/u ``dot1`` branch on actual theorem-selected ``x1``."""

    x1: ActualScheduleWideFirstPicardState

    def __post_init__(self) -> None:
        if not isinstance(self.x1, ActualScheduleWideFirstPicardState):
            raise TypeError("x1 must be ActualScheduleWideFirstPicardState")
        if not self.x1.picard_x1_materialized:
            raise ValueError("the genuine first Picard iterate must be materialized")
        if self.x1.fixed_point_materialized:
            raise ValueError("x1 must not be mislabeled as the final fixed point")
        if self.x1.Lambda != self.x1.remainder.axial.wide_pressure.amplitude.Lambda:
            raise ValueError("x1 Lambda must match the actual theorem-selected amplitude scale")

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
    def slow2_average_dot_branch_materialized(self) -> bool:
        return True

    @property
    def slow2_complete(self) -> bool:
        return False

    @property
    def natural_remainder_x1_materialized(self) -> bool:
        return False

    @property
    def fixed_point_materialized(self) -> bool:
        return False

    def _D_decimal(self) -> Decimal:
        data = actual_schedule_axis_coefficient_data(self.x1.reference)
        if data.epsilon != self.epsilon:
            raise ValueError("AxisData epsilon must match the first Picard state")
        return Decimal.from_float(data.D)

    def _u_ordinary(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[Decimal, Decimal, Decimal]:
        _, axial = self.x1.jet_pair(n, m, eta)
        return (
            axial.reference,
            axial.inverse_lambda_numerator,
            axial.inverse_lambda_squared_numerator,
        )

    def _u_pressure_numerator(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        """Return ``P`` where the pressure part of ``u1`` is ``a^2 P/Lambda``."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(
                self.x1.remainder.axial.pressure_normalized_factor(n, m, eta)
                / _TWO
            )

    def _average_ordinary(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[Decimal, Decimal, Decimal]:
        divisor = Decimal(n + 1)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return tuple(+(value / divisor) for value in self._u_ordinary(n, m, eta))

    def _average_pressure_numerator(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        divisor = Decimal(n + 1)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(self._u_pressure_numerator(n, m, eta) / divisor)

    def _left_ordinary(
        self,
        n: int,
        m: int,
        eta: float,
        D: Decimal,
    ) -> tuple[Decimal, Decimal, Decimal]:
        """Return scale numerators of ``(2*D*eta)*average(u1)``."""

        current = self._average_ordinary(n, m, eta)
        previous = self._average_ordinary(n, m - 1, eta) if m > 0 else None
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            q0 = +(Decimal(2) * D * Decimal.from_float(eta))
            q1_weight = +(Decimal(m) * Decimal(2) * D)
            values = []
            for power in range(3):
                value = q0 * current[power]
                if previous is not None:
                    value += q1_weight * previous[power]
                values.append(+value)
        return tuple(values)

    def _left_pressure_numerator(
        self,
        n: int,
        m: int,
        eta: float,
        D: Decimal,
    ) -> Decimal:
        current = self._average_pressure_numerator(n, m, eta)
        previous = self._average_pressure_numerator(n, m - 1, eta) if m > 0 else None
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            value = Decimal(2) * D * Decimal.from_float(eta) * current
            if previous is not None:
                value += Decimal(m) * Decimal(2) * D * previous
            return +value

    @staticmethod
    def _zero(
        *,
        Lambda: Decimal,
        amplitude_log: Decimal,
    ) -> MixedScaleFirstPicardSlow2AverageDotCoefficientJet:
        return MixedScaleFirstPicardSlow2AverageDotCoefficientJet(
            **{name: Decimal(0) for name in _FIELDS},
            Lambda=Lambda,
            amplitude_log=amplitude_log,
        )

    def jet(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> MixedScaleFirstPicardSlow2AverageDotCoefficientJet:
        """Evaluate the exact pinned ``dot1((2*D*eta)*average(u1),u1)`` jet."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)

        # Exercise the genuine theorem-data chain before applying dot1's exact
        # row-zero rule as well as on positive rows.
        self.x1.jet_pair(0 if n == 0 else n - 1, m, eta)
        amplitude_log = self.x1.remainder.axial.wide_pressure.amplitude.log_amplitude(eta)
        if n == 0:
            return self._zero(Lambda=self.Lambda, amplitude_log=amplitude_log)

        D = self._D_decimal()
        ordinary = [Decimal(0) for _ in range(_MAX_ORDINARY_POWER + 1)]
        pressure_linear = [Decimal(0), Decimal(0), Decimal(0), Decimal(0)]
        pressure_square_power2 = Decimal(0)
        source_degree = n - 1
        divisor = Decimal(n * n)

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for i in range(source_degree + 1):
                j = source_degree - i
                if j == 0:
                    continue
                radial_weight = Decimal(j) / divisor
                for k in range(m + 1):
                    l = m - k
                    weight = Decimal(math.comb(m, k)) * radial_weight
                    left_ordinary = self._left_ordinary(i, k, eta, D)
                    right_ordinary = self._u_ordinary(j, l, eta)

                    for left_power, left_value in enumerate(left_ordinary):
                        for right_power, right_value in enumerate(right_ordinary):
                            ordinary[left_power + right_power] += (
                                weight * left_value * right_value
                            )

                    left_pressure = self._left_pressure_numerator(i, k, eta, D)
                    right_pressure = self._u_pressure_numerator(j, l, eta)
                    for left_power, left_value in enumerate(left_ordinary):
                        pressure_linear[left_power + 1] += (
                            weight * left_value * right_pressure
                        )
                    for right_power, right_value in enumerate(right_ordinary):
                        pressure_linear[right_power + 1] += (
                            weight * left_pressure * right_value
                        )
                    pressure_square_power2 += weight * left_pressure * right_pressure

            ordinary = [+value for value in ordinary]
            pressure_linear = [+value for value in pressure_linear]
            pressure_square_power2 = +pressure_square_power2

        return MixedScaleFirstPicardSlow2AverageDotCoefficientJet(
            ordinary_reference=ordinary[0],
            ordinary_inverse_lambda_numerator=ordinary[1],
            ordinary_inverse_lambda_squared_numerator=ordinary[2],
            ordinary_inverse_lambda_cubed_numerator=ordinary[3],
            ordinary_inverse_lambda_fourth_numerator=ordinary[4],
            pressure_linear_inverse_lambda_numerator=pressure_linear[1],
            pressure_linear_inverse_lambda_squared_numerator=pressure_linear[2],
            pressure_linear_inverse_lambda_cubed_numerator=pressure_linear[3],
            pressure_square_inverse_lambda_squared_numerator=pressure_square_power2,
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )


def wide_first_picard_slow2_average_dot_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardSlow2AverageDotState:
    """Lift the pinned average/u ``dot1`` branch without theorem-scale collapse."""

    return ActualScheduleWideFirstPicardSlow2AverageDotState(x1=x1)
