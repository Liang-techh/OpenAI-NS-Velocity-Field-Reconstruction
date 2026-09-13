"""Pinned ``-param1(u1, d*u1)`` branch of ``slow2`` at first Picard state.

The pinned ``AxisContraction.naturalRemainder.slow2`` ends with

    - param1 u (product d u)

where ``param1 = inverseParamProduct 1``. At the genuine theorem-selected
first Picard state

    u1 = U0 + U1/Lambda + U2/Lambda^2 + a^2 P/Lambda,

this module propagates that split through multiplication by the actual
``AxisData.d = 1-eta^2`` and then through the pinned inverse parameter product.
Every theorem scale remains explicit: ordinary ``Lambda^0..Lambda^-4``,
pressure-linear ``a^2 Lambda^-1..Lambda^-3``, and pressure-square
``a^4 Lambda^-2``. The final minus sign is part of the materialized branch.

For output row ``n>=1`` the exact pinned coefficient identity is

    - 1/n^2 * sum_{i+j=n-1} sum_{k+l=m}
        binom(m,k) * d_eta^(k+1) u1[i] * (d*u1)[j,l].

``d*u1`` itself uses the exact eta-Leibniz product and the theorem's
radial-degree-zero ``d`` field. No caller-supplied d, Lambda, C, amplitude,
pressure/coefficient tables, surrogate state, or arbitrary cutoff is accepted.

This is one mixed-scale ``slow2(x1)`` branch only. It does not by itself claim
complete ``slow2(x1)``, ``naturalRemainder(x1)``, ``x2``, a fixed point, global
weighted ``AxisSpace`` membership, ``NaturalProfileAssembly``, or paper-exact
velocity.
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
        log_factor = +(abs(numerator).ln() - Decimal(inverse_lambda_power) * Lambda.ln())
    return SignedLogCoefficientJet(
        sign=1 if numerator > 0 else -1,
        log_scale=log_scale,
        log_factor=log_factor,
    )


@dataclass(frozen=True)
class MixedScaleFirstPicardSlow2ParamCoefficientJet:
    """One coefficient jet of the pinned ``-param1(u1,d*u1)`` branch."""

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

    def ordinary_correction_terms_decimal(self) -> tuple[Decimal, Decimal, Decimal, Decimal]:
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

    def pressure_linear_terms_log(self) -> tuple[SignedLogCoefficientJet, SignedLogCoefficientJet, SignedLogCoefficientJet]:
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
class ActualScheduleWideFirstPicardSlow2ParamState:
    """Pinned ``-param1(u1,d*u1)`` branch on theorem-selected ``x1``."""

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
    def slow2_param_branch_materialized(self) -> bool:
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

    def _data(self):
        data = actual_schedule_axis_coefficient_data(self.x1.reference)
        if data.epsilon != self.epsilon:
            raise ValueError("AxisData epsilon must match the first Picard state")
        return data

    def _u_ordinary(self, n: int, m: int, eta: float) -> tuple[Decimal, Decimal, Decimal]:
        _, axial = self.x1.jet_pair(n, m, eta)
        return (
            axial.reference,
            axial.inverse_lambda_numerator,
            axial.inverse_lambda_squared_numerator,
        )

    def _u_pressure_numerator(self, n: int, m: int, eta: float) -> Decimal:
        """Return P where the pressure part of one u1 jet is ``a^2 P/Lambda``."""
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(self.x1.remainder.axial.pressure_normalized_factor(n, m, eta) / _TWO)

    def _d_times_u(
        self,
        data,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[tuple[Decimal, Decimal, Decimal], Decimal]:
        """Return split ordinary and normalized-pressure jets of ``d*u1``."""
        ordinary = [Decimal(0), Decimal(0), Decimal(0)]
        pressure = Decimal(0)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for s in range(m + 1):
                u_order = m - s
                d_jet = Decimal.from_float(data.d.jet(0, s, eta))
                weight = Decimal(math.comb(m, s)) * d_jet
                u_ordinary = self._u_ordinary(n, u_order, eta)
                for power, value in enumerate(u_ordinary):
                    ordinary[power] += weight * value
                pressure += weight * self._u_pressure_numerator(n, u_order, eta)
            return tuple(+value for value in ordinary), +pressure

    @staticmethod
    def _zero(*, Lambda: Decimal, amplitude_log: Decimal) -> MixedScaleFirstPicardSlow2ParamCoefficientJet:
        return MixedScaleFirstPicardSlow2ParamCoefficientJet(
            **{name: Decimal(0) for name in _FIELDS},
            Lambda=Lambda,
            amplitude_log=amplitude_log,
        )

    def jet(self, n: int, m: int, eta: float) -> MixedScaleFirstPicardSlow2ParamCoefficientJet:
        """Evaluate the exact pinned ``-param1(u1,d*u1)`` coefficient jet."""
        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)

        self.x1.jet_pair(0 if n == 0 else n - 1, m + 1, eta)
        amplitude_log = self.x1.remainder.axial.wide_pressure.amplitude.log_amplitude(eta)
        data = self._data()
        if n == 0:
            return self._zero(Lambda=self.Lambda, amplitude_log=amplitude_log)

        ordinary = [Decimal(0) for _ in range(_MAX_ORDINARY_POWER + 1)]
        pressure_linear = [Decimal(0), Decimal(0), Decimal(0), Decimal(0)]
        pressure_square = Decimal(0)
        source_degree = n - 1

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            divisor = Decimal(n * n)
            for i in range(source_degree + 1):
                j = source_degree - i
                for k in range(m + 1):
                    l = m - k
                    weight = -Decimal(math.comb(m, k)) / divisor
                    left_ordinary = self._u_ordinary(i, k + 1, eta)
                    left_pressure = self._u_pressure_numerator(i, k + 1, eta)
                    right_ordinary, right_pressure = self._d_times_u(data, j, l, eta)

                    for left_power, left_value in enumerate(left_ordinary):
                        for right_power, right_value in enumerate(right_ordinary):
                            ordinary[left_power + right_power] += weight * left_value * right_value
                    for left_power, left_value in enumerate(left_ordinary):
                        pressure_linear[left_power + 1] += weight * left_value * right_pressure
                    for right_power, right_value in enumerate(right_ordinary):
                        pressure_linear[right_power + 1] += weight * left_pressure * right_value
                    pressure_square += weight * left_pressure * right_pressure

            ordinary = [+value for value in ordinary]
            pressure_linear = [+value for value in pressure_linear]
            pressure_square = +pressure_square

        return MixedScaleFirstPicardSlow2ParamCoefficientJet(
            ordinary_reference=ordinary[0],
            ordinary_inverse_lambda_numerator=ordinary[1],
            ordinary_inverse_lambda_squared_numerator=ordinary[2],
            ordinary_inverse_lambda_cubed_numerator=ordinary[3],
            ordinary_inverse_lambda_fourth_numerator=ordinary[4],
            pressure_linear_inverse_lambda_numerator=pressure_linear[1],
            pressure_linear_inverse_lambda_squared_numerator=pressure_linear[2],
            pressure_linear_inverse_lambda_cubed_numerator=pressure_linear[3],
            pressure_square_inverse_lambda_squared_numerator=pressure_square,
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )


def wide_first_picard_slow2_param_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardSlow2ParamState:
    """Lift the pinned final ``slow2`` param branch without theorem-scale collapse."""
    return ActualScheduleWideFirstPicardSlow2ParamState(x1=x1)
