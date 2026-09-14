"""Pinned genuine ``slow1(x1)`` at the actual-schedule first Picard state.

The pinned formal source defines ``bu = average(u)`` and

    slow1 =
      j2(product(product(averageCoefficient(d),bu)
                 + product(angularSlowCoefficient(d),u), phi))
      + param2(bu, product(d.d,phi))
      + dot2(product(averageCoefficient(d),bu), phi)
      + product(d.d, mixed2(bu,phi))
      - param2(phi, product(d.d,u)),

with ``averageCoefficient(d)=(2*D)*eta`` and
``angularSlowCoefficient(d)=(2*h)*eta``.  This implementation consumes one
actual ``ActualScheduleWideFirstPicardState`` and preserves ordinary
``Lambda^0..Lambda^-4`` plus pressure-linear ``a^2 Lambda^-1..-3`` families.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState

_PRECISION = 96
_TWO = Decimal(2)
_MAX_ORDINARY = 4
_MAX_PRESSURE = 3


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _eta(value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or not WINDOW_LEFT <= value <= WINDOW_RIGHT:
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]")
    return value


def _dec(value: float, name: str) -> Decimal:
    value = float(value)
    if not math.isfinite(value):
        raise ArithmeticError(f"{name} must be finite")
    return Decimal.from_float(value)


def _finite(value: Decimal, name: str) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError(f"{name} must be a finite Decimal")
    return value


def _signed_log(numerator: Decimal, *, amplitude_log: Decimal, Lambda: Decimal,
                power: int) -> SignedLogCoefficientJet:
    _finite(numerator, "numerator")
    _finite(amplitude_log, "amplitude_log")
    _finite(Lambda, "Lambda")
    if Lambda <= 0 or power < 1:
        raise ValueError("invalid pressure scale")
    if numerator == 0:
        return SignedLogCoefficientJet.zero()
    with localcontext() as ctx:
        ctx.prec = _PRECISION
        return SignedLogCoefficientJet(
            sign=1 if numerator > 0 else -1,
            log_scale=+(Decimal(2) * amplitude_log),
            log_factor=+(abs(numerator).ln() - Decimal(power) * Lambda.ln()),
        )


@dataclass(frozen=True)
class MixedScaleFirstPicardSlow1CoefficientJet:
    ordinary_reference: Decimal
    ordinary_inverse_lambda_numerator: Decimal
    ordinary_inverse_lambda_squared_numerator: Decimal
    ordinary_inverse_lambda_cubed_numerator: Decimal
    ordinary_inverse_lambda_fourth_numerator: Decimal
    pressure_linear_inverse_lambda_numerator: Decimal
    pressure_linear_inverse_lambda_squared_numerator: Decimal
    pressure_linear_inverse_lambda_cubed_numerator: Decimal
    Lambda: Decimal
    amplitude_log: Decimal

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            _finite(value, name)
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    def ordinary_factors_decimal(self) -> tuple[Decimal, ...]:
        return (
            self.ordinary_reference,
            self.ordinary_inverse_lambda_numerator,
            self.ordinary_inverse_lambda_squared_numerator,
            self.ordinary_inverse_lambda_cubed_numerator,
            self.ordinary_inverse_lambda_fourth_numerator,
        )

    def pressure_factors_decimal(self) -> tuple[Decimal, ...]:
        return (
            self.pressure_linear_inverse_lambda_numerator,
            self.pressure_linear_inverse_lambda_squared_numerator,
            self.pressure_linear_inverse_lambda_cubed_numerator,
        )

    def pressure_terms_log(self) -> tuple[SignedLogCoefficientJet, ...]:
        return tuple(
            _signed_log(value, amplitude_log=self.amplitude_log,
                        Lambda=self.Lambda, power=power)
            for power, value in enumerate(self.pressure_factors_decimal(), start=1)
        )


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardSlow1State:
    """Genuine pinned ``slow1(x1)`` bound to one actual-schedule ``x1``."""

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

    @property
    def axis_data(self):
        return self.x1.remainder.axial.axis_data

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
    def slow1_x1_materialized(self) -> bool:
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

    def _phi(self, n: int, m: int, eta: float) -> tuple[Decimal, Decimal, Decimal]:
        phi, _ = self.x1.jet_pair(n, m, eta)
        return (phi.reference, phi.inverse_lambda_numerator,
                phi.inverse_lambda_squared_numerator)

    def _u(self, n: int, m: int, eta: float) -> tuple[Decimal, Decimal, Decimal]:
        _, u = self.x1.jet_pair(n, m, eta)
        return (u.reference, u.inverse_lambda_numerator,
                u.inverse_lambda_squared_numerator)

    def _up(self, n: int, m: int, eta: float) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = _PRECISION
            return +(self.x1.remainder.axial.pressure_normalized_factor(n, m, eta) / _TWO)

    def _bu(self, n: int, m: int, eta: float) -> tuple[Decimal, Decimal, Decimal]:
        divisor = Decimal(n + 1)
        with localcontext() as ctx:
            ctx.prec = _PRECISION
            return tuple(+(value / divisor) for value in self._u(n, m, eta))

    def _bup(self, n: int, m: int, eta: float) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = _PRECISION
            return +(self._up(n, m, eta) / Decimal(n + 1))

    def _fixed(self, field, m: int, eta: float, scalar: float = 1.0) -> Decimal:
        return _dec(scalar * field.jet(0, m, eta), "fixed coefficient jet")

    def _fixed_times_state(self, fixed, ordinary, pressure,
                           n: int, m: int, eta: float):
        out_o = [Decimal(0), Decimal(0), Decimal(0)]
        out_p = Decimal(0)
        with localcontext() as ctx:
            ctx.prec = _PRECISION
            for k in range(m + 1):
                l = m - k
                weight = Decimal(math.comb(m, k)) * fixed(k, eta)
                for power, value in enumerate(ordinary(n, l, eta)):
                    out_o[power] += weight * value
                if pressure is not None:
                    out_p += weight * pressure(n, l, eta)
        return tuple(out_o), +out_p

    def _d_times_phi(self, n: int, m: int, eta: float):
        return self._fixed_times_state(
            lambda k, e: self._fixed(self.axis_data.d, k, e),
            self._phi, None, n, m, eta,
        )[0]

    def _d_times_u(self, n: int, m: int, eta: float):
        return self._fixed_times_state(
            lambda k, e: self._fixed(self.axis_data.d, k, e),
            self._u, self._up, n, m, eta,
        )

    def _avgcoef_times_bu(self, n: int, m: int, eta: float):
        scalar = 2.0 * self.axis_data.D
        return self._fixed_times_state(
            lambda k, e: self._fixed(self.axis_data.eta, k, e, scalar),
            self._bu, self._bup, n, m, eta,
        )

    def _slowcoef_times_u(self, n: int, m: int, eta: float):
        scalar = 2.0 * self.axis_data.h
        return self._fixed_times_state(
            lambda k, e: self._fixed(self.axis_data.eta, k, e, scalar),
            self._u, self._up, n, m, eta,
        )

    def _left_sum(self, n: int, m: int, eta: float):
        ao, ap = self._avgcoef_times_bu(n, m, eta)
        so, sp = self._slowcoef_times_u(n, m, eta)
        return tuple(a + b for a, b in zip(ao, so)), +(ap + sp)

    def _scaled_product(self, left_ord, left_p, right_ord, right_p,
                        n: int, m: int, eta: float, *,
                        left_eta_shift: int = 0,
                        right_radial_factor: bool = False):
        ordinary = [Decimal(0) for _ in range(_MAX_ORDINARY + 1)]
        pressure = [Decimal(0) for _ in range(_MAX_PRESSURE + 1)]
        with localcontext() as ctx:
            ctx.prec = _PRECISION
            for i in range(n + 1):
                j = n - i
                if right_radial_factor and j == 0:
                    continue
                radial = Decimal(j if right_radial_factor else 1)
                for k in range(m + 1):
                    l = m - k
                    weight = Decimal(math.comb(m, k)) * radial
                    lo = left_ord(i, k + left_eta_shift, eta)
                    ro = right_ord(j, l, eta)
                    for a, av in enumerate(lo):
                        for b, bv in enumerate(ro):
                            ordinary[a + b] += weight * av * bv
                    lp = left_p(i, k + left_eta_shift, eta) if left_p else None
                    rp = right_p(j, l, eta) if right_p else None
                    if rp is not None:
                        for a, av in enumerate(lo):
                            pressure[a + 1] += weight * av * rp
                    if lp is not None:
                        for b, bv in enumerate(ro):
                            pressure[b + 1] += weight * lp * bv
                    if lp is not None and rp is not None and lp != 0 and rp != 0:
                        raise AssertionError("slow1 must not create a pressure-square family")
        return ordinary, pressure

    def jet(self, n: int, m: int, eta: float) -> MixedScaleFirstPicardSlow1CoefficientJet:
        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta(eta)
        amplitude_log = self.x1.remainder.axial.wide_pressure.amplitude.log_amplitude(eta)
        ordinary = [Decimal(0) for _ in range(_MAX_ORDINARY + 1)]
        pressure = [Decimal(0) for _ in range(_MAX_PRESSURE + 1)]
        if n == 0:
            return MixedScaleFirstPicardSlow1CoefficientJet(
                ordinary_reference=ordinary[0],
                ordinary_inverse_lambda_numerator=ordinary[1],
                ordinary_inverse_lambda_squared_numerator=ordinary[2],
                ordinary_inverse_lambda_cubed_numerator=ordinary[3],
                ordinary_inverse_lambda_fourth_numerator=ordinary[4],
                pressure_linear_inverse_lambda_numerator=pressure[1],
                pressure_linear_inverse_lambda_squared_numerator=pressure[2],
                pressure_linear_inverse_lambda_cubed_numerator=pressure[3],
                Lambda=self.Lambda,
                amplitude_log=amplitude_log,
            )

        source_degree = n - 1
        divisor = Decimal(n) * Decimal(n + 1)
        with localcontext() as ctx:
            ctx.prec = _PRECISION

            # j2(product(product(avgCoeff,bu)+product(slowCoeff,u),phi))
            o, p = self._scaled_product(
                lambda i, k, e: self._left_sum(i, k, e)[0],
                lambda i, k, e: self._left_sum(i, k, e)[1],
                self._phi, None, source_degree, m, eta,
            )
            for q in range(5):
                ordinary[q] += o[q]
            for q in range(1, 4):
                pressure[q] += p[q]

            # + param2(bu, product(d,phi))
            o, p = self._scaled_product(
                self._bu, self._bup, self._d_times_phi, None,
                source_degree, m, eta, left_eta_shift=1,
            )
            for q in range(5):
                ordinary[q] += o[q]
            for q in range(1, 4):
                pressure[q] += p[q]

            # + dot2(product(avgCoeff,bu),phi)
            o, p = self._scaled_product(
                lambda i, k, e: self._avgcoef_times_bu(i, k, e)[0],
                lambda i, k, e: self._avgcoef_times_bu(i, k, e)[1],
                self._phi, None, source_degree, m, eta,
                right_radial_factor=True,
            )
            for q in range(5):
                ordinary[q] += o[q]
            for q in range(1, 4):
                pressure[q] += p[q]

            # + product(d,mixed2(bu,phi)); mixed2 has the same n(n+1) divisor.
            for k in range(m + 1):
                l = m - k
                d_weight = Decimal(math.comb(m, k)) * self._fixed(self.axis_data.d, k, eta)
                o, p = self._scaled_product(
                    self._bu, self._bup, self._phi, None,
                    source_degree, l, eta,
                    left_eta_shift=1, right_radial_factor=True,
                )
                for q in range(5):
                    ordinary[q] += d_weight * o[q]
                for q in range(1, 4):
                    pressure[q] += d_weight * p[q]

            # - param2(phi, product(d,u))
            o, p = self._scaled_product(
                self._phi, None,
                lambda i, k, e: self._d_times_u(i, k, e)[0],
                lambda i, k, e: self._d_times_u(i, k, e)[1],
                source_degree, m, eta, left_eta_shift=1,
            )
            for q in range(5):
                ordinary[q] -= o[q]
            for q in range(1, 4):
                pressure[q] -= p[q]

            ordinary = [+(value / divisor) for value in ordinary]
            pressure = [+(value / divisor) for value in pressure]

        return MixedScaleFirstPicardSlow1CoefficientJet(
            ordinary_reference=ordinary[0],
            ordinary_inverse_lambda_numerator=ordinary[1],
            ordinary_inverse_lambda_squared_numerator=ordinary[2],
            ordinary_inverse_lambda_cubed_numerator=ordinary[3],
            ordinary_inverse_lambda_fourth_numerator=ordinary[4],
            pressure_linear_inverse_lambda_numerator=pressure[1],
            pressure_linear_inverse_lambda_squared_numerator=pressure[2],
            pressure_linear_inverse_lambda_cubed_numerator=pressure[3],
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )


def wide_first_picard_slow1_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardSlow1State:
    return ActualScheduleWideFirstPicardSlow1State(x1=x1)
