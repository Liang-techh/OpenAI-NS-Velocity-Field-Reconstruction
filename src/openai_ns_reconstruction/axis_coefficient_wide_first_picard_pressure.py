"""Wide mixed-scale pressure source and pressure chain at the genuine ``x1``.

The pinned natural-axis pressure source at a first Picard state is

    source = a^2 * (phi1 * phi1),

and its pressure contribution is

    j1(-4*A*eta*primitive(source)
       + d*parameterPrimitive(source)
       - 2*eta*mulY(source)).

The angular first-Picard square is already split into five ordinary
``Lambda^-p`` channels.  This adapter applies the exact amplitude-square Bell
factors to the full eta derivatives, propagates every channel through the
pinned radial maps, and exposes the result as normalized Decimal five-tuples:

    d_eta^m pressure_n = a(eta)^2 * sum_{p=0}^4 C_p[n,m](eta)/Lambda^p.

The common ``a^2`` magnitude is never formed in binary64.  ``pressure_jet_log``
provides signed-log views when a caller needs the full scale.  This module
does not apply the later outer ``inverseL`` multiplication.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, localcontext
import math
from typing import Callable

from .axis_amplitude_log_scale import AmplitudeLogSource
from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_data import (
    ActualScheduleAxisCoefficientData,
    actual_schedule_axis_coefficient_data,
)
from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState
from .axis_coefficient_wide_first_picard_angular_square import (
    ActualScheduleWideFirstPicardAngularSquareState,
    wide_first_picard_angular_square_state,
)
from .axis_coefficient_wide_natural_source import (
    ActualScheduleReferenceNaturalSourceLogState,
)


_DECIMAL_PRECISION = 96
_CHANNELS = 5
_SQUARE_FIELDS = (
    "reference",
    "inverse_lambda_numerator",
    "inverse_lambda_squared_numerator",
    "inverse_lambda_cubed_numerator",
    "inverse_lambda_fourth_numerator",
)


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
    amplitude_source: AmplitudeLogSource,
    Lambda: Decimal,
    inverse_lambda_power: int,
) -> SignedLogCoefficientJet:
    _finite_decimal(numerator, "normalized pressure numerator")
    if not isinstance(amplitude_source, AmplitudeLogSource):
        raise TypeError("amplitude_source must be AmplitudeLogSource")
    _finite_decimal(Lambda, "Lambda")
    if Lambda <= 0:
        raise ValueError("Lambda must be positive")
    if Lambda != amplitude_source.enclosure.Lambda:
        raise ValueError("Lambda must match amplitude_source enclosure Lambda")
    if inverse_lambda_power < 0:
        raise ValueError("inverse_lambda_power must be nonnegative")
    if numerator == 0:
        return SignedLogCoefficientJet.zero()
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        log_scale = +(Decimal(2) * amplitude_source.midpoint)
        return SignedLogCoefficientJet(
            sign=1 if numerator > 0 else -1,
            log_scale=log_scale,
            log_factor=+(
                abs(numerator).ln()
                - Decimal(inverse_lambda_power) * Lambda.ln()
            ),
            amplitude_log_scale=amplitude_source.power(2, log_scale),
        )


_WideVectorProvider = Callable[[int, int], tuple[Decimal, ...]]


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardPressureState:
    """The five-channel ``a^2`` source and pressure chain on one actual ``x1``."""

    x1: ActualScheduleWideFirstPicardState
    angular_square: ActualScheduleWideFirstPicardAngularSquareState = field(
        init=False,
        repr=False,
    )
    axis_data: ActualScheduleAxisCoefficientData = field(init=False, repr=False)
    wide_source: ActualScheduleReferenceNaturalSourceLogState = field(
        init=False,
        repr=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.x1, ActualScheduleWideFirstPicardState):
            raise TypeError("x1 must be ActualScheduleWideFirstPicardState")
        if not self.x1.picard_x1_materialized:
            raise ValueError("the genuine first Picard iterate must be materialized")
        if self.x1.fixed_point_materialized:
            raise ValueError("x1 must not be mislabeled as the final fixed point")

        axis_data = actual_schedule_axis_coefficient_data(self.x1.reference)
        wide_pressure = self.x1.remainder.axial.wide_pressure
        wide_source = wide_pressure.source
        if axis_data.epsilon != self.x1.epsilon:
            raise ValueError("AxisData epsilon must match the first Picard state")
        if not isinstance(wide_source, ActualScheduleReferenceNaturalSourceLogState):
            raise TypeError("x1 wide pressure must carry the pinned wide source")
        if wide_source.epsilon != self.x1.epsilon:
            raise ValueError("wide source epsilon must match the first Picard state")
        if wide_pressure.amplitude.Lambda != self.x1.Lambda:
            raise ValueError("x1 Lambda must match the actual pressure amplitude scale")
        if wide_source.amplitude.Lambda != self.x1.Lambda:
            raise ValueError("wide source Lambda must match the first Picard state")
        if not isinstance(self.x1.Lambda, Decimal) or not self.x1.Lambda.is_finite():
            raise ValueError("x1 Lambda must be a finite Decimal")
        if self.x1.Lambda <= 0:
            raise ValueError("x1 Lambda must be positive")

        object.__setattr__(self, "axis_data", axis_data)
        object.__setattr__(self, "wide_source", wide_source)
        object.__setattr__(
            self,
            "angular_square",
            wide_first_picard_angular_square_state(self.x1),
        )

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
    def wide_first_picard_pressure_materialized(self) -> bool:
        return True

    @property
    def pressure_x1_materialized(self) -> bool:
        return True

    @property
    def natural_remainder_x1_materialized(self) -> bool:
        return False

    @property
    def fixed_point_materialized(self) -> bool:
        return False

    @property
    def amplitude(self):
        return self.wide_source.amplitude

    def amplitude_log(self, eta: float) -> Decimal:
        """Return the actual theorem-selected ``log(a(eta))``."""

        return self.amplitude.log_amplitude(_eta_in_window(eta))

    def _source_normalized_factors(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[Decimal, ...]:
        """Build ``C_p`` for ``source = a^2 sum_p Q_p/Lambda^p``."""

        squares = tuple(
            self.angular_square.jet(n, m - k, eta)
            for k in range(m + 1)
        )
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            totals = [Decimal(0) for _ in range(_CHANNELS)]
            for k, square in enumerate(squares):
                bell = self.wide_source._amplitude_square_bell_factor(k, eta)
                weight = Decimal(math.comb(m, k)) * bell
                for channel, name in enumerate(_SQUARE_FIELDS):
                    totals[channel] += weight * getattr(square, name)
            return tuple(+value for value in totals)

    def normalized_source_factors(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[Decimal, ...]:
        """Return normalized source numerators for ``Lambda^0`` through ``Lambda^-4``."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        return self._source_normalized_factors(n, m, eta)

    # Descriptive aliases make the tuple API convenient to downstream assembly.
    source_normalized_factors = normalized_source_factors
    source_jet = normalized_source_factors

    @staticmethod
    def _primitive_vector(
        source: _WideVectorProvider,
        n: int,
        m: int,
    ) -> tuple[Decimal, ...]:
        if n == 0:
            source(0, m)
            return (Decimal(0),) * _CHANNELS
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return tuple(value / Decimal(n) for value in source(n - 1, m))

    @staticmethod
    def _parameter_primitive_vector(
        source: _WideVectorProvider,
        n: int,
        m: int,
    ) -> tuple[Decimal, ...]:
        if n == 0:
            source(0, m + 1)
            return (Decimal(0),) * _CHANNELS
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return tuple(value / Decimal(n) for value in source(n - 1, m + 1))

    @staticmethod
    def _multiply_y_vector(
        source: _WideVectorProvider,
        n: int,
        m: int,
    ) -> tuple[Decimal, ...]:
        if n == 0:
            source(0, m)
            return (Decimal(0),) * _CHANNELS
        return source(n - 1, m)

    def _ordinary_product_vector(
        self,
        ordinary,
        wide: _WideVectorProvider,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[Decimal, ...]:
        """Apply exact radial convolution and eta Leibniz product to a wide vector."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            totals = [Decimal(0) for _ in range(_CHANNELS)]
            for i in range(n + 1):
                j = n - i
                for k in range(m + 1):
                    weight = Decimal(math.comb(m, k))
                    ordinary_value = _decimal_from_float(
                        ordinary.jet(j, m - k, eta),
                        "ordinary coefficient jet",
                    )
                    source_values = wide(i, k)
                    for channel in range(_CHANNELS):
                        totals[channel] += weight * source_values[channel] * ordinary_value
            return tuple(+value for value in totals)

    def _pressure_input_vectors(
        self,
        source: _WideVectorProvider,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[Decimal, ...]:
        data = self.axis_data
        primitive = lambda i, k: self._primitive_vector(source, i, k)
        parameter_primitive = lambda i, k: self._parameter_primitive_vector(source, i, k)
        multiply_y = lambda i, k: self._multiply_y_vector(source, i, k)
        eta_primitive = self._ordinary_product_vector(data.eta, primitive, n, m, eta)
        d_parameter_primitive = self._ordinary_product_vector(
            data.d,
            parameter_primitive,
            n,
            m,
            eta,
        )
        eta_multiply_y = self._ordinary_product_vector(
            data.eta,
            multiply_y,
            n,
            m,
            eta,
        )
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            A = _decimal_from_float(data.A, "AxisData.A")
            return tuple(
                +(-Decimal(4) * A * eta_primitive[channel]
                  + d_parameter_primitive[channel]
                  - Decimal(2) * eta_multiply_y[channel])
                for channel in range(_CHANNELS)
            )

    def normalized_pressure_factors(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[Decimal, ...]:
        """Return normalized pressure numerators ``C_0,...,C_4``."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)

        cache: dict[tuple[int, int], tuple[Decimal, ...]] = {}

        def source(i: int, k: int) -> tuple[Decimal, ...]:
            key = (i, k)
            if key not in cache:
                cache[key] = self._source_normalized_factors(i, k, eta)
            return cache[key]

        if n == 0:
            self._pressure_input_vectors(source, 0, m, eta)
            return (Decimal(0),) * _CHANNELS

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            predecessor = self._pressure_input_vectors(source, n - 1, m, eta)
            divisor = Decimal(n * n)
            return tuple(+(value / divisor) for value in predecessor)

    pressure_normalized_factors = normalized_pressure_factors
    pressure_jet = normalized_pressure_factors

    def source_jet_log(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[SignedLogCoefficientJet, ...]:
        """Return the five source channels with common ``a^2`` signed-log scale."""

        eta = _eta_in_window(eta)
        values = self.normalized_source_factors(n, m, eta)
        amplitude_source = self.amplitude.log_amplitude_source(eta)
        return tuple(
            _signed_log_term(
                value,
                amplitude_source=amplitude_source,
                Lambda=self.Lambda,
                inverse_lambda_power=power,
            )
            for power, value in enumerate(values)
        )

    def pressure_jet_log(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[SignedLogCoefficientJet, ...]:
        """Return the five pressure channels with common ``a^2`` signed-log scale."""

        eta = _eta_in_window(eta)
        values = self.normalized_pressure_factors(n, m, eta)
        amplitude_source = self.amplitude.log_amplitude_source(eta)
        return tuple(
            _signed_log_term(
                value,
                amplitude_source=amplitude_source,
                Lambda=self.Lambda,
                inverse_lambda_power=power,
            )
            for power, value in enumerate(values)
        )


def wide_first_picard_pressure_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardPressureState:
    """Bind the actual five-channel pressure source/chain to one genuine ``x1``."""

    return ActualScheduleWideFirstPicardPressureState(x1=x1)


__all__ = [
    "ActualScheduleWideFirstPicardPressureState",
    "wide_first_picard_pressure_state",
]
