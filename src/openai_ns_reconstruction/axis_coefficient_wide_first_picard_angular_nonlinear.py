"""Mixed-scale angular nonlinear terms at the first Picard state.

The pinned coefficient remainder contains two angular nonlinear pieces at a
first Picard input ``x1``:

``quad1 = j2((d * normalizedGradient) * u1 * phi1)``

and

``slow1 = j2((averageCoefficient * average(u1)
              + angularSlowCoefficient * u1) * phi1)
          + param2(average(u1), d * phi1)
          + dot2(averageCoefficient * average(u1), phi1)
          + d * mixed2(average(u1), phi1)
          - param2(phi1, d * u1)``.

This module evaluates both expressions with a compact mixed-scale channel
vector.  The ordinary channels are the five numerators of
``Lambda^0`` through ``Lambda^-4``.  The pressure-linear channels are the
three normalized numerators of ``a^2 Lambda^-1`` through ``a^2 Lambda^-3``.
Products have no pressure-pressure channel, so no ``a^4`` term is admitted.
Every radial product uses the finite coefficient convolution and eta Leibniz
rule, while ``j2``, ``param2``, ``dot2``, and ``mixed2`` use the exact divisor
``n * (n + 1)``.  If an input branch would produce a nonzero term outside the
declared channel widths, evaluation fails closed instead of silently
truncating it.

The input is restricted to the actual theorem-selected
``ActualScheduleWideFirstPicardState``.  Fixed coefficient fields are rebuilt
from its actual ``AxisData``.  The layer stops before the outer ``inverseL``
and natural resolvent used by the full remainder; it remains formal structure
and does not claim a complete ``naturalRemainder(x1)`` or a paper-exact field.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, localcontext
from functools import lru_cache
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


_DECIMAL_PRECISION = 96
_ORDINARY_WIDTH = 5
_PRESSURE_WIDTH = 3
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


@dataclass(frozen=True)
class _ChannelVector:
    """Finite mixed-scale channels with no pressure-square channel."""

    ordinary: tuple[Decimal, ...]
    pressure_linear: tuple[Decimal, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.ordinary, tuple) or len(self.ordinary) != _ORDINARY_WIDTH:
            raise ValueError("ordinary channel vector must have five entries")
        if not isinstance(self.pressure_linear, tuple) or len(self.pressure_linear) != _PRESSURE_WIDTH:
            raise ValueError("pressure-linear channel vector must have three entries")
        for index, value in enumerate(self.ordinary):
            _finite_decimal(value, f"ordinary[{index}]")
        for index, value in enumerate(self.pressure_linear):
            _finite_decimal(value, f"pressure_linear[{index}]")

    @classmethod
    def zero(cls) -> "_ChannelVector":
        return cls(
            ordinary=(Decimal(0),) * _ORDINARY_WIDTH,
            pressure_linear=(Decimal(0),) * _PRESSURE_WIDTH,
        )


Family = Callable[[int, int, float], _ChannelVector]


def _add_vectors(*vectors: _ChannelVector) -> _ChannelVector:
    if not vectors:
        return _ChannelVector.zero()
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ordinary = tuple(
            +sum(vector.ordinary[index] for vector in vectors)
            for index in range(_ORDINARY_WIDTH)
        )
        pressure_linear = tuple(
            +sum(vector.pressure_linear[index] for vector in vectors)
            for index in range(_PRESSURE_WIDTH)
        )
    return _ChannelVector(ordinary=ordinary, pressure_linear=pressure_linear)


def _scale_vector(vector: _ChannelVector, scalar: Decimal) -> _ChannelVector:
    scalar = _finite_decimal(scalar, "channel scalar")
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        return _ChannelVector(
            ordinary=tuple(+(scalar * value) for value in vector.ordinary),
            pressure_linear=tuple(
                +(scalar * value) for value in vector.pressure_linear
            ),
        )


def _product_family(left: Family, right: Family) -> Family:
    """Return the exact radial/eta product in the declared channels."""

    def product(n: int, m: int, eta: float) -> _ChannelVector:
        ordinary = [Decimal(0) for _ in range(_ORDINARY_WIDTH)]
        pressure_linear = [Decimal(0) for _ in range(_PRESSURE_WIDTH)]
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for i in range(n + 1):
                j = n - i
                for k in range(m + 1):
                    l = m - k
                    weight = Decimal(math.comb(m, k))
                    left_value = left(i, k, eta)
                    right_value = right(j, l, eta)

                    _accumulate_bilinear(
                        ordinary,
                        pressure_linear,
                        left_value,
                        right_value,
                        weight,
                    )

            return _ChannelVector(
                ordinary=tuple(+value for value in ordinary),
                pressure_linear=tuple(+value for value in pressure_linear),
            )

    return product


def _average_family(source: Family) -> Family:
    def average(n: int, m: int, eta: float) -> _ChannelVector:
        value = source(n, m, eta)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            divisor = Decimal(n + 1)
            return _scale_vector(value, +(Decimal(1) / divisor))

    return average


def _regular_inverse_two_family(source: Family) -> Family:
    """Return the zero-datum ``J2`` family with divisor ``n * (n + 1)``."""

    def inverse(n: int, m: int, eta: float) -> _ChannelVector:
        if n == 0:
            source(0, m, eta)
            return _ChannelVector.zero()
        value = source(n - 1, m, eta)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            divisor = Decimal(n * (n + 1))
            return _scale_vector(value, +(Decimal(1) / divisor))

    return inverse


def _inverse_param_two_family(left: Family, right: Family) -> Family:
    """Return ``param2(left,right)`` with the left eta derivative."""

    def param(n: int, m: int, eta: float) -> _ChannelVector:
        if n == 0:
            left(0, m + 1, eta)
            right(0, m, eta)
            return _ChannelVector.zero()
        ordinary = [Decimal(0) for _ in range(_ORDINARY_WIDTH)]
        pressure_linear = [Decimal(0) for _ in range(_PRESSURE_WIDTH)]
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            divisor = Decimal(n * (n + 1))
            for i in range(n):
                j = n - 1 - i
                for k in range(m + 1):
                    l = m - k
                    weight = Decimal(math.comb(m, k)) / divisor
                    left_value = left(i, k + 1, eta)
                    right_value = right(j, l, eta)
                    _accumulate_bilinear(
                        ordinary,
                        pressure_linear,
                        left_value,
                        right_value,
                        weight,
                    )
        return _ChannelVector(
            ordinary=tuple(+value for value in ordinary),
            pressure_linear=tuple(+value for value in pressure_linear),
        )

    return param


def _inverse_dot_two_family(left: Family, right: Family) -> Family:
    """Return ``dot2(left,right)`` with the right Euler factor."""

    def dot(n: int, m: int, eta: float) -> _ChannelVector:
        if n == 0:
            left(0, m, eta)
            right(0, m, eta)
            return _ChannelVector.zero()
        ordinary = [Decimal(0) for _ in range(_ORDINARY_WIDTH)]
        pressure_linear = [Decimal(0) for _ in range(_PRESSURE_WIDTH)]
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            divisor = Decimal(n * (n + 1))
            for i in range(n):
                j = n - 1 - i
                for k in range(m + 1):
                    l = m - k
                    weight = Decimal(math.comb(m, k)) * Decimal(j) / divisor
                    left_value = left(i, k, eta)
                    right_value = right(j, l, eta)
                    _accumulate_bilinear(
                        ordinary,
                        pressure_linear,
                        left_value,
                        right_value,
                        weight,
                    )
        return _ChannelVector(
            ordinary=tuple(+value for value in ordinary),
            pressure_linear=tuple(+value for value in pressure_linear),
        )

    return dot


def _inverse_mixed_two_family(left: Family, right: Family) -> Family:
    """Return ``mixed2(left,right)`` with left eta and right Euler factors."""

    def mixed(n: int, m: int, eta: float) -> _ChannelVector:
        if n == 0:
            left(0, m + 1, eta)
            right(0, m, eta)
            return _ChannelVector.zero()
        ordinary = [Decimal(0) for _ in range(_ORDINARY_WIDTH)]
        pressure_linear = [Decimal(0) for _ in range(_PRESSURE_WIDTH)]
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            divisor = Decimal(n * (n + 1))
            for i in range(n):
                j = n - 1 - i
                for k in range(m + 1):
                    l = m - k
                    weight = Decimal(math.comb(m, k)) * Decimal(j) / divisor
                    left_value = left(i, k + 1, eta)
                    right_value = right(j, l, eta)
                    _accumulate_bilinear(
                        ordinary,
                        pressure_linear,
                        left_value,
                        right_value,
                        weight,
                    )
        return _ChannelVector(
            ordinary=tuple(+value for value in ordinary),
            pressure_linear=tuple(+value for value in pressure_linear),
        )

    return mixed


def _accumulate_bilinear(
    ordinary: list[Decimal],
    pressure_linear: list[Decimal],
    left: _ChannelVector,
    right: _ChannelVector,
    weight: Decimal,
) -> None:
    """Accumulate one weighted product, rejecting absent channels."""

    def accumulate(
        target: list[Decimal],
        power: int,
        left_item: Decimal,
        right_item: Decimal,
        width: int,
        label: str,
    ) -> None:
        if left_item == 0 or right_item == 0 or weight == 0:
            return
        if power >= width:
            raise ArithmeticError(
                f"nonzero {label} product exceeds the declared channel width"
            )
        target[power] += weight * left_item * right_item

    for left_power, left_item in enumerate(left.ordinary):
        for right_power, right_item in enumerate(right.ordinary):
            power = left_power + right_power
            accumulate(
                ordinary,
                power,
                left_item,
                right_item,
                _ORDINARY_WIDTH,
                "ordinary",
            )
    for left_power, left_item in enumerate(left.ordinary):
        for right_power, right_item in enumerate(right.pressure_linear):
            power = left_power + right_power
            accumulate(
                pressure_linear,
                power,
                left_item,
                right_item,
                _PRESSURE_WIDTH,
                "ordinary-pressure",
            )
    for left_power, left_item in enumerate(left.pressure_linear):
        for right_power, right_item in enumerate(right.ordinary):
            power = left_power + right_power
            accumulate(
                pressure_linear,
                power,
                left_item,
                right_item,
                _PRESSURE_WIDTH,
                "pressure-ordinary",
            )
    for left_item in left.pressure_linear:
        for right_item in right.pressure_linear:
            if left_item != 0 and right_item != 0 and weight != 0:
                raise ArithmeticError(
                    "nonzero pressure-pressure product has no declared a^4 channel"
                )


def _fixed_field_family(
    data: ActualScheduleAxisCoefficientData,
    name: str,
    scalar: Decimal = Decimal(1),
) -> Family:
    field_value = getattr(data, name)

    def fixed(n: int, m: int, eta: float) -> _ChannelVector:
        value = _decimal_from_float(field_value.jet(n, m, eta), f"AxisData.{name} jet")
        if n != 0:
            if value != 0:
                raise ValueError(f"AxisData.{name} must have radial degree zero")
            return _ChannelVector.zero()
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return _ChannelVector(
                ordinary=(+(scalar * value), Decimal(0), Decimal(0), Decimal(0), Decimal(0)),
                pressure_linear=(Decimal(0), Decimal(0), Decimal(0)),
            )

    return fixed


def _signed_log_term(
    numerator: Decimal,
    *,
    amplitude_log: Decimal,
    Lambda: Decimal,
    inverse_lambda_power: int,
    amplitude_log_source: AmplitudeLogSource | None = None,
) -> SignedLogCoefficientJet:
    _finite_decimal(numerator, "pressure-linear numerator")
    _finite_decimal(amplitude_log, "amplitude_log")
    _finite_decimal(Lambda, "Lambda")
    if Lambda <= 0:
        raise ValueError("Lambda must be positive")
    if amplitude_log_source is not None:
        if not isinstance(amplitude_log_source, AmplitudeLogSource):
            raise TypeError("amplitude_log_source must be AmplitudeLogSource")
        if amplitude_log_source.midpoint != amplitude_log:
            raise ValueError("amplitude_log_source midpoint must match amplitude_log")
        if amplitude_log_source.enclosure.Lambda != Lambda:
            raise ValueError("Lambda must match amplitude_log_source enclosure Lambda")
    if inverse_lambda_power < 1:
        raise ValueError("pressure-linear inverse-Lambda power must be positive")
    if numerator == 0:
        return SignedLogCoefficientJet.zero()
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        log_scale = +(_TWO * amplitude_log)
        log_factor = +(
            abs(numerator).ln() - Decimal(inverse_lambda_power) * Lambda.ln()
        )
    return SignedLogCoefficientJet(
        sign=1 if numerator > 0 else -1,
        log_scale=log_scale,
        log_factor=log_factor,
        amplitude_log_scale=(
            None
            if amplitude_log_source is None
            else amplitude_log_source.power(2, log_scale)
        ),
    )


@dataclass(frozen=True)
class MixedScaleFirstPicardAngularNonlinearCoefficientJet:
    """One angular ``quad1`` or ``slow1`` jet with explicit scale channels."""

    ordinary: tuple[Decimal, ...]
    pressure_linear: tuple[Decimal, ...]
    Lambda: Decimal
    amplitude_log: Decimal
    amplitude_log_source: AmplitudeLogSource | None = None

    def __post_init__(self) -> None:
        _ChannelVector(self.ordinary, self.pressure_linear)
        _finite_decimal(self.Lambda, "Lambda")
        _finite_decimal(self.amplitude_log, "amplitude_log")
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")
        if self.amplitude_log_source is not None:
            if not isinstance(self.amplitude_log_source, AmplitudeLogSource):
                raise TypeError("amplitude_log_source must be AmplitudeLogSource")
            if self.amplitude_log_source.midpoint != self.amplitude_log:
                raise ValueError(
                    "amplitude_log_source midpoint must match amplitude_log"
                )
            if self.amplitude_log_source.enclosure.Lambda != self.Lambda:
                raise ValueError(
                    "Lambda must match amplitude_log_source enclosure Lambda"
                )

    def ordinary_correction_terms_decimal(self) -> tuple[Decimal, Decimal, Decimal, Decimal]:
        """Return ordinary ``Lambda^-1`` through ``Lambda^-4`` terms."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            inverse = +(Decimal(1) / self.Lambda)
            inverse2 = +(inverse * inverse)
            return (
                +(self.ordinary[1] * inverse),
                +(self.ordinary[2] * inverse2),
                +(self.ordinary[3] * inverse2 * inverse),
                +(self.ordinary[4] * inverse2 * inverse2),
            )

    def correction_terms_decimal(self) -> tuple[Decimal, Decimal, Decimal, Decimal]:
        """Alias for the ordinary inverse-Lambda correction terms."""

        return self.ordinary_correction_terms_decimal()

    def ordinary_terms_decimal(self) -> tuple[Decimal, Decimal, Decimal, Decimal, Decimal]:
        """Return all five ordinary terms, including the reference channel."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            inverse = +(Decimal(1) / self.Lambda)
            return tuple(
                +(
                    value
                    * (inverse**power if power else Decimal(1))
                )
                for power, value in enumerate(self.ordinary)
            )

    def pressure_linear_terms_log(
        self,
    ) -> tuple[SignedLogCoefficientJet, SignedLogCoefficientJet, SignedLogCoefficientJet]:
        """Return signed-log ``a^2 Lambda^-1..Lambda^-3`` terms."""

        return tuple(
            _signed_log_term(
                numerator,
                amplitude_log=self.amplitude_log,
                Lambda=self.Lambda,
                inverse_lambda_power=power,
                amplitude_log_source=self.amplitude_log_source,
            )
            for power, numerator in enumerate(self.pressure_linear, start=1)
        )


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardAngularNonlinearState:
    """Actual-schedule angular ``quad1`` and ``slow1`` at ``x1``."""

    x1: ActualScheduleWideFirstPicardState
    data: ActualScheduleAxisCoefficientData = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.x1, ActualScheduleWideFirstPicardState):
            raise TypeError("x1 must be ActualScheduleWideFirstPicardState")
        if not self.x1.picard_x1_materialized:
            raise ValueError("the genuine first Picard iterate must be materialized")
        if self.x1.fixed_point_materialized:
            raise ValueError("x1 must not be mislabeled as the final fixed point")
        amplitude = self.x1.remainder.axial.wide_pressure.amplitude
        if self.x1.Lambda != amplitude.Lambda:
            raise ValueError("x1 Lambda must match the actual theorem-selected amplitude scale")
        if amplitude.epsilon != self.x1.epsilon:
            raise ValueError("amplitude epsilon must match the first Picard state")
        _finite_decimal(self.x1.Lambda, "x1 Lambda")
        if self.x1.Lambda <= 0:
            raise ValueError("x1 Lambda must be positive")
        _finite_decimal(
            amplitude.log_amplitude_source(0.0).midpoint,
            "actual amplitude log at eta=0",
        )
        data = actual_schedule_axis_coefficient_data(self.x1.reference)
        if data.epsilon != self.x1.epsilon:
            raise ValueError("AxisData epsilon must match the first Picard state")
        if data.j != self.x1.reference.reference.j:
            raise ValueError("AxisData j must match the first Picard state")
        if data.sigma != self.x1.reference.reference.sigma:
            raise ValueError("AxisData sigma must match the first Picard state")
        object.__setattr__(self, "data", data)

    def _validate_identity(self) -> None:
        if self.data.epsilon != self.x1.epsilon:
            raise ValueError("AxisData epsilon must match the first Picard state")
        if self.data.j != self.x1.reference.reference.j:
            raise ValueError("AxisData j must match the first Picard state")
        if self.data.sigma != self.x1.reference.reference.sigma:
            raise ValueError("AxisData sigma must match the first Picard state")
        if self.x1.Lambda != self.x1.remainder.axial.wide_pressure.amplitude.Lambda:
            raise ValueError("x1 Lambda must match the actual theorem-selected amplitude scale")
        amplitude = self.x1.remainder.axial.wide_pressure.amplitude
        if amplitude.epsilon != self.x1.epsilon:
            raise ValueError("amplitude epsilon must match the first Picard state")

    def _actual_amplitude_log_source(self, eta: float) -> AmplitudeLogSource:
        return self.x1.remainder.axial.wide_pressure.amplitude.log_amplitude_source(
            eta
        )

    def _validate_jet_pair_identity(
        self,
        angular: object,
        axial: object,
        eta: float,
    ) -> Decimal:
        if angular.Lambda != self.Lambda:
            raise ValueError("phi1 Lambda mismatch")
        if axial.Lambda != self.Lambda:
            raise ValueError("u1 Lambda mismatch")
        amplitude_log_source = self._actual_amplitude_log_source(eta)
        amplitude_log = _finite_decimal(
            amplitude_log_source.midpoint,
            "actual amplitude log",
        )
        pressure = axial.pressure_over_two_lambda
        if not isinstance(pressure, SignedLogCoefficientJet):
            raise TypeError("u1 pressure channel must be a SignedLogCoefficientJet")
        if pressure.sign != 0:
            assert pressure.log_scale is not None
            with localcontext() as ctx:
                ctx.prec = _DECIMAL_PRECISION
                expected_log_scale = +(_TWO * amplitude_log)
            if pressure.log_scale != expected_log_scale:
                raise ValueError("u1 pressure amplitude log mismatch")
            source = pressure.amplitude_log_scale
            if source is None or source.q != 2:
                raise ValueError("u1 pressure amplitude source metadata is missing")
            source.source.assert_compatible(amplitude_log_source)
        return amplitude_log

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
    def quad1_materialized(self) -> bool:
        return True

    @property
    def slow1_materialized(self) -> bool:
        return True

    @property
    def angular_nonlinear_materialized(self) -> bool:
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

    def _phi_family(self) -> Family:
        @lru_cache(maxsize=None)
        def phi(n: int, m: int, eta: float) -> _ChannelVector:
            angular, axial = self.x1.jet_pair(n, m, eta)
            self._validate_jet_pair_identity(angular, axial, eta)
            return _ChannelVector(
                ordinary=(
                    angular.reference,
                    angular.inverse_lambda_numerator,
                    angular.inverse_lambda_squared_numerator,
                    Decimal(0),
                    Decimal(0),
                ),
                pressure_linear=(Decimal(0), Decimal(0), Decimal(0)),
            )

        return phi

    def _u_family(self) -> Family:
        @lru_cache(maxsize=None)
        def axial(n: int, m: int, eta: float) -> _ChannelVector:
            angular, value = self.x1.jet_pair(n, m, eta)
            self._validate_jet_pair_identity(angular, value, eta)
            with localcontext() as ctx:
                ctx.prec = _DECIMAL_PRECISION
                pressure = +(
                    self.x1.remainder.axial.pressure_normalized_factor(n, m, eta)
                    / _TWO
                )
            return _ChannelVector(
                ordinary=(
                    value.reference,
                    value.inverse_lambda_numerator,
                    value.inverse_lambda_squared_numerator,
                    Decimal(0),
                    Decimal(0),
                ),
                pressure_linear=(pressure, Decimal(0), Decimal(0)),
            )

        return axial

    def _fixed_families(self) -> tuple[Family, Family, Family, Family]:
        data = self.data
        d = _fixed_field_family(data, "d")
        normalized_gradient = _fixed_field_family(data, "normalizedGradient")
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            average_scalar = +(Decimal(2) * Decimal.from_float(data.D))
            angular_slow_scalar = +(Decimal(2) * Decimal.from_float(data.h))
        average_coefficient = _fixed_field_family(data, "eta", average_scalar)
        angular_slow = _fixed_field_family(data, "eta", angular_slow_scalar)
        angular_quadratic = _product_family(d, normalized_gradient)
        return d, angular_quadratic, average_coefficient, angular_slow

    def _result(
        self,
        vector: _ChannelVector,
        eta: float,
    ) -> MixedScaleFirstPicardAngularNonlinearCoefficientJet:
        amplitude_log_source = self._actual_amplitude_log_source(eta)
        return MixedScaleFirstPicardAngularNonlinearCoefficientJet(
            ordinary=vector.ordinary,
            pressure_linear=vector.pressure_linear,
            Lambda=self.Lambda,
            amplitude_log=amplitude_log_source.midpoint,
            amplitude_log_source=amplitude_log_source,
        )

    def quad1_jet(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> MixedScaleFirstPicardAngularNonlinearCoefficientJet:
        """Evaluate ``j2((d*normalizedGradient)*u1*phi1)``."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        self._validate_identity()
        u = self._u_family()
        phi = self._phi_family()
        _, angular_quadratic, _, _ = self._fixed_families()
        source = _product_family(_product_family(angular_quadratic, u), phi)
        result = _regular_inverse_two_family(source)(n, m, eta)
        return self._result(result, eta)

    def slow1_jet(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> MixedScaleFirstPicardAngularNonlinearCoefficientJet:
        """Evaluate the complete pinned angular ``slow1(x1)`` expression."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        self._validate_identity()

        u = self._u_family()
        phi = self._phi_family()
        d, _, average_coefficient, angular_slow = self._fixed_families()
        bu = _average_family(u)

        average_transport = _product_family(average_coefficient, bu)
        angular_transport = _product_family(angular_slow, u)
        transport = lambda i, q, z: _add_vectors(
            average_transport(i, q, z),
            angular_transport(i, q, z),
        )

        transport_term = _regular_inverse_two_family(
            _product_family(transport, phi)
        )
        d_phi = _product_family(d, phi)
        d_u = _product_family(d, u)
        param_average_term = _inverse_param_two_family(bu, d_phi)
        dot_term = _inverse_dot_two_family(average_transport, phi)
        mixed_term = _product_family(
            d,
            _inverse_mixed_two_family(bu, phi),
        )
        param_phi_term = _inverse_param_two_family(phi, d_u)

        source = lambda i, q, z: _add_vectors(
            transport_term(i, q, z),
            param_average_term(i, q, z),
            dot_term(i, q, z),
            mixed_term(i, q, z),
            _scale_vector(param_phi_term(i, q, z), Decimal(-1)),
        )
        result = source(n, m, eta)
        return self._result(result, eta)


def wide_first_picard_angular_nonlinear_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardAngularNonlinearState:
    """Build the actual-schedule angular nonlinear layer at ``x1``."""

    return ActualScheduleWideFirstPicardAngularNonlinearState(x1=x1)


# Short alias for callers that do not need the first-Picard prefix in the type
# name.  The long name remains the canonical public API.
MixedScaleAngularNonlinearCoefficientJet = MixedScaleFirstPicardAngularNonlinearCoefficientJet


__all__ = [
    "ActualScheduleWideFirstPicardAngularNonlinearState",
    "MixedScaleAngularNonlinearCoefficientJet",
    "MixedScaleFirstPicardAngularNonlinearCoefficientJet",
    "wide_first_picard_angular_nonlinear_state",
]
