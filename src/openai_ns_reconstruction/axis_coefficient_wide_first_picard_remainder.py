"""Complete coefficientwise ``naturalRemainder(x1)`` on the actual schedule.

The first Picard state is already represented by the split jets in
``axis_coefficient_wide_first_picard``.  This module composes the landed raw
``lin1``, ``quad1``, ``slow1``, ``lin2``, ``slow2``, and pressure layers at one
genuine ``x1`` and applies the outer ``inverseL`` and natural resolvent.

The scale channels stay symbolic in the coefficient numerators.  The
ordinary channels are ``Lambda^0`` through ``Lambda^-5``; normalized
pressure-linear channels are ``a^2 Lambda^0`` through ``a^2 Lambda^-4``; and
the one pressure-square channel is normalized ``a^4 Lambda^-3``.  The natural
resolvent is evaluated coefficientwise from

    R[n,m] = A[n,m] - sum_k binom(m,k) chi[0,k] R[n-1,m-k]
                         / (2*n*(n+1)),

with the exact row-zero source.  This finite radial recurrence is a local
coefficient construction.  It is not a global weighted-space convergence
certificate, a fixed point, or a paper-exact velocity construction.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, localcontext
import math
from typing import Callable, Literal

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_data import (
    ActualScheduleAxisCoefficientData,
    actual_schedule_axis_coefficient_data,
)
from .axis_coefficient_natural_operator import actual_schedule_chi_coefficient_state
from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_wide_first_picard import (
    ActualScheduleWideFirstPicardState,
)
from .axis_coefficient_wide_first_picard_angular_nonlinear import (
    ActualScheduleWideFirstPicardAngularNonlinearState,
    wide_first_picard_angular_nonlinear_state,
)
from .axis_coefficient_wide_first_picard_lin1 import (
    ActualScheduleWideFirstPicardLin1State,
    wide_first_picard_lin1_state,
)
from .axis_coefficient_wide_first_picard_lin2 import (
    ActualScheduleWideFirstPicardLin2State,
    wide_first_picard_lin2_state,
)
from .axis_coefficient_wide_first_picard_pressure import (
    ActualScheduleWideFirstPicardPressureState,
    wide_first_picard_pressure_state,
)
from .axis_coefficient_wide_first_picard_slow2 import (
    ActualScheduleWideFirstPicardSlow2State,
    wide_first_picard_slow2_state,
)


_DECIMAL_PRECISION = 96
_ORDINARY_WIDTH = 6
_PRESSURE_LINEAR_WIDTH = 5
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
    amplitude_power: int,
    Lambda: Decimal,
    inverse_lambda_power: int,
) -> SignedLogCoefficientJet:
    """Encode one normalized pressure channel without forming its magnitude."""

    _finite_decimal(numerator, "pressure numerator")
    _finite_decimal(amplitude_log, "amplitude_log")
    _finite_decimal(Lambda, "Lambda")
    if amplitude_power <= 0 or inverse_lambda_power < 0:
        raise ValueError("invalid signed-log scale powers")
    if Lambda <= 0:
        raise ValueError("Lambda must be positive")
    if numerator == 0:
        return SignedLogCoefficientJet.zero()
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        return SignedLogCoefficientJet(
            sign=1 if numerator > 0 else -1,
            log_scale=+(Decimal(amplitude_power) * amplitude_log),
            log_factor=+(
                abs(numerator).ln()
                - Decimal(inverse_lambda_power) * Lambda.ln()
            ),
        )


@dataclass(frozen=True)
class _RemainderVector:
    """One normalized mixed-scale source/resolvent coefficient."""

    ordinary: tuple[Decimal, ...]
    pressure_linear: tuple[Decimal, ...]
    pressure_square_inverse_lambda_cubed: Decimal

    def __post_init__(self) -> None:
        if not isinstance(self.ordinary, tuple) or len(self.ordinary) != _ORDINARY_WIDTH:
            raise ValueError("ordinary channels must have six entries")
        if not isinstance(self.pressure_linear, tuple) or len(self.pressure_linear) != _PRESSURE_LINEAR_WIDTH:
            raise ValueError("pressure-linear channels must have five entries")
        for index, value in enumerate(self.ordinary):
            _finite_decimal(value, f"ordinary[{index}]")
        for index, value in enumerate(self.pressure_linear):
            _finite_decimal(value, f"pressure_linear[{index}]")
        _finite_decimal(
            self.pressure_square_inverse_lambda_cubed,
            "pressure_square_inverse_lambda_cubed",
        )

    @classmethod
    def zero(cls) -> "_RemainderVector":
        return cls(
            ordinary=(Decimal(0),) * _ORDINARY_WIDTH,
            pressure_linear=(Decimal(0),) * _PRESSURE_LINEAR_WIDTH,
            pressure_square_inverse_lambda_cubed=Decimal(0),
        )


def _subtract_scaled(left: _RemainderVector, right: _RemainderVector, scalar: Decimal) -> _RemainderVector:
    """Return ``left - scalar * right`` in one 96-digit context."""

    _finite_decimal(scalar, "resolvent scalar")
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        return _RemainderVector(
            ordinary=tuple(
                +(a - scalar * b)
                for a, b in zip(left.ordinary, right.ordinary)
            ),
            pressure_linear=tuple(
                +(a - scalar * b)
                for a, b in zip(left.pressure_linear, right.pressure_linear)
            ),
            pressure_square_inverse_lambda_cubed=+(
                left.pressure_square_inverse_lambda_cubed
                - scalar * right.pressure_square_inverse_lambda_cubed
            ),
        )


def _add_shifted(
    destination: list[Decimal],
    values: tuple[Decimal, ...],
    *,
    shift: int = 0,
    sign: int = 1,
) -> None:
    """Add a finite channel tuple after an explicit scale-power shift."""

    if sign not in (-1, 1):
        raise ValueError("channel sign must be plus or minus one")
    for index, value in enumerate(values):
        target = index + shift
        if not 0 <= target < len(destination):
            if value != 0:
                raise ValueError("nonzero channel would fall outside declared scale width")
            continue
        destination[target] += sign * value


@dataclass(frozen=True)
class MixedScaleFirstPicardNaturalRemainderCoefficientJet:
    """One angular or axial coefficient of ``naturalRemainder(x1)``.

    ``ordinary[p]`` is the numerator of ``Lambda^-p``.  The normalized
    ``pressure_linear[p]`` is the numerator of ``a(eta)^2 Lambda^-p``.  The
    pressure-square field is the numerator of ``a(eta)^4 Lambda^-3``.
    """

    ordinary: tuple[Decimal, ...]
    pressure_linear: tuple[Decimal, ...]
    pressure_square_inverse_lambda_cubed: Decimal
    Lambda: Decimal
    amplitude_log: Decimal

    def __post_init__(self) -> None:
        _RemainderVector(
            self.ordinary,
            self.pressure_linear,
            self.pressure_square_inverse_lambda_cubed,
        )
        _finite_decimal(self.Lambda, "Lambda")
        _finite_decimal(self.amplitude_log, "amplitude_log")
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    @property
    def ordinary_reference(self) -> Decimal:
        return self.ordinary[0]

    @property
    def ordinary_inverse_lambda_numerator(self) -> Decimal:
        return self.ordinary[1]

    @property
    def ordinary_inverse_lambda_squared_numerator(self) -> Decimal:
        return self.ordinary[2]

    @property
    def ordinary_inverse_lambda_cubed_numerator(self) -> Decimal:
        return self.ordinary[3]

    @property
    def ordinary_inverse_lambda_fourth_numerator(self) -> Decimal:
        return self.ordinary[4]

    @property
    def ordinary_inverse_lambda_fifth_numerator(self) -> Decimal:
        return self.ordinary[5]

    def ordinary_terms_decimal(self) -> tuple[Decimal, ...]:
        """Return the six ordinary terms, retaining local 96-digit arithmetic."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            inverse = +(Decimal(1) / self.Lambda)
            return tuple(
                +(value * (inverse**power if power else Decimal(1)))
                for power, value in enumerate(self.ordinary)
            )

    def ordinary_correction_terms_decimal(self) -> tuple[Decimal, ...]:
        """Return ordinary channels after the reference channel."""

        return self.ordinary_terms_decimal()[1:]

    def pressure_linear_terms_log(self) -> tuple[SignedLogCoefficientJet, ...]:
        """Return normalized ``a^2 Lambda^0..Lambda^-4`` signed-log terms."""

        return tuple(
            _signed_log_term(
                numerator,
                amplitude_log=self.amplitude_log,
                amplitude_power=2,
                Lambda=self.Lambda,
                inverse_lambda_power=power,
            )
            for power, numerator in enumerate(self.pressure_linear)
        )

    def pressure_square_term_log(self) -> SignedLogCoefficientJet:
        """Return the normalized ``a^4 Lambda^-3`` signed-log term."""

        return _signed_log_term(
            self.pressure_square_inverse_lambda_cubed,
            amplitude_log=self.amplitude_log,
            amplitude_power=4,
            Lambda=self.Lambda,
            inverse_lambda_power=3,
        )


@dataclass(frozen=True)
class SecondPicardCoefficientJet:
    """One coefficient of ``x2 = referencePair + naturalRemainder(x1)/(2 Lambda)``.

    The zeroth ordinary channel is the genuine reference coefficient.  Every
    remainder ordinary and normalized pressure numerator is shifted one scale
    power and divided by two.  Hence ``pressure_linear[0]`` is structurally
    zero and the pressure-square channel is normalized ``a^4 Lambda^-4``.
    """

    ordinary: tuple[Decimal, ...]
    pressure_linear: tuple[Decimal, ...]
    pressure_square_inverse_lambda_fourth: Decimal
    Lambda: Decimal
    amplitude_log: Decimal

    def __post_init__(self) -> None:
        if not isinstance(self.ordinary, tuple) or len(self.ordinary) != 7:
            raise ValueError("second-Picard ordinary channels must have seven entries")
        if not isinstance(self.pressure_linear, tuple) or len(self.pressure_linear) != 6:
            raise ValueError("second-Picard pressure-linear channels must have six entries")
        for index, value in enumerate(self.ordinary):
            _finite_decimal(value, f"ordinary[{index}]")
        for index, value in enumerate(self.pressure_linear):
            _finite_decimal(value, f"pressure_linear[{index}]")
        _finite_decimal(
            self.pressure_square_inverse_lambda_fourth,
            "pressure_square_inverse_lambda_fourth",
        )
        _finite_decimal(self.Lambda, "Lambda")
        _finite_decimal(self.amplitude_log, "amplitude_log")
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    def ordinary_terms_decimal(self) -> tuple[Decimal, ...]:
        """Return all seven ordinary terms at local 96-digit precision."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            inverse = +(Decimal(1) / self.Lambda)
            return tuple(
                +(value * (inverse**power if power else Decimal(1)))
                for power, value in enumerate(self.ordinary)
            )

    def ordinary_correction_terms_decimal(self) -> tuple[Decimal, ...]:
        return self.ordinary_terms_decimal()[1:]

    def pressure_linear_terms_log(self) -> tuple[SignedLogCoefficientJet, ...]:
        return tuple(
            _signed_log_term(
                numerator,
                amplitude_log=self.amplitude_log,
                amplitude_power=2,
                Lambda=self.Lambda,
                inverse_lambda_power=power,
            )
            for power, numerator in enumerate(self.pressure_linear)
        )

    def pressure_square_term_log(self) -> SignedLogCoefficientJet:
        return _signed_log_term(
            self.pressure_square_inverse_lambda_fourth,
            amplitude_log=self.amplitude_log,
            amplitude_power=4,
            Lambda=self.Lambda,
            inverse_lambda_power=4,
        )


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardRemainderState:
    """Complete coefficientwise ``naturalRemainder(x1)`` for one actual ``x1``."""

    x1: ActualScheduleWideFirstPicardState
    axis_data: ActualScheduleAxisCoefficientData = field(init=False, repr=False)
    chi: object = field(init=False, repr=False)
    lin1: ActualScheduleWideFirstPicardLin1State = field(init=False, repr=False)
    angular_nonlinear: ActualScheduleWideFirstPicardAngularNonlinearState = field(
        init=False,
        repr=False,
    )
    lin2: ActualScheduleWideFirstPicardLin2State = field(init=False, repr=False)
    slow2: ActualScheduleWideFirstPicardSlow2State = field(init=False, repr=False)
    pressure: ActualScheduleWideFirstPicardPressureState = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.x1, ActualScheduleWideFirstPicardState):
            raise TypeError("x1 must be ActualScheduleWideFirstPicardState")
        if not self.x1.picard_x1_materialized:
            raise ValueError("the genuine first Picard iterate must be materialized")
        if self.x1.fixed_point_materialized:
            raise ValueError("x1 must not be mislabeled as the final fixed point")
        _finite_decimal(self.x1.Lambda, "x1 Lambda")
        if self.x1.Lambda <= 0:
            raise ValueError("x1 Lambda must be positive")

        wide_pressure = self.x1.remainder.axial.wide_pressure
        if self.x1.Lambda != wide_pressure.amplitude.Lambda:
            raise ValueError("x1 Lambda must match the actual pressure amplitude scale")

        data = actual_schedule_axis_coefficient_data(self.x1.reference)
        if data.epsilon != self.x1.epsilon:
            raise ValueError("AxisData epsilon must match the first Picard state")
        if data.j != self.x1.reference.reference.j:
            raise ValueError("AxisData j must match the first Picard state")
        if data.sigma != self.x1.reference.reference.sigma:
            raise ValueError("AxisData sigma must match the first Picard state")

        # Construct every branch from this exact x1 object.  The constituent
        # factories repeat their own schedule and amplitude checks.
        branches = {
            "lin1": wide_first_picard_lin1_state(self.x1),
            "angular_nonlinear": wide_first_picard_angular_nonlinear_state(self.x1),
            "lin2": wide_first_picard_lin2_state(self.x1),
            "slow2": wide_first_picard_slow2_state(self.x1),
            "pressure": wide_first_picard_pressure_state(self.x1),
        }
        for name, branch in branches.items():
            if getattr(branch, "x1", None) is not self.x1:
                raise ValueError(f"{name} branch must use the exact same x1 object")
            if getattr(branch, "epsilon", None) != self.x1.epsilon:
                raise ValueError(f"{name} branch epsilon mismatch")
            if getattr(branch, "Lambda", self.x1.Lambda) != self.x1.Lambda:
                raise ValueError(f"{name} branch Lambda mismatch")

        object.__setattr__(self, "axis_data", data)
        object.__setattr__(self, "chi", actual_schedule_chi_coefficient_state(self.x1.reference))
        for name, branch in branches.items():
            object.__setattr__(self, name, branch)

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
    def natural_remainder_x1_materialized(self) -> bool:
        return True

    @property
    def naturalRemainder_x1_materialized(self) -> bool:
        return True

    @property
    def natural_remainder_complete(self) -> bool:
        return True

    @property
    def naturalRemainder_complete(self) -> bool:
        return True

    @property
    def second_picard_materialized(self) -> bool:
        return True

    @property
    def picard_x2_materialized(self) -> bool:
        return True

    @property
    def fixed_point_materialized(self) -> bool:
        return False

    @property
    def fixed_point_convergence_certified(self) -> bool:
        return False

    @property
    def amplitude(self):
        return self.x1.remainder.axial.wide_pressure.amplitude

    def amplitude_log(self, eta: float) -> Decimal:
        return _finite_decimal(
            self.amplitude.log_amplitude(_eta_in_window(eta)),
            "actual amplitude log",
        )

    @staticmethod
    def _validate_metadata(
        value: object,
        *,
        name: str,
        Lambda: Decimal,
        amplitude_log: Decimal,
    ) -> None:
        jet_Lambda = _finite_decimal(getattr(value, "Lambda", None), f"{name}.Lambda")
        if jet_Lambda != Lambda:
            raise ValueError(f"{name} Lambda mismatch")
        jet_amplitude_log = _finite_decimal(
            getattr(value, "amplitude_log", None),
            f"{name}.amplitude_log",
        )
        if jet_amplitude_log != amplitude_log:
            raise ValueError(f"{name} amplitude log mismatch")

    def _inverse_l_jet(self, m: int, eta: float) -> Decimal:
        value = _decimal_from_float(
            self.axis_data.inverseL.jet(0, m, eta),
            "AxisData.inverseL jet",
        )
        # inverseL is pinned radially constant; force a failure if a future
        # provider silently violates that structural fact.
        if self.axis_data.inverseL.jet(1, 0, eta) != 0.0:
            raise ValueError("AxisData.inverseL must have radial degree zero")
        return value

    def _multiply_inverse_l(
        self,
        source: Callable[[int, int], _RemainderVector],
        n: int,
        m: int,
        eta: float,
    ) -> _RemainderVector:
        ordinary = [Decimal(0)] * _ORDINARY_WIDTH
        pressure_linear = [Decimal(0)] * _PRESSURE_LINEAR_WIDTH
        pressure_square = Decimal(0)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for k in range(m + 1):
                weight = Decimal(math.comb(m, k)) * self._inverse_l_jet(k, eta)
                value = source(n, m - k)
                for index, item in enumerate(value.ordinary):
                    ordinary[index] += weight * item
                for index, item in enumerate(value.pressure_linear):
                    pressure_linear[index] += weight * item
                pressure_square += weight * value.pressure_square_inverse_lambda_cubed
            return _RemainderVector(
                ordinary=tuple(+value for value in ordinary),
                pressure_linear=tuple(+value for value in pressure_linear),
                pressure_square_inverse_lambda_cubed=+pressure_square,
            )

    def _raw_angular_source(self, n: int, m: int, eta: float) -> _RemainderVector:
        lin = self.lin1.jet(n, m, eta)
        quad = self.angular_nonlinear.quad1_jet(n, m, eta)
        slow = self.angular_nonlinear.slow1_jet(n, m, eta)
        amplitude_log = self.amplitude_log(eta)
        if _finite_decimal(lin.Lambda, "lin1.Lambda") != self.Lambda:
            raise ValueError("lin1 Lambda mismatch")
        self._validate_metadata(
            quad,
            name="quad1",
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )
        self._validate_metadata(
            slow,
            name="slow1",
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )

        ordinary = [Decimal(0)] * _ORDINARY_WIDTH
        pressure_linear = [Decimal(0)] * _PRESSURE_LINEAR_WIDTH
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            _add_shifted(
                ordinary,
                (
                    lin.reference,
                    lin.inverse_lambda_numerator,
                    lin.inverse_lambda_squared_numerator,
                ),
            )
            _add_shifted(ordinary, quad.ordinary)
            _add_shifted(ordinary, slow.ordinary, shift=1, sign=-1)
            # Nonlinear pressure tuples are indexed by their first physical
            # power a^2/Lambda, so retain that leading shift explicitly.
            _add_shifted(pressure_linear, quad.pressure_linear, shift=1)
            _add_shifted(pressure_linear, slow.pressure_linear, shift=2, sign=-1)
            value = _RemainderVector(
                ordinary=tuple(+item for item in ordinary),
                pressure_linear=tuple(+item for item in pressure_linear),
                pressure_square_inverse_lambda_cubed=Decimal(0),
            )
        return value

    def _raw_axial_source(self, n: int, m: int, eta: float) -> _RemainderVector:
        lin = self.lin2.jet(n, m, eta)
        slow = self.slow2.jet(n, m, eta)
        pressure = self.pressure.normalized_pressure_factors(n, m, eta)
        amplitude_log = self.amplitude_log(eta)
        self._validate_metadata(
            lin,
            name="lin2",
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )
        self._validate_metadata(
            slow,
            name="slow2",
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )
        if not isinstance(pressure, tuple) or len(pressure) != 5:
            raise ValueError("pressure normalized factors must have five channels")
        pressure = tuple(
            _finite_decimal(item, f"pressure.normalized_pressure_factors[{index}]")
            for index, item in enumerate(pressure)
        )

        ordinary = [Decimal(0)] * _ORDINARY_WIDTH
        pressure_linear = [Decimal(0)] * _PRESSURE_LINEAR_WIDTH
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            _add_shifted(
                ordinary,
                (
                    lin.ordinary_reference,
                    lin.ordinary_inverse_lambda_numerator,
                    lin.ordinary_inverse_lambda_squared_numerator,
                ),
            )
            _add_shifted(
                ordinary,
                (
                    slow.ordinary_reference,
                    slow.ordinary_inverse_lambda_numerator,
                    slow.ordinary_inverse_lambda_squared_numerator,
                    slow.ordinary_inverse_lambda_cubed_numerator,
                    slow.ordinary_inverse_lambda_fourth_numerator,
                ),
                shift=1,
                sign=-1,
            )
            # The pressure state supplies the a^2 Lambda^0..Lambda^-4
            # channels.  The raw lin2 pressure starts at Lambda^-1 and the
            # slow2 pressure is divided by one more Lambda.
            _add_shifted(pressure_linear, pressure)
            pressure_linear[1] += lin.pressure_linear_inverse_lambda_numerator
            _add_shifted(
                pressure_linear,
                (
                    slow.pressure_linear_inverse_lambda_numerator,
                    slow.pressure_linear_inverse_lambda_squared_numerator,
                    slow.pressure_linear_inverse_lambda_cubed_numerator,
                ),
                shift=2,
                sign=-1,
            )
            value = _RemainderVector(
                ordinary=tuple(+item for item in ordinary),
                pressure_linear=tuple(+item for item in pressure_linear),
                pressure_square_inverse_lambda_cubed=+(
                    -slow.pressure_square_inverse_lambda_squared_numerator
                ),
            )
        return value

    def _source(
        self,
        side: Literal["angular", "axial"],
        n: int,
        m: int,
        eta: float,
        cache: dict[tuple[str, int, int], _RemainderVector],
        raw_cache: dict[tuple[str, int, int], _RemainderVector],
    ) -> _RemainderVector:
        key = (side, n, m)
        if key not in cache:
            # The raw formulas are supplied at one coefficient coordinate;
            # inverseL is the common outer fixed field multiplication.
            def raw_provider(i: int, k: int) -> _RemainderVector:
                return self._raw_cached(side, i, k, eta, raw_cache)

            cache[key] = self._multiply_inverse_l(raw_provider, n, m, eta)
        return cache[key]

    def _raw_cached(
        self,
        side: Literal["angular", "axial"],
        n: int,
        m: int,
        eta: float,
        raw_cache: dict[tuple[str, int, int], _RemainderVector],
    ) -> _RemainderVector:
        key = (side, n, m)
        if key not in raw_cache:
            raw_cache[key] = (
                self._raw_angular_source(n, m, eta)
                if side == "angular"
                else self._raw_axial_source(n, m, eta)
            )
        return raw_cache[key]

    def _resolve(
        self,
        side: Literal["angular", "axial"],
        n: int,
        m: int,
        eta: float,
        source_cache: dict[tuple[str, int, int], _RemainderVector],
        raw_cache: dict[tuple[str, int, int], _RemainderVector],
        resolve_cache: dict[tuple[str, int, int], _RemainderVector],
        chi_cache: dict[int, Decimal],
    ) -> _RemainderVector:
        key = (side, n, m)
        if key in resolve_cache:
            return resolve_cache[key]
        source = self._source(side, n, m, eta, source_cache, raw_cache)
        if n == 0:
            resolve_cache[key] = source
            return source

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            result = source
            denominator = Decimal(2 * n * (n + 1))
            for k in range(m + 1):
                if k not in chi_cache:
                    chi_cache[k] = _decimal_from_float(
                        self.chi.jet(0, k, eta),
                        "chi coefficient jet",
                    )
                previous = self._resolve(
                    side,
                    n - 1,
                    m - k,
                    eta,
                    source_cache,
                    raw_cache,
                    resolve_cache,
                    chi_cache,
                )
                scalar = +(Decimal(math.comb(m, k)) * chi_cache[k] / denominator)
                result = _subtract_scaled(result, previous, scalar)
            resolve_cache[key] = result
            return result

    def _resolved_pair(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[_RemainderVector, _RemainderVector]:
        source_cache: dict[tuple[str, int, int], _RemainderVector] = {}
        raw_cache: dict[tuple[str, int, int], _RemainderVector] = {}
        resolve_cache: dict[tuple[str, int, int], _RemainderVector] = {}
        chi_cache: dict[int, Decimal] = {}
        return (
            self._resolve(
                "angular",
                n,
                m,
                eta,
                source_cache,
                raw_cache,
                resolve_cache,
                chi_cache,
            ),
            self._source("axial", n, m, eta, source_cache, raw_cache),
        )

    def _jet_from_vector(
        self,
        value: _RemainderVector,
        eta: float,
    ) -> MixedScaleFirstPicardNaturalRemainderCoefficientJet:
        return MixedScaleFirstPicardNaturalRemainderCoefficientJet(
            ordinary=value.ordinary,
            pressure_linear=value.pressure_linear,
            pressure_square_inverse_lambda_cubed=value.pressure_square_inverse_lambda_cubed,
            Lambda=self.Lambda,
            amplitude_log=self.amplitude_log(eta),
        )

    def jet_pair(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[
        MixedScaleFirstPicardNaturalRemainderCoefficientJet,
        MixedScaleFirstPicardNaturalRemainderCoefficientJet,
    ]:
        """Return angular and axial ``naturalRemainder(x1)`` coefficient jets."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        angular, axial = self._resolved_pair(n, m, eta)
        return self._jet_from_vector(angular, eta), self._jet_from_vector(axial, eta)

    def second_picard_jet_pair(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[SecondPicardCoefficientJet, SecondPicardCoefficientJet]:
        """Return coefficientwise ``x2 = referencePair + R(x1)/(2 Lambda)``."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        remainder_angular, remainder_axial = self._resolved_pair(n, m, eta)
        phi0, u0 = self.x1.reference.jet_pair(n, m, eta)
        amplitude_log = self.amplitude_log(eta)

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            half = Decimal(1) / _TWO
            angular = SecondPicardCoefficientJet(
                ordinary=(
                    _decimal_from_float(phi0, "reference phi coefficient jet"),
                )
                + tuple(+(half * value) for value in remainder_angular.ordinary),
                pressure_linear=(Decimal(0),)
                + tuple(+(half * value) for value in remainder_angular.pressure_linear),
                pressure_square_inverse_lambda_fourth=+(
                    half * remainder_angular.pressure_square_inverse_lambda_cubed
                ),
                Lambda=self.Lambda,
                amplitude_log=amplitude_log,
            )
            axial = SecondPicardCoefficientJet(
                ordinary=(
                    _decimal_from_float(u0, "reference u coefficient jet"),
                )
                + tuple(+(half * value) for value in remainder_axial.ordinary),
                pressure_linear=(Decimal(0),)
                + tuple(+(half * value) for value in remainder_axial.pressure_linear),
                pressure_square_inverse_lambda_fourth=+(
                    half * remainder_axial.pressure_square_inverse_lambda_cubed
                ),
                Lambda=self.Lambda,
                amplitude_log=amplitude_log,
            )
        return angular, axial


# The long names document the actual layer; short aliases keep downstream
# coefficient tests readable without introducing a second representation.
MixedScaleFirstPicardRemainderCoefficientJet = (
    MixedScaleFirstPicardNaturalRemainderCoefficientJet
)
MixedScaleNaturalRemainderCoefficientJet = MixedScaleFirstPicardNaturalRemainderCoefficientJet
MixedScaleSecondPicardCoefficientJet = SecondPicardCoefficientJet
ActualScheduleWideFirstPicardNaturalRemainderState = ActualScheduleWideFirstPicardRemainderState


def wide_first_picard_remainder_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardRemainderState:
    """Bind the complete coefficientwise natural remainder to one actual x1."""

    return ActualScheduleWideFirstPicardRemainderState(x1=x1)


wide_first_picard_natural_remainder_state = wide_first_picard_remainder_state


__all__ = [
    "ActualScheduleWideFirstPicardNaturalRemainderState",
    "ActualScheduleWideFirstPicardRemainderState",
    "MixedScaleFirstPicardNaturalRemainderCoefficientJet",
    "MixedScaleFirstPicardRemainderCoefficientJet",
    "MixedScaleNaturalRemainderCoefficientJet",
    "MixedScaleSecondPicardCoefficientJet",
    "SecondPicardCoefficientJet",
    "wide_first_picard_natural_remainder_state",
    "wide_first_picard_remainder_state",
]
