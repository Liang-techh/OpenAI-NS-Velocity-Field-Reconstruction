"""Pinned first-Picard pressure chain for the genuine actual-schedule x1 source.

PR #323 materializes

    source(x1) = a^2 * phi1^2

as five separate inverse-Lambda coefficient families with the common nonzero
``a^2`` scale retained in signed-log form.  The pressure constituent in the
pinned ``AxisContraction.naturalRemainder`` is

    pressure = j1(
        -(4*A*eta) * primitive(source)
        + d * parameterPrimitive(source)
        -(2*eta) * mulY(source)
    ).

This module advances exactly that linear pressure chain at x1.  Each of the
five inverse-Lambda powers is propagated independently through the exact
coefficient-row maps from ``AxisOperators``:

* ``primitive``: output row n>0 is input row n-1 divided by n;
* ``parameterPrimitive``: the same radial map applied to the next eta jet;
* ``mulY``: output row n>0 is input row n-1;
* ``j1``: output row n>0 is input row n-1 divided by n^2.

Ordinary multiplication by ``eta`` and ``d`` uses the full finite radial
convolution and eta Leibniz sum.  The source amplitude is never collapsed into
binary64, no missing coefficient is filled with zero except where the pinned
row maps themselves define row zero to be zero, and callers cannot replace the
source, amplitude, Lambda, coefficient data, or pressure rows.

This materializes only ``pressure(x1)``.  It does not yet recombine the
remaining angular/axial terms into complete ``naturalRemainder(x1)``, produce
x2, certify a global AxisCoefficientSpace norm, or prove fixed-point closure.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math
from typing import Callable

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_reference_state import (
    AxisCoefficientJetState,
    WINDOW_LEFT,
    WINDOW_RIGHT,
)
from .axis_coefficient_wide_first_picard_source import (
    ActualScheduleWideFirstPicardNaturalSourceState,
)


_DECIMAL_PRECISION = 96
_FACTOR_FIELDS = (
    "reference_factor",
    "inverse_lambda_numerator",
    "inverse_lambda_squared_numerator",
    "inverse_lambda_cubed_numerator",
    "inverse_lambda_fourth_numerator",
)
_WideProvider = Callable[[int, int, float], Decimal]


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _power(value: int) -> int:
    value = _index(value, "inverse_lambda_power")
    if value >= len(_FACTOR_FIELDS):
        raise ValueError("inverse_lambda_power must lie in 0..4")
    return value


def _eta_in_window(value: float) -> float:
    eta = float(value)
    if not math.isfinite(eta) or not WINDOW_LEFT <= eta <= WINDOW_RIGHT:
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]")
    return eta


def _decimal_from_float(value: float, name: str) -> Decimal:
    value = float(value)
    if not math.isfinite(value):
        raise ArithmeticError(f"{name} must be finite")
    return Decimal.from_float(value)


def _finite_decimal(value: Decimal, name: str) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError(f"{name} must be a finite Decimal")
    return value


def _signed_log_pressure_term(
    numerator: Decimal,
    *,
    amplitude_log: Decimal,
    Lambda: Decimal,
    inverse_lambda_power: int,
) -> SignedLogCoefficientJet:
    """Encode ``a^2 * numerator / Lambda^p`` without scale collapse."""

    _finite_decimal(numerator, "numerator")
    _finite_decimal(amplitude_log, "amplitude_log")
    _finite_decimal(Lambda, "Lambda")
    inverse_lambda_power = _power(inverse_lambda_power)
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
class MixedScaleFirstPicardNaturalPressureCoefficientJet:
    """One pressure coefficient jet with all five x1 inverse-Lambda powers split."""

    reference_factor: Decimal
    inverse_lambda_numerator: Decimal
    inverse_lambda_squared_numerator: Decimal
    inverse_lambda_cubed_numerator: Decimal
    inverse_lambda_fourth_numerator: Decimal
    Lambda: Decimal
    amplitude_log: Decimal

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            _finite_decimal(value, name)
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    def normalized_factors_decimal(
        self,
    ) -> tuple[Decimal, Decimal, Decimal, Decimal, Decimal]:
        """Return factors before their explicit inverse-Lambda denominators."""

        return tuple(getattr(self, name) for name in _FACTOR_FIELDS)

    def pressure_terms_log(self) -> tuple[SignedLogCoefficientJet, ...]:
        """Return the five pressure scales as signed-log coefficient jets."""

        return tuple(
            _signed_log_pressure_term(
                numerator,
                amplitude_log=self.amplitude_log,
                Lambda=self.Lambda,
                inverse_lambda_power=power,
            )
            for power, numerator in enumerate(self.normalized_factors_decimal())
        )


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardNaturalPressureState:
    """The pinned pressure constituent evaluated on genuine first Picard x1."""

    source: ActualScheduleWideFirstPicardNaturalSourceState

    def __post_init__(self) -> None:
        if not isinstance(self.source, ActualScheduleWideFirstPicardNaturalSourceState):
            raise TypeError(
                "source must be ActualScheduleWideFirstPicardNaturalSourceState"
            )
        if not self.source.natural_source_x1_materialized:
            raise ValueError("the genuine source(x1) must be materialized")
        if self.source.x1.epsilon != self.axis_data.epsilon:
            raise ValueError("source(x1) and AxisData must use exactly matching epsilon")
        if self.source.Lambda != self.source.amplitude.Lambda:
            raise ValueError("source(x1) must retain the theorem-selected Lambda")

    @property
    def x1(self):
        return self.source.x1

    @property
    def axis_data(self):
        return self.source.x1.remainder.axial.axis_data

    @property
    def amplitude(self):
        return self.source.amplitude

    @property
    def epsilon(self) -> float:
        return self.source.epsilon

    @property
    def Lambda(self) -> Decimal:
        return self.source.Lambda

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def natural_source_x1_materialized(self) -> bool:
        return True

    @property
    def natural_pressure_x1_materialized(self) -> bool:
        return True

    @property
    def natural_remainder_x1_materialized(self) -> bool:
        return False

    @property
    def picard_x2_materialized(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    def source_factor(
        self,
        inverse_lambda_power: int,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        """Return one of the five normalized source factors ``s_p[n,m]``."""

        power = _power(inverse_lambda_power)
        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        return self.source.jet(n, m, eta).normalized_factors_decimal()[power]

    def primitive_factor(
        self,
        inverse_lambda_power: int,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        """Apply the pinned radial primitive to one source scale."""

        power = _power(inverse_lambda_power)
        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        if n == 0:
            self.source_factor(power, 0, m, eta)
            return Decimal(0)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(self.source_factor(power, n - 1, m, eta) / Decimal(n))

    def parameter_primitive_factor(
        self,
        inverse_lambda_power: int,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        """Apply ``primitive(partial_eta ·)`` to one source scale."""

        power = _power(inverse_lambda_power)
        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        if n == 0:
            self.source_factor(power, 0, m + 1, eta)
            return Decimal(0)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(
                self.source_factor(power, n - 1, m + 1, eta) / Decimal(n)
            )

    def multiply_y_factor(
        self,
        inverse_lambda_power: int,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        """Apply the pinned predecessor-row ``mulY`` map to one source scale."""

        power = _power(inverse_lambda_power)
        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        if n == 0:
            self.source_factor(power, 0, m, eta)
            return Decimal(0)
        return self.source_factor(power, n - 1, m, eta)

    def _ordinary_product_factor(
        self,
        ordinary: AxisCoefficientJetState,
        wide: _WideProvider,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        """Exact finite radial convolution and eta Leibniz rule."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        if not isinstance(ordinary, AxisCoefficientJetState):
            raise TypeError("ordinary factor must be AxisCoefficientJetState")
        if ordinary.epsilon != self.epsilon:
            raise ValueError("ordinary factor must use the actual source epsilon")

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            total = Decimal(0)
            for i in range(n + 1):
                j = n - i
                for k in range(m + 1):
                    total += (
                        Decimal(math.comb(m, k))
                        * wide(i, k, eta)
                        * _decimal_from_float(
                            ordinary.jet(j, m - k, eta),
                            "ordinary coefficient jet",
                        )
                    )
            return +total

    def pressure_input_factor(
        self,
        inverse_lambda_power: int,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        """Return the factor before the final pinned ``j1`` regular inverse."""

        power = _power(inverse_lambda_power)
        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        d = self.axis_data

        primitive = lambda i, k, x: self.primitive_factor(power, i, k, x)
        parameter_primitive = lambda i, k, x: self.parameter_primitive_factor(
            power, i, k, x
        )
        multiply_y = lambda i, k, x: self.multiply_y_factor(power, i, k, x)

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            eta_primitive = self._ordinary_product_factor(
                d.eta, primitive, n, m, eta
            )
            d_parameter_primitive = self._ordinary_product_factor(
                d.d, parameter_primitive, n, m, eta
            )
            eta_multiply_y = self._ordinary_product_factor(
                d.eta, multiply_y, n, m, eta
            )
            return +(
                -Decimal(4)
                * _decimal_from_float(d.A, "A")
                * eta_primitive
                + d_parameter_primitive
                - Decimal(2) * eta_multiply_y
            )

    def pressure_factor(
        self,
        inverse_lambda_power: int,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        """Apply final ``j1``: predecessor row divided by ``n^2``."""

        power = _power(inverse_lambda_power)
        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        if n == 0:
            self.pressure_input_factor(power, 0, m, eta)
            return Decimal(0)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(
                self.pressure_input_factor(power, n - 1, m, eta)
                / Decimal(n * n)
            )

    def jet(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> MixedScaleFirstPicardNaturalPressureCoefficientJet:
        """Return one five-scale pressure coefficient jet at genuine x1."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        factors = tuple(self.pressure_factor(p, n, m, eta) for p in range(5))
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            amplitude_log = +self.amplitude.log_amplitude(eta)
        return MixedScaleFirstPicardNaturalPressureCoefficientJet(
            reference_factor=factors[0],
            inverse_lambda_numerator=factors[1],
            inverse_lambda_squared_numerator=factors[2],
            inverse_lambda_cubed_numerator=factors[3],
            inverse_lambda_fourth_numerator=factors[4],
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )


def wide_first_picard_natural_pressure_state(
    source: ActualScheduleWideFirstPicardNaturalSourceState,
) -> ActualScheduleWideFirstPicardNaturalPressureState:
    """Bind the pinned pressure chain to the genuine x1 source only."""

    return ActualScheduleWideFirstPicardNaturalPressureState(source=source)
