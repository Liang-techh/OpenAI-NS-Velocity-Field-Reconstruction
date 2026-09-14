"""The theorem-selected natural source ``a^2 * phi1^2`` at genuine first Picard x1.

The landed ``ActualScheduleWideFirstPicardState`` keeps the angular field as

    phi1 = r + b / Lambda + c / Lambda^2,

and ``axis_coefficient_wide_first_picard_angular_square`` already applies the
pinned radial convolution and eta-Leibniz rule to obtain ``phi1^2`` with powers
``Lambda^0`` through ``Lambda^-4`` kept separate.  The source term in
``AxisContraction.naturalRemainder`` is

    source = a^2 * phi^2,

where the actual theorem-selected amplitude is
``a(eta) = exp(Lambda * realPhase(eta)) / C``.  This module performs exactly
that missing x1 multiplication without narrowing ``a`` to binary64.

Because ``a`` has radial degree zero, only eta Leibniz factors are added.  Each
requested eta derivative of ``a^2`` is evaluated with the complete Bell
recurrence and the same actual normalized-gradient jets used by the landed
amplitude state.  The common ``a^2`` scale remains in signed-log form while the
five inverse-Lambda coefficient factors remain Decimal values.

This closes the source(x1) constituent only.  It does not materialize the x1
pressure chain, the angular/axial recombination, complete naturalRemainder(x1),
x2, a global AxisCoefficientSpace norm certificate, or a fixed point.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState
from .axis_coefficient_wide_first_picard_angular_square import (
    wide_first_picard_angular_square_state,
)


_DECIMAL_PRECISION = 96
_FACTOR_FIELDS = (
    "reference_factor",
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


def _decimal_from_float(value: float, name: str) -> Decimal:
    value = float(value)
    if not math.isfinite(value):
        raise ArithmeticError(f"{name} must be finite")
    return Decimal.from_float(value)


def _finite_decimal(value: Decimal, name: str) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError(f"{name} must be a finite Decimal")
    return value


def _signed_log_source_term(
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
    if Lambda <= 0:
        raise ValueError("Lambda must be positive")
    if not isinstance(inverse_lambda_power, int) or inverse_lambda_power < 0:
        raise ValueError("inverse_lambda_power must be a nonnegative integer")
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
class MixedScaleFirstPicardNaturalSourceCoefficientJet:
    """One coefficient jet of ``a^2 * phi1^2`` after extracting ``a^2``.

    The represented source coefficient is

      ``a^2 * (s0 + s1/Lambda + ... + s4/Lambda^4)``.

    The five ``s`` values are stored below; ``source_terms_log`` restores the
    common amplitude scale without ever forming ``a`` or ``a^2`` in binary64.
    """

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

    def normalized_factors_decimal(self) -> tuple[Decimal, Decimal, Decimal, Decimal, Decimal]:
        """Return the five factors before inverse-Lambda denominators are applied."""

        return tuple(getattr(self, name) for name in _FACTOR_FIELDS)

    def source_terms_log(self) -> tuple[SignedLogCoefficientJet, ...]:
        """Return all five source scales as signed-log coefficient jets."""

        return tuple(
            _signed_log_source_term(
                numerator,
                amplitude_log=self.amplitude_log,
                Lambda=self.Lambda,
                inverse_lambda_power=power,
            )
            for power, numerator in enumerate(self.normalized_factors_decimal())
        )


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardNaturalSourceState:
    """Typed genuine ``source(x1)=a^2*phi1^2`` on the actual schedule chain."""

    x1: ActualScheduleWideFirstPicardState

    def __post_init__(self) -> None:
        if not isinstance(self.x1, ActualScheduleWideFirstPicardState):
            raise TypeError("x1 must be ActualScheduleWideFirstPicardState")
        if not self.x1.picard_x1_materialized:
            raise ValueError("the genuine first Picard iterate must be materialized")
        if self.x1.fixed_point_materialized:
            raise ValueError("x1 must not be mislabeled as the final fixed point")
        if self.x1.Lambda != self.amplitude.Lambda:
            raise ValueError("x1 Lambda must match the theorem-selected amplitude")
        if self.x1.epsilon != self.amplitude.epsilon:
            raise ValueError("x1 epsilon must match the theorem-selected amplitude")

    @property
    def amplitude(self):
        return self.x1.remainder.axial.wide_pressure.amplitude

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
    def natural_source_x1_materialized(self) -> bool:
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

    def _amplitude_square_bell_factor(self, order: int, eta: float) -> Decimal:
        """Return ``d^order(a^2)/a^2`` from the pinned realGradient jets."""

        order = _index(order, "order")
        eta = _eta_in_window(eta)
        if order == 0:
            return Decimal(1)

        q: list[Decimal] = [Decimal(0)]
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for derivative_order in range(order):
                gradient_jet = self.amplitude.data.normalizedGradient.jet(
                    0, derivative_order, eta
                )
                q.append(
                    +(
                        Decimal(2)
                        * self.Lambda
                        * _decimal_from_float(gradient_jet, "realGradient jet")
                    )
                )

            bell = [Decimal(1)]
            for n in range(order):
                total = Decimal(0)
                for k in range(n + 1):
                    total += (
                        Decimal(math.comb(n, k))
                        * q[k + 1]
                        * bell[n - k]
                    )
                bell.append(+total)
            return bell[order]

    def jet(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> MixedScaleFirstPicardNaturalSourceCoefficientJet:
        """Evaluate one exact mixed-scale source jet at genuine first Picard x1."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        angular_square = wide_first_picard_angular_square_state(self.x1)

        totals = [Decimal(0) for _ in _FACTOR_FIELDS]
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for k in range(m + 1):
                square = angular_square.jet(n, m - k, eta)
                square_factors = (
                    square.reference,
                    square.inverse_lambda_numerator,
                    square.inverse_lambda_squared_numerator,
                    square.inverse_lambda_cubed_numerator,
                    square.inverse_lambda_fourth_numerator,
                )
                weight = Decimal(math.comb(m, k))
                amplitude_factor = self._amplitude_square_bell_factor(k, eta)
                for power, factor in enumerate(square_factors):
                    totals[power] += weight * amplitude_factor * factor
            totals = [+value for value in totals]
            amplitude_log = +self.amplitude.log_amplitude(eta)

        return MixedScaleFirstPicardNaturalSourceCoefficientJet(
            reference_factor=totals[0],
            inverse_lambda_numerator=totals[1],
            inverse_lambda_squared_numerator=totals[2],
            inverse_lambda_cubed_numerator=totals[3],
            inverse_lambda_fourth_numerator=totals[4],
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )


def wide_first_picard_natural_source_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardNaturalSourceState:
    """Bind the pinned natural source to one genuine first Picard state."""

    return ActualScheduleWideFirstPicardNaturalSourceState(x1=x1)
