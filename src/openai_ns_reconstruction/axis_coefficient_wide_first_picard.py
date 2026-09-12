"""Wide/mixed-scale first Picard iterate for the actual SchedulePressure datum.

The pinned natural-axis contraction uses the fixed-point map

    F(x) = x0 + (1 / (2 * Lambda)) * naturalRemainder(x),

with ``x0 = referencePair``.  The complete theorem-selected mixed-scale
``naturalRemainder(x0)`` is already materialized without narrowing its
``1/Lambda`` or signed-log ``a^2`` pieces.  This module closes exactly the next
seam: it applies the outer ``1/(2*Lambda)`` factor and adds ``x0`` while
preserving every scale as an explicit decomposition.

For one angular coefficient jet the represented value is

    phi0 + base/(2*Lambda) + slow_numerator/(2*Lambda^2),

and for one axial coefficient jet it is

    u0 + base/(2*Lambda) + slow_numerator/(2*Lambda^2)
       + pressure/(2*Lambda).

The reference value is deliberately *not* numerically added to the tiny
corrections at fixed Decimal precision: doing so would erase a mathematically
nonzero theorem-scale increment.  Likewise the pressure correction remains in
split signed-log form.  No caller-supplied sigma, epsilon, Lambda, C, amplitude,
coefficient table, or cutoff is accepted.

This materializes the first Picard iterate only.  It does not certify Picard
convergence, a final fixed point, global weighted ``AxisSpace`` membership, or a
paper-exact velocity profile.  The underlying phase value is still evaluated by
the landed numerical quadrature, so the truth status remains ``formal-structure``.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_reference_state import ActualScheduleReferenceAxisState
from .axis_coefficient_wide_natural_remainder import (
    ActualScheduleReferenceWideNaturalRemainderState,
    actual_schedule_reference_wide_natural_remainder_state,
)
from .outgoing_tail import TailData


_DECIMAL_PRECISION = 96
_TWO = Decimal(2)


def _finite_decimal(value: Decimal, name: str) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError(f"{name} must be a finite Decimal")
    return value


def _decimal_from_float(value: float, name: str) -> Decimal:
    value = float(value)
    if not math.isfinite(value):
        raise ArithmeticError(f"{name} must be finite")
    return Decimal.from_float(value)


def _half(value: Decimal) -> Decimal:
    _finite_decimal(value, "value")
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        return +(value / _TWO)


def _pressure_over_two_lambda(
    pressure: SignedLogCoefficientJet,
    Lambda: Decimal,
) -> SignedLogCoefficientJet:
    """Scale one signed-log pressure jet by ``1/(2*Lambda)`` losslessly.

    The enormous common amplitude log is left untouched.  The moderate
    ``log(2*Lambda)`` correction is applied to ``log_factor`` instead, avoiding
    cancellation against the current O(10^784) common log scale.
    """

    if not isinstance(pressure, SignedLogCoefficientJet):
        raise TypeError("pressure must be SignedLogCoefficientJet")
    _finite_decimal(Lambda, "Lambda")
    if Lambda <= 0:
        raise ValueError("Lambda must be positive")
    if pressure.sign == 0:
        return SignedLogCoefficientJet.zero()

    assert pressure.log_scale is not None
    assert pressure.log_factor is not None
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        log_divisor = +(_TWO * Lambda).ln()
        log_factor = +(pressure.log_factor - log_divisor)
    return SignedLogCoefficientJet(
        sign=pressure.sign,
        log_scale=pressure.log_scale,
        log_factor=log_factor,
    )


@dataclass(frozen=True)
class MixedScaleFirstPicardAngularCoefficientJet:
    """One angular coefficient jet of genuine theorem-scale ``x1``.

    The represented value is

    ``reference + inverse_lambda_numerator/Lambda
      + inverse_lambda_squared_numerator/Lambda^2``.

    The three pieces remain separate so the nonzero correction cannot disappear
    when the O(1) reference value is present.
    """

    reference: Decimal
    inverse_lambda_numerator: Decimal
    inverse_lambda_squared_numerator: Decimal
    Lambda: Decimal

    def __post_init__(self) -> None:
        _finite_decimal(self.reference, "reference")
        _finite_decimal(self.inverse_lambda_numerator, "inverse_lambda_numerator")
        _finite_decimal(
            self.inverse_lambda_squared_numerator,
            "inverse_lambda_squared_numerator",
        )
        _finite_decimal(self.Lambda, "Lambda")
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    def correction_terms_decimal(self) -> tuple[Decimal, Decimal]:
        """Return the two non-pressure corrections without adding to reference."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            first = +(self.inverse_lambda_numerator / self.Lambda)
            second = +(
                self.inverse_lambda_squared_numerator
                / self.Lambda
                / self.Lambda
            )
        return first, second


@dataclass(frozen=True)
class MixedScaleFirstPicardAxialCoefficientJet:
    """One axial coefficient jet of genuine theorem-scale ``x1``.

    The represented value is the reference, two explicit inverse-Lambda
    corrections, and the signed-log pressure correction.  No total-value
    method is exposed because a fixed-precision sum would erase the tiny terms.
    """

    reference: Decimal
    inverse_lambda_numerator: Decimal
    inverse_lambda_squared_numerator: Decimal
    Lambda: Decimal
    pressure_over_two_lambda: SignedLogCoefficientJet

    def __post_init__(self) -> None:
        _finite_decimal(self.reference, "reference")
        _finite_decimal(self.inverse_lambda_numerator, "inverse_lambda_numerator")
        _finite_decimal(
            self.inverse_lambda_squared_numerator,
            "inverse_lambda_squared_numerator",
        )
        _finite_decimal(self.Lambda, "Lambda")
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")
        if not isinstance(self.pressure_over_two_lambda, SignedLogCoefficientJet):
            raise TypeError("pressure_over_two_lambda must be SignedLogCoefficientJet")

    def ordinary_correction_terms_decimal(self) -> tuple[Decimal, Decimal]:
        """Return the two ordinary correction terms without scale collapse."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            first = +(self.inverse_lambda_numerator / self.Lambda)
            second = +(
                self.inverse_lambda_squared_numerator
                / self.Lambda
                / self.Lambda
            )
        return first, second


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardState:
    """The first Picard iterate ``x1 = x0 + R(x0)/(2*Lambda)``.

    ``remainder`` is the already-validated complete mixed-scale actual-schedule
    ``naturalRemainder(x0)`` pair, so this type inherits its exact agreement on
    TailData, ``j``, certified ``sigma``, epsilon, and theorem-selected Lambda.
    """

    remainder: ActualScheduleReferenceWideNaturalRemainderState

    def __post_init__(self) -> None:
        if not isinstance(
            self.remainder,
            ActualScheduleReferenceWideNaturalRemainderState,
        ):
            raise TypeError(
                "remainder must be ActualScheduleReferenceWideNaturalRemainderState"
            )
        if not self.remainder.mixed_scale_natural_remainder_x0_complete:
            raise ValueError("naturalRemainder(x0) must be complete")
        _finite_decimal(self.remainder.Lambda, "Lambda")
        if self.remainder.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    @property
    def reference(self) -> ActualScheduleReferenceAxisState:
        return self.remainder.reference

    @property
    def epsilon(self) -> float:
        return self.remainder.epsilon

    @property
    def Lambda(self) -> Decimal:
        return self.remainder.Lambda

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def picard_x1_materialized(self) -> bool:
        return True

    @property
    def fixed_point_materialized(self) -> bool:
        return False

    @property
    def fixed_point_convergence_certified(self) -> bool:
        return False

    def jet_pair(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[
        MixedScaleFirstPicardAngularCoefficientJet,
        MixedScaleFirstPicardAxialCoefficientJet,
    ]:
        """Return one mixed-scale coefficient-jet pair of the first Picard iterate."""

        angular_remainder, axial_remainder = self.remainder.jet_pair(n, m, eta)
        phi0, u0 = self.reference.jet_pair(n, m, eta)
        Lambda = self.Lambda

        angular = MixedScaleFirstPicardAngularCoefficientJet(
            reference=_decimal_from_float(phi0, "phi0 coefficient jet"),
            inverse_lambda_numerator=_half(angular_remainder.ordinary_base),
            inverse_lambda_squared_numerator=_half(
                angular_remainder.inverse_lambda_numerator
            ),
            Lambda=Lambda,
        )
        axial = MixedScaleFirstPicardAxialCoefficientJet(
            reference=_decimal_from_float(u0, "u0 coefficient jet"),
            inverse_lambda_numerator=_half(axial_remainder.ordinary_base),
            inverse_lambda_squared_numerator=_half(
                axial_remainder.inverse_lambda_numerator
            ),
            Lambda=Lambda,
            pressure_over_two_lambda=_pressure_over_two_lambda(
                axial_remainder.pressure,
                Lambda,
            ),
        )
        return angular, axial


def actual_schedule_wide_first_picard_state(
    data: TailData,
    j: float,
    *,
    phase_samples: int = 4001,
) -> ActualScheduleWideFirstPicardState:
    """Materialize the theorem-selected first Picard iterate on the actual schedule.

    The only caller inputs are the schedule datum and ``j`` already required by
    the landed construction.  All sigma/epsilon/Lambda/C/amplitude choices are
    inherited from the actual SchedulePressure theorem chain.
    """

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    remainder = actual_schedule_reference_wide_natural_remainder_state(
        data,
        j,
        phase_samples=phase_samples,
    )
    return ActualScheduleWideFirstPicardState(remainder=remainder)
