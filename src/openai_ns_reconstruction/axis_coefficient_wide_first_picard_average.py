"""Pinned radial average of the theorem-scale first Picard axial field.

The next evaluation of ``AxisContraction.naturalRemainder`` at the genuine
first Picard iterate ``x1`` starts with

    bu = coefficientOperators.average(u1).

The landed ``u1`` is intentionally represented as a mixed-scale decomposition,
not as one binary64 coefficient: an O(1) reference term, explicit
``1/Lambda`` and ``1/Lambda^2`` corrections, and a signed-log pressure
correction.  This module lifts the pinned radial average

    J_avg[n,m] = J[n,m] / (n + 1)

componentwise to that decomposition without collapsing any theorem-scale term.
No caller-supplied scale, Lambda, pressure value, coefficient table, or
surrogate state is accepted.

This is only the first mixed-scale operator needed for ``naturalRemainder(x1)``.
It does not yet materialize the products/bilinear terms, ``naturalRemainder(x1)``,
``x2``, or a certified fixed point.  The Stage-1 truth boundary therefore
remains ``formal-structure`` / not paper-exact.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState


_DECIMAL_PRECISION = 96


def _finite_decimal(value: Decimal, name: str) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError(f"{name} must be a finite Decimal")
    return value


def _average_signed_log(
    value: SignedLogCoefficientJet,
    radial_degree: int,
) -> SignedLogCoefficientJet:
    """Apply the exact positive factor ``1/(n+1)`` in split log form."""

    if not isinstance(value, SignedLogCoefficientJet):
        raise TypeError("value must be SignedLogCoefficientJet")
    if not isinstance(radial_degree, int) or radial_degree < 0:
        raise ValueError("radial_degree must be a nonnegative integer")
    if value.sign == 0:
        return SignedLogCoefficientJet.zero()

    assert value.log_scale is not None
    assert value.log_factor is not None
    divisor = Decimal(radial_degree + 1)
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        log_factor = +(value.log_factor - divisor.ln())
    return SignedLogCoefficientJet(
        sign=value.sign,
        log_scale=value.log_scale,
        log_factor=log_factor,
    )


@dataclass(frozen=True)
class MixedScaleFirstPicardAxialAverageCoefficientJet:
    """One coefficient jet of ``average(u1)`` with all scales preserved."""

    reference: Decimal
    inverse_lambda_numerator: Decimal
    inverse_lambda_squared_numerator: Decimal
    Lambda: Decimal
    pressure_over_two_lambda_average: SignedLogCoefficientJet

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
        if not isinstance(
            self.pressure_over_two_lambda_average,
            SignedLogCoefficientJet,
        ):
            raise TypeError(
                "pressure_over_two_lambda_average must be SignedLogCoefficientJet"
            )

    def ordinary_correction_terms_decimal(self) -> tuple[Decimal, Decimal]:
        """Return the two ordinary correction terms without summing into reference."""

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
class ActualScheduleWideFirstPicardAverageState:
    """Exact mixed-scale radial average ``b u1`` on the actual schedule datum."""

    x1: ActualScheduleWideFirstPicardState

    def __post_init__(self) -> None:
        if not isinstance(self.x1, ActualScheduleWideFirstPicardState):
            raise TypeError("x1 must be ActualScheduleWideFirstPicardState")
        if not self.x1.picard_x1_materialized:
            raise ValueError("the genuine first Picard iterate must be materialized")
        if self.x1.fixed_point_materialized:
            raise ValueError("x1 must not be mislabeled as the final fixed point")

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
    def mixed_scale_average_u1_materialized(self) -> bool:
        return True

    @property
    def natural_remainder_x1_materialized(self) -> bool:
        return False

    @property
    def fixed_point_materialized(self) -> bool:
        return False

    def jet(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> MixedScaleFirstPicardAxialAverageCoefficientJet:
        """Return one exact mixed-scale coefficient jet of pinned ``average(u1)``."""

        if not isinstance(n, int) or n < 0:
            raise ValueError("n must be a nonnegative integer")
        _, axial = self.x1.jet_pair(n, m, eta)
        divisor = Decimal(n + 1)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            reference = +(axial.reference / divisor)
            first_numerator = +(axial.inverse_lambda_numerator / divisor)
            second_numerator = +(
                axial.inverse_lambda_squared_numerator / divisor
            )

        return MixedScaleFirstPicardAxialAverageCoefficientJet(
            reference=reference,
            inverse_lambda_numerator=first_numerator,
            inverse_lambda_squared_numerator=second_numerator,
            Lambda=axial.Lambda,
            pressure_over_two_lambda_average=_average_signed_log(
                axial.pressure_over_two_lambda,
                n,
            ),
        )


def wide_first_picard_average_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardAverageState:
    """Lift the pinned radial average to the genuine theorem-scale ``u1`` state."""

    return ActualScheduleWideFirstPicardAverageState(x1=x1)
