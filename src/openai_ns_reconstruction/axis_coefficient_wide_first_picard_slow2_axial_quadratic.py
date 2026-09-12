"""First pressure-aware quadratic branch of pinned ``slow2`` at ``x1``.

The first term of ``AxisContraction.naturalRemainder.slow2`` is

    j1 (axialQuadraticCoefficient * (u * u)),

with ``axialQuadraticCoefficient = 2 * A * eta``.  The genuine theorem-scale
``u1 * u1`` product is already materialized as a split mixed-scale state.  This
module propagates that state through exactly this ordinary coefficient
multiplication and the zero-datum ``j1`` regular inverse without collapsing any
inverse-Lambda or signed-log pressure scale.

Because ``2*A*eta`` has radial degree zero and only its zeroth and first eta
derivatives are nonzero, the coefficient product is exact at each source row:

    d_eta^m[(2 A eta) F]
      = (2 A eta) d_eta^m F + m (2 A) d_eta^(m-1) F.

The pinned ``j1`` then maps source row ``n-1`` to output row ``n`` by division
by ``radialDivisor(1,n-1)=n^2`` and keeps row zero exactly zero.

No caller-supplied A, Lambda, amplitude, pressure/coefficient/derivative table,
or radial cutoff is accepted.  This is only the first quadratic branch of
``slow2(x1)``; it is not complete ``slow2``, ``naturalRemainder(x1)``, ``x2``,
or a fixed-point/paper-exact certificate.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_data import actual_schedule_axis_coefficient_data
from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState
from .axis_coefficient_wide_first_picard_axial_square import (
    MixedScaleFirstPicardAxialSquareCoefficientJet,
    wide_first_picard_axial_square_state,
)


_DECIMAL_PRECISION = 96


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
class MixedScaleFirstPicardSlow2AxialQuadraticCoefficientJet:
    """One jet of ``j1((2*A*eta) * (u1*u1))`` with theorem scales split."""

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
class ActualScheduleWideFirstPicardSlow2AxialQuadraticState:
    """The first ``slow2`` branch on the genuine theorem-scale first Picard state."""

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
    def slow2_axial_quadratic_branch_materialized(self) -> bool:
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

    def _A_decimal(self) -> Decimal:
        data = actual_schedule_axis_coefficient_data(self.x1.reference)
        if data.epsilon != self.epsilon:
            raise ValueError("AxisData epsilon must match the first Picard state")
        return Decimal.from_float(data.A)

    @staticmethod
    def _zero_like(
        source: MixedScaleFirstPicardAxialSquareCoefficientJet,
    ) -> MixedScaleFirstPicardSlow2AxialQuadraticCoefficientJet:
        return MixedScaleFirstPicardSlow2AxialQuadraticCoefficientJet(
            **{name: Decimal(0) for name in _FIELDS},
            Lambda=source.Lambda,
            amplitude_log=source.amplitude_log,
        )

    def jet(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> MixedScaleFirstPicardSlow2AxialQuadraticCoefficientJet:
        """Evaluate ``j1((2*A*eta)*(u1*u1))`` at one coefficient eta-jet."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        square_state = wide_first_picard_axial_square_state(self.x1)

        # Route row-zero through the real square state so all upstream
        # theorem-data/coordinate guards are still exercised.
        if n == 0:
            return self._zero_like(square_state.jet(0, m, eta))

        current = square_state.jet(n - 1, m, eta)
        previous = square_state.jet(n - 1, m - 1, eta) if m > 0 else None
        A = self._A_decimal()

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            q0 = +(Decimal(2) * A * Decimal.from_float(eta))
            q1_weight = +(Decimal(m) * Decimal(2) * A)
            divisor = Decimal(n * n)

            transformed: dict[str, Decimal] = {}
            for name in _FIELDS:
                value = q0 * getattr(current, name)
                if previous is not None:
                    value += q1_weight * getattr(previous, name)
                transformed[name] = +(value / divisor)

        return MixedScaleFirstPicardSlow2AxialQuadraticCoefficientJet(
            **transformed,
            Lambda=current.Lambda,
            amplitude_log=current.amplitude_log,
        )


def wide_first_picard_slow2_axial_quadratic_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardSlow2AxialQuadraticState:
    """Lift the first pinned ``slow2`` branch without theorem-scale collapse."""

    return ActualScheduleWideFirstPicardSlow2AxialQuadraticState(x1=x1)
