"""Complete pinned ``slow2(x1)`` at the genuine first Picard state.

The four theorem-scale constituents of ``AxisContraction.naturalRemainder.slow2``
are already materialized separately on the actual SchedulePressure first Picard
state.  This module performs the missing typed composition only:

    j1 ((2*A*eta) * (u1*u1))
  + dot1 ((2*D*eta) * average(u1)) u1
  + d * mixed1 (average(u1)) u1
  - param1 u1 (d*u1).

Every constituent uses the same genuine ``ActualScheduleWideFirstPicardState``.
The composition keeps the existing theorem-scale representation: ordinary
``Lambda^0..Lambda^-4`` numerators, pressure-linear ``a^2 Lambda^-1..-3``
numerators, and the pressure-square ``a^4 Lambda^-2`` numerator.  No branch is
recomputed here, no caller may replace coefficient data, and no signed-log
pressure term is collapsed to binary64.

This closes ``slow2(x1)`` only.  It does not materialize the remaining angular
or source/pressure/resolvent pieces of ``naturalRemainder(x1)``, ``x2``, a
fixed point, ``NaturalProfileAssembly``, or paper-exact velocity.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState
from .axis_coefficient_wide_first_picard_slow2_axial_quadratic import (
    wide_first_picard_slow2_axial_quadratic_state,
)
from .axis_coefficient_wide_first_picard_slow2_average_dot import (
    wide_first_picard_slow2_average_dot_state,
)
from .axis_coefficient_wide_first_picard_slow2_average_mixed import (
    wide_first_picard_slow2_average_mixed_state,
)
from .axis_coefficient_wide_first_picard_slow2_param import (
    wide_first_picard_slow2_param_state,
)


_DECIMAL_PRECISION = 96
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
class MixedScaleFirstPicardSlow2CoefficientJet:
    """One coefficient eta-jet of the complete four-term ``slow2(x1)`` sum."""

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


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardSlow2State:
    """Typed complete ``slow2(x1)`` bound to one genuine theorem-selected ``x1``."""

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
    def slow2_all_constituent_branches_materialized(self) -> bool:
        return True

    @property
    def slow2_complete(self) -> bool:
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

    def _branch_states(self):
        branches = (
            wide_first_picard_slow2_axial_quadratic_state(self.x1),
            wide_first_picard_slow2_average_dot_state(self.x1),
            wide_first_picard_slow2_average_mixed_state(self.x1),
            wide_first_picard_slow2_param_state(self.x1),
        )
        if any(branch.x1 is not self.x1 for branch in branches):
            raise ValueError("all slow2 branches must share the identical first Picard state")
        if any(branch.Lambda != self.Lambda for branch in branches):
            raise ValueError("all slow2 branches must share the first Picard Lambda")
        if any(branch.epsilon != self.epsilon for branch in branches):
            raise ValueError("all slow2 branches must share the first Picard epsilon")
        return branches

    def jet(self, n: int, m: int, eta: float) -> MixedScaleFirstPicardSlow2CoefficientJet:
        """Sum the four already-materialized pinned ``slow2(x1)`` branches."""

        jets = tuple(branch.jet(n, m, eta) for branch in self._branch_states())
        amplitude_log = jets[0].amplitude_log
        if any(jet.Lambda != self.Lambda for jet in jets):
            raise ValueError("slow2 branch Lambda mismatch")
        if any(jet.amplitude_log != amplitude_log for jet in jets[1:]):
            raise ValueError("slow2 branches must share the same theorem amplitude")

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            values = {
                name: +sum((getattr(jet, name) for jet in jets), Decimal(0))
                for name in _FIELDS
            }
        return MixedScaleFirstPicardSlow2CoefficientJet(
            **values,
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )


def wide_first_picard_slow2_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardSlow2State:
    """Compose all four pinned first-Picard ``slow2`` branches."""

    return ActualScheduleWideFirstPicardSlow2State(x1=x1)
