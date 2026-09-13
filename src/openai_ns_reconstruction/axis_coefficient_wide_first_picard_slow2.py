"""Complete mixed-scale ``slow2(x1)`` for the actual first Picard state.

The pinned natural-axis remainder has the axial slow term

    j1((2*A*eta) * (u*u))
      + dot1((2*D*eta) * average(u), u)
      + d * mixed1(average(u), u)
      - param1(u, d*u).

The four theorem-scale branches are materialized by the companion modules
``axis_coefficient_wide_first_picard_slow2_axial_quadratic``,
``..._average_dot``, ``..._average_mixed``, and ``..._param``.  This module
does only their exact typed assembly.  Every branch is evaluated on one genuine
``ActualScheduleWideFirstPicardState`` and its nine Decimal numerators are added
at the pinned 96-digit precision.  In particular, the parameter branch already
contains the displayed minus sign and is added without another sign change.

The ordinary terms remain split as ``Lambda^0`` through ``Lambda^-4``.  The
pressure-linear terms remain normalized ``a^2`` numerators for
``Lambda^-1`` through ``Lambda^-3`` and the pressure-square term remains the
normalized ``a^4 Lambda^-2`` numerator.  Signed-log views are rebuilt only
after these normalized numerators have been summed, so no nonzero theorem-scale
pressure contribution is narrowed through binary64.

This is a complete ``slow2(x1)`` value only.  It does not materialize the
remaining angular/axial pieces of ``naturalRemainder(x1)``, a later Picard
iterate, a fixed point, a global weighted ``AxisSpace`` certificate, or a
paper-exact velocity profile.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, localcontext
import math

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_reference_state import WINDOW_LEFT, WINDOW_RIGHT
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState
from .axis_coefficient_wide_first_picard_slow2_axial_quadratic import (
    ActualScheduleWideFirstPicardSlow2AxialQuadraticState,
    wide_first_picard_slow2_axial_quadratic_state,
)
from .axis_coefficient_wide_first_picard_slow2_average_dot import (
    ActualScheduleWideFirstPicardSlow2AverageDotState,
    wide_first_picard_slow2_average_dot_state,
)
from .axis_coefficient_wide_first_picard_slow2_average_mixed import (
    ActualScheduleWideFirstPicardSlow2AverageMixedState,
    wide_first_picard_slow2_average_mixed_state,
)
from .axis_coefficient_wide_first_picard_slow2_param import (
    ActualScheduleWideFirstPicardSlow2ParamState,
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
    """Encode one normalized numerator in the companion signed-log format."""

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
class MixedScaleFirstPicardSlow2CoefficientJet:
    """One complete coefficient jet of the pinned ``slow2(x1)`` sum."""

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
        """Return the four ordinary inverse-Lambda correction terms."""

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
        """Return the normalized ``a^2 Lambda^-1..Lambda^-3`` terms."""

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
        """Return the normalized ``a^4 Lambda^-2`` term."""

        return _signed_log_term(
            self.pressure_square_inverse_lambda_squared_numerator,
            amplitude_log=self.amplitude_log,
            amplitude_power=4,
            Lambda=self.Lambda,
            inverse_lambda_power=2,
        )


@dataclass(frozen=True)
class ActualScheduleWideFirstPicardSlow2State:
    """Complete mixed-scale ``slow2(x1)`` on one actual first Picard state."""

    x1: ActualScheduleWideFirstPicardState
    axial_quadratic: ActualScheduleWideFirstPicardSlow2AxialQuadraticState = field(
        init=False,
        repr=False,
    )
    average_dot: ActualScheduleWideFirstPicardSlow2AverageDotState = field(
        init=False,
        repr=False,
    )
    average_mixed: ActualScheduleWideFirstPicardSlow2AverageMixedState = field(
        init=False,
        repr=False,
    )
    param: ActualScheduleWideFirstPicardSlow2ParamState = field(
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
        if self.x1.Lambda != self.x1.remainder.axial.wide_pressure.amplitude.Lambda:
            raise ValueError("x1 Lambda must match the actual theorem-selected amplitude scale")
        _finite_decimal(self.x1.Lambda, "x1 Lambda")
        if self.x1.Lambda <= 0:
            raise ValueError("x1 Lambda must be positive")

        # All four factories receive this exact x1 object.  This prevents a
        # branch from a different schedule datum from entering the aggregate.
        object.__setattr__(
            self,
            "axial_quadratic",
            wide_first_picard_slow2_axial_quadratic_state(self.x1),
        )
        object.__setattr__(
            self,
            "average_dot",
            wide_first_picard_slow2_average_dot_state(self.x1),
        )
        object.__setattr__(
            self,
            "average_mixed",
            wide_first_picard_slow2_average_mixed_state(self.x1),
        )
        object.__setattr__(self, "param", wide_first_picard_slow2_param_state(self.x1))
        self._validate_branches()

    def _validate_branches(self) -> None:
        branches = (self.axial_quadratic, self.average_dot, self.average_mixed, self.param)
        for branch in branches:
            if branch.x1 is not self.x1:
                raise ValueError("all slow2 branches must use the exact same x1 object")
            if branch.epsilon != self.epsilon:
                raise ValueError("slow2 branch epsilon mismatch")
            if branch.Lambda != self.Lambda:
                raise ValueError("slow2 branch Lambda mismatch")
        if not self.axial_quadratic.slow2_axial_quadratic_branch_materialized:
            raise ValueError("axial quadratic slow2 branch is incomplete")
        if not self.average_dot.slow2_average_dot_branch_materialized:
            raise ValueError("average dot slow2 branch is incomplete")
        if not self.average_mixed.slow2_average_mixed_branch_materialized:
            raise ValueError("average mixed slow2 branch is incomplete")
        if not self.param.slow2_param_branch_materialized:
            raise ValueError("parameter slow2 branch is incomplete")

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
    def slow2_average_dot_branch_materialized(self) -> bool:
        return True

    @property
    def slow2_average_mixed_branch_materialized(self) -> bool:
        return True

    @property
    def slow2_param_branch_materialized(self) -> bool:
        return True

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
    def fixed_point_materialized(self) -> bool:
        return False

    @property
    def fixed_point_convergence_certified(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @staticmethod
    def _validate_jet_metadata(
        jet: object,
        *,
        branch_name: str,
        Lambda: Decimal,
        amplitude_log: Decimal,
    ) -> None:
        for name in _FIELDS:
            value = getattr(jet, name, None)
            _finite_decimal(value, f"{branch_name}.{name}")
        jet_Lambda = _finite_decimal(getattr(jet, "Lambda", None), f"{branch_name}.Lambda")
        if jet_Lambda != Lambda:
            raise ValueError(f"{branch_name} jet Lambda mismatch")
        jet_amplitude_log = _finite_decimal(
            getattr(jet, "amplitude_log", None),
            f"{branch_name}.amplitude_log",
        )
        if jet_amplitude_log != amplitude_log:
            raise ValueError(f"{branch_name} jet amplitude log mismatch")

    def jet(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> MixedScaleFirstPicardSlow2CoefficientJet:
        """Return one complete theorem-scale ``slow2(x1)`` coefficient jet."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        self._validate_branches()

        branch_jets = (
            (
                "axial quadratic",
                self.axial_quadratic.jet(n, m, eta),
            ),
            (
                "average dot",
                self.average_dot.jet(n, m, eta),
            ),
            (
                "average mixed",
                self.average_mixed.jet(n, m, eta),
            ),
            (
                "parameter",
                self.param.jet(n, m, eta),
            ),
        )
        amplitude_log = _finite_decimal(
            self.x1.remainder.axial.wide_pressure.amplitude.log_amplitude(eta),
            "actual amplitude log",
        )
        for branch_name, branch_jet in branch_jets:
            self._validate_jet_metadata(
                branch_jet,
                branch_name=branch_name,
                Lambda=self.Lambda,
                amplitude_log=amplitude_log,
            )

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            sums = {
                name: +sum(getattr(branch_jet, name) for _, branch_jet in branch_jets)
                for name in _FIELDS
            }

        return MixedScaleFirstPicardSlow2CoefficientJet(
            **sums,
            Lambda=self.Lambda,
            amplitude_log=amplitude_log,
        )


def wide_first_picard_slow2_state(
    x1: ActualScheduleWideFirstPicardState,
) -> ActualScheduleWideFirstPicardSlow2State:
    """Assemble all four pinned ``slow2(x1)`` branches from one genuine x1."""

    return ActualScheduleWideFirstPicardSlow2State(x1=x1)


__all__ = [
    "ActualScheduleWideFirstPicardSlow2State",
    "MixedScaleFirstPicardSlow2CoefficientJet",
    "wide_first_picard_slow2_state",
]
