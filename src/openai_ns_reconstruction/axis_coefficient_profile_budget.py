"""Conditional actual-schedule norm budget for formal profile coefficients.

This module binds the existing actual SchedulePressure scale chain to the
formal coefficient solver.  It carries the theorem-shaped reference-pair
majorant and one-step fixed-point radius into a conditional profile norm
budget.  The remaining identification of the executable coefficient family
with the theorem's compatible ``AxisSpace`` fixed point is recorded explicitly
and is never promoted to a global norm certificate here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, ROUND_CEILING, localcontext

from .axis_coefficient_formal_solver import FormalAxisCoefficientSolverState
from .axis_coefficient_radial_tail import (
    AxisCoefficientRadialTailBound,
    axis_coefficient_radial_tail_bound,
)
from .axis_fixed_point_picard import NaturalPicardContractionCertificate
from .natural_scale_selection_wide import WideNaturalScaleSelection
from .stage1_scale_chain_wide import diagnose_actual_schedule_scale_chain_wide


_DECIMAL_PRECISION = 96
_PRODUCT_NORM = Decimal(64)
_PRIMITIVE_NORM = Decimal(80)


def _finite_decimal(value: Decimal, name: str) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError(f"{name} must be a finite Decimal")
    return value


def _add_up(left: Decimal, right: Decimal) -> Decimal:
    _finite_decimal(left, "left bound")
    _finite_decimal(right, "right bound")
    if left < 0 or right < 0:
        raise ValueError("profile bounds must be nonnegative")
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.rounding = ROUND_CEILING
        return +(left + right)


def _mul_up(left: Decimal, right: Decimal, name: str) -> Decimal:
    """Multiply nonnegative Decimal bounds with directed upward rounding."""

    _finite_decimal(left, f"{name} left bound")
    _finite_decimal(right, f"{name} right bound")
    if left < 0 or right < 0:
        raise ValueError(f"{name} bounds must be nonnegative")
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.rounding = ROUND_CEILING
        return +(left * right)


def _float_bound(value: float, name: str) -> Decimal:
    """Import an actual analytic float bound without decimal under-rounding."""

    if isinstance(value, bool) or not isinstance(value, float):
        raise TypeError(f"{name} must be a float")
    if not value == value or value in (float("inf"), float("-inf")):
        raise ValueError(f"{name} must be finite")
    result = Decimal.from_float(value)
    _finite_decimal(result, name)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


@dataclass(frozen=True)
class ActualScheduleProfileBudget:
    """Actual-parameter-derived conditional profile norm budget.

    Only ``solver`` is accepted at construction.  All scalar bounds are
    recomputed from its genuine SchedulePressure anchor and the landed wide
    theorem chain; caller-fitted norm inputs cannot be injected through this
    state.
    """

    solver: FormalAxisCoefficientSolverState
    scale: WideNaturalScaleSelection = field(init=False, repr=False)
    certificate: NaturalPicardContractionCertificate = field(init=False, repr=False)
    Lambda: Decimal = field(init=False)
    epsilon: Decimal = field(init=False)
    C_exponent: Decimal = field(init=False)
    reference_norm_upper: Decimal = field(init=False)
    fixed_point_distance_upper: Decimal = field(init=False)
    profile_norm_upper: Decimal = field(init=False)
    amplitude_norm_upper: Decimal = field(init=False)
    angular_ratio_norm_upper: Decimal = field(init=False)
    pressure_norm_upper: Decimal = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.solver, FormalAxisCoefficientSolverState):
            raise TypeError("solver must be FormalAxisCoefficientSolverState")

        x1 = self.solver.x1
        reference_pair = x1.reference.reference
        diagnostic = diagnose_actual_schedule_scale_chain_wide(
            reference_pair.data,
            reference_pair.j,
        )
        if diagnostic.upstream_obstruction is not None:
            raise RuntimeError(diagnostic.upstream_obstruction)
        if diagnostic.narrow_prefix.analytic_norms is None:
            raise RuntimeError("actual-schedule chain did not produce analytic norms")
        if diagnostic.remainder is None:
            raise RuntimeError("actual-schedule chain did not produce a wide remainder")
        if diagnostic.scale is None:
            raise RuntimeError("actual-schedule chain did not produce a wide scale")

        recomputed_scale = diagnostic.scale
        amplitude = x1.remainder.axial.wide_pressure.amplitude
        source_scale = amplitude.scale
        if recomputed_scale.Lambda != source_scale.Lambda:
            raise RuntimeError("recomputed Lambda does not match the actual amplitude scale")
        if recomputed_scale.C.exponent_upper != source_scale.C.exponent_upper:
            raise RuntimeError("recomputed symbolic C does not match the actual amplitude scale")
        if recomputed_scale.Lambda != self.solver.Lambda:
            raise RuntimeError("recomputed Lambda does not match the formal solver")

        analytic_epsilon = diagnostic.narrow_prefix.schedule_inputs.neighborhood.epsilon
        source_epsilon = amplitude.epsilon
        solver_epsilon = self.solver.epsilon
        if (
            analytic_epsilon != source_epsilon
            or analytic_epsilon != solver_epsilon
            or analytic_epsilon != x1.reference.epsilon
        ):
            raise RuntimeError("actual analytic epsilon does not match the formal solver")

        certificate = NaturalPicardContractionCertificate.from_scale(recomputed_scale)
        reference_norm = diagnostic.remainder.reference_norm_upper
        distance = certificate.one_step_radius_upper
        profile_norm = _add_up(reference_norm, distance)
        amplitude_norm = _float_bound(
            diagnostic.narrow_prefix.analytic_norms.amplitude_norm_upper,
            "amplitude norm",
        )
        if amplitude_norm <= 0:
            raise ValueError("amplitude norm must be positive")
        angular_ratio_norm = _mul_up(
            _mul_up(_PRODUCT_NORM, amplitude_norm, "amplitude-product norm"),
            profile_norm,
            "angular-ratio norm",
        )
        pressure_product_norm = _mul_up(
            angular_ratio_norm,
            angular_ratio_norm,
            "pressure source-product norm",
        )
        pressure_norm = _mul_up(
            _mul_up(_PRIMITIVE_NORM, _PRODUCT_NORM, "pressure primitive-product norm"),
            pressure_product_norm,
            "pressure norm",
        )
        _finite_decimal(recomputed_scale.Lambda, "Lambda")
        _finite_decimal(recomputed_scale.C.exponent_upper, "C exponent")
        _finite_decimal(reference_norm, "reference norm")
        _finite_decimal(distance, "fixed-point distance")
        object.__setattr__(self, "scale", recomputed_scale)
        object.__setattr__(self, "certificate", certificate)
        object.__setattr__(self, "Lambda", recomputed_scale.Lambda)
        object.__setattr__(self, "epsilon", Decimal.from_float(solver_epsilon))
        object.__setattr__(self, "C_exponent", recomputed_scale.C.exponent_upper)
        object.__setattr__(self, "reference_norm_upper", reference_norm)
        object.__setattr__(self, "fixed_point_distance_upper", distance)
        object.__setattr__(self, "profile_norm_upper", profile_norm)
        object.__setattr__(self, "amplitude_norm_upper", amplitude_norm)
        object.__setattr__(self, "angular_ratio_norm_upper", angular_ratio_norm)
        object.__setattr__(self, "pressure_norm_upper", pressure_norm)

    @property
    def conditional_on_axis_space_identification(self) -> bool:
        return True

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def includes_coefficient_roundoff(self) -> bool:
        return False

    @property
    def paper_exact(self) -> bool:
        return False

    def conditional_tail_bound(
        self,
        Y: Decimal,
        max_n: int,
        eta_order: int = 0,
        radial_order: int = 0,
        average: bool = False,
    ) -> AxisCoefficientRadialTailBound:
        """Delegate an omitted-row bound using this actual profile budget."""

        return axis_coefficient_radial_tail_bound(
            self.profile_norm_upper,
            self.epsilon,
            Y,
            max_n,
            eta_order=eta_order,
            radial_order=radial_order,
            average=average,
        )

    def conditional_angular_ratio_tail_bound(
        self,
        Y: Decimal,
        max_n: int,
        eta_order: int = 0,
        radial_order: int = 0,
        average: bool = False,
    ) -> AxisCoefficientRadialTailBound:
        """Bound omitted rows of the final smooth ratio ``F=a*phi``.

        The physical ``E=sqrt(2*X)*F`` factor is applied by the profile
        assembly and is intentionally outside this coefficient-space bound.
        """

        return axis_coefficient_radial_tail_bound(
            self.angular_ratio_norm_upper,
            self.epsilon,
            Y,
            max_n,
            eta_order=eta_order,
            radial_order=radial_order,
            average=average,
        )

    def conditional_pressure_tail_bound(
        self,
        Y: Decimal,
        max_n: int,
        eta_order: int = 0,
        radial_order: int = 0,
        average: bool = False,
    ) -> AxisCoefficientRadialTailBound:
        """Bound omitted rows of ``P=primitive(F^2)`` conditionally."""

        return axis_coefficient_radial_tail_bound(
            self.pressure_norm_upper,
            self.epsilon,
            Y,
            max_n,
            eta_order=eta_order,
            radial_order=radial_order,
            average=average,
        )


def actual_schedule_profile_budget(
    solver: FormalAxisCoefficientSolverState,
) -> ActualScheduleProfileBudget:
    """Build the actual-parameter-derived conditional profile budget."""

    return ActualScheduleProfileBudget(solver=solver)


__all__ = [
    "ActualScheduleProfileBudget",
    "actual_schedule_profile_budget",
]
