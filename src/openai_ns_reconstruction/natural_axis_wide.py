"""Sparse physical wide-prefix assembly for the formal profile.

The formal coefficient solver works in the radial coordinate ``Y = Lambda X``
and stores coefficient channels as ``(q, p)`` for ``a**q Lambda**(-p)``.
This module assembles the physical quantities while retaining those channels;
it deliberately does not project the result to binary64 or form amplitude
powers.  The physical chart is ``X >= 0`` and ``|eta| <= 1``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from decimal import Decimal, ROUND_CEILING, localcontext
from fractions import Fraction
import math
from types import MappingProxyType

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_formal_solver import FormalAxisCoefficientSolverState
from .axis_coefficient_mixed_scale import Channel, MixedScaleCoefficient
from .axis_coefficient_profile_budget import ActualScheduleProfileBudget
from .axis_coefficient_profile_prefix import formal_axis_profile_prefix
from .schedule_axis_pressure import axis_pressure


_DECIMAL_PRECISION = 96


def _finite_decimal(value: Decimal, name: str) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError(f"{name} must be a finite Decimal")
    return value


def _physical_eta(eta: float) -> float:
    if isinstance(eta, bool) or not isinstance(eta, (float, int)):
        raise ValueError("eta must be a finite real number")
    eta_float = float(eta)
    if not math.isfinite(eta_float) or abs(eta_float) > 1.0:
        raise ValueError("eta must lie in the physical chart |eta| <= 1")
    return eta_float


def _exact_decimal_product(left: Decimal, right: Decimal) -> Decimal:
    """Multiply finite Decimals with enough precision for an exact product."""

    left = _finite_decimal(left, "left")
    right = _finite_decimal(right, "right")
    precision = max(
        _DECIMAL_PRECISION,
        len(left.as_tuple().digits) + len(right.as_tuple().digits) + 2,
    )
    with localcontext() as context:
        context.prec = precision
        return +(left * right)


def _sqrt_two_x(X: Decimal) -> Decimal:
    radicand = _exact_decimal_product(Decimal(2), X)
    with localcontext() as context:
        context.prec = _DECIMAL_PRECISION
        return +radicand.sqrt()


def _ordinary(value: Decimal) -> MixedScaleCoefficient:
    return MixedScaleCoefficient.channel(0, 0, value)


def _decimal_from_float(value: float, name: str) -> Decimal:
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return Decimal.from_float(value)


def _transport_geometry(
    solver: FormalAxisCoefficientSolverState,
    eta: float,
) -> tuple[Decimal, Decimal, Decimal, Decimal]:
    """Return actual ``D``, ``h``, ``eta`` and positive ``L`` at Decimal96."""

    if not isinstance(solver, FormalAxisCoefficientSolverState):
        raise TypeError("solver must be a FormalAxisCoefficientSolverState")
    eta_float = _physical_eta(eta)
    D = _decimal_from_float(solver.data.D, "AxisData.D")
    h = _decimal_from_float(solver.data.h, "AxisData.h")
    eta_decimal = _decimal_from_float(eta_float, "eta")
    with localcontext() as context:
        context.prec = _DECIMAL_PRECISION
        L = +(Decimal(1) - Decimal(2) * h * eta_decimal * eta_decimal)
    if L <= 0:
        raise ValueError("physical chart requires L = 1 - 2*h*eta^2 > 0")
    return D, h, eta_decimal, L


def _fraction_to_decimal_ceiling(value: Fraction) -> Decimal:
    if value < 0:
        raise ValueError("bound must be nonnegative")
    with localcontext() as context:
        context.prec = _DECIMAL_PRECISION
        context.rounding = ROUND_CEILING
        return +(Decimal(value.numerator) / Decimal(value.denominator))


def _sqrt_two_x_upper(X: Decimal) -> Decimal:
    """Return a Decimal square-root upper enclosure verified by Fraction."""

    radicand = Fraction(2) * Fraction(X)
    candidate = _sqrt_two_x(X)
    for _ in range(4):
        if Fraction(candidate) * Fraction(candidate) >= radicand:
            return candidate
        with localcontext() as context:
            context.prec = _DECIMAL_PRECISION
            candidate = candidate.next_plus()
    raise ArithmeticError("could not enclose sqrt(2*X) upward")


@dataclass(frozen=True)
class ConditionalTruncationBounds(Mapping[str, Decimal]):
    """Immutable conditional omitted-tail bounds for the physical maps."""

    _values: Mapping[str, Decimal]
    conditional_on_axis_space_identification: bool = True
    includes_coefficient_roundoff: bool = False
    includes_axis_pressure_quadrature_roundoff: bool = False
    paper_exact: bool = False

    def __post_init__(self) -> None:
        required = {
            "F",
            "E",
            "U",
            "dU_deta",
            "average_U",
            "d_average_U_deta",
            "V0",
            "Pi",
        }
        values = dict(self._values)
        if set(values) != required:
            raise ValueError("conditional bounds must contain all eight physical fields")
        for name, value in values.items():
            if not isinstance(value, Decimal) or not value.is_finite() or value < 0:
                raise ValueError(f"conditional bound {name} must be nonnegative Decimal")
        if self.conditional_on_axis_space_identification is not True:
            raise ValueError("bounds are conditional on AxisSpace identification")
        if self.includes_coefficient_roundoff is not False:
            raise ValueError("coefficient roundoff is not included")
        if self.includes_axis_pressure_quadrature_roundoff is not False:
            raise ValueError("axis pressure quadrature roundoff is not included")
        if self.paper_exact is not False:
            raise ValueError("bounds are not paper-exact")
        object.__setattr__(self, "_values", MappingProxyType(values))

    def __getitem__(self, key: str) -> Decimal:
        return self._values[key]

    def __iter__(self):
        return iter(self._values)

    def __len__(self) -> int:
        return len(self._values)


@dataclass(frozen=True)
class WideNaturalProfilePrefix:
    """Formal physical profile fields retained as sparse mixed-scale maps."""

    F: MixedScaleCoefficient
    E: MixedScaleCoefficient
    U: MixedScaleCoefficient
    dU_deta: MixedScaleCoefficient
    average_U: MixedScaleCoefficient
    d_average_U_deta: MixedScaleCoefficient
    V0: MixedScaleCoefficient
    Pi: MixedScaleCoefficient
    X: Decimal
    Y: Decimal
    eta: float
    max_n: int
    Lambda: Decimal
    amplitude_log: Decimal
    _solver: FormalAxisCoefficientSolverState = field(repr=False, compare=False)

    def __post_init__(self) -> None:
        X = _finite_decimal(self.X, "X")
        Y = _finite_decimal(self.Y, "Y")
        Lambda = _finite_decimal(self.Lambda, "Lambda")
        _finite_decimal(self.amplitude_log, "amplitude_log")
        if X < 0:
            raise ValueError("X must be nonnegative")
        if Y < 0 or Y >= Decimal(20):
            raise ValueError("Y must satisfy 0 <= Y < 20")
        if Lambda <= 0:
            raise ValueError("Lambda must be positive")
        if isinstance(self.max_n, bool) or not isinstance(self.max_n, int) or self.max_n < 0:
            raise ValueError("max_n must be a nonnegative integer")
        _physical_eta(self.eta)
        _transport_geometry(self._solver, self.eta)
        for name in (
            "F",
            "E",
            "U",
            "dU_deta",
            "average_U",
            "d_average_U_deta",
            "V0",
            "Pi",
        ):
            if not isinstance(getattr(self, name), MixedScaleCoefficient):
                raise TypeError(f"{name} must be a MixedScaleCoefficient")

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def fixed_point_materialized(self) -> bool:
        return False

    @property
    def truncation_certified(self) -> bool:
        return False

    @property
    def coefficient_roundoff_certified(self) -> bool:
        return False

    @property
    def fields(self) -> Mapping[str, MixedScaleCoefficient]:
        return MappingProxyType(
            {
                "F": self.F,
                "E": self.E,
                "U": self.U,
                "dU_deta": self.dU_deta,
                "average_U": self.average_U,
                "d_average_U_deta": self.d_average_U_deta,
                "V0": self.V0,
                "Pi": self.Pi,
            }
        )

    @property
    def solver(self) -> FormalAxisCoefficientSolverState:
        return self._solver

    def conditional_truncation_bounds(self, budget) -> ConditionalTruncationBounds:
        """Return conditional omitted-row bounds from this prefix's solver.

        The budget must be the actual-parameter budget built from this exact
        solver.  Every conversion uses exact Fraction scaling and rounds only
        upward to the 96-digit Decimal representation.
        """

        if not isinstance(budget, ActualScheduleProfileBudget):
            raise TypeError("budget must be an ActualScheduleProfileBudget")
        if budget.solver is not self._solver:
            raise ValueError("budget must be built from this exact solver")
        Y = self.Y
        max_n = self.max_n
        f_tail = budget.conditional_angular_ratio_tail_bound(Y, max_n)
        u_tail = budget.conditional_tail_bound(Y, max_n)
        u_eta_tail = budget.conditional_tail_bound(Y, max_n, eta_order=1)
        average_tail = budget.conditional_tail_bound(Y, max_n, average=True)
        average_eta_tail = budget.conditional_tail_bound(
            Y,
            max_n,
            eta_order=1,
            average=True,
        )
        pressure_tail = budget.conditional_pressure_tail_bound(Y, max_n)

        def divide_lambda(bound: Decimal) -> Decimal:
            return _fraction_to_decimal_ceiling(
                Fraction(bound) / Fraction(self.Lambda)
            )

        sqrt_upper = _sqrt_two_x_upper(self.X)
        E_bound = _fraction_to_decimal_ceiling(
            Fraction(sqrt_upper) * Fraction(f_tail.upper_bound)
        )
        d_average_bound = divide_lambda(average_eta_tail.upper_bound)
        D, h, eta_decimal, L = _transport_geometry(self._solver, self.eta)
        if self.X == 0:
            V0_bound = Decimal(0)
        else:
            x_fraction = Fraction(self.X)
            eta_fraction = Fraction(eta_decimal)
            D_fraction = Fraction(D)
            h_fraction = Fraction(h)
            L_fraction = (
                Fraction(1)
                - Fraction(2) * h_fraction * eta_fraction * eta_fraction
            )
            if L_fraction <= 0:
                raise ValueError("physical chart requires rational L > 0")
            u_fraction = Fraction(u_tail.upper_bound) / Fraction(self.Lambda)
            average_fraction = (
                Fraction(average_tail.upper_bound) / Fraction(self.Lambda)
            )
            d_average_fraction = (
                Fraction(average_eta_tail.upper_bound) / Fraction(self.Lambda)
            )
            bracket = (
                abs(Fraction(2) * eta_fraction) * u_fraction
                + abs(Fraction(2) * D_fraction * eta_fraction) * average_fraction
                + abs(Fraction(1) - eta_fraction * eta_fraction) * d_average_fraction
            )
            V0_bound = _fraction_to_decimal_ceiling(
                abs(x_fraction) / L_fraction * bracket
            )
        return ConditionalTruncationBounds(
            {
                "F": f_tail.upper_bound,
                "E": E_bound,
                "U": divide_lambda(u_tail.upper_bound),
                "dU_deta": divide_lambda(u_eta_tail.upper_bound),
                "average_U": divide_lambda(average_tail.upper_bound),
                "d_average_U_deta": d_average_bound,
                "V0": V0_bound,
                "Pi": divide_lambda(pressure_tail.upper_bound),
            }
        )

    def terms_log(self, field_name: str) -> Mapping[Channel, SignedLogCoefficientJet]:
        """Return signed-log terms without collapsing sparse scale channels."""

        try:
            coefficient = self.fields[field_name]
        except KeyError as error:
            raise ValueError(f"unknown wide-profile field {field_name!r}") from error

        result: dict[Channel, SignedLogCoefficientJet] = {}
        for (q, p), numerator in coefficient.items():
            if numerator == 0:
                continue
            scale_precision = max(
                _DECIMAL_PRECISION,
                len(self.amplitude_log.as_tuple().digits) + len(str(q)) + 2,
            )
            with localcontext() as context:
                context.prec = scale_precision
                log_scale = +(Decimal(q) * self.amplitude_log)
            with localcontext() as context:
                context.prec = _DECIMAL_PRECISION
                log_factor = +(
                    abs(numerator).ln() - Decimal(p) * self.Lambda.ln()
                )
            result[(q, p)] = SignedLogCoefficientJet(
                sign=1 if numerator > 0 else -1,
                log_scale=log_scale,
                log_factor=log_factor,
            )
        return MappingProxyType(result)


def wide_natural_profile_prefix(
    solver: FormalAxisCoefficientSolverState,
    max_n: int,
    X: Decimal,
    eta: float,
) -> WideNaturalProfilePrefix:
    """Assemble a formal physical profile prefix in sparse mixed-scale form.

    ``formal_axis_profile_prefix`` supplies the value and eta-derivative jets
    in ``Y``.  The physical factor ``F = a phi`` is represented by shifting
    every angular channel by one amplitude power; no amplitude power is formed.
    """

    if not isinstance(solver, FormalAxisCoefficientSolverState):
        raise TypeError("solver must be a FormalAxisCoefficientSolverState")
    if isinstance(max_n, bool) or not isinstance(max_n, int) or max_n < 0:
        raise ValueError("max_n must be a nonnegative integer")
    X = _finite_decimal(X, "X")
    if X < 0:
        raise ValueError("X must be nonnegative")
    eta = _physical_eta(eta)

    Lambda = _finite_decimal(solver.Lambda, "solver.Lambda")
    if Lambda <= 0:
        raise ValueError("solver.Lambda must be positive")
    Y = _exact_decimal_product(Lambda, X)
    if Y >= Decimal(20):
        raise ValueError("physical radial coordinate must satisfy Y < 20")

    D, h, eta_decimal, L = _transport_geometry(solver, eta)

    profile = formal_axis_profile_prefix(solver, max_n, Y, eta)
    eta_profile = formal_axis_profile_prefix(
        solver,
        max_n,
        Y,
        eta,
        eta_order=1,
    )

    # The angular coefficient is the value jet only.  Multiplication by the
    # physical amplitude is represented by the odd q=1 channel shift.
    F = profile.angular.multiply(MixedScaleCoefficient.channel(1, 0, Decimal(1)))
    E = F.scale(_sqrt_two_x(X))

    reference = solver.x1.reference.reference
    j = _decimal_from_float(reference.j, "reference.j")
    with localcontext() as context:
        context.prec = _DECIMAL_PRECISION
        ordinary_u = _ordinary(+(Decimal(4) * eta_decimal + j))

    U = ordinary_u.add(profile.axial.shift_lambda(1))
    dU_deta = _ordinary(Decimal(4)).add(eta_profile.axial.shift_lambda(1))
    average_U = ordinary_u.add(profile.axial_average.shift_lambda(1))
    d_average_U_deta = _ordinary(Decimal(4)).add(
        eta_profile.axial_average.shift_lambda(1)
    )

    # ``V0`` is the regular radial-velocity numerator from the physical
    # profile formula.  The explicit axis branch keeps its structural zero.
    if X == 0:
        V0 = MixedScaleCoefficient.zero()
    else:
        with localcontext() as context:
            context.prec = _DECIMAL_PRECISION
            two_eta = +(Decimal(2) * eta_decimal)
            two_D_eta = +(Decimal(2) * D * eta_decimal)
            d_factor = +(Decimal(1) - eta_decimal * eta_decimal)
            x_over_L = +(X / L)
        V0 = (
            U.scale(two_eta)
            .add(average_U.scale(two_D_eta.copy_negate()))
            .add(d_average_U_deta.scale(d_factor.copy_negate()))
            .scale(x_over_L)
        )

    axis_pressure_value = _decimal_from_float(
        axis_pressure(reference.data, eta),
        "axis pressure",
    )
    Pi = _ordinary(axis_pressure_value).add(profile.pressure.shift_lambda(1))

    return WideNaturalProfilePrefix(
        F=F,
        E=E,
        U=U,
        dU_deta=dU_deta,
        average_U=average_U,
        d_average_U_deta=d_average_U_deta,
        V0=V0,
        Pi=Pi,
        X=X,
        Y=Y,
        eta=eta,
        max_n=max_n,
        Lambda=Lambda,
        amplitude_log=profile.amplitude_log,
        _solver=solver,
    )


__all__ = [
    "ConditionalTruncationBounds",
    "WideNaturalProfilePrefix",
    "wide_natural_profile_prefix",
]
