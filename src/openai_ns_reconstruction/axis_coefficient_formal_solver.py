"""Finite-row formal coefficient solver on the generic mixed-scale algebra.

This module composes the pinned natural-axis coefficient formulas with the
sparse channel operations in :mod:`axis_coefficient_mixed_scale`.  A requested
row is evaluated by a local triangular recursion: the remainder row at radial
degree ``n`` uses only solved coefficient rows below ``n`` because every
``J_r`` shifts its source to ``n - 1``.  The local caches are discarded after
each ``jet_pair`` call.

The genuine first-Picard state supplies only the actual reference anchor,
epsilon, Lambda, amplitude logarithm, and the normalized amplitude-square Bell
factors.  Its finite iterate coefficient values are deliberately not used.
This is an executable formal coefficient composition, not a global
fixed-point, norm, convergence, or paper-exact certificate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from fractions import Fraction
import math
from functools import lru_cache
from typing import Callable

from .axis_coefficient_data import (
    ActualScheduleAxisCoefficientData,
    actual_schedule_axis_coefficient_data,
)
from .axis_coefficient_mixed_scale import (
    MixedScaleCoefficient,
    MixedScaleFamily,
    mixed_scale_average,
    mixed_scale_eta_derivative,
    mixed_scale_euler,
    mixed_scale_j_r,
    mixed_scale_primitive,
    mixed_scale_product,
)
from .axis_coefficient_reference_state import (
    WINDOW_LEFT,
    WINDOW_RIGHT,
)
from .axis_coefficient_rational_data import RationalAxisCoefficientData
from .axis_coefficient_wide_first_picard import ActualScheduleWideFirstPicardState


_DECIMAL_PRECISION = 96
_ONE = Decimal(1)
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


def _fraction_to_decimal96(value: Fraction, name: str) -> Decimal:
    """Round one exact rational scalar to nearest Decimal96 independently."""

    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be a Fraction")
    with localcontext() as context:
        context.prec = _DECIMAL_PRECISION
        context.rounding = ROUND_HALF_EVEN
        result = +(Decimal(value.numerator) / Decimal(value.denominator))
    return _finite_decimal(result, name)


def _family_add(*families: MixedScaleFamily) -> MixedScaleFamily:
    def added(n: int, m: int, eta: float) -> MixedScaleCoefficient:
        total = MixedScaleCoefficient.zero()
        for family in families:
            total = total.add(family(n, m, eta))
        return total

    return added


def _family_scale(source: MixedScaleFamily, scalar: Decimal) -> MixedScaleFamily:
    scalar = _finite_decimal(scalar, "family scalar")

    def scaled(n: int, m: int, eta: float) -> MixedScaleCoefficient:
        return source(n, m, eta).scale(scalar)

    return scaled


def _family_shift_lambda(source: MixedScaleFamily, shift: int) -> MixedScaleFamily:
    def shifted(n: int, m: int, eta: float) -> MixedScaleCoefficient:
        return source(n, m, eta).shift_lambda(shift)

    return shifted


def _param_r(left: MixedScaleFamily, right: MixedScaleFamily, r: int) -> MixedScaleFamily:
    return mixed_scale_j_r(
        mixed_scale_product(mixed_scale_eta_derivative(left), right),
        r,
    )


def _dot_r(left: MixedScaleFamily, right: MixedScaleFamily, r: int) -> MixedScaleFamily:
    return mixed_scale_j_r(
        mixed_scale_product(left, mixed_scale_euler(right)),
        r,
    )


def _mixed_r(left: MixedScaleFamily, right: MixedScaleFamily, r: int) -> MixedScaleFamily:
    return mixed_scale_j_r(
        mixed_scale_product(
            mixed_scale_eta_derivative(left),
            mixed_scale_euler(right),
        ),
        r,
    )


def _mul_y(source: MixedScaleFamily) -> MixedScaleFamily:
    """Build the pinned predecessor-row multiplication by ``Y``."""

    def multiply(n: int, m: int, eta: float) -> MixedScaleCoefficient:
        n = _index(n, "n")
        _index(m, "m")
        if n == 0:
            return MixedScaleCoefficient.zero()
        return source(n - 1, m, eta)

    return multiply


@dataclass(frozen=True)
class FormalAxisCoefficientSolverState:
    """Locally recursive formal coefficient state anchored by genuine ``x1``."""

    x1: ActualScheduleWideFirstPicardState
    data: ActualScheduleAxisCoefficientData = field(init=False, repr=False)
    rational_data: RationalAxisCoefficientData = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.x1, ActualScheduleWideFirstPicardState):
            raise TypeError("x1 must be ActualScheduleWideFirstPicardState")
        if not self.x1.picard_x1_materialized:
            raise ValueError("the genuine first Picard iterate must be materialized")
        if self.x1.fixed_point_materialized:
            raise ValueError("x1 must not be mislabeled as the final fixed point")

        amplitude = self.x1.remainder.axial.wide_pressure.amplitude
        source = self.x1.remainder.axial.wide_pressure.source
        _finite_decimal(self.x1.Lambda, "x1 Lambda")
        if self.x1.Lambda <= 0:
            raise ValueError("x1 Lambda must be positive")
        if amplitude.Lambda != self.x1.Lambda or source.amplitude.Lambda != self.x1.Lambda:
            raise ValueError("x1 Lambda must match the actual amplitude source")
        if amplitude.epsilon != self.x1.epsilon or source.epsilon != self.x1.epsilon:
            raise ValueError("amplitude source epsilon must match x1")
        _finite_decimal(amplitude.log_amplitude(0.0), "anchor amplitude log")

        data = actual_schedule_axis_coefficient_data(self.x1.reference)
        if data.epsilon != self.x1.epsilon:
            raise ValueError("AxisData epsilon must match x1")
        if data.j != self.x1.reference.reference.j:
            raise ValueError("AxisData j must match x1")
        if data.sigma != self.x1.reference.reference.sigma:
            raise ValueError("AxisData sigma must match x1")
        object.__setattr__(self, "data", data)
        object.__setattr__(self, "rational_data", RationalAxisCoefficientData(data))

    @property
    def Lambda(self) -> Decimal:
        return self.x1.Lambda

    @property
    def epsilon(self) -> float:
        return self.x1.epsilon

    def amplitude_log(self, eta: float) -> Decimal:
        eta = _eta_in_window(eta)
        amplitude_log = self.x1.remainder.axial.wide_pressure.amplitude.log_amplitude(eta)
        return _finite_decimal(amplitude_log, "anchor amplitude log")

    @property
    def formal_coefficients_materialized(self) -> bool:
        return True

    @property
    def fixed_point_materialized(self) -> bool:
        return False

    @property
    def fixed_point_convergence_certified(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def natural_remainder_x1_materialized(self) -> bool:
        return False

    def reference_jet_pair(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[Decimal, Decimal]:
        """Return Decimal96 angular reference and converted axial reference jets."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        with localcontext() as context:
            context.prec = _DECIMAL_PRECISION
            angular = self.rational_data.angular_reference_jet_decimal(n, m, eta)
        axial = _decimal_from_float(
            self.x1.reference.u.jet(n, m, eta),
            "reference u jet",
        )
        return _finite_decimal(angular, "reference phi jet"), axial

    def _fixed_family(
        self,
        name: str,
        scalar: Decimal = _ONE,
    ) -> MixedScaleFamily:
        rational_data = self.rational_data
        scalar = _finite_decimal(scalar, f"AxisData.{name} scalar")

        def fixed(n: int, m: int, eta: float) -> MixedScaleCoefficient:
            n = _index(n, "n")
            m = _index(m, "m")
            eta = _eta_in_window(eta)
            value = rational_data.jet_decimal(name, n, m, eta)
            if n != 0:
                if value != 0:
                    raise ValueError(f"AxisData.{name} must have radial degree zero")
                return MixedScaleCoefficient.zero()
            with localcontext() as ctx:
                ctx.prec = _DECIMAL_PRECISION
                return MixedScaleCoefficient.channel(0, 0, +(scalar * value))

        return fixed

    def _build_families(
        self,
        eta: float,
        input_pair: tuple[MixedScaleFamily, MixedScaleFamily] | None = None,
    ) -> tuple[
        MixedScaleFamily,
        MixedScaleFamily,
        Callable[[int, int], MixedScaleCoefficient],
        Callable[[int, int], MixedScaleCoefficient],
        MixedScaleFamily,
    ]:
        """Build one eta-local graph for ``T(input) = reference + R(input)/(2*Lambda)``.

        With ``input_pair`` omitted, the input families are the output
        families themselves and the historical triangular recursion is
        preserved.  Supplying a pair routes every nonlinear remainder operand
        through that pair while the returned families remain the updated
        reference-plus-remainder output.
        """

        eta = _eta_in_window(eta)
        if input_pair is not None:
            if not isinstance(input_pair, tuple) or len(input_pair) != 2:
                raise TypeError("input_pair must be a (phi_family, u_family) tuple")
            external_phi, external_u = input_pair
            if not callable(external_phi) or not callable(external_u):
                raise TypeError("input_pair entries must be callable mixed-scale families")
        else:
            external_phi = None
            external_u = None
        half = +(Decimal(1) / _TWO)

        d = self._fixed_family("d")
        normalized_gradient = self._fixed_family("normalizedGradient")
        eta_field = self._fixed_family("eta")
        one = self._fixed_family("one")
        u_star = self._fixed_family("uStar")
        u_star_eta = self._fixed_family("uStarEta")
        w_star = self._fixed_family("wStar")
        h_star = self._fixed_family("hStar")
        inverse_l = self._fixed_family("inverseL")

        A = _fraction_to_decimal96(self.rational_data.A, "RationalAxisData.A")
        D = _fraction_to_decimal96(self.rational_data.D, "RationalAxisData.D")
        h = _fraction_to_decimal96(self.rational_data.h, "RationalAxisData.h")
        average_coefficient = _family_scale(eta_field, +(Decimal(2) * D))
        angular_slow = _family_scale(eta_field, +(Decimal(2) * h))
        axial_quadratic = _family_scale(eta_field, +(Decimal(2) * A))
        angular_linear = _family_add(
            w_star,
            _family_scale(one, h),
            _family_scale(
                mixed_scale_product(eta_field, u_star),
                -(Decimal(2) * h),
            ),
        )
        axial_linear = _family_add(
            _family_scale(one, A),
            _family_scale(
                mixed_scale_product(eta_field, u_star),
                -(Decimal(4) * A),
            ),
            mixed_scale_product(d, u_star_eta),
        )
        angular_quadratic = mixed_scale_product(d, normalized_gradient)

        def amplitude_square(n: int, m: int, z: float) -> MixedScaleCoefficient:
            n = _index(n, "n")
            m = _index(m, "m")
            z = _eta_in_window(z)
            if n != 0:
                return MixedScaleCoefficient.zero()
            bell = self.rational_data.amplitude_power_bell_decimal(
                2,
                m,
                z,
                self.Lambda,
            )
            return MixedScaleCoefficient.channel(
                2,
                0,
                _finite_decimal(bell, "a^2 Bell factor"),
            )

        def _local_eta(z: float) -> None:
            if _eta_in_window(z) != eta:
                raise ValueError("formal solver family graph requires one fixed eta")

        @lru_cache(maxsize=None)
        def external_phi_value(n: int, m: int) -> MixedScaleCoefficient:
            if external_phi is None:
                return get_phi(n, m)
            value = external_phi(n, m, eta)
            if not isinstance(value, MixedScaleCoefficient):
                raise TypeError("external phi family must return MixedScaleCoefficient")
            return value

        @lru_cache(maxsize=None)
        def external_u_value(n: int, m: int) -> MixedScaleCoefficient:
            if external_u is None:
                return get_u(n, m)
            value = external_u(n, m, eta)
            if not isinstance(value, MixedScaleCoefficient):
                raise TypeError("external u family must return MixedScaleCoefficient")
            return value

        def input_phi_family(n: int, m: int, z: float) -> MixedScaleCoefficient:
            _local_eta(z)
            return external_phi_value(_index(n, "n"), _index(m, "m"))

        def input_u_family(n: int, m: int, z: float) -> MixedScaleCoefficient:
            _local_eta(z)
            return external_u_value(_index(n, "n"), _index(m, "m"))

        bu = mixed_scale_average(input_u_family)
        average_transport = mixed_scale_product(average_coefficient, bu)
        angular_transport = mixed_scale_product(angular_slow, input_u_family)
        transport = _family_add(average_transport, angular_transport)

        lin1 = _family_add(
            mixed_scale_j_r(mixed_scale_product(angular_linear, input_phi_family), 2),
            _dot_r(w_star, input_phi_family, 2),
            _param_r(input_phi_family, h_star, 2),
        )
        quad1 = mixed_scale_j_r(
            mixed_scale_product(
                mixed_scale_product(angular_quadratic, input_u_family),
                input_phi_family,
            ),
            2,
        )
        slow1 = _family_add(
            mixed_scale_j_r(mixed_scale_product(transport, input_phi_family), 2),
            _param_r(bu, mixed_scale_product(d, input_phi_family), 2),
            _dot_r(average_transport, input_phi_family, 2),
            mixed_scale_product(d, _mixed_r(bu, input_phi_family, 2)),
            _family_scale(_param_r(input_phi_family, mixed_scale_product(d, input_u_family), 2), -_ONE),
        )

        lin2 = _family_add(
            mixed_scale_j_r(mixed_scale_product(axial_linear, input_u_family), 1),
            _dot_r(w_star, input_u_family, 1),
            _param_r(input_u_family, h_star, 1),
        )
        slow2 = _family_add(
            mixed_scale_j_r(
                mixed_scale_product(
                    axial_quadratic,
                    mixed_scale_product(input_u_family, input_u_family),
                ),
                1,
            ),
            _dot_r(average_transport, input_u_family, 1),
            mixed_scale_product(d, _mixed_r(bu, input_u_family, 1)),
            _family_scale(_param_r(input_u_family, mixed_scale_product(d, input_u_family), 1), -_ONE),
        )

        source = mixed_scale_product(
            amplitude_square,
            mixed_scale_product(input_phi_family, input_phi_family),
        )
        primitive_source = mixed_scale_primitive(source)
        eta_primitive_source = mixed_scale_eta_derivative(primitive_source)
        pressure_input = _family_add(
            _family_scale(
                mixed_scale_product(eta_field, primitive_source),
                -(Decimal(4) * A),
            ),
            mixed_scale_product(d, eta_primitive_source),
            _family_scale(mixed_scale_product(eta_field, _mul_y(source)), -_TWO),
        )
        pressure = mixed_scale_j_r(pressure_input, 1)

        angular_raw = mixed_scale_product(
            inverse_l,
            _family_add(
                lin1,
                quad1,
                _family_scale(_family_shift_lambda(slow1, 1), -_ONE),
            ),
        )
        axial_raw = mixed_scale_product(
            inverse_l,
            _family_add(
                lin2,
                _family_scale(_family_shift_lambda(slow2, 1), -_ONE),
                pressure,
            ),
        )
        @lru_cache(maxsize=None)
        def get_angular_r(n: int, m: int) -> MixedScaleCoefficient:
            n = _index(n, "n")
            m = _index(m, "m")
            if n == 0:
                return MixedScaleCoefficient.zero()
            with localcontext() as ctx:
                ctx.prec = _DECIMAL_PRECISION
                divisor = Decimal(2 * n * (n + 1))
                result = angular_raw(n, m, eta)
                for k in range(m + 1):
                    chi = self.rational_data.jet_decimal("chi", 0, k, eta)
                    correction = get_angular_r(n - 1, m - k).scale(
                        +(Decimal(math.comb(m, k)) * chi / divisor)
                    )
                    result = result.add(correction.scale(-_ONE))
                return result

        @lru_cache(maxsize=None)
        def get_axial_r(n: int, m: int) -> MixedScaleCoefficient:
            n = _index(n, "n")
            m = _index(m, "m")
            if n == 0:
                return MixedScaleCoefficient.zero()
            return axial_raw(n, m, eta)

        @lru_cache(maxsize=None)
        def get_phi(n: int, m: int) -> MixedScaleCoefficient:
            n = _index(n, "n")
            m = _index(m, "m")
            phi0 = self.rational_data.angular_reference_jet_decimal(n, m, eta)
            reference = MixedScaleCoefficient.channel(
                0,
                0,
                _finite_decimal(phi0, "reference phi jet"),
            )
            return reference.add(get_angular_r(n, m).shift_lambda(1).scale(half))

        @lru_cache(maxsize=None)
        def get_u(n: int, m: int) -> MixedScaleCoefficient:
            n = _index(n, "n")
            m = _index(m, "m")
            u0 = _decimal_from_float(
                self.x1.reference.u.jet(n, m, eta),
                "reference u jet",
            )
            reference = MixedScaleCoefficient.channel(
                0,
                0,
                _finite_decimal(u0, "reference u jet"),
            )
            return reference.add(get_axial_r(n, m).shift_lambda(1).scale(half))

        # This is the scaled pressure primitive
        # P = primitive((a * phi)^2) = primitive(a^2 * phi^2).
        # Keep it separate from the axial remainder forcing above: with an
        # external input pair, this fifth family is P(input), without J1[d
        # P_eta - 4 A eta P - 2 eta Y P_Y] and without an outer Lambda factor.
        # With the default self-recursive input, it follows that same local
        # recursive input family rather than exposing an updated iterate's P.
        def output_phi_family(n: int, m: int, z: float) -> MixedScaleCoefficient:
            _local_eta(z)
            return get_phi(_index(n, "n"), _index(m, "m"))

        def output_u_family(n: int, m: int, z: float) -> MixedScaleCoefficient:
            _local_eta(z)
            return get_u(_index(n, "n"), _index(m, "m"))

        return output_phi_family, output_u_family, get_angular_r, get_axial_r, primitive_source

    def picard_map_families(
        self,
        input_pair: tuple[MixedScaleFamily, MixedScaleFamily],
        eta: float,
    ) -> tuple[MixedScaleFamily, MixedScaleFamily]:
        """Return output families for one finite Picard map application.

        The supplied pair is used only as the input to the pinned remainder;
        the returned pair is the updated ``reference + R(input)/(2*Lambda)``
        output.  All row and eta caches belong to this local graph and are
        released when the returned families become unreachable.
        """

        if not isinstance(input_pair, tuple) or len(input_pair) != 2:
            raise TypeError("input_pair must be a (phi_family, u_family) tuple")
        if not callable(input_pair[0]) or not callable(input_pair[1]):
            raise TypeError("input_pair entries must be callable mixed-scale families")
        eta = _eta_in_window(eta)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            output_phi, output_u, _, _, _ = self._build_families(eta, input_pair)
        return output_phi, output_u

    def picard_map_jet_prefix(
        self,
        input_pair: tuple[MixedScaleFamily, MixedScaleFamily],
        max_n: int,
        m: int,
        eta: float,
    ) -> tuple[tuple[MixedScaleCoefficient, MixedScaleCoefficient], ...]:
        """Evaluate one Picard map application on rows ``0..max_n``."""

        max_n = _index(max_n, "max_n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            output_phi, output_u = self.picard_map_families(input_pair, eta)
            return tuple(
                (output_phi(n, m, eta), output_u(n, m, eta))
                for n in range(max_n + 1)
            )

    def apply_map_prefix(
        self,
        input_pair: tuple[MixedScaleFamily, MixedScaleFamily],
        max_n: int,
        m: int,
        eta: float,
    ) -> tuple[tuple[MixedScaleCoefficient, MixedScaleCoefficient], ...]:
        """Compatibility name for :meth:`picard_map_jet_prefix`."""

        return self.picard_map_jet_prefix(input_pair, max_n, m, eta)

    def jet_pair(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[MixedScaleCoefficient, MixedScaleCoefficient]:
        """Return formal angular and axial coefficient channels at one row."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            phi_family, u_family, _, _, _ = self._build_families(eta)
            return phi_family(n, m, eta), u_family(n, m, eta)

    def jet_prefix(
        self,
        max_n: int,
        m: int,
        eta: float,
    ) -> tuple[tuple[MixedScaleCoefficient, MixedScaleCoefficient], ...]:
        """Return rows ``0..max_n`` from one eta-local recursive family graph."""

        max_n = _index(max_n, "max_n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            phi_family, u_family, _, _, _ = self._build_families(eta)
            return tuple(
                (phi_family(n, m, eta), u_family(n, m, eta))
                for n in range(max_n + 1)
            )

    def profile_jet_prefix(
        self,
        max_n: int,
        m: int,
        eta: float,
    ) -> tuple[
        tuple[MixedScaleCoefficient, MixedScaleCoefficient, MixedScaleCoefficient],
        ...,
    ]:
        """Return finite rows of ``(phi, u, P)`` from one local family graph.

        ``P`` is the scaled pressure primitive
        ``primitive(a^2 * phi^2)``.  It is kept independent of the axial
        remainder pressure forcing and of any outer inverse-Lambda factor.
        """

        max_n = _index(max_n, "max_n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            phi_family, u_family, _, _, pressure_family = self._build_families(eta)
            return tuple(
                (
                    phi_family(n, m, eta),
                    u_family(n, m, eta),
                    pressure_family(n, m, eta),
                )
                for n in range(max_n + 1)
            )


def formal_axis_coefficient_solver(
    x1: ActualScheduleWideFirstPicardState,
) -> FormalAxisCoefficientSolverState:
    """Bind the formal triangular solver to one genuine first-Picard anchor."""

    return FormalAxisCoefficientSolverState(x1=x1)


__all__ = [
    "FormalAxisCoefficientSolverState",
    "formal_axis_coefficient_solver",
]
