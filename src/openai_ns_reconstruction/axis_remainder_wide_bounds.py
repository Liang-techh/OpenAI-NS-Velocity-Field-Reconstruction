"""Wide-decimal propagation for the natural-axis fixed-point remainder.

The binary64 implementation in :mod:`axis_remainder_bounds` mirrors the pinned
``AxisContraction.Controlled`` bookkeeping, but the actual SchedulePressure
ledger can overflow while multiplying already-certified upper bounds.  This
module evaluates the *same positive scalar majorant algebra* in a wider Decimal
representation so that a representation overflow is not mistaken for a
mathematical obstruction.

Inputs are still the independently supplied/certified ``NaturalOperatorNormBounds``
and ``AxisDataNormBounds`` objects.  No field values are sampled and no
fixed-point data are manufactured.  Binary floats are converted exactly with
``Decimal.from_float`` and every positive addition/multiplication performed in
this module is rounded toward ``+infinity`` at fixed high precision.  This
preserves the direction of the existing conservative upper-bound ledger at this
arithmetic layer; it is not an interval/Lean proof of the upstream float inputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, ROUND_CEILING, localcontext
import math
from typing import Union

from .axis_remainder_bounds import AxisDataNormBounds, NaturalOperatorNormBounds

DecimalLike = Union[Decimal, float, int]
_DECIMAL_PRECISION = 96


def _decimal_finite(value: DecimalLike, name: str) -> Decimal:
    if isinstance(value, Decimal):
        result = value
    elif isinstance(value, int):
        result = Decimal(value)
    else:
        f = float(value)
        if not math.isfinite(f):
            raise ValueError(f"{name} must be finite")
        result = Decimal.from_float(f)
    if not result.is_finite():
        raise ValueError(f"{name} must be finite")
    return result


def _decimal_nonnegative(value: DecimalLike, name: str) -> Decimal:
    result = _decimal_finite(value, name)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _add_up(*values: DecimalLike) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.rounding = ROUND_CEILING
        total = Decimal(0)
        for value in values:
            total += _decimal_nonnegative(value, "summand")
        return +total


def _mul_up(*values: DecimalLike) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.rounding = ROUND_CEILING
        product = Decimal(1)
        for value in values:
            product *= _decimal_nonnegative(value, "factor")
        return +product


@dataclass(frozen=True)
class WideControlledEstimate:
    """Decimal image of the nonnegative ``Controlled.bound/lip`` algebra."""

    bound: Decimal
    lip: Decimal

    def __post_init__(self) -> None:
        _decimal_nonnegative(self.bound, "bound")
        _decimal_nonnegative(self.lip, "lip")

    @classmethod
    def const(cls, bound: DecimalLike) -> "WideControlledEstimate":
        return cls(_decimal_nonnegative(bound, "constant bound"), Decimal(0))

    @classmethod
    def coordinate(cls, radius: DecimalLike) -> "WideControlledEstimate":
        return cls(_decimal_nonnegative(radius, "radius"), Decimal(1))

    def add(self, other: "WideControlledEstimate") -> "WideControlledEstimate":
        return WideControlledEstimate(
            _add_up(self.bound, other.bound),
            _add_up(self.lip, other.lip),
        )

    def sub(self, other: "WideControlledEstimate") -> "WideControlledEstimate":
        # Norm bookkeeping for ``Controlled.sub`` is addition after negation.
        return self.add(other)

    def linear(self, operator_norm: DecimalLike) -> "WideControlledEstimate":
        norm = _decimal_nonnegative(operator_norm, "operator norm")
        return WideControlledEstimate(
            _mul_up(norm, self.bound),
            _mul_up(norm, self.lip),
        )

    def bilinear(
        self,
        other: "WideControlledEstimate",
        operator_norm: DecimalLike,
    ) -> "WideControlledEstimate":
        norm = _decimal_nonnegative(operator_norm, "bilinear operator norm")
        lip_sum = _add_up(
            _mul_up(self.lip, other.bound),
            _mul_up(self.bound, other.lip),
        )
        return WideControlledEstimate(
            _mul_up(norm, self.bound, other.bound),
            _mul_up(norm, lip_sum),
        )

    def unit_smul(self) -> "WideControlledEstimate":
        return self


@dataclass(frozen=True)
class WideDerivedAxisNormBounds:
    angular_linear: Decimal
    angular_quadratic: Decimal
    average_coefficient: Decimal
    angular_slow: Decimal
    axial_linear: Decimal
    axial_quadratic: Decimal
    four_A_eta: Decimal
    two_eta: Decimal

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            _decimal_nonnegative(value, name)


def derive_axis_coefficient_norm_bounds_wide(
    operators: NaturalOperatorNormBounds,
    data: AxisDataNormBounds,
) -> WideDerivedAxisNormBounds:
    """Re-evaluate the six pinned coefficient majorants in wide arithmetic."""

    p = _decimal_nonnegative(operators.product, "product")
    abs_A = abs(_decimal_finite(data.A, "A"))
    abs_D = abs(_decimal_finite(data.D, "D"))
    abs_h = abs(_decimal_finite(data.h, "h"))
    one = _decimal_nonnegative(data.one, "one")
    eta = _decimal_nonnegative(data.eta, "eta")
    d = _decimal_nonnegative(data.d, "d")
    u_star = _decimal_nonnegative(data.u_star, "u_star")
    u_star_eta = _decimal_nonnegative(data.u_star_eta, "u_star_eta")
    w_star = _decimal_nonnegative(data.w_star, "w_star")
    gradient = _decimal_nonnegative(data.normalized_gradient, "normalized_gradient")

    return WideDerivedAxisNormBounds(
        angular_linear=_add_up(
            w_star,
            _mul_up(abs_h, one),
            _mul_up(2, abs_h, p, eta, u_star),
        ),
        angular_quadratic=_mul_up(p, d, gradient),
        average_coefficient=_mul_up(2, abs_D, eta),
        angular_slow=_mul_up(2, abs_h, eta),
        axial_linear=_add_up(
            _mul_up(abs_A, one),
            _mul_up(4, abs_A, p, eta, u_star),
            _mul_up(p, d, u_star_eta),
        ),
        axial_quadratic=_mul_up(2, abs_A, eta),
        four_A_eta=_mul_up(4, abs_A, eta),
        two_eta=_mul_up(2, eta),
    )


@dataclass(frozen=True)
class NaturalRemainderWideCertificate:
    """Wide representation of the same conservative remainder ledger."""

    reference_norm_upper: Decimal
    radius_upper: Decimal
    remainder_bound_upper: Decimal
    remainder_lipschitz_upper: Decimal
    derived: WideDerivedAxisNormBounds

    def __post_init__(self) -> None:
        _decimal_nonnegative(self.reference_norm_upper, "reference_norm_upper")
        radius = _decimal_nonnegative(self.radius_upper, "radius_upper")
        if radius <= 0:
            raise ValueError("radius_upper must be positive")
        _decimal_nonnegative(self.remainder_bound_upper, "remainder_bound_upper")
        _decimal_nonnegative(
            self.remainder_lipschitz_upper,
            "remainder_lipschitz_upper",
        )

    @property
    def decimal_order_bound(self) -> int | None:
        if self.remainder_bound_upper == 0:
            return None
        return int(self.remainder_bound_upper.adjusted())

    @property
    def decimal_order_lipschitz(self) -> int | None:
        if self.remainder_lipschitz_upper == 0:
            return None
        return int(self.remainder_lipschitz_upper.adjusted())


def reference_pair_norm_upper_wide(
    operators: NaturalOperatorNormBounds,
    data: AxisDataNormBounds,
) -> Decimal:
    """Bound the product-space reference pair without binary64 products."""

    phi0 = _mul_up(operators.resolvent, data.one)
    u0 = _mul_up(
        Decimal("0.5"),
        operators.j1,
        operators.product,
        data.inverse_l,
        data.z_star,
    )
    return max(phi0, u0)


def natural_remainder_bound_certificate_wide(
    *,
    operators: NaturalOperatorNormBounds,
    data: AxisDataNormBounds,
    amplitude_norm_upper: float,
) -> NaturalRemainderWideCertificate:
    """Mirror ``controlledRemainder`` using upward-rounded Decimal arithmetic.

    This function changes only the scalar representation.  It uses the same
    operator/data norm inputs and the same triangle/product/Lipschitz rules as
    :func:`axis_remainder_bounds.natural_remainder_bound_certificate`.
    """

    M = _decimal_nonnegative(amplitude_norm_upper, "amplitude_norm_upper")
    derived = derive_axis_coefficient_norm_bounds_wide(operators, data)
    reference_upper = reference_pair_norm_upper_wide(operators, data)
    radius = _add_up(reference_upper, 1)

    phi = WideControlledEstimate.coordinate(radius)
    u = WideControlledEstimate.coordinate(radius)
    bu = u.linear(operators.average)

    def c(value: DecimalLike) -> WideControlledEstimate:
        return WideControlledEstimate.const(value)

    def mul(
        left: WideControlledEstimate,
        right: WideControlledEstimate,
    ) -> WideControlledEstimate:
        return left.bilinear(right, operators.product)

    lin1 = (
        mul(c(derived.angular_linear), phi)
        .linear(operators.j2)
        .add(c(data.w_star).bilinear(phi, operators.dot2))
        .add(phi.bilinear(c(data.h_star), operators.param2))
    )
    quad1 = mul(mul(c(derived.angular_quadratic), u), phi).linear(operators.j2)

    slow1 = (
        mul(
            mul(c(derived.average_coefficient), bu).add(
                mul(c(derived.angular_slow), u)
            ),
            phi,
        )
        .linear(operators.j2)
        .add(bu.bilinear(mul(c(data.d), phi), operators.param2))
        .add(
            mul(c(derived.average_coefficient), bu).bilinear(
                phi,
                operators.dot2,
            )
        )
        .add(mul(c(data.d), bu.bilinear(phi, operators.mixed2)))
        .sub(phi.bilinear(mul(c(data.d), u), operators.param2))
    )

    lin2 = (
        mul(c(derived.axial_linear), u)
        .linear(operators.j1)
        .add(c(data.w_star).bilinear(u, operators.dot1))
        .add(u.bilinear(c(data.h_star), operators.param1))
    )
    slow2 = (
        mul(c(derived.axial_quadratic), mul(u, u))
        .linear(operators.j1)
        .add(
            mul(c(derived.average_coefficient), bu).bilinear(
                u,
                operators.dot1,
            )
        )
        .add(mul(c(data.d), bu.bilinear(u, operators.mixed1)))
        .sub(u.bilinear(mul(c(data.d), u), operators.param1))
    )

    ac = WideControlledEstimate(M, Decimal(0))
    source = mul(mul(ac, ac), mul(phi, phi))
    pressure_inside = (
        mul(c(derived.four_A_eta), source.linear(operators.primitive))
        .add(mul(c(data.d), source.linear(operators.parameter_primitive)))
        .sub(mul(c(derived.two_eta), source.linear(operators.mul_y)))
    )
    pressure = pressure_inside.linear(operators.j1)

    first = mul(
        c(data.inverse_l),
        lin1.add(quad1).sub(slow1.unit_smul()),
    ).linear(operators.resolvent)
    second = mul(
        c(data.inverse_l),
        lin2.sub(slow2.unit_smul()).add(pressure),
    )

    total = first.add(second)
    return NaturalRemainderWideCertificate(
        reference_norm_upper=reference_upper,
        radius_upper=radius,
        remainder_bound_upper=total.bound,
        remainder_lipschitz_upper=total.lip,
        derived=derived,
    )
