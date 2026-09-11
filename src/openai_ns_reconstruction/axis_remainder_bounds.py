"""Certified upper-bound propagation for the natural-axis fixed-point remainder.

This module mirrors the scalar ``bound``/``lip`` bookkeeping in the pinned
OpenAI Lean definitions
``AxisContraction.Controlled`` and ``AxisContraction.controlledRemainder``.
It does **not** infer coefficient norms from samples and it does not manufacture
fixed-point data.  Instead, it turns independently certified upper bounds for
operator norms and the fixed ``AxisData`` coefficients into conservative
upper bounds for ``remainderBound`` and ``remainderLip``.

The distinction matters: Lean defines ``remainderBound``/``remainderLip``
using the actual Banach-space norms.  Replacing those norms by certified upper
bounds is monotone and therefore sufficient for choosing a safe ``Lambda``, but
it is not claimed to be definitionally equal to the opaque Lean quantities.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


def _finite_nonnegative(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


def _finite_positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


@dataclass(frozen=True)
class ControlledEstimate:
    """Scalar image of ``AxisContraction.Controlled``'s two numeric fields."""

    bound: float
    lip: float

    def __post_init__(self) -> None:
        _finite_nonnegative(self.bound, "bound")
        _finite_nonnegative(self.lip, "lip")

    @classmethod
    def const(cls, bound: float) -> "ControlledEstimate":
        return cls(_finite_nonnegative(bound, "constant bound"), 0.0)

    @classmethod
    def coordinate(cls, radius: float) -> "ControlledEstimate":
        return cls(_finite_nonnegative(radius, "radius"), 1.0)

    def add(self, other: "ControlledEstimate") -> "ControlledEstimate":
        return ControlledEstimate(self.bound + other.bound, self.lip + other.lip)

    def sub(self, other: "ControlledEstimate") -> "ControlledEstimate":
        # ``Controlled.sub`` is ``add`` with ``neg``; negation preserves both fields.
        return self.add(other)

    def linear(self, operator_norm: float) -> "ControlledEstimate":
        norm = _finite_nonnegative(operator_norm, "operator norm")
        return ControlledEstimate(norm * self.bound, norm * self.lip)

    def bilinear(
        self,
        other: "ControlledEstimate",
        operator_norm: float,
    ) -> "ControlledEstimate":
        norm = _finite_nonnegative(operator_norm, "bilinear operator norm")
        return ControlledEstimate(
            norm * self.bound * other.bound,
            norm * (self.lip * other.bound + self.bound * other.lip),
        )

    def unit_smul(self) -> "ControlledEstimate":
        """Mirror ``Controlled.unitSmul`` for a scalar with ``|c| <= 1``."""

        return self


@dataclass(frozen=True)
class NaturalOperatorNormBounds:
    """Certified upper bounds for the operators in ``NaturalOperators``.

    ``pinned_axis_operators`` fills every coefficient-space operator using the
    explicit inequalities proved in ``AxisOperators.lean``.  The angular
    resolvent is intentionally supplied separately because the pinned theorem
    bounds it by a nontrivial factorial series depending on ``||chi||``.
    """

    product: float
    average: float
    primitive: float
    parameter_primitive: float
    mul_y: float
    j1: float
    j2: float
    param1: float
    param2: float
    dot1: float
    dot2: float
    mixed1: float
    mixed2: float
    resolvent: float

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            _finite_nonnegative(value, name)

    @classmethod
    def pinned_axis_operators(
        cls,
        *,
        epsilon: float,
        resolvent_norm_upper: float,
    ) -> "NaturalOperatorNormBounds":
        """Use the explicit operator-norm upper bounds from pinned Lean.

        The constants are the proved bounds in ``AxisOperators.lean``:
        product ``64``; average ``1``; primitive, regular inverses and ``mulY``
        ``80``; parameter primitive ``80/epsilon``; inverse parameter/mixed
        bilinear maps ``5120/epsilon``; inverse-dot maps ``5120``.
        """

        eps = _finite_positive(epsilon, "epsilon")
        resolvent = _finite_nonnegative(resolvent_norm_upper, "resolvent_norm_upper")
        return cls(
            product=64.0,
            average=1.0,
            primitive=80.0,
            parameter_primitive=80.0 / eps,
            mul_y=80.0,
            j1=80.0,
            j2=80.0,
            param1=5120.0 / eps,
            param2=5120.0 / eps,
            dot1=5120.0,
            dot2=5120.0,
            mixed1=5120.0 / eps,
            mixed2=5120.0 / eps,
            resolvent=resolvent,
        )


@dataclass(frozen=True)
class AxisDataNormBounds:
    """Certified norm upper bounds for the fixed ``AxisData`` coefficient fields."""

    A: float
    D: float
    h: float
    one: float
    eta: float
    d: float
    inverse_l: float
    u_star: float
    u_star_eta: float
    w_star: float
    h_star: float
    normalized_gradient: float
    z_star: float

    def __post_init__(self) -> None:
        for name in (
            "one",
            "eta",
            "d",
            "inverse_l",
            "u_star",
            "u_star_eta",
            "w_star",
            "h_star",
            "normalized_gradient",
            "z_star",
        ):
            _finite_nonnegative(getattr(self, name), name)
        for name in ("A", "D", "h"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")


@dataclass(frozen=True)
class DerivedAxisNormBounds:
    angular_linear: float
    angular_quadratic: float
    average_coefficient: float
    angular_slow: float
    axial_linear: float
    axial_quadratic: float
    four_A_eta: float
    two_eta: float

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            _finite_nonnegative(value, name)


def derive_axis_coefficient_norm_bounds(
    operators: NaturalOperatorNormBounds,
    data: AxisDataNormBounds,
) -> DerivedAxisNormBounds:
    """Propagate primitive norm certificates through the six Lean coefficients."""

    p = operators.product
    abs_A = abs(float(data.A))
    abs_D = abs(float(data.D))
    abs_h = abs(float(data.h))

    angular_linear = (
        data.w_star
        + abs_h * data.one
        + 2.0 * abs_h * p * data.eta * data.u_star
    )
    angular_quadratic = p * data.d * data.normalized_gradient
    average_coefficient = 2.0 * abs_D * data.eta
    angular_slow = 2.0 * abs_h * data.eta
    axial_linear = (
        abs_A * data.one
        + 4.0 * abs_A * p * data.eta * data.u_star
        + p * data.d * data.u_star_eta
    )
    axial_quadratic = 2.0 * abs_A * data.eta
    return DerivedAxisNormBounds(
        angular_linear=angular_linear,
        angular_quadratic=angular_quadratic,
        average_coefficient=average_coefficient,
        angular_slow=angular_slow,
        axial_linear=axial_linear,
        axial_quadratic=axial_quadratic,
        four_A_eta=4.0 * abs_A * data.eta,
        two_eta=2.0 * data.eta,
    )


@dataclass(frozen=True)
class NaturalRemainderBoundCertificate:
    """Conservative inputs for the theorem's ``remainderBound/remainderLip``.

    ``reference_norm_upper`` bounds the product-space reference pair norm,
    ``radius_upper = reference_norm_upper + 1`` is therefore an admissible
    replacement for the Lean radius ``||referencePair|| + 1``.  The returned
    ``remainder_bound_upper`` and ``remainder_lipschitz_upper`` may be fed into
    the already-landed natural-scale selector as safe upper bounds once every
    primitive norm input has an independent certificate.
    """

    reference_norm_upper: float
    radius_upper: float
    remainder_bound_upper: float
    remainder_lipschitz_upper: float
    derived: DerivedAxisNormBounds

    def __post_init__(self) -> None:
        _finite_nonnegative(self.reference_norm_upper, "reference_norm_upper")
        _finite_positive(self.radius_upper, "radius_upper")
        _finite_nonnegative(self.remainder_bound_upper, "remainder_bound_upper")
        _finite_nonnegative(self.remainder_lipschitz_upper, "remainder_lipschitz_upper")


def reference_pair_norm_upper(
    operators: NaturalOperatorNormBounds,
    data: AxisDataNormBounds,
) -> float:
    """Bound ``||referencePair||`` using the product-space max norm.

    The pair is ``(S one, -1/2 * j1(product inverseL zStar))``.  Mathlib's
    product norm is the max norm, so the maximum of the two component bounds
    is a valid certificate.
    """

    phi0 = operators.resolvent * data.one
    u0 = 0.5 * operators.j1 * operators.product * data.inverse_l * data.z_star
    value = max(phi0, u0)
    if not math.isfinite(value):
        raise ArithmeticError("reference-pair norm bound overflowed")
    return value


def natural_remainder_bound_certificate(
    *,
    operators: NaturalOperatorNormBounds,
    data: AxisDataNormBounds,
    amplitude_norm_upper: float,
) -> NaturalRemainderBoundCertificate:
    """Mirror ``controlledRemainder`` using only certified scalar norm bounds.

    ``amplitude_norm_upper`` is the theorem's uniform ``M`` for normalized
    angular input ``a``.  No sampled field values enter this computation.
    """

    M = _finite_nonnegative(amplitude_norm_upper, "amplitude_norm_upper")
    derived = derive_axis_coefficient_norm_bounds(operators, data)
    reference_upper = reference_pair_norm_upper(operators, data)
    radius = reference_upper + 1.0
    if not math.isfinite(radius):
        raise ArithmeticError("reference radius bound overflowed")

    phi = ControlledEstimate.coordinate(radius)
    u = ControlledEstimate.coordinate(radius)
    bu = u.linear(operators.average)

    def c(value: float) -> ControlledEstimate:
        return ControlledEstimate.const(value)

    def mul(left: ControlledEstimate, right: ControlledEstimate) -> ControlledEstimate:
        return left.bilinear(right, operators.product)

    lin1 = (
        mul(c(derived.angular_linear), phi).linear(operators.j2)
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
        .add(mul(c(derived.average_coefficient), bu).bilinear(phi, operators.dot2))
        .add(mul(c(data.d), bu.bilinear(phi, operators.mixed2)))
        .sub(phi.bilinear(mul(c(data.d), u), operators.param2))
    )

    lin2 = (
        mul(c(derived.axial_linear), u).linear(operators.j1)
        .add(c(data.w_star).bilinear(u, operators.dot1))
        .add(u.bilinear(c(data.h_star), operators.param1))
    )
    slow2 = (
        mul(c(derived.axial_quadratic), mul(u, u)).linear(operators.j1)
        .add(mul(c(derived.average_coefficient), bu).bilinear(u, operators.dot1))
        .add(mul(c(data.d), bu.bilinear(u, operators.mixed1)))
        .sub(u.bilinear(mul(c(data.d), u), operators.param1))
    )

    ac = ControlledEstimate(M, 0.0)
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

    # ``Controlled.pair`` uses the sum of component bounds/Lipschitz constants.
    total = first.add(second)
    if not math.isfinite(total.bound) or not math.isfinite(total.lip):
        raise ArithmeticError("natural-remainder bound propagation overflowed")

    return NaturalRemainderBoundCertificate(
        reference_norm_upper=reference_upper,
        radius_upper=radius,
        remainder_bound_upper=total.bound,
        remainder_lipschitz_upper=total.lip,
        derived=derived,
    )
