import math

import pytest

from openai_ns_reconstruction.axis_remainder_bounds import (
    AxisDataNormBounds,
    NaturalOperatorNormBounds,
    derive_axis_coefficient_norm_bounds,
    natural_remainder_bound_certificate,
    reference_pair_norm_upper,
)


def _axis_data(**updates: float) -> AxisDataNormBounds:
    values = dict(
        A=0.3,
        D=0.2,
        h=0.01,
        one=1.0,
        eta=0.5,
        d=0.25,
        inverse_l=0.75,
        u_star=0.4,
        u_star_eta=0.6,
        w_star=0.7,
        h_star=0.8,
        normalized_gradient=0.9,
        z_star=1.1,
    )
    values.update(updates)
    return AxisDataNormBounds(**values)


def test_pinned_operator_norm_ledger_matches_lean_constants() -> None:
    ops = NaturalOperatorNormBounds.pinned_axis_operators(
        epsilon=0.5,
        resolvent_norm_upper=7.0,
    )
    assert ops.product == 64.0
    assert ops.average == 1.0
    assert ops.primitive == 80.0
    assert ops.parameter_primitive == 160.0
    assert ops.mul_y == 80.0
    assert ops.j1 == 80.0
    assert ops.j2 == 80.0
    assert ops.param1 == 10240.0
    assert ops.param2 == 10240.0
    assert ops.dot1 == 5120.0
    assert ops.dot2 == 5120.0
    assert ops.mixed1 == 10240.0
    assert ops.mixed2 == 10240.0
    assert ops.resolvent == 7.0


def test_derived_axis_norms_use_triangle_and_operator_bounds() -> None:
    ops = NaturalOperatorNormBounds.pinned_axis_operators(
        epsilon=1.0,
        resolvent_norm_upper=2.0,
    )
    data = _axis_data()
    derived = derive_axis_coefficient_norm_bounds(ops, data)

    assert derived.angular_linear == pytest.approx(
        data.w_star
        + abs(data.h) * data.one
        + 2.0 * abs(data.h) * ops.product * data.eta * data.u_star
    )
    assert derived.angular_quadratic == pytest.approx(
        ops.product * data.d * data.normalized_gradient
    )
    assert derived.average_coefficient == pytest.approx(2.0 * abs(data.D) * data.eta)
    assert derived.angular_slow == pytest.approx(2.0 * abs(data.h) * data.eta)
    assert derived.axial_linear == pytest.approx(
        abs(data.A) * data.one
        + 4.0 * abs(data.A) * ops.product * data.eta * data.u_star
        + ops.product * data.d * data.u_star_eta
    )
    assert derived.axial_quadratic == pytest.approx(2.0 * abs(data.A) * data.eta)
    assert derived.four_A_eta == pytest.approx(4.0 * abs(data.A) * data.eta)
    assert derived.two_eta == pytest.approx(2.0 * data.eta)


def test_reference_pair_uses_product_space_max_norm_bound() -> None:
    ops = NaturalOperatorNormBounds.pinned_axis_operators(
        epsilon=1.0,
        resolvent_norm_upper=2.0,
    )
    data = _axis_data()
    phi0 = ops.resolvent * data.one
    u0 = 0.5 * ops.j1 * ops.product * data.inverse_l * data.z_star
    assert reference_pair_norm_upper(ops, data) == max(phi0, u0)


def test_controlled_remainder_wiring_has_an_exact_simple_case() -> None:
    # This is a unit test for the scalar Controlled bookkeeping only, not a
    # Navier--Stokes profile.  j1/j2 and all derivative composites are switched
    # off so only the two dot terms survive; final multiplication by inverseL
    # and the angular resolvent both have norm one.
    ops = NaturalOperatorNormBounds(
        product=1.0,
        average=0.0,
        primitive=0.0,
        parameter_primitive=0.0,
        mul_y=0.0,
        j1=0.0,
        j2=0.0,
        param1=0.0,
        param2=0.0,
        dot1=1.0,
        dot2=1.0,
        mixed1=0.0,
        mixed2=0.0,
        resolvent=1.0,
    )
    data = AxisDataNormBounds(
        A=0.0,
        D=0.0,
        h=0.0,
        one=0.0,
        eta=0.0,
        d=0.0,
        inverse_l=1.0,
        u_star=0.0,
        u_star_eta=0.0,
        w_star=1.0,
        h_star=0.0,
        normalized_gradient=0.0,
        z_star=0.0,
    )
    cert = natural_remainder_bound_certificate(
        operators=ops,
        data=data,
        amplitude_norm_upper=0.0,
    )

    assert cert.reference_norm_upper == 0.0
    assert cert.radius_upper == 1.0
    # phi and u each have Controlled(radius=1, lip=1); dot2/dot1 contribute
    # one copy apiece, and Controlled.pair sums the two component estimates.
    assert cert.remainder_bound_upper == 2.0
    assert cert.remainder_lipschitz_upper == 2.0


def test_amplitude_bound_only_increases_the_uniform_certificate() -> None:
    ops = NaturalOperatorNormBounds.pinned_axis_operators(
        epsilon=2.0,
        resolvent_norm_upper=1.5,
    )
    data = _axis_data()
    zero = natural_remainder_bound_certificate(
        operators=ops,
        data=data,
        amplitude_norm_upper=0.0,
    )
    positive = natural_remainder_bound_certificate(
        operators=ops,
        data=data,
        amplitude_norm_upper=0.2,
    )
    assert math.isfinite(positive.remainder_bound_upper)
    assert math.isfinite(positive.remainder_lipschitz_upper)
    assert positive.remainder_bound_upper >= zero.remainder_bound_upper
    assert positive.remainder_lipschitz_upper >= zero.remainder_lipschitz_upper


def test_fail_closed_for_uncertified_norm_inputs() -> None:
    with pytest.raises(ValueError, match="epsilon"):
        NaturalOperatorNormBounds.pinned_axis_operators(
            epsilon=0.0,
            resolvent_norm_upper=1.0,
        )
    with pytest.raises(ValueError, match="resolvent_norm_upper"):
        NaturalOperatorNormBounds.pinned_axis_operators(
            epsilon=1.0,
            resolvent_norm_upper=-1.0,
        )
    with pytest.raises(ValueError, match="z_star"):
        _axis_data(z_star=-0.1)
    with pytest.raises(ValueError, match="amplitude_norm_upper"):
        natural_remainder_bound_certificate(
            operators=NaturalOperatorNormBounds.pinned_axis_operators(
                epsilon=1.0,
                resolvent_norm_upper=1.0,
            ),
            data=_axis_data(),
            amplitude_norm_upper=float("nan"),
        )
