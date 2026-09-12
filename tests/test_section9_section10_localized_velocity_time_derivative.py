import math
from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.local_field import LocalField
from openai_ns_reconstruction.section9_eq921_prefix_jet import (
    Section9FinitePrefixJetCertificate,
    required_spacetime_multiindices,
)
from openai_ns_reconstruction.section9_section10_localized_velocity_time_derivative import (
    localize_section9_prefix_time_derivative,
)
from openai_ns_reconstruction.spatial_localization import section10_localized_field


def _A(x: float, y: float, z: float, t: float) -> np.ndarray:
    f = (
        x * y
        + 0.3 * x
        + 0.2 * z
        + t * (0.4 * x - 0.1 * y + 0.05 * z)
        + 0.2 * t * t
    )
    return np.array([0.0, 0.0, f])


def _B(x: float, y: float, z: float, t: float) -> float:
    del z
    return math.hypot(x, y) * (0.25 + 0.2 * t)


def _prefix_at(x: float, y: float, z: float, t: float, *, order: int = 2, **overrides):
    required = required_spacetime_multiindices(order)
    A = {alpha: np.zeros(3) for alpha in required}
    B = {alpha: 0.0 for alpha in required}
    p = {alpha: 0.0 for alpha in required}

    A[(0, 0, 0, 0)] = _A(x, y, z, t)
    B[(0, 0, 0, 0)] = _B(x, y, z, t)
    p[(0, 0, 0, 0)] = 1.1 + 0.07 * t

    if order >= 1:
        A[(1, 0, 0, 0)] = np.array(
            [0.0, 0.0, 0.4 * x - 0.1 * y + 0.05 * z + 0.4 * t]
        )
        A[(0, 1, 0, 0)] = np.array([0.0, 0.0, y + 0.3 + 0.4 * t])
        A[(0, 0, 1, 0)] = np.array([0.0, 0.0, x - 0.1 * t])
        A[(0, 0, 0, 1)] = np.array([0.0, 0.0, 0.2 + 0.05 * t])

        r = math.hypot(x, y)
        k = 0.25 + 0.2 * t
        B[(1, 0, 0, 0)] = 0.2 * r
        B[(0, 1, 0, 0)] = k * x / r
        B[(0, 0, 1, 0)] = k * y / r
        p[(1, 0, 0, 0)] = 0.07

    if order >= 2:
        A[(2, 0, 0, 0)] = np.array([0.0, 0.0, 0.4])
        A[(1, 1, 0, 0)] = np.array([0.0, 0.0, 0.4])
        A[(1, 0, 1, 0)] = np.array([0.0, 0.0, -0.1])
        A[(1, 0, 0, 1)] = np.array([0.0, 0.0, 0.05])
        A[(0, 1, 1, 0)] = np.array([0.0, 0.0, 1.0])

        r = math.hypot(x, y)
        k = 0.25 + 0.2 * t
        r3 = r**3
        B[(1, 1, 0, 0)] = 0.2 * x / r
        B[(1, 0, 1, 0)] = 0.2 * y / r
        B[(0, 2, 0, 0)] = k * y * y / r3
        B[(0, 1, 1, 0)] = -k * x * y / r3
        B[(0, 0, 2, 0)] = k * x * x / r3

    kwargs = dict(
        q=Fraction(1, 2),
        prefix_order=1,
        stages=(1,),
        derivative_order=order,
        A_derivatives=A,
        B_derivatives=B,
        p_derivatives=p,
        provider_id="manufactured-time-prefix",
        provider_revision="test-r1",
        provider_provenance="analytic polynomial/radial fixture",
    )
    kwargs.update(overrides)
    return Section9FinitePrefixJetCertificate(**kwargs)


def test_localized_time_derivative_matches_independent_nested_finite_difference():
    x, y, z, t = 0.19, 0.06, 0.02, 0.6
    prefix = _prefix_at(x, y, z, t)
    certificate = localize_section9_prefix_time_derivative(
        prefix, x=x, y=y, z=z, t=t
    )

    # Hand differentiation of curl(A):
    # curl(A)=(x-0.1t, -(y+0.3+0.4t), 0).
    np.testing.assert_allclose(
        certificate.curl_A_time_derivative,
        np.array([-0.1, -0.4, 0.0]),
        rtol=0.0,
        atol=1.0e-15,
    )

    # Independent path: finite-difference the complete localized velocity in
    # time.  Its poloidal part finite-differences curl(cA) as a whole because
    # LocalField carries no analytic curl, so it does not reuse the production
    # time-derivative/product-rule implementation under test.
    field = section10_localized_field(LocalField(_A, _B))
    h = 4.0e-6
    plus = field.velocity(x, y, z, t + h, eps=2.0e-6)
    minus = field.velocity(x, y, z, t - h, eps=2.0e-6)
    numeric_time_derivative = (plus - minus) / (2.0 * h)
    np.testing.assert_allclose(
        certificate.localized_velocity_time_derivative,
        numeric_time_derivative,
        rtol=2.0e-4,
        atol=2.0e-5,
    )

    assert certificate.formal_localized_time_derivative_ready
    assert not certificate.residual_artifact_ready
    assert not certificate.paper_exact_velocity_available


def test_time_derivative_is_exactly_zero_outside_fixed_spatial_support():
    x, y, z, t = 0.30, 0.02, 0.01, 0.6
    certificate = localize_section9_prefix_time_derivative(
        _prefix_at(x, y, z, t), x=x, y=y, z=z, t=t
    )
    assert certificate.base.cutoff_value == 0.0
    np.testing.assert_array_equal(certificate.base.cutoff_gradient, np.zeros(3))
    np.testing.assert_array_equal(
        certificate.localized_velocity_time_derivative,
        np.zeros(3),
    )


def test_requires_second_order_prefix_jet_and_preserves_fail_closed_truth_gate():
    x, y, z, t = 0.19, 0.06, 0.02, 0.6
    with pytest.raises(ValueError, match="at least two"):
        localize_section9_prefix_time_derivative(
            _prefix_at(x, y, z, t, order=1), x=x, y=y, z=z, t=t
        )

    untrusted = _prefix_at(
        x, y, z, t, paper_exact_velocity_available=True
    )
    with pytest.raises(ValueError, match="truth gate"):
        localize_section9_prefix_time_derivative(
            untrusted, x=x, y=y, z=z, t=t
        )
