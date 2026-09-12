import math
from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.local_field import LocalField
from openai_ns_reconstruction.section9_eq921_prefix_jet import (
    Section9FinitePrefixJetCertificate,
    required_spacetime_multiindices,
)
from openai_ns_reconstruction.section9_section10_localized_velocity_spatial_jacobian import (
    localize_section9_prefix_spatial_jacobian,
)
from openai_ns_reconstruction.spatial_localization import section10_localized_field
from openai_ns_reconstruction.verify import jacobian_numeric


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


def _prefix_at(
    x: float,
    y: float,
    z: float,
    t: float,
    *,
    order: int = 2,
    **overrides,
) -> Section9FinitePrefixJetCertificate:
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
        provider_id="manufactured-spatial-prefix",
        provider_revision="test-r1",
        provider_provenance="analytic polynomial/radial fixture",
    )
    kwargs.update(overrides)
    return Section9FinitePrefixJetCertificate(**kwargs)


def test_localized_spatial_jacobian_matches_independent_nested_finite_difference():
    # 16 r^2 and 4 z are both strictly in the executable cutoff transition
    # collar, so the comparison exercises gradient, Hessian, and mixed terms.
    x, y, z, t = 0.20, 0.02, 0.16, 0.6
    certificate = localize_section9_prefix_spatial_jacobian(
        _prefix_at(x, y, z, t), x=x, y=y, z=z, t=t
    )

    # Independent hand derivatives of curl(A) and of B e_theta.  Here
    # curl(A)=(x-0.1t, -(y+0.3+0.4t), 0) and
    # B e_theta=k(-y,x,0), k=0.25+0.2t.
    np.testing.assert_allclose(
        certificate.curl_A_spatial_jacobian,
        np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 0.0],
            ]
        ),
        rtol=0.0,
        atol=1.0e-15,
    )
    k = 0.25 + 0.2 * t
    np.testing.assert_allclose(
        certificate.swirl_spatial_jacobian,
        np.array(
            [
                [0.0, -k, 0.0],
                [k, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ]
        ),
        rtol=2.0e-14,
        atol=2.0e-14,
    )

    # Independent path: finite-difference the complete LocalizedField velocity.
    # The nested velocity evaluator has no analytic curl, so its poloidal piece
    # finite-differences curl(cA) as a whole instead of reusing the production
    # spatial-Jacobian/product-rule implementation under test.
    field = section10_localized_field(LocalField(_A, _B))

    def independent_velocity(x0: float, y0: float, z0: float, t0: float) -> np.ndarray:
        return field.velocity(x0, y0, z0, t0, eps=2.0e-6)

    independent = jacobian_numeric(
        independent_velocity, x, y, z, t, eps=4.0e-5
    )
    np.testing.assert_allclose(
        certificate.localized_velocity_spatial_jacobian,
        independent,
        rtol=2.5e-3,
        atol=3.0e-3,
    )

    assert certificate.formal_localized_spatial_jacobian_ready
    assert not certificate.axis_spatial_derivative_covered
    assert not certificate.residual_artifact_ready
    assert not certificate.paper_exact_velocity_available


def test_localized_spatial_jacobian_is_exactly_zero_outside_fixed_support():
    x, y, z, t = 0.30, 0.02, 0.01, 0.6
    certificate = localize_section9_prefix_spatial_jacobian(
        _prefix_at(x, y, z, t), x=x, y=y, z=z, t=t
    )
    assert certificate.base.cutoff_value == 0.0
    np.testing.assert_array_equal(certificate.base.cutoff_gradient, np.zeros(3))
    np.testing.assert_array_equal(certificate.cutoff_hessian, np.zeros((3, 3)))
    np.testing.assert_array_equal(
        certificate.localized_velocity_spatial_jacobian,
        np.zeros((3, 3)),
    )


def test_spatial_jacobian_requires_second_order_prefix_and_fail_closed_truth():
    x, y, z, t = 0.20, 0.02, 0.16, 0.6
    with pytest.raises(ValueError, match="at least two"):
        localize_section9_prefix_spatial_jacobian(
            _prefix_at(x, y, z, t, order=1), x=x, y=y, z=z, t=t
        )

    untrusted = _prefix_at(
        x, y, z, t, paper_exact_velocity_available=True
    )
    with pytest.raises(ValueError, match="truth gate"):
        localize_section9_prefix_spatial_jacobian(
            untrusted, x=x, y=y, z=z, t=t
        )
