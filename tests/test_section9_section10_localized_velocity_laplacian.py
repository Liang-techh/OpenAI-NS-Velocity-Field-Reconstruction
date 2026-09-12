import math
from dataclasses import replace
from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.local_field import LocalField
from openai_ns_reconstruction.section9_eq921_prefix_jet import (
    Section9FinitePrefixJetCertificate,
    required_spacetime_multiindices,
)
from openai_ns_reconstruction.section9_section10_localized_velocity_laplacian import (
    localize_section9_prefix_laplacian,
)
from openai_ns_reconstruction.spatial_localization import section10_localized_field
from openai_ns_reconstruction.verify import laplacian_vector_numeric

_K = 0.4


def _A(x: float, y: float, z: float, t: float) -> np.ndarray:
    del t
    return np.array([0.0, 0.0, x**3 + y**3 + z**3])


def _B(x: float, y: float, z: float, t: float) -> float:
    del z, t
    return _K * math.hypot(x, y) ** 3


def _analytic_curl(x: float, y: float, z: float, t: float) -> np.ndarray:
    del z, t
    return np.array([3.0 * y * y, -3.0 * x * x, 0.0])


def _prefix_at(
    x: float,
    y: float,
    z: float,
    t: float,
    *,
    order: int = 3,
    **overrides,
) -> Section9FinitePrefixJetCertificate:
    del t
    required = required_spacetime_multiindices(order)
    A = {alpha: np.zeros(3) for alpha in required}
    B = {alpha: 0.0 for alpha in required}
    p = {alpha: 0.0 for alpha in required}

    A[(0, 0, 0, 0)] = _A(x, y, z, 0.0)
    B[(0, 0, 0, 0)] = _B(x, y, z, 0.0)

    if order >= 1:
        A[(0, 1, 0, 0)] = np.array([0.0, 0.0, 3.0 * x * x])
        A[(0, 0, 1, 0)] = np.array([0.0, 0.0, 3.0 * y * y])
        A[(0, 0, 0, 1)] = np.array([0.0, 0.0, 3.0 * z * z])

        r = math.hypot(x, y)
        B[(0, 1, 0, 0)] = 3.0 * _K * r * x
        B[(0, 0, 1, 0)] = 3.0 * _K * r * y

    if order >= 2:
        A[(0, 2, 0, 0)] = np.array([0.0, 0.0, 6.0 * x])
        A[(0, 0, 2, 0)] = np.array([0.0, 0.0, 6.0 * y])
        A[(0, 0, 0, 2)] = np.array([0.0, 0.0, 6.0 * z])

        r = math.hypot(x, y)
        B[(0, 2, 0, 0)] = 3.0 * _K * (r + x * x / r)
        B[(0, 1, 1, 0)] = 3.0 * _K * x * y / r
        B[(0, 0, 2, 0)] = 3.0 * _K * (r + y * y / r)

    if order >= 3:
        A[(0, 3, 0, 0)] = np.array([0.0, 0.0, 6.0])
        A[(0, 0, 3, 0)] = np.array([0.0, 0.0, 6.0])
        A[(0, 0, 0, 3)] = np.array([0.0, 0.0, 6.0])

        r = math.hypot(x, y)
        r3 = r**3
        B[(0, 3, 0, 0)] = _K * (9.0 * x / r - 3.0 * x**3 / r3)
        B[(0, 2, 1, 0)] = _K * (3.0 * y / r - 3.0 * x * x * y / r3)
        B[(0, 1, 2, 0)] = _K * (3.0 * x / r - 3.0 * x * y * y / r3)
        B[(0, 0, 3, 0)] = _K * (9.0 * y / r - 3.0 * y**3 / r3)

    kwargs = dict(
        q=Fraction(1, 2),
        prefix_order=1,
        stages=(1,),
        derivative_order=order,
        A_derivatives=A,
        B_derivatives=B,
        p_derivatives=p,
        provider_id="manufactured-laplacian-prefix",
        provider_revision="test-r1",
        provider_provenance="analytic cubic/radial fixture",
    )
    kwargs.update(overrides)
    return Section9FinitePrefixJetCertificate(**kwargs)


def test_localized_laplacian_matches_independent_final_velocity_difference():
    # Both 16 r^2 and 4 z lie strictly in the cutoff transition collar.  The
    # fixture also has nonzero Delta curl(A) and nonzero Delta(B e_theta), so
    # the comparison exercises the new third-order A path and swirl Laplacian.
    x, y, z, t = 0.20, 0.02, 0.16, 0.6
    certificate = localize_section9_prefix_laplacian(
        _prefix_at(x, y, z, t), x=x, y=y, z=z, t=t
    )

    np.testing.assert_allclose(
        certificate.curl_A_laplacian,
        np.array([6.0, -6.0, 0.0]),
        rtol=0.0,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        certificate.swirl_laplacian,
        8.0 * _K * np.array([-y, x, 0.0]),
        rtol=2.0e-13,
        atol=2.0e-13,
    )

    # Independent oracle: finite-difference the complete already-landed
    # LocalizedField velocity.  It receives only the hand-derived analytic
    # curl(A); it does not call the new Laplacian adapter or consume its third
    # derivatives/Hessian contractions.
    field = section10_localized_field(
        LocalField(_A, _B, analytic_curl=_analytic_curl)
    )
    independent = laplacian_vector_numeric(
        field.velocity, x, y, z, t, eps=2.0e-5
    )
    np.testing.assert_allclose(
        certificate.localized_velocity_laplacian,
        independent,
        rtol=1.0e-5,
        atol=3.0e-3,
    )

    assert certificate.formal_localized_laplacian_ready
    assert not certificate.axis_laplacian_covered
    assert not certificate.residual_artifact_ready
    assert not certificate.paper_exact_velocity_available


def test_localized_laplacian_is_exactly_zero_outside_fixed_support():
    x, y, z, t = 0.30, 0.02, 0.01, 0.6
    certificate = localize_section9_prefix_laplacian(
        _prefix_at(x, y, z, t), x=x, y=y, z=z, t=t
    )
    assert certificate.spatial.base.cutoff_value == 0.0
    assert certificate.cutoff_laplacian == 0.0
    np.testing.assert_array_equal(
        certificate.cutoff_laplacian_gradient, np.zeros(3)
    )
    np.testing.assert_array_equal(
        certificate.localized_velocity_laplacian, np.zeros(3)
    )


def test_localized_laplacian_requires_third_order_prefix_and_truth_gate():
    x, y, z, t = 0.20, 0.02, 0.16, 0.6
    with pytest.raises(ValueError, match="at least three"):
        localize_section9_prefix_laplacian(
            _prefix_at(x, y, z, t, order=2), x=x, y=y, z=z, t=t
        )

    untrusted = _prefix_at(
        x, y, z, t, paper_exact_velocity_available=True
    )
    with pytest.raises(ValueError, match="truth gate"):
        localize_section9_prefix_laplacian(
            untrusted, x=x, y=y, z=z, t=t
        )
