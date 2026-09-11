import numpy as np
import pytest

from openai_ns_reconstruction.curl_realization_algebra import (
    TangentCurlRealization,
    curl_remainder_from_coefficient_curl,
    inverse_carrier,
    normal_coefficient,
    normal_cross,
    normal_dot,
    phase_factor,
    principal_curl_coefficient,
    tangential_projection,
)


def test_normal_coefficient_and_triple_product_against_numpy():
    n = np.array([1.0, 2.0, 2.0])
    a = np.array([3.0 + 2.0j, -4.0 + 1.0j, 5.0 - 3.0j])

    expected_cross = np.cross(n.astype(complex), a)
    expected_b = expected_cross / np.dot(n, n)
    assert np.allclose(normal_cross(n, a), expected_cross, rtol=0.0, atol=0.0)
    assert np.allclose(normal_coefficient(n, a), expected_b, rtol=0.0, atol=0.0)
    assert normal_dot(n, a) == np.dot(n.astype(complex), a)

    # Independent vector-triple-product evaluation of -n x B.
    expected_principal = a - n * (np.dot(n.astype(complex), a) / np.dot(n, n))
    assert np.allclose(principal_curl_coefficient(n, a), expected_principal, rtol=1e-15, atol=1e-15)
    assert np.allclose(tangential_projection(n, a), expected_principal, rtol=1e-15, atol=1e-15)


def test_exact_tangent_specialization_recovers_raw_amplitude():
    # Both real and imaginary parts are exactly orthogonal to n in binary arithmetic.
    n = (1.0, 2.0, 2.0)
    a = (2.0 + 0.0j, -1.0 + 2.0j, 0.0 - 2.0j)
    cert = TangentCurlRealization(frequency=8.0, normal=n, amplitude=a)

    assert normal_dot(n, a) == 0j
    assert np.allclose(cert.principal, np.asarray(a, dtype=complex), rtol=0.0, atol=0.0)
    assert cert.carrier_identity == -1.0 + 0.0j


def test_displayed_derivative_remainder_is_independently_reconstructed():
    n = (1.0, 2.0, 2.0)
    a = (2.0 + 0.0j, -1.0 + 2.0j, 0.0 - 2.0j)
    K = -4.0
    curl_b = np.array([1.5 - 2.0j, -3.0 + 0.5j, 2.25 + 1.0j])
    cert = TangentCurlRealization(frequency=K, normal=n, amplitude=a)

    expected_remainder = (1j / K) * curl_b
    expected_realized = np.asarray(a, dtype=complex) + expected_remainder
    assert np.allclose(
        curl_remainder_from_coefficient_curl(K, curl_b),
        expected_remainder,
        rtol=0.0,
        atol=0.0,
    )
    assert np.allclose(cert.realized_coefficient(curl_b), expected_realized, rtol=0.0, atol=0.0)
    assert inverse_carrier(K) * phase_factor(K) == -1.0 + 0.0j


def test_non_tangent_data_keep_projection_defect_explicit_and_fail_certificate():
    n = (0.0, 0.0, 2.0)
    a = (1.0 + 0.0j, 2.0 + 0.0j, 3.0 + 0.0j)
    projected = principal_curl_coefficient(n, a)
    assert np.allclose(projected, np.array([1.0, 2.0, 0.0], dtype=complex))
    assert not np.allclose(projected, np.asarray(a, dtype=complex))
    with pytest.raises(ValueError, match="exact tangency"):
        TangentCurlRealization(frequency=3.0, normal=n, amplitude=a)


def test_fail_closed_on_degenerate_or_nonfinite_inputs():
    with pytest.raises(ValueError, match="phase normal"):
        normal_coefficient((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    with pytest.raises(ValueError, match="frequency"):
        inverse_carrier(0.0)
    with pytest.raises(ValueError, match="frequency"):
        phase_factor(float("inf"))
    with pytest.raises(ValueError, match="finite complex"):
        normal_cross((1.0, 0.0, 0.0), (complex(float("nan"), 0.0), 0.0, 0.0))
