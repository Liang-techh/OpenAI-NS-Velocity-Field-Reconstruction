import math

import numpy as np
import pytest

from openai_ns_reconstruction.cylindrical_curl_jet import (
    TangentCylindricalCurlJet,
    cylindrical_coefficient_curl,
)


def _cyl_field(r: float, theta: float, z: float) -> np.ndarray:
    """Explicit smooth complex coefficient used only as an independent oracle."""
    return np.array(
        [
            r * r * z + 1j * math.sin(theta),
            r * z + 1j * math.cos(theta),
            z * z + r * math.cos(theta) + 1j * r * math.sin(theta),
        ],
        dtype=complex,
    )


def _cyl_to_cart(v: np.ndarray, theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array(
        [v[0] * c - v[1] * s, v[0] * s + v[1] * c, v[2]],
        dtype=complex,
    )


def _cart_field(x: float, y: float, z: float) -> np.ndarray:
    r = math.hypot(x, y)
    theta = math.atan2(y, x)
    return _cyl_to_cart(_cyl_field(r, theta, z), theta)


def test_cylindrical_curl_against_independent_cartesian_finite_difference():
    r, theta, z = 1.7, 0.4, -0.3
    coefficient = _cyl_field(r, theta, z)

    d_r = np.array(
        [2.0 * r * z, z, math.cos(theta) + 1j * math.sin(theta)], dtype=complex
    )
    d_theta = np.array(
        [
            1j * math.cos(theta),
            -1j * math.sin(theta),
            -r * math.sin(theta) + 1j * r * math.cos(theta),
        ],
        dtype=complex,
    )
    d_z = np.array([r * r, r, 2.0 * z], dtype=complex)

    cylindrical = cylindrical_coefficient_curl(
        radius=r,
        coefficient=coefficient,
        radial_derivative=d_r,
        angular_derivative=d_theta,
        axial_derivative=d_z,
    )

    # Independent oracle: convert the explicit coefficient field to Cartesian
    # components first, then differentiate those Cartesian components directly.
    # This path never uses the production cylindrical curl formula.
    x = r * math.cos(theta)
    y = r * math.sin(theta)
    step = 1e-6

    def derivative(axis: int) -> np.ndarray:
        plus = [x, y, z]
        minus = [x, y, z]
        plus[axis] += step
        minus[axis] -= step
        return (_cart_field(*plus) - _cart_field(*minus)) / (2.0 * step)

    d_x = derivative(0)
    d_y = derivative(1)
    d_z_cart = derivative(2)
    cartesian_curl = np.array(
        [
            d_y[2] - d_z_cart[1],
            d_z_cart[0] - d_x[2],
            d_x[1] - d_y[0],
        ],
        dtype=complex,
    )

    assert np.allclose(
        _cyl_to_cart(cylindrical, theta),
        cartesian_curl,
        rtol=2e-9,
        atol=2e-9,
    )


def test_tangent_jet_uses_pinned_coefficient_and_structured_curl():
    # Both real and imaginary amplitude parts are exactly tangent to n.
    normal = (1.0, 2.0, 2.0)
    amplitude = (2.0 + 0.0j, -1.0 + 2.0j, 0.0 - 2.0j)
    cert = TangentCylindricalCurlJet(
        frequency=8.0,
        normal=normal,
        amplitude=amplitude,
        radius=2.0,
        radial_derivative=(1.0 + 1.0j, 2.0 - 1.0j, 3.0 + 0.0j),
        angular_derivative=(4.0 + 0.0j, 5.0 + 2.0j, 6.0 - 1.0j),
        axial_derivative=(7.0 - 1.0j, 8.0 + 0.0j, 9.0 + 2.0j),
    )

    b = cert.coefficient
    d_r = np.asarray(cert.radial_derivative, dtype=complex)
    d_theta = np.asarray(cert.angular_derivative, dtype=complex)
    d_z = np.asarray(cert.axial_derivative, dtype=complex)
    independent_curl = np.array(
        [
            d_theta[2] / 2.0 - d_z[1],
            d_z[0] - d_r[2],
            d_r[1] + b[1] / 2.0 - d_theta[0] / 2.0,
        ],
        dtype=complex,
    )

    assert np.allclose(cert.coefficient_curl, independent_curl, rtol=0.0, atol=0.0)
    assert np.allclose(
        cert.realized_coefficient,
        np.asarray(amplitude, dtype=complex) + (1j / 8.0) * independent_curl,
        rtol=0.0,
        atol=0.0,
    )


def test_cylindrical_jet_fails_closed_on_invalid_geometry_or_tangency():
    with pytest.raises(ValueError, match="radius>0"):
        cylindrical_coefficient_curl(
            radius=0.0,
            coefficient=(0.0, 0.0, 0.0),
            radial_derivative=(0.0, 0.0, 0.0),
            angular_derivative=(0.0, 0.0, 0.0),
            axial_derivative=(0.0, 0.0, 0.0),
        )

    with pytest.raises(ValueError, match="tangency"):
        TangentCylindricalCurlJet(
            frequency=8.0,
            normal=(1.0, 0.0, 0.0),
            amplitude=(1.0 + 0.0j, 0.0, 0.0),
            radius=1.0,
            radial_derivative=(0.0, 0.0, 0.0),
            angular_derivative=(0.0, 0.0, 0.0),
            axial_derivative=(0.0, 0.0, 0.0),
        )
