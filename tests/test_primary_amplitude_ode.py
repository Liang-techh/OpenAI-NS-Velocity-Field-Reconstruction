import numpy as np
import pytest

from openai_ns_reconstruction.primary_amplitude_ode import (
    PrimaryModalDatum,
    harmonic_damping,
    modal_derivative_to_physical,
    modal_errors,
    modal_forcing,
    modal_operator,
    modal_rhs,
    modal_to_physical,
    moving_coefficients,
    physical_forcing_components,
    physical_rhs,
    physical_to_modal,
)


def _datum():
    return PrimaryModalDatum(
        rho=0.3,
        rho_dot=-0.2,
        g_K=0.7,
        F=-0.8,
        N_theta=0.6,
        g_N=-0.5,
        rotation=0.15,
        eigenvalue=1.2,
        eigenvector=-0.9,
        eigen_rate=0.11,
        viscosity=0.04,
    )


def test_coefficients_and_modal_entries_match_pinned_formulas():
    data = _datum()
    denom = 1.0 + data.rho**2
    a = data.rho * (data.g_K - data.rho_dot) / denom
    b = (2.0 * data.F * data.N_theta - data.rho * data.rotation) / denom
    c = -(2.0 * data.F * data.N_theta + data.g_N) + data.rho * data.rotation
    coeff = moving_coefficients(data)
    assert (coeff.a, coeff.b, coeff.c) == pytest.approx((a, b, c), rel=0.0, abs=0.0)

    B = b - data.eigenvalue / data.eigenvector
    C = c - data.eigenvalue * data.eigenvector
    expected_errors = np.array(
        [
            (a + data.eigenvector * B + C / data.eigenvector - data.eigen_rate) / 2.0,
            (a - data.eigenvector * B + C / data.eigenvector + data.eigen_rate) / 2.0,
            (a + data.eigenvector * B - C / data.eigenvector + data.eigen_rate) / 2.0,
            (a - data.eigenvector * B - C / data.eigenvector - data.eigen_rate) / 2.0,
        ]
    )
    err = modal_errors(data)
    actual_errors = np.array([err.e11, err.e12, err.e21, err.e22])
    assert np.allclose(actual_errors, expected_errors, rtol=0.0, atol=0.0)

    d = 9.0 * data.viscosity
    expected_matrix = np.array(
        [
            [data.eigenvalue - d + expected_errors[0], expected_errors[1]],
            [expected_errors[2], -data.eigenvalue - d + expected_errors[3]],
        ]
    )
    assert harmonic_damping(data, 3) == d
    assert np.allclose(modal_operator(data, 3), expected_matrix, rtol=0.0, atol=0.0)


def test_modal_ode_reconstructs_independent_physical_rhs():
    data = _datum()
    state = np.array([0.7, -0.2])
    forcing = dict(f_radial=0.2, f_K=-0.4, f_N=0.1)

    # Modal path uses the PrimaryODE coefficient/forcing transform.
    z_dot = modal_rhs(data, 3, state, **forcing)
    reconstructed = modal_derivative_to_physical(data, state, z_dot)

    # Independent path evaluates MovingFrameODE.rhsX/rhsY directly in (x,y).
    x_y = modal_to_physical(data, state)
    direct = physical_rhs(data, 3, x_y, **forcing)
    assert np.allclose(reconstructed, direct, rtol=5e-15, atol=5e-15)

    # The basis conversion itself is invertible for the required h != 0.
    assert np.allclose(physical_to_modal(data, x_y), state, rtol=0.0, atol=1e-16)


def test_forcing_transform_is_checked_from_independent_physical_components():
    data = _datum()
    force_xy = physical_forcing_components(
        rho=data.rho, f_radial=0.35, f_K=-0.2, f_N=0.45
    )
    expected_modal = np.array(
        [
            (force_xy[0] + force_xy[1] / data.eigenvector) / 2.0,
            (force_xy[0] - force_xy[1] / data.eigenvector) / 2.0,
        ]
    )
    actual = modal_forcing(data, f_radial=0.35, f_K=-0.2, f_N=0.45)
    assert np.allclose(actual, expected_modal, rtol=0.0, atol=0.0)

    zero = modal_forcing(data, f_radial=0.0, f_K=0.0, f_N=0.0)
    assert np.array_equal(zero, np.zeros(2))


def test_fail_closed_on_invalid_basis_or_inputs():
    kwargs = dict(
        rho=0.0,
        rho_dot=0.0,
        g_K=0.0,
        F=0.0,
        N_theta=1.0,
        g_N=0.0,
        rotation=0.0,
        eigenvalue=1.0,
        eigenvector=1.0,
        eigen_rate=0.0,
        viscosity=0.0,
    )
    with pytest.raises(ValueError, match="eigenvector"):
        PrimaryModalDatum(**{**kwargs, "eigenvector": 0.0})
    with pytest.raises(ValueError, match="viscosity"):
        PrimaryModalDatum(**{**kwargs, "viscosity": float("nan")})

    data = PrimaryModalDatum(**kwargs)
    with pytest.raises(ValueError, match="harmonic"):
        modal_operator(data, True)
    with pytest.raises(ValueError, match="finite two-vector"):
        modal_rhs(data, 1, (1.0, float("inf")))
    with pytest.raises(ValueError, match="f_N"):
        modal_forcing(data, f_radial=0.0, f_K=0.0, f_N=float("nan"))
