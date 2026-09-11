import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_lower_history_source import (
    ProfileSecondJet,
    positive_axis_point_data_from_lower_history,
    preceding_diffusion_from_second_jet,
)
from openai_ns_reconstruction.background_recurrence import omega_over_x_eq_5_6


def _poly(coeffs, X, eta):
    c0, c1, c2, c3, c4, c5 = coeffs
    return c0 + c1 * X + c2 * eta + c3 * X * eta + c4 * X**2 + c5 * eta**2


def _jet(coeffs, X, eta):
    c0, c1, c2, c3, c4, c5 = coeffs
    return ProfileSecondJet(
        value=_poly(coeffs, X, eta),
        radial=c1 + c3 * eta + 2 * c4 * X,
        radial2=2 * c4,
        parameter=c2 + c3 * X + 2 * c5 * eta,
        radial_parameter=c3,
        parameter2=2 * c5,
    )


def _direct_z(h, power, coeffs, X, eta):
    jet = _jet(coeffs, X, eta)
    d = 1 - eta**2
    ell = 1 - 2 * h * eta**2
    return (
        2 * eta * power * jet.value
        + d * jet.parameter
        - 2 * eta * X * jet.radial
    ) / ell


def _independent_preceding_diffusion(h, power, order, coeffs, X, eta):
    """Direct Z_(b-D) Z_b oracle using finite differences only for outer derivatives."""
    b = power + 2 * (order - 1) * h
    D = 0.5 - h
    eps = 2e-6

    def first(x, e):
        return _direct_z(h, b, coeffs, x, e)

    value = first(X, eta)
    dX = (first(X + eps, eta) - first(X - eps, eta)) / (2 * eps)
    dEta = (first(X, eta + eps) - first(X, eta - eps)) / (2 * eps)
    d = 1 - eta**2
    ell = 1 - 2 * h * eta**2
    return (2 * eta * (b - D) * value + d * dEta - 2 * eta * X * dX) / ell


def test_preceding_diffusion_matches_direct_composed_Z_oracle():
    h, order, X, eta = 0.005, 3, 0.19, 0.23
    coeffs = (1.1, -0.3, 0.2, 0.7, 0.4, -0.25)
    power = -1.0 - h
    got = preceding_diffusion_from_second_jet(
        h, power, order, X, eta, _jet(coeffs, X, eta)
    )
    expected = _independent_preceding_diffusion(h, power, order, coeffs, X, eta)
    assert got == pytest.approx(expected, rel=2e-8, abs=2e-8)


def test_point_data_builds_pinned_strict_lower_source_at_order_two():
    h, order, X, eta = 0.005, 2, 0.16, -0.21
    phi_coeffs = [
        (0.9, -0.1, 0.05, 0.2, 0.03, -0.04),
        (-0.3, 0.25, -0.1, 0.12, -0.06, 0.08),
    ]
    axial_coeffs = [
        (0.2, 0.3, -0.07, -0.11, 0.04, 0.05),
        (-0.15, -0.2, 0.09, 0.17, 0.02, -0.03),
    ]
    beta_coeffs = [
        (0.1, -0.08, 0.03, 0.06, 0.01, -0.02),
        (-0.05, 0.12, -0.04, -0.09, 0.025, 0.015),
    ]
    phi = [_jet(c, X, eta) for c in phi_coeffs]
    axial = [_jet(c, X, eta) for c in axial_coeffs]
    beta = [_jet(c, X, eta) for c in beta_coeffs]

    got = positive_axis_point_data_from_lower_history(
        h, order, X, eta, phi, axial, beta
    )

    angular_power = -1.0 - h
    axial_power = -0.5 - h
    lam1 = 2 * h
    d = 1 - eta**2
    ell = 1 - 2 * h * eta**2

    def axial_value(power, jet):
        return (
            2 * eta * power * jet.value
            + d * jet.parameter
            - 2 * eta * X * jet.radial
        ) / ell

    angular_conv = beta[1].value * (X * phi[1].radial + phi[1].value)
    angular_conv += axial[1].value * axial_value(angular_power + lam1, phi[1])
    axial_conv = beta[1].value * X * axial[1].radial
    axial_conv += axial[1].value * axial_value(axial_power + lam1, axial[1])

    expected_angular = angular_conv - _independent_preceding_diffusion(
        h, angular_power, order, phi_coeffs[1], X, eta
    )
    expected_axial = axial_conv - _independent_preceding_diffusion(
        h, axial_power, order, axial_coeffs[1], X, eta
    )
    expected_omega = omega_over_x_eq_5_6(
        1,
        X,
        eta,
        [j.regular_flux_jet() for j in beta],
        [j.value for j in axial],
        h=h,
    )

    assert got.base.phi == phi[0].first_jet()
    assert got.base.axial == axial[0].first_jet()
    assert got.base.beta == beta[0].value
    assert got.source.angular == pytest.approx(expected_angular, rel=2e-8, abs=2e-8)
    assert got.source.axial == pytest.approx(expected_axial, rel=2e-8, abs=2e-8)
    assert got.source.pressure_product == pytest.approx(phi[1].value**2, rel=1e-14)
    assert got.source.omega_quotient == pytest.approx(expected_omega, rel=1e-14)


def test_order_one_uses_only_order_zero_history_and_empty_lower_convolution():
    h, X, eta = 0.005, 0.07, 0.18
    phi = [_jet((0.7, 0.1, -0.05, 0.03, 0.02, 0.01), X, eta)]
    axial = [_jet((-0.2, 0.04, 0.02, -0.06, -0.01, 0.03), X, eta)]
    beta = [_jet((0.08, -0.02, 0.01, 0.04, 0.005, -0.015), X, eta)]
    got = positive_axis_point_data_from_lower_history(h, 1, X, eta, phi, axial, beta)
    assert got.source.pressure_product == 0.0
    assert math.isfinite(got.source.angular)
    assert math.isfinite(got.source.axial)
    assert math.isfinite(got.source.omega_quotient)


def test_lower_history_adapter_fails_closed_on_current_or_malformed_data():
    jet = ProfileSecondJet(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    with pytest.raises(ValueError, match="orders 0 through 1"):
        positive_axis_point_data_from_lower_history(
            0.005, 2, 0.1, 0.0, [jet], [jet, jet], [jet, jet]
        )
    with pytest.raises(ValueError):
        positive_axis_point_data_from_lower_history(
            0.005, 0, 0.1, 0.0, [jet], [jet], [jet]
        )
    with pytest.raises(TypeError):
        positive_axis_point_data_from_lower_history(
            0.005, 1, 0.1, 0.0, [jet], [jet], [(0.0,)]
        )
