import math

import numpy as np
import pytest
from scipy.integrate import quad

from openai_ns_reconstruction.background_lower_history_solver import (
    first_positive_order_term_from_lower_history,
    lower_history_positive_axis_providers,
    positive_axis_eq_5_7_fields_from_lower_history,
)
from openai_ns_reconstruction.background_lower_history_source import ProfileSecondJet


def _x_polynomial_jet(coeffs, X, eta):
    """Second jet of c0+c1 X+c2 X^2, independent of eta."""
    del eta
    c0, c1, c2 = coeffs
    return ProfileSecondJet(
        value=c0 + c1 * X + c2 * X * X,
        radial=c1 + 2.0 * c2 * X,
        radial2=2.0 * c2,
        parameter=0.0,
        radial_parameter=0.0,
        parameter2=0.0,
    )


def test_first_picard_term_from_order_one_history_matches_independent_closed_form_integrals():
    """Cross-check G f_1 without using the production source/Volterra evaluator."""
    h = 0.005
    C = 3.0
    eta = 0.0
    xi = 0.73
    phi_coeffs = (1.2, -0.4, 0.15)
    axial_coeffs = (-0.3, 0.55, -0.08)
    beta_coeffs = (0.4, -0.2, 0.07)

    def phi_provider(j, X, e):
        assert j == 0
        return _x_polynomial_jet(phi_coeffs, X, e)

    def axial_provider(j, X, e):
        assert j == 0
        return _x_polynomial_jet(axial_coeffs, X, e)

    def beta_provider(j, X, e):
        assert j == 0
        return _x_polynomial_jet(beta_coeffs, X, e)

    got = first_positive_order_term_from_lower_history(
        1,
        xi,
        eta,
        h=h,
        C=C,
        phi_provider=phi_provider,
        axial_provider=axial_provider,
        beta_provider=beta_provider,
        quadrature_points=96,
    )

    angular_power = -1.0 - h
    axial_power = -0.5 - h

    def poly(c, X):
        return c[0] + c[1] * X + c[2] * X * X

    def dpoly(c, X):
        return c[1] + 2.0 * c[2] * X

    def d2poly(c, X):
        return 2.0 * c[2]

    # At eta=0 and n=1, the strict-lower convolutions are empty. Directly
    # composing the paper Z operators gives
    #   actualLowerAngular = 2 X phi_0' - 2 b_phi phi_0,
    # and the analogous formula for U_0.  Eq. (5.6) at k=0 reduces to the
    # explicit regular polynomial expression below.
    def independent_forcing(s):
        X = s * s
        phi = poly(phi_coeffs, X)
        phi_x = dpoly(phi_coeffs, X)
        axial = poly(axial_coeffs, X)
        axial_x = dpoly(axial_coeffs, X)
        beta = poly(beta_coeffs, X)
        beta_x = dpoly(beta_coeffs, X)
        beta_xx = d2poly(beta_coeffs, X)

        lower_angular = 2.0 * X * phi_x - 2.0 * angular_power * phi
        lower_axial = 2.0 * X * axial_x - 2.0 * axial_power * axial
        omega_over_x = (
            beta
            + X * beta_x
            + beta * (0.5 * beta + X * beta_x)
            - 2.0 * (2.0 * beta_x + X * beta_xx)
        )
        pressure_source = -0.5 * omega_over_x
        return np.array(
            [
                0.0,
                0.0,
                0.0,
                2.0 * s * pressure_source,
                2.0 * lower_angular,
                2.0 * lower_axial,
            ]
        )

    weights = (0, 0, 2, 0, 3, 1)
    expected = np.array(
        [
            quad(
                lambda s, i=i, weight=weight: (s / xi) ** weight
                * independent_forcing(s)[i],
                0.0,
                xi,
                epsabs=1e-12,
                epsrel=1e-12,
            )[0]
            for i, weight in enumerate(weights)
        ]
    )
    assert got == pytest.approx(expected, rel=3e-12, abs=3e-12)
    assert np.all(got[:3] == 0.0)
    assert np.linalg.norm(got[3:]) > 0.0


def test_fields_query_only_strict_lower_orders_and_use_X_equal_xi_squared():
    h = 0.005
    calls = {"phi": [], "axial": [], "beta": []}

    def make_provider(name, base):
        def provider(j, X, eta):
            assert j < 2
            calls[name].append((j, X, eta))
            return _x_polynomial_jet((base + j, 0.1 * (j + 1), 0.0), X, eta)

        return provider

    fields = positive_axis_eq_5_7_fields_from_lower_history(
        h,
        2,
        2.5,
        make_provider("phi", 1.0),
        make_provider("axial", -0.2),
        make_provider("beta", 0.3),
    )
    xi, eta = 0.41, -0.17
    forcing = fields.forcing(xi, eta)
    assert forcing.shape == (6,)
    assert np.all(np.isfinite(forcing))
    for name in calls:
        assert [entry[0] for entry in calls[name]] == [0, 1]
        assert all(entry[1] == pytest.approx(xi * xi) for entry in calls[name])
        assert all(entry[2] == eta for entry in calls[name])


def test_provider_adapter_fails_closed_on_wrong_type_or_invalid_geometry():
    good = lambda j, X, eta: _x_polynomial_jet((1.0, 0.0, 0.0), X, eta)
    bad = lambda j, X, eta: (1.0, 0.0)

    providers = lower_history_positive_axis_providers(0.005, 1, good, good, bad)
    with pytest.raises(TypeError, match="beta_provider must return ProfileSecondJet"):
        providers.source(0.2, 0.0)

    with pytest.raises(ValueError, match="order must be a positive integer"):
        lower_history_positive_axis_providers(0.005, 0, good, good, good)
    with pytest.raises(TypeError, match="phi_provider must be callable"):
        lower_history_positive_axis_providers(0.005, 1, None, good, good)

    providers = lower_history_positive_axis_providers(0.005, 1, good, good, good)
    with pytest.raises(ValueError, match="xi must be finite and nonnegative"):
        providers.source(-0.1, 0.0)
    with pytest.raises(ValueError, match="eta must be finite"):
        providers.base(0.1, 1.2)
