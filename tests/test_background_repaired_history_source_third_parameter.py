import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair_phi_fifth_mixed_jets import (
    ProfileFifthMixedJet,
)
from openai_ns_reconstruction.background_positive_axis import positive_axis_forcing
from openai_ns_reconstruction.background_regular_flux_fifth_mixed_jets import (
    AxialSixthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history_forcing_third_parameter import (
    hierarchy_owned_positive_axis_forcing_third_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_sixth_mixed import (
    Section5LowerHistorySixthMixedHierarchy,
    Section5SixthMixedCoefficientJetSource,
)
from openai_ns_reconstruction.background_repaired_history_source_second_parameter import (
    hierarchy_owned_lower_source_second_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_source_third_parameter import (
    hierarchy_owned_lower_source_third_parameter_jet,
)


H = 0.005
C = 2.75


def _falling(power, derivative):
    if derivative > power:
        return 0
    out = 1
    for j in range(derivative):
        out *= power - j
    return out


def _poly_derivative(coefficients, eta, derivative):
    return math.fsum(
        float(coefficient)
        * _falling(power, derivative)
        * float(eta) ** (power - derivative)
        for power, coefficient in enumerate(coefficients)
        if power >= derivative
    )


def _x_jet(X, scale):
    X = float(X)
    scale = float(scale)
    return (
        scale * (1.0 + 0.23 * X + 0.07 * X * X),
        scale * (0.23 + 0.14 * X),
        scale * 0.14,
    )


def _phi_fifth(scale, X, eta):
    coefficients = (0.31, -0.18, 0.07, 0.025, -0.006, 0.0008)
    x0, x1, x2 = _x_jet(X, scale)
    q = [_poly_derivative(coefficients, eta, derivative) for derivative in range(6)]
    return ProfileFifthMixedJet(
        value=x0 * q[0],
        radial=x1 * q[0],
        radial2=x2 * q[0],
        parameter=x0 * q[1],
        radial_parameter=x1 * q[1],
        parameter2=x0 * q[2],
        radial2_parameter=x2 * q[1],
        radial_parameter2=x1 * q[2],
        parameter3=x0 * q[3],
        radial2_parameter2=x2 * q[2],
        radial_parameter3=x1 * q[3],
        parameter4=x0 * q[4],
        radial2_parameter3=x2 * q[3],
        radial_parameter4=x1 * q[4],
        parameter5=x0 * q[5],
    )


def _axial_sixth(scale, X, eta):
    coefficients = (0.22, 0.13, -0.055, 0.019, 0.004, -0.0009, 0.00012)
    x0, x1, x2 = _x_jet(X, scale)
    q = [_poly_derivative(coefficients, eta, derivative) for derivative in range(7)]
    return AxialSixthMixedJet(
        value=x0 * q[0],
        radial=x1 * q[0],
        radial2=x2 * q[0],
        parameter=x0 * q[1],
        radial_parameter=x1 * q[1],
        parameter2=x0 * q[2],
        radial2_parameter=x2 * q[1],
        radial_parameter2=x1 * q[2],
        parameter3=x0 * q[3],
        radial2_parameter2=x2 * q[2],
        radial_parameter3=x1 * q[3],
        parameter4=x0 * q[4],
        radial2_parameter3=x2 * q[3],
        radial_parameter4=x1 * q[4],
        parameter5=x0 * q[5],
        radial2_parameter4=x2 * q[4],
        radial_parameter5=x1 * q[5],
        parameter6=x0 * q[6],
    )


def _source(order, scale):
    def phi(X, eta):
        return _phi_fifth(scale, X, eta)

    def axial(X, eta):
        return _axial_sixth(scale, X, eta)

    return Section5SixthMixedCoefficientJetSource(
        order=order,
        phi_second_jet_provider=lambda X, eta: phi(X, eta).fourth().third().second(),
        axial_third_mixed_jet_provider=lambda X, eta: axial(X, eta).fifth().fourth().third(),
        provenance=(
            "analytic hierarchy-owned regression fixture; not a paper coefficient"
        ),
        normalization_C=None if order == 0 else C,
        axial_fourth_mixed_jet_provider=lambda X, eta: axial(X, eta).fifth().fourth(),
        phi_third_mixed_jet_provider=lambda X, eta: phi(X, eta).fourth().third(),
        axial_fifth_mixed_jet_provider=lambda X, eta: axial(X, eta).fifth(),
        phi_fourth_mixed_jet_provider=lambda X, eta: phi(X, eta).fourth(),
        phi_fifth_mixed_jet_provider=phi,
        axial_sixth_mixed_jet_provider=axial,
    )


def _hierarchy():
    return Section5LowerHistorySixthMixedHierarchy(
        H,
        C,
        (_source(0, 1.0), _source(1, 0.37)),
        quadrature_points=32,
    )


def _components(source):
    return np.array(
        [
            source.angular,
            source.axial,
            source.pressure_product,
            source.omega_quotient,
        ],
        dtype=float,
    )


def test_third_eta_lower_source_delegates_lower_rows_exactly():
    hierarchy = _hierarchy()
    actual = hierarchy_owned_lower_source_third_parameter_jet(
        hierarchy, 2, 0.73, 0.19
    )
    lower = hierarchy_owned_lower_source_second_parameter_jet(
        hierarchy, 2, 0.73, 0.19
    )
    assert actual.second() == lower


def test_third_eta_lower_source_matches_derivative_of_landed_second_eta_path():
    hierarchy = _hierarchy()
    X = 0.73
    eta = -0.21
    eps = 2.0e-5
    actual = hierarchy_owned_lower_source_third_parameter_jet(
        hierarchy, 2, X, eta
    )
    plus = hierarchy_owned_lower_source_second_parameter_jet(
        hierarchy, 2, X, eta + eps
    ).parameter2
    minus = hierarchy_owned_lower_source_second_parameter_jet(
        hierarchy, 2, X, eta - eps
    ).parameter2
    oracle = (_components(plus) - _components(minus)) / (2.0 * eps)
    np.testing.assert_allclose(
        _components(actual.parameter3),
        oracle,
        rtol=5.0e-5,
        atol=5.0e-7,
    )


def test_third_eta_lower_source_is_regular_on_axis():
    hierarchy = _hierarchy()
    actual = hierarchy_owned_lower_source_third_parameter_jet(
        hierarchy, 1, 0.0, 0.17
    )
    assert np.isfinite(_components(actual.parameter3)).all()
    assert actual.second() == hierarchy_owned_lower_source_second_parameter_jet(
        hierarchy, 1, 0.0, 0.17
    )


def test_third_eta_lower_source_fails_closed_on_wrong_hierarchy_type():
    with pytest.raises(TypeError, match="Section5LowerHistorySixthMixedHierarchy"):
        hierarchy_owned_lower_source_third_parameter_jet(object(), 1, 0.0, 0.0)


def test_hierarchy_owned_positive_axis_forcing_value_delegates_exactly():
    hierarchy = _hierarchy()
    xi = 0.81
    eta = 0.16
    actual = hierarchy_owned_positive_axis_forcing_third_parameter_jet(
        hierarchy, 2, xi, eta
    )
    source = hierarchy_owned_lower_source_third_parameter_jet(
        hierarchy, 2, xi * xi, eta
    )
    expected = positive_axis_forcing(H, C, xi, eta, source.value)
    np.testing.assert_array_equal(actual.value, expected)


def test_hierarchy_owned_positive_axis_forcing_third_eta_matches_second_row_derivative():
    hierarchy = _hierarchy()
    xi = 0.77
    eta = -0.18
    eps = 2.0e-5
    actual = hierarchy_owned_positive_axis_forcing_third_parameter_jet(
        hierarchy, 2, xi, eta
    )
    plus = hierarchy_owned_positive_axis_forcing_third_parameter_jet(
        hierarchy, 2, xi, eta + eps
    ).parameter2
    minus = hierarchy_owned_positive_axis_forcing_third_parameter_jet(
        hierarchy, 2, xi, eta - eps
    ).parameter2
    oracle = (plus - minus) / (2.0 * eps)
    np.testing.assert_allclose(actual.parameter3, oracle, rtol=7.0e-5, atol=7.0e-7)


def test_hierarchy_owned_positive_axis_forcing_is_regular_on_axis():
    hierarchy = _hierarchy()
    actual = hierarchy_owned_positive_axis_forcing_third_parameter_jet(
        hierarchy, 1, 0.0, 0.17
    )
    assert np.isfinite(actual.value).all()
    assert np.isfinite(actual.parameter).all()
    assert np.isfinite(actual.parameter2).all()
    assert np.isfinite(actual.parameter3).all()
    assert actual.parameter3[3] == 0.0


def test_hierarchy_owned_positive_axis_forcing_fails_closed_on_wrong_hierarchy():
    with pytest.raises(TypeError, match="Section5LowerHistorySixthMixedHierarchy"):
        hierarchy_owned_positive_axis_forcing_third_parameter_jet(
            object(), 1, 0.0, 0.0
        )
