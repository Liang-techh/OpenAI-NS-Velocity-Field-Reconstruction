from types import SimpleNamespace

import numpy as np
import pytest

from openai_ns_reconstruction.background_inner_solver import (
    picard_map_eq_5_7,
    singular_inverse_eq_5_7,
)
from openai_ns_reconstruction.background_picard_second_parameter_jet import (
    picard_second_parameter_jet_eq_5_7,
)


M0 = np.diag([1.0, 1.2, 0.8, 1.1, 0.9, 1.3])
N0 = np.array(
    [
        [0.0, 0.2, 0.0, 0.0, 0.0, 0.0],
        [-0.1, 0.0, 0.1, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.15, 0.0, 0.0],
        [0.05, 0.0, 0.0, 0.0, 0.1, 0.0],
        [0.0, 0.0, -0.08, 0.0, 0.0, 0.12],
        [0.0, 0.04, 0.0, -0.06, 0.0, 0.0],
    ]
)
P0 = np.diag([0.04, -0.03, 0.02, 0.01, -0.02, 0.03])
M1 = np.diag([0.3, -0.2, 0.25, 0.15, -0.1, 0.2])
N1 = np.array(
    [
        [0.0, 0.0, 0.05, 0.0, 0.0, 0.0],
        [0.03, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, -0.04, 0.0, 0.0, 0.02, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, -0.03],
        [0.02, 0.0, 0.0, 0.01, 0.0, 0.0],
        [0.0, 0.0, 0.04, 0.0, -0.02, 0.0],
    ]
)
P1 = np.diag([-0.02, 0.03, 0.01, -0.01, 0.025, -0.015])

V0 = np.array([0.4, -0.3, 0.2, 0.5, -0.1, 0.25])
V1 = np.array([-0.15, 0.2, 0.1, -0.05, 0.3, -0.2])
V2 = np.array([0.05, -0.04, 0.08, 0.02, -0.06, 0.03])
V3 = np.array([0.02, 0.01, -0.015, 0.025, 0.005, -0.02])
F0 = np.array([0.2, 0.1, -0.2, 0.3, -0.1, 0.05])
F1 = np.array([-0.04, 0.06, 0.02, -0.03, 0.01, 0.05])
F2 = np.array([0.015, -0.01, 0.02, 0.005, -0.012, 0.008])


def A0(s, eta):
    return (1.0 + 0.2 * s) * M0 + eta * N0 + eta * eta * P0


def A1(s, eta):
    return (0.7 + 0.1 * s) * M1 + eta * N1 + eta * eta * P1


def matrix_jets(s, eta):
    return SimpleNamespace(
        A0=SimpleNamespace(
            value=A0(s, eta),
            parameter=N0 + 2.0 * eta * P0,
            parameter2=2.0 * P0,
        ),
        A1=SimpleNamespace(
            value=A1(s, eta),
            parameter=N1 + 2.0 * eta * P1,
            parameter2=2.0 * P1,
        ),
    )


def previous(s, eta):
    scale = 1.0 + 0.25 * s
    return scale * (V0 + eta * V1 + eta**2 * V2 + eta**3 * V3)


def previous_parameter(s, eta):
    scale = 1.0 + 0.25 * s
    return scale * (V1 + 2.0 * eta * V2 + 3.0 * eta**2 * V3)


def previous_jet(s, eta):
    scale = 1.0 + 0.25 * s
    return SimpleNamespace(
        value=previous(s, eta),
        parameter=previous_parameter(s, eta),
        parameter2=scale * (2.0 * V2 + 6.0 * eta * V3),
        parameter3=scale * (6.0 * V3),
    )


def forcing(s, eta):
    scale = 1.0 + 0.1 * s * s
    return scale * (F0 + eta * F1 + eta * eta * F2)


def forcing_jet(s, eta):
    scale = 1.0 + 0.1 * s * s
    return SimpleNamespace(
        value=forcing(s, eta),
        parameter=scale * (F1 + 2.0 * eta * F2),
        parameter2=scale * (2.0 * F2),
    )


def direct_picard(xi, eta):
    return picard_map_eq_5_7(
        2,
        xi,
        eta,
        previous,
        previous_parameter,
        A0,
        A1,
        forcing,
        quadrature_points=20,
    )


def test_second_eta_picard_jet_matches_landed_value_and_test_only_difference_oracle():
    xi = 0.43
    eta = 0.17
    step = 1.0e-4

    result = picard_second_parameter_jet_eq_5_7(
        2,
        xi,
        eta,
        matrix_jets,
        previous_jet,
        forcing_jet,
        quadrature_points=20,
    )
    center = direct_picard(xi, eta)
    plus = direct_picard(xi, eta + step)
    minus = direct_picard(xi, eta - step)

    np.testing.assert_allclose(result.value, center, rtol=2.0e-13, atol=2.0e-13)
    np.testing.assert_allclose(
        result.parameter,
        (plus - minus) / (2.0 * step),
        rtol=2.0e-8,
        atol=2.0e-9,
    )
    np.testing.assert_allclose(
        result.parameter2,
        (plus - 2.0 * center + minus) / (step * step),
        rtol=2.0e-6,
        atol=2.0e-7,
    )


def test_previous_third_eta_derivative_enters_only_second_output_derivative():
    xi = 0.39
    eta = -0.21
    delta = np.array([0.07, -0.04, 0.02, 0.03, -0.01, 0.05])

    baseline = picard_second_parameter_jet_eq_5_7(
        2, xi, eta, matrix_jets, previous_jet, forcing_jet, quadrature_points=20
    )

    def perturbed_previous_jet(s, eta_value):
        jet = previous_jet(s, eta_value)
        return SimpleNamespace(
            value=jet.value,
            parameter=jet.parameter,
            parameter2=jet.parameter2,
            parameter3=jet.parameter3 + (1.0 + 0.25 * s) * delta,
        )

    perturbed = picard_second_parameter_jet_eq_5_7(
        2,
        xi,
        eta,
        matrix_jets,
        perturbed_previous_jet,
        forcing_jet,
        quadrature_points=20,
    )

    expected_delta = singular_inverse_eq_5_7(
        xi,
        eta,
        lambda s, e: A1(s, e) @ ((1.0 + 0.25 * s) * delta),
        quadrature_points=20,
    )
    np.testing.assert_array_equal(perturbed.value, baseline.value)
    np.testing.assert_array_equal(perturbed.parameter, baseline.parameter)
    np.testing.assert_allclose(
        perturbed.parameter2 - baseline.parameter2,
        expected_delta,
        rtol=3.0e-13,
        atol=3.0e-13,
    )


def test_axis_preflight_does_not_hide_missing_third_eta_data():
    def malformed_previous_jet(s, eta):
        jet = previous_jet(s, eta)
        return SimpleNamespace(
            value=jet.value,
            parameter=jet.parameter,
            parameter2=jet.parameter2,
            parameter3=np.zeros(5),
        )

    with pytest.raises(ValueError, match="previous.parameter3"):
        picard_second_parameter_jet_eq_5_7(
            2,
            0.0,
            0.1,
            matrix_jets,
            malformed_previous_jet,
            forcing_jet,
        )
