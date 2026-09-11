import numpy as np
import pytest

from openai_ns_reconstruction.background_first_picard_parameter_jet import (
    FirstPicardParameterJet,
    first_picard_parameter_jet_eq_5_7,
)
from openai_ns_reconstruction.background_inner_solver import (
    EQ_5_7_SINGULAR_WEIGHTS,
    first_picard_term_eq_5_7,
)


def _analytic_source_data():
    amplitudes = np.array([0.8, -0.3, 1.2, 0.5, -0.7, 0.9])
    slopes = np.array([0.2, -0.4, 0.1, 0.6, -0.3, 0.5])
    curvatures = np.array([0.05, 0.02, -0.03, 0.04, 0.01, -0.02])
    powers = np.array([0, 1, 2, 3, 1, 2], dtype=int)

    def source(xi: float, eta: float):
        factor = 1.0 + slopes * eta + curvatures * eta * eta
        return amplitudes * xi**powers * factor

    def source_parameter(xi: float, eta: float):
        return amplitudes * xi**powers * (slopes + 2.0 * curvatures * eta)

    return amplitudes, slopes, curvatures, powers, source, source_parameter


def test_parameter_jet_matches_independent_closed_form_integral() -> None:
    amplitudes, slopes, curvatures, powers, source, source_parameter = _analytic_source_data()
    xi, eta = 0.73, -0.21

    got = first_picard_parameter_jet_eq_5_7(
        2, xi, eta, source, source_parameter, quadrature_points=24
    )

    weights = np.asarray(EQ_5_7_SINGULAR_WEIGHTS, dtype=float)
    denominators = powers + weights + 1.0
    radial = amplitudes * xi ** (powers + 1) / denominators
    expected_value = radial * (1.0 + slopes * eta + curvatures * eta * eta)
    expected_parameter = radial * (slopes + 2.0 * curvatures * eta)

    assert isinstance(got, FirstPicardParameterJet)
    assert np.allclose(got.value, expected_value, rtol=2e-13, atol=2e-14)
    assert np.allclose(got.parameter, expected_parameter, rtol=2e-13, atol=2e-14)
    assert not got.value.flags.writeable
    assert not got.parameter.flags.writeable


def test_parameter_jet_agrees_with_independent_eta_finite_difference_of_existing_integral() -> None:
    _, _, _, _, source, source_parameter = _analytic_source_data()
    xi, eta = 0.61, 0.17
    got = first_picard_parameter_jet_eq_5_7(
        1, xi, eta, source, source_parameter, quadrature_points=32
    )

    eps = 2.0e-6
    plus = first_picard_term_eq_5_7(
        1, xi, eta + eps, source, quadrature_points=32
    )
    minus = first_picard_term_eq_5_7(
        1, xi, eta - eps, source, quadrature_points=32
    )
    independent_fd = (plus - minus) / (2.0 * eps)
    assert np.allclose(got.parameter, independent_fd, rtol=2e-9, atol=2e-10)


def test_axis_value_and_parameter_derivative_are_exactly_zero() -> None:
    _, _, _, _, source, source_parameter = _analytic_source_data()
    got = first_picard_parameter_jet_eq_5_7(1, 0.0, 0.3, source, source_parameter)
    assert np.array_equal(got.value, np.zeros(6))
    assert np.array_equal(got.parameter, np.zeros(6))


def test_derivative_provider_is_required_and_shape_checked() -> None:
    _, _, _, _, source, _ = _analytic_source_data()
    with pytest.raises(TypeError, match="lower_order_source_parameter"):
        first_picard_parameter_jet_eq_5_7(1, 0.4, 0.1, source, None)  # type: ignore[arg-type]

    def bad_parameter(_xi: float, _eta: float):
        return [1.0, 2.0]

    with pytest.raises(ValueError, match=r"shape \(6,\)"):
        first_picard_parameter_jet_eq_5_7(1, 0.4, 0.1, source, bad_parameter)


def test_domain_and_positive_order_gates_remain_fail_closed() -> None:
    _, _, _, _, source, source_parameter = _analytic_source_data()
    with pytest.raises(ValueError, match="positive integer"):
        first_picard_parameter_jet_eq_5_7(0, 0.4, 0.1, source, source_parameter)
    with pytest.raises(ValueError, match=r"\|eta\| <= 1"):
        first_picard_parameter_jet_eq_5_7(1, 0.4, 1.01, source, source_parameter)
