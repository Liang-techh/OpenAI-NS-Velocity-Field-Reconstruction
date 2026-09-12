import numpy as np
import pytest

from openai_ns_reconstruction.background_first_picard_parameter_jet import (
    first_picard_parameter_jet_eq_5_7,
)
from openai_ns_reconstruction.background_first_picard_second_parameter_jet import (
    FirstPicardSecondParameterJet,
    first_picard_second_parameter_jet_eq_5_7,
)
from openai_ns_reconstruction.background_inner_solver import EQ_5_7_SINGULAR_WEIGHTS


def _analytic_source_data():
    amplitudes = np.array([0.7, -0.4, 1.1, 0.6, -0.8, 0.9])
    slopes = np.array([0.2, -0.3, 0.1, 0.5, -0.2, 0.4])
    curvatures = np.array([0.05, 0.03, -0.02, 0.04, 0.01, -0.03])
    cubics = np.array([-0.01, 0.02, 0.03, -0.02, 0.015, 0.01])
    powers = np.array([0, 1, 2, 3, 1, 2], dtype=int)

    def source(xi: float, eta: float):
        factor = 1.0 + slopes * eta + curvatures * eta**2 + cubics * eta**3
        return amplitudes * xi**powers * factor

    def source_parameter(xi: float, eta: float):
        factor = slopes + 2.0 * curvatures * eta + 3.0 * cubics * eta**2
        return amplitudes * xi**powers * factor

    def source_second_parameter(xi: float, eta: float):
        factor = 2.0 * curvatures + 6.0 * cubics * eta
        return amplitudes * xi**powers * factor

    return (
        amplitudes,
        slopes,
        curvatures,
        cubics,
        powers,
        source,
        source_parameter,
        source_second_parameter,
    )


def test_second_parameter_jet_matches_independent_closed_form_integral() -> None:
    (
        amplitudes,
        slopes,
        curvatures,
        cubics,
        powers,
        source,
        source_parameter,
        source_second_parameter,
    ) = _analytic_source_data()
    xi, eta = 0.71, -0.19

    got = first_picard_second_parameter_jet_eq_5_7(
        2,
        xi,
        eta,
        source,
        source_parameter,
        source_second_parameter,
        quadrature_points=24,
    )

    weights = np.asarray(EQ_5_7_SINGULAR_WEIGHTS, dtype=float)
    radial = amplitudes * xi ** (powers + 1) / (powers + weights + 1.0)
    expected_value = radial * (
        1.0 + slopes * eta + curvatures * eta**2 + cubics * eta**3
    )
    expected_parameter = radial * (
        slopes + 2.0 * curvatures * eta + 3.0 * cubics * eta**2
    )
    expected_second = radial * (2.0 * curvatures + 6.0 * cubics * eta)

    assert isinstance(got, FirstPicardSecondParameterJet)
    assert np.allclose(got.value, expected_value, rtol=2e-13, atol=2e-14)
    assert np.allclose(got.parameter, expected_parameter, rtol=2e-13, atol=2e-14)
    assert np.allclose(got.second_parameter, expected_second, rtol=2e-13, atol=2e-14)
    assert not got.value.flags.writeable
    assert not got.parameter.flags.writeable
    assert not got.second_parameter.flags.writeable


def test_second_parameter_matches_eta_difference_of_landed_first_parameter_path() -> None:
    *_, source, source_parameter, source_second_parameter = _analytic_source_data()
    xi, eta = 0.63, 0.16
    got = first_picard_second_parameter_jet_eq_5_7(
        1,
        xi,
        eta,
        source,
        source_parameter,
        source_second_parameter,
        quadrature_points=32,
    )

    eps = 2.0e-6
    plus = first_picard_parameter_jet_eq_5_7(
        1, xi, eta + eps, source, source_parameter, quadrature_points=32
    ).parameter
    minus = first_picard_parameter_jet_eq_5_7(
        1, xi, eta - eps, source, source_parameter, quadrature_points=32
    ).parameter
    independent_fd = (plus - minus) / (2.0 * eps)

    assert np.allclose(got.second_parameter, independent_fd, rtol=3e-9, atol=3e-10)


def test_axis_value_and_both_parameter_derivatives_are_exactly_zero() -> None:
    *_, source, source_parameter, source_second_parameter = _analytic_source_data()
    got = first_picard_second_parameter_jet_eq_5_7(
        1, 0.0, 0.2, source, source_parameter, source_second_parameter
    )
    assert np.array_equal(got.value, np.zeros(6))
    assert np.array_equal(got.parameter, np.zeros(6))
    assert np.array_equal(got.second_parameter, np.zeros(6))


def test_second_derivative_provider_is_required_and_shape_checked() -> None:
    *_, source, source_parameter, _ = _analytic_source_data()
    with pytest.raises(TypeError, match="lower_order_source_second_parameter"):
        first_picard_second_parameter_jet_eq_5_7(
            1, 0.4, 0.1, source, source_parameter, None  # type: ignore[arg-type]
        )

    def bad_second_parameter(_xi: float, _eta: float):
        return [1.0, 2.0]

    with pytest.raises(ValueError, match=r"shape \(6,\)"):
        first_picard_second_parameter_jet_eq_5_7(
            1, 0.4, 0.1, source, source_parameter, bad_second_parameter
        )


def test_domain_and_positive_order_gates_remain_fail_closed() -> None:
    *_, source, source_parameter, source_second_parameter = _analytic_source_data()
    with pytest.raises(ValueError, match="positive integer"):
        first_picard_second_parameter_jet_eq_5_7(
            0, 0.4, 0.1, source, source_parameter, source_second_parameter
        )
    with pytest.raises(ValueError, match=r"\|eta\| <= 1"):
        first_picard_second_parameter_jet_eq_5_7(
            1, 0.4, 1.01, source, source_parameter, source_second_parameter
        )
