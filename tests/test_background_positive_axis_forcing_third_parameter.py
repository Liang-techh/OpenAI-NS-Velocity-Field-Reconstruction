from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.background_positive_axis import (
    PositiveAxisSourceJet,
    pressure_source,
)
from openai_ns_reconstruction.background_positive_axis_forcing_third_parameter import (
    PositiveAxisSourceThirdParameterJet,
    positive_axis_forcing_third_parameter_from_source_jet,
)


_COEFFS = {
    "angular": np.array([0.4, -0.7, 0.2, 0.5]),
    "axial": np.array([-0.3, 0.6, -0.4, 0.25]),
    "pressure_product": np.array([1.2, -0.5, 0.3, -0.2]),
    "omega_quotient": np.array([0.1, 0.4, -0.1, 0.15]),
}


def _poly_derivative(coefficients: np.ndarray, eta: float, order: int) -> float:
    coefficients = np.asarray(coefficients, dtype=float)
    for _ in range(order):
        coefficients = np.array(
            [(j + 1) * coefficients[j + 1] for j in range(len(coefficients) - 1)],
            dtype=float,
        )
    return float(sum(c * eta**j for j, c in enumerate(coefficients)))


def _source_at(eta: float, order: int) -> PositiveAxisSourceJet:
    return PositiveAxisSourceJet(
        angular=_poly_derivative(_COEFFS["angular"], eta, order),
        axial=_poly_derivative(_COEFFS["axial"], eta, order),
        pressure_product=_poly_derivative(_COEFFS["pressure_product"], eta, order),
        omega_quotient=_poly_derivative(_COEFFS["omega_quotient"], eta, order),
    )


def _source_third_jet(eta: float) -> PositiveAxisSourceThirdParameterJet:
    return PositiveAxisSourceThirdParameterJet(
        value=_source_at(eta, 0),
        parameter=_source_at(eta, 1),
        parameter2=_source_at(eta, 2),
        parameter3=_source_at(eta, 3),
    )


def _forcing_second(h: float, C: float, xi: float, eta: float) -> np.ndarray:
    source0 = _source_at(eta, 0)
    source1 = _source_at(eta, 1)
    source2 = _source_at(eta, 2)
    p0 = pressure_source(C, source0)
    p1 = pressure_source(C, source1)
    p2 = pressure_source(C, source2)
    ell = 1.0 - 2.0 * h * eta * eta
    g = eta / ell
    g1 = 1.0 / ell + 4.0 * h * eta * eta / (ell * ell)
    g2 = 12.0 * h * eta / (ell * ell) + 32.0 * h * h * eta**3 / ell**3
    X = xi * xi
    out = np.zeros(6)
    out[3] = 2.0 * xi * p2
    out[4] = 2.0 * source2.angular
    out[5] = 2.0 * source2.axial - 4.0 * X * (g2 * p0 + 2.0 * g1 * p1 + g * p2)
    return out


def test_third_forcing_matches_test_only_derivative_of_exact_second_formula() -> None:
    h = 0.007
    C = 1.6
    xi = 0.73
    eta = 0.23
    actual = positive_axis_forcing_third_parameter_from_source_jet(
        h, C, xi, eta, _source_third_jet(eta)
    )

    eps = 1.0e-5
    expected = (
        _forcing_second(h, C, xi, eta + eps)
        - _forcing_second(h, C, xi, eta - eps)
    ) / (2.0 * eps)
    np.testing.assert_allclose(actual, expected, rtol=2.0e-7, atol=2.0e-8)
    assert not actual.flags.writeable


def test_pressure_third_jet_enters_only_the_exact_expected_rows() -> None:
    h = 0.009
    C = 1.9
    xi = 0.61
    eta = -0.31
    zero = PositiveAxisSourceJet(0.0, 0.0, 0.0, 0.0)
    delta = 0.47
    source = PositiveAxisSourceThirdParameterJet(
        value=zero,
        parameter=zero,
        parameter2=zero,
        parameter3=PositiveAxisSourceJet(0.0, 0.0, delta, 0.0),
    )
    actual = positive_axis_forcing_third_parameter_from_source_jet(
        h, C, xi, eta, source
    )

    ell = 1.0 - 2.0 * h * eta * eta
    g = eta / ell
    dp3 = delta / (C * C)
    expected = np.zeros(6)
    expected[3] = 2.0 * xi * dp3
    expected[5] = -4.0 * xi * xi * g * dp3
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=2.0e-15)


def test_axis_does_not_hide_nonpressure_third_source_data() -> None:
    zero = PositiveAxisSourceJet(0.0, 0.0, 0.0, 0.0)
    source = PositiveAxisSourceThirdParameterJet(
        value=zero,
        parameter=zero,
        parameter2=zero,
        parameter3=PositiveAxisSourceJet(0.35, -0.2, 0.0, 0.0),
    )
    actual = positive_axis_forcing_third_parameter_from_source_jet(
        0.006, 1.4, 0.0, 0.17, source
    )
    np.testing.assert_allclose(actual, np.array([0.0, 0.0, 0.0, 0.0, 0.7, -0.4]))


def test_third_forcing_primitive_fails_closed_on_malformed_inputs() -> None:
    zero = PositiveAxisSourceJet(0.0, 0.0, 0.0, 0.0)
    with pytest.raises(TypeError):
        PositiveAxisSourceThirdParameterJet(zero, zero, zero, object())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        positive_axis_forcing_third_parameter_from_source_jet(  # type: ignore[arg-type]
            0.006, 1.4, 0.5, 0.1, object()
        )
    source = PositiveAxisSourceThirdParameterJet(zero, zero, zero, zero)
    with pytest.raises(ValueError):
        positive_axis_forcing_third_parameter_from_source_jet(
            0.006, 1.4, -0.1, 0.1, source
        )
    with pytest.raises(ValueError):
        positive_axis_forcing_third_parameter_from_source_jet(
            0.006, 1.4, 0.1, 1.1, source
        )
