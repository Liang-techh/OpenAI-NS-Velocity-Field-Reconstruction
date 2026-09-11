import numpy as np
import pytest

from openai_ns_reconstruction.background_lower_history_source import ProfileSecondJet
from openai_ns_reconstruction.background_omega_parameter_jet import (
    RegularFluxThirdMixedJet,
    omega_over_x_parameter_jet_eq_5_6,
)
from openai_ns_reconstruction.background_recurrence import omega_over_x_eq_5_6


def _third_polynomial(scale: float, X: float, eta: float) -> RegularFluxThirdMixedJet:
    # Deliberately includes X^2 eta, X eta^2 and eta^3 so every additional
    # third-mixed field enters the production derivative.
    value = scale * (
        1.0
        + 0.2 * X
        + 0.3 * X * X
        + 0.4 * eta
        + 0.5 * X * eta
        + 0.6 * eta * eta
        + 0.7 * X * X * eta
        + 0.8 * X * eta * eta
        + 0.9 * eta**3
    )
    radial = scale * (
        0.2 + 0.6 * X + 0.5 * eta + 1.4 * X * eta + 0.8 * eta * eta
    )
    radial2 = scale * (0.6 + 1.4 * eta)
    parameter = scale * (
        0.4
        + 0.5 * X
        + 1.2 * eta
        + 0.7 * X * X
        + 1.6 * X * eta
        + 2.7 * eta * eta
    )
    radial_parameter = scale * (0.5 + 1.4 * X + 1.6 * eta)
    parameter2 = scale * (1.2 + 1.6 * X + 5.4 * eta)
    return RegularFluxThirdMixedJet(
        value=value,
        radial=radial,
        radial2=radial2,
        parameter=parameter,
        radial_parameter=radial_parameter,
        parameter2=parameter2,
        radial2_parameter=1.4 * scale,
        radial_parameter2=1.6 * scale,
        parameter3=5.4 * scale,
    )


def _axial_polynomial(scale: float, X: float, eta: float) -> ProfileSecondJet:
    beta = _third_polynomial(scale, X, eta)
    return ProfileSecondJet(
        value=beta.value,
        radial=beta.radial,
        radial2=beta.radial2,
        parameter=beta.parameter,
        radial_parameter=beta.radial_parameter,
        parameter2=beta.parameter2,
    )


def _history(order: int, X: float, eta: float):
    beta = [
        _third_polynomial(0.7 + 0.11 * j, X, eta) for j in range(order + 1)
    ]
    axial = [
        _axial_polynomial(1.1 - 0.07 * j, X, eta) for j in range(order + 1)
    ]
    return beta, axial


def _value_only(order: int, X: float, eta: float, h: float) -> float:
    beta, axial = _history(order, X, eta)
    return omega_over_x_eq_5_6(
        order,
        X,
        eta,
        [jet.second() for jet in beta],
        [jet.value for jet in axial],
        h=h,
    )


@pytest.mark.parametrize("order", [0, 1, 2])
@pytest.mark.parametrize("X", [0.0, 0.37, 1.1])
@pytest.mark.parametrize("eta", [-0.31, 0.21])
def test_omega_parameter_jet_matches_independent_value_finite_difference(
    order: int, X: float, eta: float
):
    h = 0.005
    beta, axial = _history(order, X, eta)
    got = omega_over_x_parameter_jet_eq_5_6(
        order, X, eta, beta, axial, h=h
    )

    assert got.value == pytest.approx(_value_only(order, X, eta, h), rel=2e-13, abs=2e-13)

    eps = 2.0e-6
    finite_difference = (
        _value_only(order, X, eta + eps, h)
        - _value_only(order, X, eta - eps, h)
    ) / (2.0 * eps)
    assert got.parameter == pytest.approx(finite_difference, rel=3e-7, abs=2e-8)


def test_third_mixed_terms_are_not_silently_ignored():
    h = 0.005
    order, X, eta = 1, 0.43, 0.19
    beta, axial = _history(order, X, eta)
    reference = omega_over_x_parameter_jet_eq_5_6(
        order, X, eta, beta, axial, h=h
    ).parameter

    changed = list(beta)
    j = changed[0]
    changed[0] = RegularFluxThirdMixedJet(
        value=j.value,
        radial=j.radial,
        radial2=j.radial2,
        parameter=j.parameter,
        radial_parameter=j.radial_parameter,
        parameter2=j.parameter2,
        radial2_parameter=j.radial2_parameter + 0.37,
        radial_parameter2=j.radial_parameter2 - 0.21,
        parameter3=j.parameter3 + 0.44,
    )
    perturbed = omega_over_x_parameter_jet_eq_5_6(
        order, X, eta, changed, axial, h=h
    ).parameter
    assert not np.isclose(reference, perturbed, rtol=0.0, atol=1e-8)


def test_omega_parameter_jet_fails_closed_on_incomplete_or_wrong_history():
    beta, axial = _history(1, 0.2, 0.1)
    with pytest.raises(ValueError, match="orders 0 through 1"):
        omega_over_x_parameter_jet_eq_5_6(1, 0.2, 0.1, beta[:1], axial, h=0.005)
    with pytest.raises(TypeError, match="RegularFluxThirdMixedJet"):
        omega_over_x_parameter_jet_eq_5_6(1, 0.2, 0.1, [object(), beta[1]], axial, h=0.005)
    with pytest.raises(ValueError, match="nonnegative integer"):
        omega_over_x_parameter_jet_eq_5_6(-1, 0.2, 0.1, beta, axial, h=0.005)
