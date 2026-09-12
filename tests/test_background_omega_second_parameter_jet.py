from dataclasses import replace

import numpy as np
import pytest

from openai_ns_reconstruction.background_lower_history_source import ProfileSecondJet
from openai_ns_reconstruction.background_omega_parameter_jet import (
    omega_over_x_parameter_jet_eq_5_6,
)
from openai_ns_reconstruction.background_omega_second_parameter_jet import (
    RegularFluxFourthMixedJet,
    omega_over_x_second_parameter_jet_eq_5_6,
)


def _fourth_polynomial(
    scale: float, X: float, eta: float
) -> RegularFluxFourthMixedJet:
    # Degree two in X and degree four in eta.  The final three monomials make
    # beta_XXetaeta, beta_Xetaetaeta and beta_etaetaetaeta all nonzero.
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
        + 1.1 * X * X * eta * eta
        + 1.2 * X * eta**3
        + 1.3 * eta**4
    )
    radial = scale * (
        0.2
        + 0.6 * X
        + 0.5 * eta
        + 1.4 * X * eta
        + 0.8 * eta * eta
        + 2.2 * X * eta * eta
        + 1.2 * eta**3
    )
    radial2 = scale * (0.6 + 1.4 * eta + 2.2 * eta * eta)
    parameter = scale * (
        0.4
        + 0.5 * X
        + 1.2 * eta
        + 0.7 * X * X
        + 1.6 * X * eta
        + 2.7 * eta * eta
        + 2.2 * X * X * eta
        + 3.6 * X * eta * eta
        + 5.2 * eta**3
    )
    radial_parameter = scale * (
        0.5 + 1.4 * X + 1.6 * eta + 4.4 * X * eta + 3.6 * eta * eta
    )
    parameter2 = scale * (
        1.2
        + 1.6 * X
        + 5.4 * eta
        + 2.2 * X * X
        + 7.2 * X * eta
        + 15.6 * eta * eta
    )
    return RegularFluxFourthMixedJet(
        value=value,
        radial=radial,
        radial2=radial2,
        parameter=parameter,
        radial_parameter=radial_parameter,
        parameter2=parameter2,
        radial2_parameter=scale * (1.4 + 4.4 * eta),
        radial_parameter2=scale * (1.6 + 4.4 * X + 7.2 * eta),
        parameter3=scale * (5.4 + 7.2 * X + 31.2 * eta),
        radial2_parameter2=4.4 * scale,
        radial_parameter3=7.2 * scale,
        parameter4=31.2 * scale,
    )


def _axial_polynomial(scale: float, X: float, eta: float) -> ProfileSecondJet:
    beta = _fourth_polynomial(scale, X, eta)
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
        _fourth_polynomial(0.7 + 0.11 * j, X, eta) for j in range(order + 1)
    ]
    axial = [
        _axial_polynomial(1.1 - 0.07 * j, X, eta) for j in range(order + 1)
    ]
    return beta, axial


def _first_parameter(order: int, X: float, eta: float, h: float):
    beta, axial = _history(order, X, eta)
    return omega_over_x_parameter_jet_eq_5_6(
        order, X, eta, [jet.third() for jet in beta], axial, h=h
    )


@pytest.mark.parametrize("order", [0, 1, 2])
@pytest.mark.parametrize("X", [0.0, 0.37, 1.1])
@pytest.mark.parametrize("eta", [-0.31, 0.21])
def test_omega_second_parameter_jet_matches_landed_first_path_and_difference(
    order: int, X: float, eta: float
):
    h = 0.005
    beta, axial = _history(order, X, eta)
    got = omega_over_x_second_parameter_jet_eq_5_6(
        order, X, eta, beta, axial, h=h
    )
    first = omega_over_x_parameter_jet_eq_5_6(
        order, X, eta, [jet.third() for jet in beta], axial, h=h
    )

    assert got.value == first.value
    assert got.parameter == first.parameter

    eps = 1.0e-5
    finite_difference = (
        _first_parameter(order, X, eta + eps, h).parameter
        - _first_parameter(order, X, eta - eps, h).parameter
    ) / (2.0 * eps)
    assert got.parameter2 == pytest.approx(
        finite_difference, rel=4.0e-7, abs=5.0e-7
    )


def test_every_new_fourth_mixed_beta_entry_is_structurally_active():
    h = 0.005
    order, X, eta = 1, 0.43, 0.19
    beta, axial = _history(order, X, eta)
    reference = omega_over_x_second_parameter_jet_eq_5_6(
        order, X, eta, beta, axial, h=h
    ).parameter2

    perturbations = (
        (0, "radial2_parameter2", 0.37),
        (0, "radial_parameter3", -0.21),
        (0, "parameter4", 0.44),
    )
    for index, field, delta in perturbations:
        changed = list(beta)
        changed[index] = replace(
            changed[index], **{field: getattr(changed[index], field) + delta}
        )
        perturbed = omega_over_x_second_parameter_jet_eq_5_6(
            order, X, eta, changed, axial, h=h
        ).parameter2
        assert not np.isclose(reference, perturbed, rtol=0.0, atol=1.0e-8)


def test_current_order_radial2_parameter2_enters_radial_viscosity():
    h = 0.005
    order, X, eta = 1, 0.43, 0.19
    beta, axial = _history(order, X, eta)
    reference = omega_over_x_second_parameter_jet_eq_5_6(
        order, X, eta, beta, axial, h=h
    ).parameter2
    changed = list(beta)
    changed[order] = replace(
        changed[order],
        radial2_parameter2=changed[order].radial2_parameter2 + 0.53,
    )
    perturbed = omega_over_x_second_parameter_jet_eq_5_6(
        order, X, eta, changed, axial, h=h
    ).parameter2
    assert not np.isclose(reference, perturbed, rtol=0.0, atol=1.0e-8)


def test_omega_second_parameter_jet_fails_closed_on_bad_stronger_history():
    beta, axial = _history(1, 0.2, 0.1)
    with pytest.raises(ValueError, match="orders 0 through 1"):
        omega_over_x_second_parameter_jet_eq_5_6(
            1, 0.2, 0.1, beta[:1], axial, h=0.005
        )
    with pytest.raises(TypeError, match="RegularFluxFourthMixedJet"):
        omega_over_x_second_parameter_jet_eq_5_6(
            1, 0.2, 0.1, [object(), beta[1]], axial, h=0.005
        )
    with pytest.raises(ValueError, match="nonnegative integer"):
        omega_over_x_second_parameter_jet_eq_5_6(
            -1, 0.2, 0.1, beta, axial, h=0.005
        )
    with pytest.raises(ValueError, match="parameter4 must be finite"):
        replace(beta[0], parameter4=np.inf)
