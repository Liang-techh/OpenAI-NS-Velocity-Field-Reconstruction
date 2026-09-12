import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_omega_second_parameter_jet import (
    RegularFluxFourthMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_fourth_mixed_jets import (
    AxialFifthMixedJet,
    regular_flux_fourth_mixed_jet_eq_5_2,
)
from openai_ns_reconstruction.background_regular_flux_third_mixed_jets import (
    AxialFourthMixedJet,
    regular_flux_third_mixed_jet_eq_5_2,
)


# Independent X-degree-two / eta-degree-five polynomial.  The final three
# monomials make U_XXetaetaeta, U_Xetaetaetaeta and U_eta^5 genuinely nonzero.
_MONOMIALS = (
    (1.10, 0, 0),
    (0.20, 1, 0),
    (-0.05, 2, 0),
    (0.17, 0, 1),
    (0.08, 1, 1),
    (-0.04, 2, 1),
    (0.03, 0, 2),
    (0.025, 1, 2),
    (0.020, 2, 2),
    (-0.015, 0, 3),
    (0.012, 1, 3),
    (0.009, 2, 3),
    (0.008, 0, 4),
    (-0.006, 1, 4),
    (0.004, 0, 5),
)


def _falling(power: int, derivative: int) -> int:
    if derivative > power:
        return 0
    out = 1
    for j in range(derivative):
        out *= power - j
    return out


def _poly_derivative(X: float, eta: float, dx: int, de: int, scale: float) -> float:
    total = 0.0
    for coefficient, px, pe in _MONOMIALS:
        fx = _falling(px, dx)
        fe = _falling(pe, de)
        if fx == 0 or fe == 0:
            continue
        total += (
            scale
            * coefficient
            * fx
            * fe
            * X ** (px - dx)
            * eta ** (pe - de)
        )
    return total


def _fifth_provider(scale: float):
    def provider(X: float, eta: float) -> AxialFifthMixedJet:
        d = lambda dx, de: _poly_derivative(X, eta, dx, de, scale)
        return AxialFifthMixedJet(
            value=d(0, 0),
            radial=d(1, 0),
            radial2=d(2, 0),
            parameter=d(0, 1),
            radial_parameter=d(1, 1),
            parameter2=d(0, 2),
            radial2_parameter=d(2, 1),
            radial_parameter2=d(1, 2),
            parameter3=d(0, 3),
            radial2_parameter2=d(2, 2),
            radial_parameter3=d(1, 3),
            parameter4=d(0, 4),
            radial2_parameter3=d(2, 3),
            radial_parameter4=d(1, 4),
            parameter5=d(0, 5),
        )

    return provider


def _fourth_provider(provider):
    return lambda X, eta: provider(X, eta).fourth()


@pytest.mark.parametrize("order", [0, 1, 2])
@pytest.mark.parametrize("X", [0.0, 0.37, 1.1])
def test_fourth_mixed_beta_matches_independent_eta_difference_of_landed_third_jet(
    order, X
):
    h = 0.005
    eta = 0.21
    eps = 2.0e-5
    provider = _fifth_provider(1.0 + 0.17 * order)

    actual = regular_flux_fourth_mixed_jet_eq_5_2(
        h, order, X, eta, provider, quadrature_points=48
    )

    plus = regular_flux_third_mixed_jet_eq_5_2(
        h,
        order,
        X,
        eta + eps,
        _fourth_provider(provider),
        quadrature_points=48,
    )
    minus = regular_flux_third_mixed_jet_eq_5_2(
        h,
        order,
        X,
        eta - eps,
        _fourth_provider(provider),
        quadrature_points=48,
    )

    oracle = np.array(
        [
            (plus.radial2_parameter - minus.radial2_parameter) / (2.0 * eps),
            (plus.radial_parameter2 - minus.radial_parameter2) / (2.0 * eps),
            (plus.parameter3 - minus.parameter3) / (2.0 * eps),
        ]
    )
    observed = np.array(
        [
            actual.radial2_parameter2,
            actual.radial_parameter3,
            actual.parameter4,
        ]
    )
    np.testing.assert_allclose(observed, oracle, rtol=4.0e-6, atol=4.0e-8)

    # Every lower field remains exactly the landed third-mixed Eq. (5.2) jet.
    center = regular_flux_third_mixed_jet_eq_5_2(
        h, order, X, eta, _fourth_provider(provider), quadrature_points=48
    )
    assert actual.third() == center


def test_all_three_new_fifth_u_inputs_are_structurally_active():
    provider = _fifth_provider(1.0)
    h, order, X, eta = 0.005, 2, 0.53, 0.18
    baseline = provider(X, eta)

    def varied(field: str):
        def inner(x: float, e: float) -> AxialFifthMixedJet:
            jet = provider(x, e)
            values = jet.__dict__.copy()
            values[field] += 0.25
            return AxialFifthMixedJet(**values)

        return inner

    base = regular_flux_fourth_mixed_jet_eq_5_2(h, order, X, eta, provider)
    assert baseline.radial2_parameter3 != 0.0
    assert baseline.radial_parameter4 != 0.0
    assert baseline.parameter5 != 0.0

    xxeee = regular_flux_fourth_mixed_jet_eq_5_2(
        h, order, X, eta, varied("radial2_parameter3")
    )
    xeeee = regular_flux_fourth_mixed_jet_eq_5_2(
        h, order, X, eta, varied("radial_parameter4")
    )
    eeeee = regular_flux_fourth_mixed_jet_eq_5_2(
        h, order, X, eta, varied("parameter5")
    )
    assert not math.isclose(xxeee.radial2_parameter2, base.radial2_parameter2)
    assert not math.isclose(xeeee.radial_parameter3, base.radial_parameter3)
    assert not math.isclose(eeeee.parameter4, base.parameter4)


def test_fail_closed_on_missing_fifth_mixed_u_data_and_bad_domain():
    provider = _fifth_provider(1.0)

    def only_fourth(X, eta):
        return AxialFourthMixedJet(**provider(X, eta).fourth().__dict__)

    with pytest.raises(TypeError, match="AxialFifthMixedJet"):
        regular_flux_fourth_mixed_jet_eq_5_2(0.005, 1, 0.2, 0.1, only_fourth)
    with pytest.raises(ValueError, match="order"):
        regular_flux_fourth_mixed_jet_eq_5_2(0.005, -1, 0.2, 0.1, provider)
    with pytest.raises(ValueError, match="nonnegative"):
        regular_flux_fourth_mixed_jet_eq_5_2(0.005, 1, -0.1, 0.1, provider)
    with pytest.raises(ValueError, match="eta"):
        regular_flux_fourth_mixed_jet_eq_5_2(0.005, 1, 0.2, 1.1, provider)

    values = provider(0.2, 0.1).__dict__.copy()
    values["parameter5"] = float("nan")
    with pytest.raises(ValueError, match="parameter5"):
        AxialFifthMixedJet(**values)


def test_result_is_the_eq_56_second_parameter_jet_input_type():
    result = regular_flux_fourth_mixed_jet_eq_5_2(
        0.005, 1, 0.4, 0.2, _fifth_provider(1.0)
    )
    assert isinstance(result, RegularFluxFourthMixedJet)
