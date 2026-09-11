import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_omega_parameter_jet import (
    RegularFluxThirdMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_second_jets import (
    AxialThirdMixedJet,
    regular_flux_second_jet_eq_5_2,
)
from openai_ns_reconstruction.background_regular_flux_third_mixed_jets import (
    AxialFourthMixedJet,
    regular_flux_third_mixed_jet_eq_5_2,
)


# A polynomial with independent X^2 eta^2, X eta^3 and eta^4 terms, so each
# newly required fourth-mixed U derivative is genuinely nonzero.
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
    (0.008, 0, 4),
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


def _fourth_provider(scale: float):
    def provider(X: float, eta: float) -> AxialFourthMixedJet:
        d = lambda dx, de: _poly_derivative(X, eta, dx, de, scale)
        return AxialFourthMixedJet(
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
        )

    return provider


def _third_provider(provider):
    return lambda X, eta: provider(X, eta).third()


@pytest.mark.parametrize("order", [0, 1, 2])
@pytest.mark.parametrize("X", [0.0, 0.37, 1.1])
def test_third_mixed_beta_matches_independent_eta_difference_of_landed_second_jet(
    order, X
):
    h = 0.005
    eta = 0.21
    eps = 2.0e-5
    provider = _fourth_provider(1.0 + 0.17 * order)

    actual = regular_flux_third_mixed_jet_eq_5_2(
        h, order, X, eta, provider, quadrature_points=48
    )

    plus = regular_flux_second_jet_eq_5_2(
        h,
        order,
        X,
        eta + eps,
        _third_provider(provider),
        quadrature_points=48,
    )
    minus = regular_flux_second_jet_eq_5_2(
        h,
        order,
        X,
        eta - eps,
        _third_provider(provider),
        quadrature_points=48,
    )

    oracle = np.array(
        [
            (plus.radial2 - minus.radial2) / (2.0 * eps),
            (plus.radial_parameter - minus.radial_parameter) / (2.0 * eps),
            (plus.parameter2 - minus.parameter2) / (2.0 * eps),
        ]
    )
    observed = np.array(
        [
            actual.radial2_parameter,
            actual.radial_parameter2,
            actual.parameter3,
        ]
    )
    np.testing.assert_allclose(observed, oracle, rtol=3.0e-6, atol=3.0e-8)

    # The first six fields remain exactly the already-landed Eq. (5.2) second
    # jet, so the new adapter is a strict extension rather than a competing
    # implementation of lower derivative orders.
    center = regular_flux_second_jet_eq_5_2(
        h, order, X, eta, _third_provider(provider), quadrature_points=48
    )
    np.testing.assert_allclose(
        np.array(
            [
                actual.value,
                actual.radial,
                actual.radial2,
                actual.parameter,
                actual.radial_parameter,
                actual.parameter2,
            ]
        ),
        np.array(
            [
                center.value,
                center.radial,
                center.radial2,
                center.parameter,
                center.radial_parameter,
                center.parameter2,
            ]
        ),
        rtol=0.0,
        atol=0.0,
    )


def test_all_three_new_fourth_u_inputs_are_structurally_active():
    provider = _fourth_provider(1.0)
    h, order, X, eta = 0.005, 2, 0.53, 0.18
    baseline = provider(X, eta)

    def varied(field: str):
        def inner(x: float, e: float) -> AxialFourthMixedJet:
            jet = provider(x, e)
            values = jet.__dict__.copy()
            values[field] += 0.25
            return AxialFourthMixedJet(**values)

        return inner

    base = regular_flux_third_mixed_jet_eq_5_2(h, order, X, eta, provider)
    assert baseline.radial2_parameter2 != 0.0
    assert baseline.radial_parameter3 != 0.0
    assert baseline.parameter4 != 0.0

    xxee = regular_flux_third_mixed_jet_eq_5_2(
        h, order, X, eta, varied("radial2_parameter2")
    )
    xeee = regular_flux_third_mixed_jet_eq_5_2(
        h, order, X, eta, varied("radial_parameter3")
    )
    eeee = regular_flux_third_mixed_jet_eq_5_2(
        h, order, X, eta, varied("parameter4")
    )
    assert not math.isclose(xxee.radial2_parameter, base.radial2_parameter)
    assert not math.isclose(xeee.radial_parameter2, base.radial_parameter2)
    assert not math.isclose(eeee.parameter3, base.parameter3)


def test_fail_closed_on_missing_fourth_mixed_u_data_and_bad_domain():
    provider = _fourth_provider(1.0)

    def only_third(X, eta):
        return AxialThirdMixedJet(**provider(X, eta).third().__dict__)

    with pytest.raises(TypeError, match="AxialFourthMixedJet"):
        regular_flux_third_mixed_jet_eq_5_2(0.005, 1, 0.2, 0.1, only_third)
    with pytest.raises(ValueError, match="order"):
        regular_flux_third_mixed_jet_eq_5_2(0.005, -1, 0.2, 0.1, provider)
    with pytest.raises(ValueError, match="nonnegative"):
        regular_flux_third_mixed_jet_eq_5_2(0.005, 1, -0.1, 0.1, provider)
    with pytest.raises(ValueError, match="eta"):
        regular_flux_third_mixed_jet_eq_5_2(0.005, 1, 0.2, 1.1, provider)

    with pytest.raises(ValueError, match="parameter4"):
        AxialFourthMixedJet(
            value=0.0,
            radial=0.0,
            radial2=0.0,
            parameter=0.0,
            radial_parameter=0.0,
            parameter2=0.0,
            radial2_parameter=0.0,
            radial_parameter2=0.0,
            parameter3=0.0,
            radial2_parameter2=0.0,
            radial_parameter3=0.0,
            parameter4=float("nan"),
        )


def test_result_is_the_eq_56_parameter_jet_input_type():
    result = regular_flux_third_mixed_jet_eq_5_2(
        0.005, 1, 0.4, 0.2, _fourth_provider(1.0)
    )
    assert isinstance(result, RegularFluxThirdMixedJet)
