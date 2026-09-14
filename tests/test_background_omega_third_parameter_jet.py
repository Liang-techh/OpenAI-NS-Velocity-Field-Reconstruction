import math

import pytest

from openai_ns_reconstruction.background_lower_history_source import ProfileSecondJet
from openai_ns_reconstruction.background_omega_second_parameter_jet import (
    omega_over_x_second_parameter_jet_eq_5_6,
)
from openai_ns_reconstruction.background_omega_third_parameter_jet import (
    omega_over_x_third_parameter_jet_eq_5_6,
)
from openai_ns_reconstruction.background_regular_flux_fifth_mixed_jets import (
    RegularFluxFifthMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_second_jets import AxialThirdMixedJet


def _falling(n: int, k: int) -> float:
    if k > n:
        return 0.0
    out = 1.0
    for j in range(k):
        out *= n - j
    return out


def _poly_derivative(coeff, X, eta, dx, de):
    total = 0.0
    for (px, pe), value in coeff.items():
        if px < dx or pe < de:
            continue
        total += (
            value
            * _falling(px, dx)
            * _falling(pe, de)
            * X ** (px - dx)
            * eta ** (pe - de)
        )
    return total


def _beta(coeff, X, eta):
    d = lambda dx, de: _poly_derivative(coeff, X, eta, dx, de)
    return RegularFluxFifthMixedJet(
        value=d(0, 0), radial=d(1, 0), radial2=d(2, 0),
        parameter=d(0, 1), radial_parameter=d(1, 1), parameter2=d(0, 2),
        radial2_parameter=d(2, 1), radial_parameter2=d(1, 2), parameter3=d(0, 3),
        radial2_parameter2=d(2, 2), radial_parameter3=d(1, 3), parameter4=d(0, 4),
        radial2_parameter3=d(2, 3), radial_parameter4=d(1, 4), parameter5=d(0, 5),
    )


def _axial(coeff, X, eta):
    d = lambda dx, de: _poly_derivative(coeff, X, eta, dx, de)
    return AxialThirdMixedJet(
        value=d(0, 0), radial=d(1, 0), radial2=d(2, 0),
        parameter=d(0, 1), radial_parameter=d(1, 1), parameter2=d(0, 2),
        radial2_parameter=d(2, 1), radial_parameter2=d(1, 2), parameter3=d(0, 3),
    )


def _second_axial(jet):
    return ProfileSecondJet(
        jet.value, jet.radial, jet.radial2,
        jet.parameter, jet.radial_parameter, jet.parameter2,
    )


BETA = (
    {(0, 0): 0.7, (1, 1): -0.4, (2, 0): 0.13, (0, 3): 0.2, (1, 4): -0.03, (2, 5): 0.002},
    {(0, 0): -0.25, (1, 0): 0.18, (0, 2): 0.31, (2, 2): -0.06, (1, 5): 0.004},
)
AXIAL = (
    {(0, 0): 0.4, (1, 0): -0.2, (0, 1): 0.3, (2, 1): 0.05, (1, 3): -0.02},
    {(0, 0): -0.15, (1, 1): 0.11, (0, 2): -0.08, (2, 0): 0.03, (0, 3): 0.01},
)


def _data(X, eta):
    return (
        [_beta(c, X, eta) for c in BETA],
        [_axial(c, X, eta) for c in AXIAL],
    )


def test_third_eta_omega_delegates_lower_rows_and_matches_test_only_difference():
    h = 0.23
    X = 0.41
    eta = 0.27
    beta, axial = _data(X, eta)
    result = omega_over_x_third_parameter_jet_eq_5_6(1, X, eta, beta, axial, h=h)
    lower = omega_over_x_second_parameter_jet_eq_5_6(
        1, X, eta, [jet.fourth() for jet in beta], [_second_axial(j) for j in axial], h=h
    )
    assert result.value == lower.value
    assert result.parameter == lower.parameter
    assert result.parameter2 == lower.parameter2

    # Finite differences are deliberately confined to this regression oracle.
    eps = 2.0e-5
    bp, ap = _data(X, eta + eps)
    bm, am = _data(X, eta - eps)
    plus = omega_over_x_second_parameter_jet_eq_5_6(
        1, X, eta + eps, [j.fourth() for j in bp], [_second_axial(j) for j in ap], h=h
    ).parameter2
    minus = omega_over_x_second_parameter_jet_eq_5_6(
        1, X, eta - eps, [j.fourth() for j in bm], [_second_axial(j) for j in am], h=h
    ).parameter2
    oracle = (plus - minus) / (2.0 * eps)
    assert result.parameter3 == pytest.approx(oracle, rel=3e-6, abs=3e-6)


def test_third_eta_omega_is_regular_on_axis():
    beta, axial = _data(0.0, -0.19)
    result = omega_over_x_third_parameter_jet_eq_5_6(1, 0.0, -0.19, beta, axial, h=0.23)
    assert math.isfinite(result.parameter3)
