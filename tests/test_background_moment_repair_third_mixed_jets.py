import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_moment_repair_third_mixed_jets import (
    Lemma52RepairedThirdMixedJetAdapter,
)
from openai_ns_reconstruction.background_regular_flux_second_jets import AxialThirdMixedJet
from openai_ns_reconstruction.profiles import toy_gaussian_profile


def _repair():
    return Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=128,
    )


def _moment_jets(eta):
    e = float(eta)
    values = np.array(
        [
            0.20 + 0.10 * e + 0.03 * e**2 + 0.004 * e**3,
            -0.15 + 0.02 * e**2 - 0.003 * e**3,
            0.11 - 0.03 * e + 0.002 * e**3,
            0.25 + 0.04 * e + 0.01 * e**2 - 0.005 * e**3,
            -0.12 + 0.05 * e + 0.006 * e**2 + 0.001 * e**3,
        ]
    )
    first = np.array(
        [
            0.10 + 0.06 * e + 0.012 * e**2,
            0.04 * e - 0.009 * e**2,
            -0.03 + 0.006 * e**2,
            0.04 + 0.02 * e - 0.015 * e**2,
            0.05 + 0.012 * e + 0.003 * e**2,
        ]
    )
    second = np.array(
        [
            0.06 + 0.024 * e,
            0.04 - 0.018 * e,
            0.012 * e,
            0.02 - 0.03 * e,
            0.012 + 0.006 * e,
        ]
    )
    third = np.array([0.024, -0.018, 0.012, -0.03, 0.006])
    return np.vstack((values, first, second, third))


def _patch_jets(_eta):
    # A constant nonzero patch factor keeps the independent eta oracle simple
    # while still exercising genuine third derivatives of the repair coefficients.
    return np.array([1.5, 0.0, 0.0, 0.0])


def _base_u_third_mixed_jet(X, eta):
    e = math.exp(-float(X))
    eta = float(eta)
    return AxialThirdMixedJet(
        value=eta * e,
        radial=-eta * e,
        radial2=eta * e,
        parameter=e,
        radial_parameter=-e,
        parameter2=0.0,
        radial2_parameter=e,
        radial_parameter2=0.0,
        parameter3=0.0,
    )


def _bridge(moment_provider=_moment_jets, patch_provider=_patch_jets):
    functional = Lemma52RepairedProfileAdapter(
        _repair(),
        toy_gaussian_profile(),
        moment_provider,
        patch_provider,
    )
    return Lemma52RepairedThirdMixedJetAdapter(
        functional,
        _base_u_third_mixed_jet,
    )


def _second_x(fn, X, eta, h=4.0e-4):
    return (
        -fn(X + 2.0 * h, eta)
        + 16.0 * fn(X + h, eta)
        - 30.0 * fn(X, eta)
        + 16.0 * fn(X - h, eta)
        - fn(X - 2.0 * h, eta)
    ) / (12.0 * h * h)


def _second_eta(fn, X, eta, h=4.0e-4):
    return (
        -fn(X, eta + 2.0 * h)
        + 16.0 * fn(X, eta + h)
        - 30.0 * fn(X, eta)
        + 16.0 * fn(X, eta - h)
        - fn(X, eta - 2.0 * h)
    ) / (12.0 * h * h)


def _radial2_parameter(fn, X, eta, he=8.0e-4):
    return (_second_x(fn, X, eta + he) - _second_x(fn, X, eta - he)) / (2.0 * he)


def _radial_parameter2(fn, X, eta, hx=8.0e-4):
    return (_second_eta(fn, X + hx, eta) - _second_eta(fn, X - hx, eta)) / (2.0 * hx)


def _parameter3(fn, X, eta, h=2.0e-3):
    return (
        fn(X, eta + 2.0 * h)
        - 2.0 * fn(X, eta + h)
        + 2.0 * fn(X, eta - h)
        - fn(X, eta - 2.0 * h)
    ) / (2.0 * h**3)


def test_repaired_third_mixed_u_jet_matches_independent_value_finite_differences():
    bridge = _bridge()
    X = 0.5 * 2.17**2
    eta = 0.23
    jet = bridge.U_third_mixed_jet(X, eta)
    fn = bridge.profile.U

    assert math.isclose(jet.value, fn(X, eta), rel_tol=0.0, abs_tol=2e-13)
    assert math.isclose(
        jet.radial2_parameter,
        _radial2_parameter(fn, X, eta),
        rel_tol=2e-3,
        abs_tol=2e-3,
    )
    assert math.isclose(
        jet.radial_parameter2,
        _radial_parameter2(fn, X, eta),
        rel_tol=2e-3,
        abs_tol=2e-3,
    )
    assert math.isclose(
        jet.parameter3,
        _parameter3(fn, X, eta),
        rel_tol=2e-3,
        abs_tol=2e-3,
    )


def test_repaired_third_mixed_u_jet_feeds_eq_5_2_value_consistently():
    bridge = _bridge()
    h = 0.005
    order = 2
    X = 0.5 * 2.17**2
    eta = -0.21

    beta = bridge.beta_second_jet(h, order, X, eta, quadrature_points=64)
    repaired_profile = bridge.profile.as_leading_profile()
    expected = repaired_profile.radial_flux_factor(
        X,
        eta,
        h,
        lam=2.0 * order * h,
        n=64,
    )
    assert math.isclose(beta.value, expected, rel_tol=2e-12, abs_tol=2e-12)


def test_third_mixed_bridge_is_compact_and_skips_unneeded_eta_jets_off_patch():
    bridge = _bridge(
        lambda eta: _moment_jets(eta)[:1],
        lambda eta: _patch_jets(eta)[:1],
    )
    eta = 0.17
    X_off = 0.5 * 1.0**2
    assert bridge.U_third_mixed_jet(X_off, eta) == _base_u_third_mixed_jet(X_off, eta)

    with pytest.raises(ValueError):
        bridge.U_third_mixed_jet(0.5 * 2.17**2, eta)


def test_third_mixed_bridge_fails_closed_on_non_jet_base_provider():
    functional = Lemma52RepairedProfileAdapter(
        _repair(), toy_gaussian_profile(), _moment_jets, _patch_jets
    )
    bridge = Lemma52RepairedThirdMixedJetAdapter(
        functional,
        lambda X, eta: (X, eta),
    )
    with pytest.raises(TypeError):
        bridge.U_third_mixed_jet(0.5 * 2.17**2, 0.1)
