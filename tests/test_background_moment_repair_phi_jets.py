import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_lower_history_source import ProfileSecondJet
from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_phi_jets import (
    Lemma52RepairedPhiSecondJetAdapter,
)
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.profiles import toy_gaussian_profile


def _repair():
    return Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=128,
    )


def _moment_jets(eta):
    eta = float(eta)
    values = np.array(
        [
            0.20 + 0.10 * eta,
            -0.15 + 0.02 * eta**2,
            0.11 - 0.03 * eta,
            0.25 + 0.04 * eta + 0.01 * eta**2,
            -0.12 + 0.05 * eta,
        ]
    )
    first = np.array([0.10, 0.04 * eta, -0.03, 0.04 + 0.02 * eta, 0.05])
    second = np.array([0.0, 0.04, 0.0, 0.02, 0.0])
    return np.vstack((values, first, second))


def _patch_jets(eta):
    return np.array([1.5 + 0.2 * float(eta), 0.2, 0.0])


def _base_phi_jet(X, eta):
    X = float(X)
    eta = float(eta)
    expx = math.exp(-X)
    angular = 1.0 + 0.05 * eta * eta
    value = expx * angular
    parameter = expx * 0.10 * eta
    return ProfileSecondJet(
        value=value,
        radial=-value,
        radial2=value,
        parameter=parameter,
        radial_parameter=-parameter,
        parameter2=0.10 * expx,
    )


def _bridge(C=2.75):
    functional = Lemma52RepairedProfileAdapter(
        _repair(), toy_gaussian_profile(), _moment_jets, _patch_jets
    )
    return Lemma52RepairedPhiSecondJetAdapter(functional, C, _base_phi_jet)


def _phi_value(bridge, X, eta):
    base = _base_phi_jet(X, eta).value
    radius = math.sqrt(2.0 * X)
    if radius == 0.0:
        return base
    # Independent function-level oracle: the repaired-E increment is already
    # evaluable without using the production phi radial-derivative formulas.
    delta_E = bridge.profile.E(X, eta) - bridge.profile.base_profile.E(X, eta)
    return base + bridge.C * delta_E / radius


def _first_x(fn, X, eta, h=2.0e-5):
    return (fn(X + h, eta) - fn(X - h, eta)) / (2.0 * h)


def _second_x(fn, X, eta, h=2.0e-4):
    return (
        -fn(X + 2.0 * h, eta)
        + 16.0 * fn(X + h, eta)
        - 30.0 * fn(X, eta)
        + 16.0 * fn(X - h, eta)
        - fn(X - 2.0 * h, eta)
    ) / (12.0 * h * h)


def _first_eta(fn, X, eta, h=2.0e-5):
    return (fn(X, eta + h) - fn(X, eta - h)) / (2.0 * h)


def _second_eta(fn, X, eta, h=2.0e-4):
    return (
        -fn(X, eta + 2.0 * h)
        + 16.0 * fn(X, eta + h)
        - 30.0 * fn(X, eta)
        + 16.0 * fn(X, eta - h)
        - fn(X, eta - 2.0 * h)
    ) / (12.0 * h * h)


def _mixed(fn, X, eta, hx=8.0e-5, he=8.0e-5):
    return (
        fn(X + hx, eta + he)
        - fn(X + hx, eta - he)
        - fn(X - hx, eta + he)
        + fn(X - hx, eta - he)
    ) / (4.0 * hx * he)


def test_repaired_phi_second_jet_matches_function_level_eq_5_14_oracle():
    bridge = _bridge()
    radius = 4.17
    X = 0.5 * radius * radius
    eta = -0.23
    jet = bridge.phi_second_jet(X, eta)
    fn = lambda x, e: _phi_value(bridge, x, e)

    assert math.isclose(jet.value, fn(X, eta), rel_tol=0.0, abs_tol=2e-13)
    assert math.isclose(jet.radial, _first_x(fn, X, eta), rel_tol=3e-6, abs_tol=3e-7)
    assert math.isclose(jet.radial2, _second_x(fn, X, eta), rel_tol=4e-5, abs_tol=4e-5)
    assert math.isclose(
        jet.parameter, _first_eta(fn, X, eta), rel_tol=3e-6, abs_tol=3e-7
    )
    assert math.isclose(
        jet.parameter2, _second_eta(fn, X, eta), rel_tol=4e-5, abs_tol=4e-5
    )
    assert math.isclose(
        jet.radial_parameter, _mixed(fn, X, eta), rel_tol=5e-5, abs_tol=5e-5
    )


def test_phi_bridge_is_exactly_compact_and_axis_regular_without_division():
    functional = Lemma52RepairedProfileAdapter(
        _repair(),
        toy_gaussian_profile(),
        lambda eta: _moment_jets(eta)[:1],
        lambda eta: _patch_jets(eta)[:1],
    )
    bridge = Lemma52RepairedPhiSecondJetAdapter(functional, 2.75, _base_phi_jet)

    eta = 0.17
    for X in (0.0, 0.5 * 1.0**2):
        assert bridge.phi_second_jet(X, eta) == _base_phi_jet(X, eta)

    # On an active E-repair support the second eta data are genuinely required.
    with pytest.raises(ValueError):
        bridge.phi_second_jet(0.5 * 4.17**2, eta)


def test_phi_bridge_fails_closed_on_invalid_normalization_and_base_provider():
    functional = Lemma52RepairedProfileAdapter(
        _repair(), toy_gaussian_profile(), _moment_jets, _patch_jets
    )
    with pytest.raises(ValueError):
        Lemma52RepairedPhiSecondJetAdapter(functional, 0.0, _base_phi_jet)

    bridge = Lemma52RepairedPhiSecondJetAdapter(
        functional, 2.75, lambda X, eta: (X, eta)
    )
    with pytest.raises(TypeError):
        bridge.phi_second_jet(0.5 * 4.17**2, 0.1)
