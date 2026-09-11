import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_lower_history_source import ProfileSecondJet
from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_moment_repair_second_jets import (
    Lemma52RepairedSecondJetAdapter,
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
    first = np.array(
        [0.10, 0.04 * eta, -0.03, 0.04 + 0.02 * eta, 0.05]
    )
    second = np.array([0.0, 0.04, 0.0, 0.02, 0.0])
    return np.vstack((values, first, second))


def _patch_jets(eta):
    return np.array([1.5 + 0.2 * float(eta), 0.2, 0.0])


def _base_u_jet(X, eta):
    e = math.exp(-float(X))
    eta = float(eta)
    return ProfileSecondJet(
        value=eta * e,
        radial=-eta * e,
        radial2=eta * e,
        parameter=e,
        radial_parameter=-e,
        parameter2=0.0,
    )


def _base_e_jet(X, eta):
    X = float(X)
    eta = float(eta)
    radius = math.sqrt(2.0 * X)
    expx = math.exp(-X)
    angular = 1.0 + 0.1 * eta * eta
    value = radius * expx * angular
    log_x = 0.5 / X - 1.0
    parameter = radius * expx * 0.2 * eta
    return ProfileSecondJet(
        value=value,
        radial=value * log_x,
        radial2=value * (log_x * log_x - 0.5 / (X * X)),
        parameter=parameter,
        radial_parameter=parameter * log_x,
        parameter2=radius * expx * 0.2,
    )


def _bridge():
    functional = Lemma52RepairedProfileAdapter(
        _repair(),
        toy_gaussian_profile(),
        _moment_jets,
        _patch_jets,
    )
    return Lemma52RepairedSecondJetAdapter(
        functional,
        _base_u_jet,
        _base_e_jet,
    )


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


@pytest.mark.parametrize(
    "component,radius,eta",
    [
        ("U", 2.17, 0.23),
        ("E", 4.17, -0.29),
    ],
)
def test_repaired_second_jet_matches_independent_value_finite_differences(
    component, radius, eta
):
    bridge = _bridge()
    X = 0.5 * radius * radius
    if component == "U":
        jet = bridge.U_second_jet(X, eta)
        fn = bridge.profile.U
    else:
        jet = bridge.E_second_jet(X, eta)
        fn = bridge.profile.E

    # The oracle differentiates only the independently evaluable repaired
    # function values; it does not reuse the production bump derivative formulas.
    assert math.isclose(jet.value, fn(X, eta), rel_tol=0.0, abs_tol=2e-13)
    assert math.isclose(jet.radial, _first_x(fn, X, eta), rel_tol=2e-6, abs_tol=2e-7)
    assert math.isclose(jet.radial2, _second_x(fn, X, eta), rel_tol=2e-5, abs_tol=2e-5)
    assert math.isclose(
        jet.parameter, _first_eta(fn, X, eta), rel_tol=2e-6, abs_tol=2e-7
    )
    assert math.isclose(
        jet.parameter2, _second_eta(fn, X, eta), rel_tol=2e-5, abs_tol=2e-5
    )
    assert math.isclose(
        jet.radial_parameter, _mixed(fn, X, eta), rel_tol=3e-5, abs_tol=3e-5
    )


def test_second_jet_bridge_is_exactly_compact_and_skips_unneeded_eta_jets_off_patch():
    base = toy_gaussian_profile()
    repair = _repair()

    # Only zeroth-order moment data are available.  Off the compact repair
    # supports every correction derivative is mathematically zero, so the
    # bridge can return the base jet without pretending higher eta jets exist.
    functional = Lemma52RepairedProfileAdapter(
        repair,
        base,
        lambda eta: _moment_jets(eta)[:1],
        lambda eta: _patch_jets(eta)[:1],
    )
    bridge = Lemma52RepairedSecondJetAdapter(functional, _base_u_jet, _base_e_jet)

    X_off = 0.5 * 1.0**2
    eta = 0.17
    assert bridge.U_second_jet(X_off, eta) == _base_u_jet(X_off, eta)
    assert bridge.E_second_jet(X_off, eta) == _base_e_jet(X_off, eta)

    with pytest.raises(ValueError):
        bridge.U_second_jet(0.5 * 2.17**2, eta)
    with pytest.raises(ValueError):
        bridge.E_second_jet(0.5 * 4.17**2, eta)


def test_second_jet_bridge_fails_closed_on_non_jet_base_provider():
    functional = Lemma52RepairedProfileAdapter(
        _repair(), toy_gaussian_profile(), _moment_jets, _patch_jets
    )
    bridge = Lemma52RepairedSecondJetAdapter(
        functional,
        lambda X, eta: (X, eta),
        _base_e_jet,
    )
    with pytest.raises(TypeError):
        bridge.U_second_jet(0.5 * 2.17**2, 0.1)
