import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_phi_jets import (
    Lemma52RepairedPhiSecondJetAdapter,
)
from openai_ns_reconstruction.background_moment_repair_phi_third_mixed_jets import (
    Lemma52RepairedPhiThirdMixedJetAdapter,
)
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_preceding_diffusion_parameter_jet import (
    ProfileThirdMixedJet,
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
            0.20 + 0.10 * eta + 0.010 * eta**3,
            -0.15 + 0.020 * eta**2,
            0.11 - 0.030 * eta + 0.005 * eta**3,
            0.25 + 0.040 * eta + 0.010 * eta**2 - 0.003 * eta**3,
            -0.12 + 0.050 * eta + 0.004 * eta**3,
        ]
    )
    first = np.array(
        [
            0.10 + 0.030 * eta**2,
            0.040 * eta,
            -0.030 + 0.015 * eta**2,
            0.040 + 0.020 * eta - 0.009 * eta**2,
            0.050 + 0.012 * eta**2,
        ]
    )
    second = np.array(
        [
            0.060 * eta,
            0.040,
            0.030 * eta,
            0.020 - 0.018 * eta,
            0.024 * eta,
        ]
    )
    third = np.array([0.060, 0.0, 0.030, -0.018, 0.024])
    return np.vstack((values, first, second, third))


def _patch_jets(eta):
    eta = float(eta)
    return np.array(
        [
            1.5 + 0.2 * eta + 0.03 * eta**2 + 0.01 * eta**3,
            0.2 + 0.06 * eta + 0.03 * eta**2,
            0.06 + 0.06 * eta,
            0.06,
        ]
    )


def _base_phi_third_jet(X, eta):
    X = float(X)
    eta = float(eta)
    expx = math.exp(-X)
    angular = 1.0 + 0.04 * eta + 0.05 * eta**2 + 0.01 * eta**3
    angular1 = 0.04 + 0.10 * eta + 0.03 * eta**2
    angular2 = 0.10 + 0.06 * eta
    angular3 = 0.06
    value = expx * angular
    parameter = expx * angular1
    parameter2 = expx * angular2
    parameter3 = expx * angular3
    return ProfileThirdMixedJet(
        value=value,
        radial=-value,
        radial2=value,
        parameter=parameter,
        radial_parameter=-parameter,
        parameter2=parameter2,
        radial2_parameter=parameter,
        radial_parameter2=-parameter2,
        parameter3=parameter3,
    )


def _functional(moment_provider=_moment_jets, patch_provider=_patch_jets):
    return Lemma52RepairedProfileAdapter(
        _repair(),
        toy_gaussian_profile(),
        moment_provider,
        patch_provider,
    )


def _bridge(C=2.75):
    return Lemma52RepairedPhiThirdMixedJetAdapter(
        _functional(), C, _base_phi_third_jet
    )


def _second_bridge(functional=None, C=2.75):
    if functional is None:
        functional = _functional()
    return Lemma52RepairedPhiSecondJetAdapter(
        functional,
        C,
        lambda X, eta: _base_phi_third_jet(X, eta).second(),
    )


def _first_eta(fn, eta, h=2.0e-5):
    return (fn(eta + h) - fn(eta - h)) / (2.0 * h)


def test_repaired_phi_third_mixed_jet_extends_existing_second_bridge():
    bridge = _bridge()
    second = _second_bridge(bridge.profile, bridge.C)
    X = 0.5 * 4.17**2
    eta = -0.23

    jet = bridge.phi_third_mixed_jet(X, eta)
    assert jet.second() == second.phi_second_jet(X, eta)


def test_repaired_phi_third_mixed_terms_match_eta_derivative_of_old_value_path():
    bridge = _bridge()
    second = _second_bridge(bridge.profile, bridge.C)
    X = 0.5 * 4.17**2
    eta = -0.23
    jet = bridge.phi_third_mixed_jet(X, eta)

    oracle_xxeta = _first_eta(
        lambda e: second.phi_second_jet(X, e).radial2, eta
    )
    oracle_xetaeta = _first_eta(
        lambda e: second.phi_second_jet(X, e).radial_parameter, eta
    )
    oracle_etaetaeta = _first_eta(
        lambda e: second.phi_second_jet(X, e).parameter2, eta
    )

    assert math.isclose(
        jet.radial2_parameter, oracle_xxeta, rel_tol=5e-6, abs_tol=5e-7
    )
    assert math.isclose(
        jet.radial_parameter2, oracle_xetaeta, rel_tol=5e-6, abs_tol=5e-7
    )
    assert math.isclose(
        jet.parameter3, oracle_etaetaeta, rel_tol=5e-6, abs_tol=5e-7
    )


def test_phi_third_mixed_bridge_is_compact_and_fails_closed_only_on_active_support():
    functional = _functional(
        lambda eta: _moment_jets(eta)[:1],
        lambda eta: _patch_jets(eta)[:1],
    )
    bridge = Lemma52RepairedPhiThirdMixedJetAdapter(
        functional, 2.75, _base_phi_third_jet
    )
    eta = 0.17

    for X in (0.0, 0.5 * 1.0**2):
        assert bridge.phi_third_mixed_jet(X, eta) == _base_phi_third_jet(X, eta)

    with pytest.raises(ValueError):
        bridge.phi_third_mixed_jet(0.5 * 4.17**2, eta)


def test_phi_third_mixed_bridge_rejects_invalid_normalization_and_provider():
    functional = _functional()
    with pytest.raises(ValueError):
        Lemma52RepairedPhiThirdMixedJetAdapter(
            functional, 0.0, _base_phi_third_jet
        )

    bridge = Lemma52RepairedPhiThirdMixedJetAdapter(
        functional, 2.75, lambda X, eta: (X, eta)
    )
    with pytest.raises(TypeError):
        bridge.phi_third_mixed_jet(0.5 * 4.17**2, 0.1)
