import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_phi_fourth_mixed_jets import (
    Lemma52RepairedPhiFourthMixedJetAdapter,
)
from openai_ns_reconstruction.background_moment_repair_phi_third_mixed_jets import (
    Lemma52RepairedPhiThirdMixedJetAdapter,
)
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_preceding_diffusion_parameter_jet import (
    preceding_diffusion_parameter_jet,
)
from openai_ns_reconstruction.background_preceding_diffusion_second_parameter_jet import (
    ProfileFourthMixedJet,
    preceding_diffusion_second_parameter_jet,
)
from openai_ns_reconstruction.profiles import toy_gaussian_profile


def _repair():
    return Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=128,
    )


# Ascending eta powers. Degree four keeps the newly required fourth row
# nonzero while retaining an independent polynomial derivative oracle.
_MOMENT_POLYS = np.array(
    [
        [0.20, 0.10, 0.03, 0.004, 0.0007],
        [-0.15, 0.00, 0.02, -0.003, 0.0005],
        [0.11, -0.03, 0.00, 0.002, -0.0004],
        [0.25, 0.04, 0.01, -0.005, 0.0006],
        [-0.12, 0.05, 0.006, 0.001, -0.0003],
    ],
    dtype=float,
)


def _falling(power, derivative):
    if derivative > power:
        return 0
    out = 1
    for j in range(derivative):
        out *= power - j
    return out


def _poly_derivative(coefficients, eta, derivative):
    return math.fsum(
        float(coefficient)
        * _falling(power, derivative)
        * float(eta) ** (power - derivative)
        for power, coefficient in enumerate(coefficients)
        if power >= derivative
    )


def _moment_jets(eta):
    return np.vstack(
        [
            np.array(
                [_poly_derivative(row, eta, derivative) for row in _MOMENT_POLYS],
                dtype=float,
            )
            for derivative in range(5)
        ]
    )


def _patch_jets(eta):
    eta = float(eta)
    return np.array(
        [
            1.5 + 0.2 * eta + 0.03 * eta**2 + 0.01 * eta**3 + 0.002 * eta**4,
            0.2 + 0.06 * eta + 0.03 * eta**2 + 0.008 * eta**3,
            0.06 + 0.06 * eta + 0.024 * eta**2,
            0.06 + 0.048 * eta,
            0.048,
        ]
    )


def _base_phi_fourth_jet(X, eta):
    X = float(X)
    eta = float(eta)
    expx = math.exp(-X)
    coeff = [1.0, 0.04, 0.05, 0.01, 0.002]
    derivatives = [
        expx * _poly_derivative(coeff, eta, derivative)
        for derivative in range(5)
    ]
    value, parameter, parameter2, parameter3, parameter4 = derivatives
    return ProfileFourthMixedJet(
        value=value,
        radial=-value,
        radial2=value,
        parameter=parameter,
        radial_parameter=-parameter,
        parameter2=parameter2,
        radial2_parameter=parameter,
        radial_parameter2=-parameter2,
        parameter3=parameter3,
        radial2_parameter2=parameter2,
        radial_parameter3=-parameter3,
        parameter4=parameter4,
    )


def _functional(moment_provider=_moment_jets, patch_provider=_patch_jets):
    return Lemma52RepairedProfileAdapter(
        _repair(),
        toy_gaussian_profile(),
        moment_provider,
        patch_provider,
    )


def _fourth_bridge(functional=None, C=2.75):
    if functional is None:
        functional = _functional()
    return Lemma52RepairedPhiFourthMixedJetAdapter(
        functional,
        C,
        _base_phi_fourth_jet,
    )


def _third_bridge(functional=None, C=2.75):
    if functional is None:
        functional = _functional()
    return Lemma52RepairedPhiThirdMixedJetAdapter(
        functional,
        C,
        lambda X, eta: _base_phi_fourth_jet(X, eta).third(),
    )


def test_repaired_phi_fourth_mixed_jet_extends_existing_third_bridge():
    functional = _functional()
    fourth = _fourth_bridge(functional)
    third = _third_bridge(functional)
    X = 0.5 * 4.17**2
    eta = -0.23

    actual = fourth.phi_fourth_mixed_jet(X, eta)
    assert actual.third() == third.phi_third_mixed_jet(X, eta)


def test_repaired_phi_fourth_entries_match_eta_difference_of_third_bridge():
    functional = _functional()
    fourth = _fourth_bridge(functional)
    third = _third_bridge(functional)
    X = 0.5 * 4.17**2
    eta = 0.21
    eps = 2.0e-5

    actual = fourth.phi_fourth_mixed_jet(X, eta)
    plus = third.phi_third_mixed_jet(X, eta + eps)
    minus = third.phi_third_mixed_jet(X, eta - eps)
    observed = np.array(
        [
            actual.radial2_parameter2,
            actual.radial_parameter3,
            actual.parameter4,
        ]
    )
    oracle = np.array(
        [
            (plus.radial2_parameter - minus.radial2_parameter) / (2.0 * eps),
            (plus.radial_parameter2 - minus.radial_parameter2) / (2.0 * eps),
            (plus.parameter3 - minus.parameter3) / (2.0 * eps),
        ]
    )
    np.testing.assert_allclose(observed, oracle, rtol=8.0e-6, atol=8.0e-8)


def test_repaired_phi_fourth_jet_drives_angular_preceding_diffusion_second_eta():
    functional = _functional()
    fourth = _fourth_bridge(functional)
    third = _third_bridge(functional)
    h = 0.005
    power = -1.0 - h
    order = 2
    X = 0.5 * 4.17**2
    eta = -0.19
    eps = 2.0e-5

    actual = preceding_diffusion_second_parameter_jet(
        h,
        power,
        order,
        X,
        eta,
        fourth.phi_fourth_mixed_jet(X, eta),
    )
    plus = preceding_diffusion_parameter_jet(
        h,
        power,
        order,
        X,
        eta + eps,
        third.phi_third_mixed_jet(X, eta + eps),
    )
    minus = preceding_diffusion_parameter_jet(
        h,
        power,
        order,
        X,
        eta - eps,
        third.phi_third_mixed_jet(X, eta - eps),
    )
    oracle = (plus.parameter - minus.parameter) / (2.0 * eps)
    assert math.isclose(actual.parameter2, oracle, rel_tol=2.0e-5, abs_tol=2.0e-7)


def test_phi_fourth_mixed_bridge_is_compact_and_fails_closed_only_on_support():
    functional = _functional(
        lambda eta: _moment_jets(eta)[:1],
        lambda eta: _patch_jets(eta)[:1],
    )
    bridge = _fourth_bridge(functional)
    eta = 0.17

    for X in (0.0, 0.5 * 1.0**2):
        assert bridge.phi_fourth_mixed_jet(X, eta) == _base_phi_fourth_jet(X, eta)

    with pytest.raises(ValueError):
        bridge.phi_fourth_mixed_jet(0.5 * 4.17**2, eta)


def test_phi_fourth_mixed_bridge_rejects_invalid_normalization_and_provider():
    functional = _functional()
    with pytest.raises(ValueError):
        Lemma52RepairedPhiFourthMixedJetAdapter(
            functional,
            0.0,
            _base_phi_fourth_jet,
        )

    bridge = Lemma52RepairedPhiFourthMixedJetAdapter(
        functional,
        2.75,
        lambda X, eta: (X, eta),
    )
    with pytest.raises(TypeError, match="ProfileFourthMixedJet"):
        bridge.phi_fourth_mixed_jet(0.5 * 4.17**2, 0.1)
