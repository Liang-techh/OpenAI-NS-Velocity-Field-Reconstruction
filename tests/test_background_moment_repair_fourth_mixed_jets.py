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
from openai_ns_reconstruction.background_moment_repair_fourth_mixed_jets import (
    Lemma52RepairedFourthMixedJetAdapter,
)
from openai_ns_reconstruction.background_omega_parameter_jet import (
    RegularFluxThirdMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_third_mixed_jets import (
    AxialFourthMixedJet,
)
from openai_ns_reconstruction.profiles import toy_gaussian_profile


def _repair():
    return Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=128,
    )


# Coefficients are stored in ascending eta power.  The degree-four terms make
# the newly required fourth eta row nonzero while keeping a closed polynomial
# oracle for all lower derivative rows.
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


def _patch_jets(_eta):
    # Constant nonzero p(eta) isolates the moment-jet derivative propagation and
    # keeps the finite-difference oracle independent of quotient stencil noise.
    return np.array([1.5, 0.0, 0.0, 0.0, 0.0])


def _base_u_fourth_mixed_jet(X, eta):
    e = math.exp(-float(X))
    eta = float(eta)
    return AxialFourthMixedJet(
        value=eta * e,
        radial=-eta * e,
        radial2=eta * e,
        parameter=e,
        radial_parameter=-e,
        parameter2=0.0,
        radial2_parameter=e,
        radial_parameter2=0.0,
        parameter3=0.0,
        radial2_parameter2=0.0,
        radial_parameter3=0.0,
        parameter4=0.0,
    )


def _functional(moment_provider=_moment_jets, patch_provider=_patch_jets):
    return Lemma52RepairedProfileAdapter(
        _repair(),
        toy_gaussian_profile(),
        moment_provider,
        patch_provider,
    )


def _fourth_bridge(moment_provider=_moment_jets, patch_provider=_patch_jets):
    return Lemma52RepairedFourthMixedJetAdapter(
        _functional(moment_provider, patch_provider),
        _base_u_fourth_mixed_jet,
    )


def _third_bridge():
    functional = _functional()
    return Lemma52RepairedThirdMixedJetAdapter(
        functional,
        lambda X, eta: _base_u_fourth_mixed_jet(X, eta).third(),
    )


def test_repaired_fourth_mixed_u_jet_matches_independent_eta_difference_of_third_bridge():
    fourth = _fourth_bridge()
    third = _third_bridge()
    X = 0.5 * 2.17**2
    eta = 0.23
    eps = 2.0e-5

    actual = fourth.U_fourth_mixed_jet(X, eta)
    plus = third.U_third_mixed_jet(X, eta + eps)
    minus = third.U_third_mixed_jet(X, eta - eps)

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
    np.testing.assert_allclose(observed, oracle, rtol=5.0e-6, atol=5.0e-8)

    # The first nine fields are a strict lift of the already-landed repaired
    # third-mixed bridge rather than a competing implementation.
    center = third.U_third_mixed_jet(X, eta)
    np.testing.assert_allclose(
        np.array(
            [
                actual.value,
                actual.radial,
                actual.radial2,
                actual.parameter,
                actual.radial_parameter,
                actual.parameter2,
                actual.radial2_parameter,
                actual.radial_parameter2,
                actual.parameter3,
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
                center.radial2_parameter,
                center.radial_parameter2,
                center.parameter3,
            ]
        ),
        rtol=0.0,
        atol=0.0,
    )


def test_repaired_fourth_mixed_u_jet_feeds_eq_5_2_third_mixed_adapter():
    bridge = _fourth_bridge()
    h = 0.005
    order = 2
    X = 0.5 * 2.17**2
    eta = -0.21

    beta = bridge.beta_third_mixed_jet(h, order, X, eta, quadrature_points=64)
    assert isinstance(beta, RegularFluxThirdMixedJet)

    repaired_profile = bridge.profile.as_leading_profile()
    expected = repaired_profile.radial_flux_factor(
        X,
        eta,
        h,
        lam=2.0 * order * h,
        n=64,
    )
    assert math.isclose(beta.value, expected, rel_tol=2e-12, abs_tol=2e-12)


def test_fourth_mixed_bridge_is_compact_and_skips_unneeded_fourth_eta_row_off_patch():
    bridge = _fourth_bridge(
        lambda eta: _moment_jets(eta)[:1],
        lambda eta: _patch_jets(eta)[:1],
    )
    eta = 0.17
    X_off = 0.5 * 1.0**2
    assert bridge.U_fourth_mixed_jet(X_off, eta) == _base_u_fourth_mixed_jet(X_off, eta)

    with pytest.raises(ValueError):
        bridge.U_fourth_mixed_jet(0.5 * 2.17**2, eta)


def test_fourth_mixed_bridge_fails_closed_on_non_jet_base_provider():
    functional = _functional()
    bridge = Lemma52RepairedFourthMixedJetAdapter(
        functional,
        lambda X, eta: (X, eta),
    )
    with pytest.raises(TypeError, match="AxialFourthMixedJet"):
        bridge.U_fourth_mixed_jet(0.5 * 2.17**2, 0.1)
