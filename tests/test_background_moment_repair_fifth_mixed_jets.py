import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_fifth_mixed_jets import (
    Lemma52RepairedFifthMixedJetAdapter,
)
from openai_ns_reconstruction.background_moment_repair_fourth_mixed_jets import (
    Lemma52RepairedFourthMixedJetAdapter,
)
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_regular_flux_fourth_mixed_jets import (
    AxialFifthMixedJet,
)
from openai_ns_reconstruction.profiles import toy_gaussian_profile


def _repair():
    return Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=128,
    )


# Ascending eta powers.  Degree five makes the newly required fifth eta row
# nonzero while retaining an independent closed polynomial derivative oracle.
_MOMENT_POLYS = np.array(
    [
        [0.20, 0.10, 0.03, 0.004, 0.0007, 0.00011],
        [-0.15, 0.00, 0.02, -0.003, 0.0005, -0.00009],
        [0.11, -0.03, 0.00, 0.002, -0.0004, 0.00008],
        [0.25, 0.04, 0.01, -0.005, 0.0006, -0.00007],
        [-0.12, 0.05, 0.006, 0.001, -0.0003, 0.00006],
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
            for derivative in range(6)
        ]
    )


def _patch_jets(_eta):
    return np.array([1.5, 0.0, 0.0, 0.0, 0.0, 0.0])


def _base_u_fifth_mixed_jet(X, eta):
    e = math.exp(-float(X))
    eta = float(eta)
    return AxialFifthMixedJet(
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
        radial2_parameter3=0.0,
        radial_parameter4=0.0,
        parameter5=0.0,
    )


def _functional(moment_provider=_moment_jets, patch_provider=_patch_jets):
    return Lemma52RepairedProfileAdapter(
        _repair(),
        toy_gaussian_profile(),
        moment_provider,
        patch_provider,
    )


def _fifth_bridge(moment_provider=_moment_jets, patch_provider=_patch_jets):
    return Lemma52RepairedFifthMixedJetAdapter(
        _functional(moment_provider, patch_provider),
        _base_u_fifth_mixed_jet,
    )


def _fourth_bridge():
    functional = _functional()
    return Lemma52RepairedFourthMixedJetAdapter(
        functional,
        lambda X, eta: _base_u_fifth_mixed_jet(X, eta).fourth(),
    )


def test_repaired_fifth_mixed_u_jet_matches_eta_difference_of_fourth_bridge():
    fifth = _fifth_bridge()
    fourth = _fourth_bridge()
    X = 0.5 * 2.17**2
    eta = 0.23
    eps = 2.0e-5

    actual = fifth.U_fifth_mixed_jet(X, eta)
    plus = fourth.U_fourth_mixed_jet(X, eta + eps)
    minus = fourth.U_fourth_mixed_jet(X, eta - eps)
    oracle = np.array(
        [
            (plus.radial2_parameter2 - minus.radial2_parameter2) / (2.0 * eps),
            (plus.radial_parameter3 - minus.radial_parameter3) / (2.0 * eps),
            (plus.parameter4 - minus.parameter4) / (2.0 * eps),
        ]
    )
    observed = np.array(
        [
            actual.radial2_parameter3,
            actual.radial_parameter4,
            actual.parameter5,
        ]
    )
    np.testing.assert_allclose(observed, oracle, rtol=8.0e-6, atol=8.0e-8)

    # The lower projection is owned by the already-landed repaired fourth bridge.
    assert actual.fourth() == fourth.U_fourth_mixed_jet(X, eta)


def test_repaired_fifth_mixed_u_jet_drives_eq_5_2_fourth_mixed_beta():
    fifth = _fifth_bridge()
    fourth = _fourth_bridge()
    h = 0.005
    order = 2
    X = 0.5 * 2.17**2
    eta = -0.21
    eps = 2.0e-5

    actual = fifth.beta_fourth_mixed_jet(
        h, order, X, eta, quadrature_points=64
    )
    plus = fourth.beta_third_mixed_jet(
        h, order, X, eta + eps, quadrature_points=64
    )
    minus = fourth.beta_third_mixed_jet(
        h, order, X, eta - eps, quadrature_points=64
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
    np.testing.assert_allclose(observed, oracle, rtol=2.0e-5, atol=2.0e-7)


def test_fifth_mixed_bridge_is_compact_and_skips_unneeded_fifth_eta_row_off_patch():
    bridge = _fifth_bridge(
        lambda eta: _moment_jets(eta)[:1],
        lambda eta: _patch_jets(eta)[:1],
    )
    eta = 0.17
    X_off = 0.5 * 1.0**2
    assert bridge.U_fifth_mixed_jet(X_off, eta) == _base_u_fifth_mixed_jet(
        X_off, eta
    )

    with pytest.raises(ValueError):
        bridge.U_fifth_mixed_jet(0.5 * 2.17**2, eta)


def test_fifth_mixed_bridge_fails_closed_on_non_jet_base_provider():
    bridge = Lemma52RepairedFifthMixedJetAdapter(
        _functional(),
        lambda X, eta: (X, eta),
    )
    with pytest.raises(TypeError, match="AxialFifthMixedJet"):
        bridge.U_fifth_mixed_jet(0.5 * 2.17**2, 0.1)
