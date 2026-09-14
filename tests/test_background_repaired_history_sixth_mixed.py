import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_fifth_mixed_jets import Lemma52RepairedFifthMixedJetAdapter
from openai_ns_reconstruction.background_moment_repair_phi_fifth_mixed_jets import ProfileFifthMixedJet
from openai_ns_reconstruction.background_moment_repair_profile import Lemma52RepairedProfileAdapter
from openai_ns_reconstruction.background_moment_repair_sixth_mixed_jets import Lemma52RepairedSixthMixedJetAdapter
from openai_ns_reconstruction.background_regular_flux_fifth_mixed_jets import AxialSixthMixedJet, RegularFluxFifthMixedJet
from openai_ns_reconstruction.background_repaired_history_phi_fifth_mixed import Section5PhiFifthMixedCoefficientJetSource
from openai_ns_reconstruction.background_repaired_history_sixth_mixed import Section5LowerHistorySixthMixedHierarchy, Section5SixthMixedCoefficientJetSource
from openai_ns_reconstruction.profiles import toy_gaussian_profile

C = 2.75
H = 0.005


def _repair():
    return Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=128,
    )


_MOMENT_POLYS = np.array([
    [0.20, 0.10, 0.03, 0.004, 0.0007, 0.00011, 0.000013],
    [-0.15, 0.00, 0.02, -0.003, 0.0005, -0.00009, 0.000011],
    [0.11, -0.03, 0.00, 0.002, -0.0004, 0.00008, -0.000010],
    [0.25, 0.04, 0.01, -0.005, 0.0006, -0.00007, 0.000009],
    [-0.12, 0.05, 0.006, 0.001, -0.0003, 0.00006, -0.000008],
], dtype=float)
_U_POLY = np.array([0.2, -0.15, 0.04, 0.01, -0.003, 0.0005, -0.00007])


def _falling(power, derivative):
    if derivative > power:
        return 0
    out = 1
    for j in range(derivative):
        out *= power - j
    return out


def _poly_derivative(coefficients, eta, derivative):
    return math.fsum(
        float(coefficient) * _falling(power, derivative) * float(eta) ** (power - derivative)
        for power, coefficient in enumerate(coefficients)
        if power >= derivative
    )


def _moment_jets(eta):
    return np.vstack([
        np.array([_poly_derivative(row, eta, derivative) for row in _MOMENT_POLYS], dtype=float)
        for derivative in range(7)
    ])


def _patch_jets(_eta):
    return np.array([1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])


def _functional(moment_provider=_moment_jets, patch_provider=_patch_jets):
    return Lemma52RepairedProfileAdapter(
        _repair(), toy_gaussian_profile(), moment_provider, patch_provider,
    )


def _base_u_sixth_mixed_jet(X, eta):
    e = math.exp(-float(X))
    q = [_poly_derivative(_U_POLY, eta, derivative) * e for derivative in range(7)]
    return AxialSixthMixedJet(
        value=q[0], radial=-q[0], radial2=q[0],
        parameter=q[1], radial_parameter=-q[1], parameter2=q[2],
        radial2_parameter=q[1], radial_parameter2=-q[2], parameter3=q[3],
        radial2_parameter2=q[2], radial_parameter3=-q[3], parameter4=q[4],
        radial2_parameter3=q[3], radial_parameter4=-q[4], parameter5=q[5],
        radial2_parameter4=q[4], radial_parameter5=-q[5], parameter6=q[6],
    )


def _base_phi_fifth_mixed_jet(_X, _eta):
    return ProfileFifthMixedJet(
        value=0.0, radial=0.0, radial2=0.0, parameter=0.0,
        radial_parameter=0.0, parameter2=0.0, radial2_parameter=0.0,
        radial_parameter2=0.0, parameter3=0.0, radial2_parameter2=0.0,
        radial_parameter3=0.0, parameter4=0.0, radial2_parameter3=0.0,
        radial_parameter4=0.0, parameter5=0.0,
    )


def _sixth_bridge():
    return Lemma52RepairedSixthMixedJetAdapter(_functional(), _base_u_sixth_mixed_jet)


def _fifth_bridge():
    return Lemma52RepairedFifthMixedJetAdapter(
        _functional(), lambda X, eta: _base_u_sixth_mixed_jet(X, eta).fifth(),
    )


def _leading_source():
    return Section5SixthMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _base_phi_fifth_mixed_jet(X, eta).fourth().third().second(),
        axial_third_mixed_jet_provider=lambda X, eta: _base_u_sixth_mixed_jet(X, eta).fifth().fourth().third(),
        provenance="analytic sixth-mixed leading fixture; Issue #1 remains upstream",
        axial_fourth_mixed_jet_provider=lambda X, eta: _base_u_sixth_mixed_jet(X, eta).fifth().fourth(),
        phi_third_mixed_jet_provider=lambda X, eta: _base_phi_fifth_mixed_jet(X, eta).fourth().third(),
        axial_fifth_mixed_jet_provider=lambda X, eta: _base_u_sixth_mixed_jet(X, eta).fifth(),
        phi_fourth_mixed_jet_provider=lambda X, eta: _base_phi_fifth_mixed_jet(X, eta).fourth(),
        phi_fifth_mixed_jet_provider=_base_phi_fifth_mixed_jet,
        axial_sixth_mixed_jet_provider=_base_u_sixth_mixed_jet,
    )


def _hierarchy():
    repaired = Section5SixthMixedCoefficientJetSource.from_lemma52_repair_full_sixth_mixed(
        1, C=C, profile=_functional(),
        base_phi_fifth_mixed_jet=_base_phi_fifth_mixed_jet,
        base_U_sixth_mixed_jet=_base_u_sixth_mixed_jet,
    )
    return Section5LowerHistorySixthMixedHierarchy(
        H, C, (_leading_source(), repaired), quadrature_points=64,
    )


def test_compact_repaired_sixth_u_matches_eta_difference_of_fifth_layer():
    sixth = _sixth_bridge()
    fifth = _fifth_bridge()
    X = 0.5 * 2.17**2
    eta = 0.23
    eps = 2.0e-5
    actual = sixth.U_sixth_mixed_jet(X, eta)
    plus = fifth.U_fifth_mixed_jet(X, eta + eps)
    minus = fifth.U_fifth_mixed_jet(X, eta - eps)
    oracle = np.array([
        (plus.radial2_parameter3 - minus.radial2_parameter3) / (2.0 * eps),
        (plus.radial_parameter4 - minus.radial_parameter4) / (2.0 * eps),
        (plus.parameter5 - minus.parameter5) / (2.0 * eps),
    ])
    observed = np.array([
        actual.radial2_parameter4, actual.radial_parameter5, actual.parameter6,
    ])
    np.testing.assert_allclose(observed, oracle, rtol=1.0e-5, atol=1.0e-7)
    assert actual.fifth() == fifth.U_fifth_mixed_jet(X, eta)


def test_repaired_sixth_u_drives_eq_5_2_fifth_mixed_beta():
    sixth = _sixth_bridge()
    fifth = _fifth_bridge()
    X = 0.5 * 2.17**2
    eta = -0.21
    eps = 2.0e-5
    actual = sixth.beta_fifth_mixed_jet(H, 1, X, eta, quadrature_points=64)
    plus = fifth.beta_fourth_mixed_jet(H, 1, X, eta + eps, quadrature_points=64)
    minus = fifth.beta_fourth_mixed_jet(H, 1, X, eta - eps, quadrature_points=64)
    oracle = np.array([
        (plus.radial2_parameter2 - minus.radial2_parameter2) / (2.0 * eps),
        (plus.radial_parameter3 - minus.radial_parameter3) / (2.0 * eps),
        (plus.parameter4 - minus.parameter4) / (2.0 * eps),
    ])
    observed = np.array([
        actual.radial2_parameter3, actual.radial_parameter4, actual.parameter5,
    ])
    np.testing.assert_allclose(observed, oracle, rtol=3.0e-5, atol=3.0e-7)


def test_strong_hierarchy_owns_sixth_u_and_derives_fifth_beta():
    hierarchy = _hierarchy()
    independent = _sixth_bridge()
    X = 0.5 * 2.17**2
    eta = 0.19
    actual_u = hierarchy.axial_sixth_mixed_jet(1, X, eta)
    expected_u = independent.U_sixth_mixed_jet(X, eta)
    assert actual_u == expected_u
    assert actual_u.fifth() == hierarchy.axial_fifth_mixed_jet(1, X, eta)
    actual_beta = hierarchy.beta_fifth_mixed_jet(1, X, eta)
    expected_beta = independent.beta_fifth_mixed_jet(H, 1, X, eta, quadrature_points=64)
    assert isinstance(actual_beta, RegularFluxFifthMixedJet)
    assert actual_beta == expected_beta
    assert actual_beta.fourth() == hierarchy.beta_fourth_mixed_jet(1, X, eta)


def test_sixth_hierarchy_fails_closed_on_old_phi_fifth_source():
    weak = Section5PhiFifthMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _base_phi_fifth_mixed_jet(X, eta).fourth().third().second(),
        axial_third_mixed_jet_provider=lambda X, eta: _base_u_sixth_mixed_jet(X, eta).fifth().fourth().third(),
        provenance="old strong fixture without sixth-mixed U",
        axial_fourth_mixed_jet_provider=lambda X, eta: _base_u_sixth_mixed_jet(X, eta).fifth().fourth(),
        phi_third_mixed_jet_provider=lambda X, eta: _base_phi_fifth_mixed_jet(X, eta).fourth().third(),
        axial_fifth_mixed_jet_provider=lambda X, eta: _base_u_sixth_mixed_jet(X, eta).fifth(),
        phi_fourth_mixed_jet_provider=lambda X, eta: _base_phi_fifth_mixed_jet(X, eta).fourth(),
        phi_fifth_mixed_jet_provider=_base_phi_fifth_mixed_jet,
    )
    hierarchy = Section5LowerHistorySixthMixedHierarchy(H, C, (weak,))
    with pytest.raises(ValueError, match="sixth-mixed U"):
        hierarchy.axial_sixth_mixed_jet(0, 0.5 * 2.17**2, 0.1)


def test_sixth_mixed_bridge_is_compact_off_patch():
    bridge = Lemma52RepairedSixthMixedJetAdapter(
        _functional(lambda eta: _moment_jets(eta)[:1], lambda eta: _patch_jets(eta)[:1]),
        _base_u_sixth_mixed_jet,
    )
    eta = 0.17
    X_off = 0.5 * 1.0**2
    assert bridge.U_sixth_mixed_jet(X_off, eta) == _base_u_sixth_mixed_jet(X_off, eta)
    with pytest.raises(ValueError):
        bridge.U_sixth_mixed_jet(0.5 * 2.17**2, eta)
