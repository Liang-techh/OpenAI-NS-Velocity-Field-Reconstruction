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
from openai_ns_reconstruction.background_moment_repair_third_mixed_jets import (
    Lemma52RepairedThirdMixedJetAdapter,
)
from openai_ns_reconstruction.background_regular_flux_second_jets import AxialThirdMixedJet
from openai_ns_reconstruction.background_repaired_history import (
    Section5CoefficientJetSource,
    Section5LowerHistoryJetHierarchy,
)
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
    return np.array([1.5, 0.0, 0.0, 0.0])


def _base_phi_second_jet(X, eta):
    X = float(X)
    eta = float(eta)
    e = math.exp(-X)
    angular = 1.0 + 0.1 * eta * eta
    value = C * e * angular
    parameter = C * e * 0.2 * eta
    return ProfileSecondJet(
        value=value,
        radial=-value,
        radial2=value,
        parameter=parameter,
        radial_parameter=-parameter,
        parameter2=C * e * 0.2,
    )


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


def _functional_repair():
    return Lemma52RepairedProfileAdapter(
        _repair(),
        toy_gaussian_profile(),
        _moment_jets,
        _patch_jets,
    )


def _hierarchy():
    leading = Section5CoefficientJetSource(
        order=0,
        phi_second_jet_provider=_base_phi_second_jet,
        axial_third_mixed_jet_provider=_base_u_third_mixed_jet,
        provenance="analytic leading-order test fixture; Issue #1 remains upstream",
    )
    repaired = Section5CoefficientJetSource.from_lemma52_repair(
        1,
        C=C,
        profile=_functional_repair(),
        base_phi_second_jet=_base_phi_second_jet,
        base_U_third_mixed_jet=_base_u_third_mixed_jet,
    )
    return Section5LowerHistoryJetHierarchy(
        H,
        C,
        (leading, repaired),
        quadrature_points=64,
    )


def test_positive_order_source_owns_actual_compact_repaired_phi_and_u_jets():
    functional = _functional_repair()
    hierarchy = _hierarchy()

    X_phi = 0.5 * 4.17**2
    eta = 0.23
    expected_phi = Lemma52RepairedPhiSecondJetAdapter(
        functional, C, _base_phi_second_jet
    ).phi_second_jet(X_phi, eta)
    assert hierarchy.phi_second_jet(1, X_phi, eta) == expected_phi

    X_u = 0.5 * 2.17**2
    expected_u = Lemma52RepairedThirdMixedJetAdapter(
        functional, _base_u_third_mixed_jet
    ).U_third_mixed_jet(X_u, eta)
    actual_u = hierarchy.axial_third_mixed_jet(1, X_u, eta)
    assert actual_u == expected_u
    projected = hierarchy.axial_second_jet(1, X_u, eta)
    assert projected.value == expected_u.value
    assert projected.radial2 == expected_u.radial2
    assert projected.parameter2 == expected_u.parameter2


def test_hierarchy_derives_beta_from_owned_repaired_u_not_an_independent_provider():
    hierarchy = _hierarchy()
    functional = _functional_repair()
    X = 0.5 * 2.17**2
    eta = -0.21

    beta = hierarchy.beta_second_jet(1, X, eta)
    repaired_profile = functional.as_leading_profile()
    expected = repaired_profile.radial_flux_factor(
        X,
        eta,
        H,
        lam=2.0 * H,
        n=64,
    )
    assert math.isclose(beta.value, expected, rel_tol=2e-12, abs_tol=2e-12)


def test_contiguous_owned_history_feeds_actual_lower_source_for_next_order():
    hierarchy = _hierarchy()
    providers = hierarchy.positive_axis_providers(2)
    xi = 1.13
    eta = 0.17

    base = providers.base(xi, eta)
    source = providers.source(xi, eta)
    fields = hierarchy.positive_axis_fields(2)

    assert math.isfinite(base.phi.value)
    assert math.isfinite(base.axial.value)
    assert math.isfinite(base.beta)
    assert all(
        math.isfinite(value)
        for value in (
            source.angular,
            source.axial,
            source.pressure_product,
            source.omega_quotient,
        )
    )
    assert fields.A0(xi, eta).shape == (6, 6)
    assert fields.A1(xi, eta).shape == (6, 6)
    assert fields.forcing(xi, eta).shape == (6,)
    # order 2 is deliberately absent: the recursive bridge must consume only
    # the owned strict lower history 0,1 rather than current-order data.
    with pytest.raises(ValueError):
        hierarchy.phi_second_jet(2, xi * xi, eta)


def test_hierarchy_fails_closed_on_gaps_and_mismatched_repair_normalization():
    repaired = Section5CoefficientJetSource.from_lemma52_repair(
        1,
        C=C,
        profile=_functional_repair(),
        base_phi_second_jet=_base_phi_second_jet,
        base_U_third_mixed_jet=_base_u_third_mixed_jet,
    )
    with pytest.raises(ValueError):
        Section5LowerHistoryJetHierarchy(H, C, (repaired,))
    leading = Section5CoefficientJetSource(
        0,
        _base_phi_second_jet,
        _base_u_third_mixed_jet,
        "leading fixture",
    )
    with pytest.raises(ValueError):
        Section5LowerHistoryJetHierarchy(H, C + 0.25, (leading, repaired))
    with pytest.raises(ValueError):
        Section5CoefficientJetSource.from_lemma52_repair(
            0,
            C=C,
            profile=_functional_repair(),
            base_phi_second_jet=_base_phi_second_jet,
            base_U_third_mixed_jet=_base_u_third_mixed_jet,
        )
