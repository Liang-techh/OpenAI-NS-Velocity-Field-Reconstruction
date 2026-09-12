import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_lower_history_source import ProfileSecondJet
from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_regular_flux_third_mixed_jets import (
    AxialFourthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history import (
    Section5CoefficientJetSource,
    Section5LowerHistoryJetHierarchy,
)
from openai_ns_reconstruction.background_repaired_history_omega_parameter import (
    hierarchy_owned_omega_parameter_jet,
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
    return np.array([1.5, 0.0, 0.0, 0.0, 0.0])


def _functional_repair():
    return Lemma52RepairedProfileAdapter(
        _repair(),
        toy_gaussian_profile(),
        _moment_jets,
        _patch_jets,
    )


def _phi_second_jet(X, eta):
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


def _u_fourth_mixed_jet(X, eta):
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


def _leading(*, own_fourth=True):
    return Section5CoefficientJetSource(
        order=0,
        phi_second_jet_provider=_phi_second_jet,
        axial_third_mixed_jet_provider=lambda X, eta: _u_fourth_mixed_jet(X, eta).third(),
        axial_fourth_mixed_jet_provider=_u_fourth_mixed_jet if own_fourth else None,
        provenance="analytic leading-order test fixture; Issue #1 remains upstream",
    )


def _hierarchy(*, leading_owns_fourth=True):
    repaired = Section5CoefficientJetSource.from_lemma52_repair_fourth_mixed(
        1,
        C=C,
        profile=_functional_repair(),
        base_phi_second_jet=_phi_second_jet,
        base_U_fourth_mixed_jet=_u_fourth_mixed_jet,
    )
    return Section5LowerHistoryJetHierarchy(
        H,
        C,
        (_leading(own_fourth=leading_owns_fourth), repaired),
        quadrature_points=64,
    )


def _existing_omega_value(hierarchy, omega_order, X, eta):
    providers = hierarchy.positive_axis_providers(omega_order + 1)
    return providers.source(math.sqrt(X), eta).omega_quotient


def test_owned_omega_parameter_jet_matches_existing_value_and_independent_eta_difference():
    hierarchy = _hierarchy()
    X = 0.5 * 2.17**2
    eta = 0.19

    actual = hierarchy_owned_omega_parameter_jet(hierarchy, 1, X, eta)
    existing = _existing_omega_value(hierarchy, 1, X, eta)
    assert math.isclose(actual.value, existing, rel_tol=2e-12, abs_tol=2e-12)

    step = 2e-6
    finite_difference = (
        _existing_omega_value(hierarchy, 1, X, eta + step)
        - _existing_omega_value(hierarchy, 1, X, eta - step)
    ) / (2.0 * step)
    assert math.isclose(
        actual.parameter,
        finite_difference,
        rel_tol=3e-6,
        abs_tol=3e-7,
    )


def test_owned_omega_parameter_jet_is_regular_on_axis():
    hierarchy = _hierarchy()
    actual = hierarchy_owned_omega_parameter_jet(hierarchy, 1, 0.0, -0.23)
    existing = _existing_omega_value(hierarchy, 1, 0.0, -0.23)
    assert math.isfinite(actual.value)
    assert math.isfinite(actual.parameter)
    assert math.isclose(actual.value, existing, rel_tol=2e-12, abs_tol=2e-12)


def test_owned_omega_parameter_jet_fails_closed_without_leading_fourth_mixed_data():
    hierarchy = _hierarchy(leading_owns_fourth=False)
    with pytest.raises(ValueError, match="coefficient order 0 does not own fourth-mixed U data"):
        hierarchy_owned_omega_parameter_jet(hierarchy, 1, 0.5 * 2.17**2, 0.1)


def test_owned_omega_parameter_jet_rejects_missing_or_invalid_order():
    hierarchy = _hierarchy()
    with pytest.raises(ValueError, match="missing one or more required coefficient sources"):
        hierarchy_owned_omega_parameter_jet(hierarchy, 2, 1.0, 0.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        hierarchy_owned_omega_parameter_jet(hierarchy, -1, 1.0, 0.0)
