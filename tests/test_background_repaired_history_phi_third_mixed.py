import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_phi_third_mixed_jets import (
    Lemma52RepairedPhiThirdMixedJetAdapter,
)
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_preceding_diffusion_parameter_jet import (
    ProfileThirdMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_third_mixed_jets import (
    AxialFourthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history import (
    Section5CoefficientJetSource,
)
from openai_ns_reconstruction.background_repaired_history_phi_third_mixed import (
    Section5LowerHistoryPhiThirdMixedHierarchy,
    Section5PhiThirdMixedCoefficientJetSource,
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


def _base_phi_third_mixed_jet(X, eta):
    X = float(X)
    eta = float(eta)
    e = math.exp(-X)
    p = 1.0 + 0.1 * eta + 0.2 * eta**2 + 0.03 * eta**3
    p1 = 0.1 + 0.4 * eta + 0.09 * eta**2
    p2 = 0.4 + 0.18 * eta
    p3 = 0.18
    return ProfileThirdMixedJet(
        value=C * e * p,
        radial=-C * e * p,
        radial2=C * e * p,
        parameter=C * e * p1,
        radial_parameter=-C * e * p1,
        parameter2=C * e * p2,
        radial2_parameter=C * e * p1,
        radial_parameter2=-C * e * p2,
        parameter3=C * e * p3,
    )


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


def _leading_strong_source():
    return Section5PhiThirdMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _base_phi_third_mixed_jet(X, eta).second(),
        axial_third_mixed_jet_provider=lambda X, eta: _base_u_fourth_mixed_jet(X, eta).third(),
        provenance="analytic leading-order test fixture; Issue #1 remains upstream",
        phi_third_mixed_jet_provider=_base_phi_third_mixed_jet,
    )


def test_hierarchy_owns_actual_repaired_phi_third_mixed_and_projects_second_jet():
    functional = _functional_repair()
    repaired = Section5PhiThirdMixedCoefficientJetSource.from_lemma52_repair_full_mixed(
        1,
        C=C,
        profile=functional,
        base_phi_third_mixed_jet=_base_phi_third_mixed_jet,
        base_U_fourth_mixed_jet=_base_u_fourth_mixed_jet,
    )
    hierarchy = Section5LowerHistoryPhiThirdMixedHierarchy(
        H,
        C,
        (_leading_strong_source(), repaired),
        quadrature_points=64,
    )
    X = 0.5 * 2.17**2
    eta = 0.23

    independent = Lemma52RepairedPhiThirdMixedJetAdapter(
        functional,
        C,
        _base_phi_third_mixed_jet,
    )
    expected = independent.phi_third_mixed_jet(X, eta)
    actual = hierarchy.phi_third_mixed_jet(1, X, eta)

    assert actual == expected
    assert hierarchy.phi_second_jet(1, X, eta) == expected.second()
    assert hierarchy.axial_fourth_mixed_jet(1, X, eta) == pytest.approx(
        _independent_repaired_u(functional, X, eta)
    )


def _independent_repaired_u(functional, X, eta):
    from openai_ns_reconstruction.background_moment_repair_fourth_mixed_jets import (
        Lemma52RepairedFourthMixedJetAdapter,
    )

    return Lemma52RepairedFourthMixedJetAdapter(
        functional,
        _base_u_fourth_mixed_jet,
    ).U_fourth_mixed_jet(X, eta)


def test_phi_third_mixed_path_fails_closed_for_old_fourth_mixed_source():
    old = Section5CoefficientJetSource.from_lemma52_repair_fourth_mixed(
        1,
        C=C,
        profile=_functional_repair(),
        base_phi_second_jet=lambda X, eta: _base_phi_third_mixed_jet(X, eta).second(),
        base_U_fourth_mixed_jet=_base_u_fourth_mixed_jet,
    )
    hierarchy = Section5LowerHistoryPhiThirdMixedHierarchy(
        H,
        C,
        (_leading_strong_source(), old),
    )

    with pytest.raises(ValueError, match="third-mixed phi"):
        hierarchy.phi_third_mixed_jet(1, 0.5 * 2.17**2, 0.1)


def test_hierarchy_rejects_incoherent_third_to_second_phi_projection():
    zero_second = _base_phi_third_mixed_jet(0.0, 0.0).second()
    zero_second = type(zero_second)(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    bad = Section5PhiThirdMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda _X, _eta: zero_second,
        axial_third_mixed_jet_provider=lambda X, eta: _base_u_fourth_mixed_jet(X, eta).third(),
        provenance="deliberately incoherent test fixture",
        phi_third_mixed_jet_provider=_base_phi_third_mixed_jet,
    )
    hierarchy = Section5LowerHistoryPhiThirdMixedHierarchy(H, C, (bad,))

    with pytest.raises(ValueError, match="incoherent"):
        hierarchy.phi_third_mixed_jet(0, 0.5 * 2.17**2, 0.2)


def test_lemma52_full_mixed_constructor_rejects_order_zero():
    with pytest.raises(ValueError, match="positive integer"):
        Section5PhiThirdMixedCoefficientJetSource.from_lemma52_repair_full_mixed(
            0,
            C=C,
            profile=_functional_repair(),
            base_phi_third_mixed_jet=_base_phi_third_mixed_jet,
            base_U_fourth_mixed_jet=_base_u_fourth_mixed_jet,
        )
