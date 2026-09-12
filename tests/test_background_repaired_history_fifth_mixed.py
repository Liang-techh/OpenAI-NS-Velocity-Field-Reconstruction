import math
from dataclasses import replace

import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_fifth_mixed_jets import (
    Lemma52RepairedFifthMixedJetAdapter,
)
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_omega_second_parameter_jet import (
    RegularFluxFourthMixedJet,
)
from openai_ns_reconstruction.background_preceding_diffusion_parameter_jet import (
    ProfileThirdMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_fourth_mixed_jets import (
    AxialFifthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history_fifth_mixed import (
    Section5FifthMixedCoefficientJetSource,
    Section5LowerHistoryFifthMixedHierarchy,
)
from openai_ns_reconstruction.background_repaired_history_omega_parameter import (
    hierarchy_owned_omega_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_omega_second_parameter import (
    hierarchy_owned_omega_second_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_phi_third_mixed import (
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


def _leading_strong_source():
    return Section5FifthMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _base_phi_third_mixed_jet(
            X, eta
        ).second(),
        axial_third_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth().third(),
        provenance=(
            "analytic fifth-mixed leading test fixture; Issue #1 remains upstream"
        ),
        axial_fourth_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth(),
        phi_third_mixed_jet_provider=_base_phi_third_mixed_jet,
        axial_fifth_mixed_jet_provider=_base_u_fifth_mixed_jet,
    )


def _hierarchy():
    repaired = (
        Section5FifthMixedCoefficientJetSource.from_lemma52_repair_full_fifth_mixed(
            1,
            C=C,
            profile=_functional_repair(),
            base_phi_third_mixed_jet=_base_phi_third_mixed_jet,
            base_U_fifth_mixed_jet=_base_u_fifth_mixed_jet,
        )
    )
    return Section5LowerHistoryFifthMixedHierarchy(
        H,
        C,
        (_leading_strong_source(), repaired),
        quadrature_points=64,
    )


def test_hierarchy_owns_repaired_fifth_u_and_exact_lower_projections():
    functional = _functional_repair()
    hierarchy = _hierarchy()
    X = 0.5 * 2.17**2
    eta = 0.23

    independent = Lemma52RepairedFifthMixedJetAdapter(
        functional,
        _base_u_fifth_mixed_jet,
    )
    expected_u = independent.U_fifth_mixed_jet(X, eta)
    actual_u = hierarchy.axial_fifth_mixed_jet(1, X, eta)

    assert actual_u == expected_u
    assert actual_u.fourth() == hierarchy.axial_fourth_mixed_jet(1, X, eta)
    assert actual_u.fourth().third() == hierarchy.axial_third_mixed_jet(1, X, eta)

    expected_beta = independent.beta_fourth_mixed_jet(
        H,
        1,
        X,
        eta,
        quadrature_points=64,
    )
    actual_beta = hierarchy.beta_fourth_mixed_jet(1, X, eta)
    assert isinstance(actual_beta, RegularFluxFourthMixedJet)
    assert actual_beta == expected_beta


@pytest.mark.parametrize("X", [0.0, 0.5 * 2.17**2])
def test_hierarchy_owned_omega_second_matches_eta_difference_of_owned_first_path(X):
    hierarchy = _hierarchy()
    eta = -0.19
    eps = 2.0e-5

    actual = hierarchy_owned_omega_second_parameter_jet(
        hierarchy,
        1,
        X,
        eta,
    )
    first = hierarchy_owned_omega_parameter_jet(hierarchy, 1, X, eta)
    plus = hierarchy_owned_omega_parameter_jet(hierarchy, 1, X, eta + eps)
    minus = hierarchy_owned_omega_parameter_jet(hierarchy, 1, X, eta - eps)

    assert actual.value == first.value
    assert actual.parameter == first.parameter
    oracle = (plus.parameter - minus.parameter) / (2.0 * eps)
    assert math.isclose(actual.parameter2, oracle, rel_tol=3.0e-5, abs_tol=3.0e-7)


def test_omega_second_fails_closed_when_leading_source_lacks_fifth_mixed_u():
    weak = Section5PhiThirdMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _base_phi_third_mixed_jet(
            X, eta
        ).second(),
        axial_third_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth().third(),
        provenance="weak leading fixture without fifth-mixed U",
        axial_fourth_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth(),
        phi_third_mixed_jet_provider=_base_phi_third_mixed_jet,
    )
    hierarchy = Section5LowerHistoryFifthMixedHierarchy(H, C, (weak,))

    with pytest.raises(ValueError, match="fifth-mixed U"):
        hierarchy_owned_omega_second_parameter_jet(hierarchy, 0, 0.0, 0.1)


def test_strong_hierarchy_rejects_incoherent_fifth_to_fourth_projection():
    def bad_fourth(X, eta):
        return replace(
            _base_u_fifth_mixed_jet(X, eta).fourth(),
            parameter4=1.0,
        )

    bad = Section5FifthMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _base_phi_third_mixed_jet(
            X, eta
        ).second(),
        axial_third_mixed_jet_provider=lambda X, eta: bad_fourth(X, eta).third(),
        provenance="deliberately incoherent fifth/fourth test fixture",
        axial_fourth_mixed_jet_provider=bad_fourth,
        phi_third_mixed_jet_provider=_base_phi_third_mixed_jet,
        axial_fifth_mixed_jet_provider=_base_u_fifth_mixed_jet,
    )
    hierarchy = Section5LowerHistoryFifthMixedHierarchy(H, C, (bad,))

    with pytest.raises(ValueError, match="incoherent"):
        hierarchy.axial_fourth_mixed_jet(0, 0.5 * 2.17**2, 0.2)


def test_lemma52_fifth_mixed_constructor_rejects_order_zero():
    with pytest.raises(ValueError, match="positive integer"):
        Section5FifthMixedCoefficientJetSource.from_lemma52_repair_full_fifth_mixed(
            0,
            C=C,
            profile=_functional_repair(),
            base_phi_third_mixed_jet=_base_phi_third_mixed_jet,
            base_U_fifth_mixed_jet=_base_u_fifth_mixed_jet,
        )
