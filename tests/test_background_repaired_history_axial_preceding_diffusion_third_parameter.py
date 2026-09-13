import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_preceding_diffusion_parameter_jet import (
    ProfileThirdMixedJet,
)
from openai_ns_reconstruction.background_preceding_diffusion_second_parameter_jet import (
    ProfileFourthMixedJet,
    preceding_diffusion_second_parameter_jet,
)
from openai_ns_reconstruction.background_regular_flux_fourth_mixed_jets import (
    AxialFifthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history_axial_preceding_diffusion_third_parameter import (
    hierarchy_owned_axial_preceding_diffusion_third_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_fifth_mixed import (
    Section5FifthMixedCoefficientJetSource,
    Section5LowerHistoryFifthMixedHierarchy,
)
from openai_ns_reconstruction.profiles import toy_gaussian_profile


H = 0.005
C = 2.75

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
    repair = Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=128,
    )
    return Lemma52RepairedProfileAdapter(
        repair,
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


def _leading_source():
    return Section5FifthMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _base_phi_third_mixed_jet(
            X, eta
        ).second(),
        axial_third_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth().third(),
        provenance="analytic leading fifth-mixed test fixture; Issue #1 remains upstream",
        axial_fourth_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth(),
        phi_third_mixed_jet_provider=_base_phi_third_mixed_jet,
        axial_fifth_mixed_jet_provider=_base_u_fifth_mixed_jet,
    )


def _hierarchy():
    repaired = Section5FifthMixedCoefficientJetSource.from_lemma52_repair_full_fifth_mixed(
        1,
        C=C,
        profile=_functional_repair(),
        base_phi_third_mixed_jet=_base_phi_third_mixed_jet,
        base_U_fifth_mixed_jet=_base_u_fifth_mixed_jet,
    )
    return Section5LowerHistoryFifthMixedHierarchy(
        H,
        C,
        (_leading_source(), repaired),
        quadrature_points=64,
    )


def _profile_fourth_from_axial(jet):
    fourth = jet.fourth()
    return ProfileFourthMixedJet(
        value=fourth.value,
        radial=fourth.radial,
        radial2=fourth.radial2,
        parameter=fourth.parameter,
        radial_parameter=fourth.radial_parameter,
        parameter2=fourth.parameter2,
        radial2_parameter=fourth.radial2_parameter,
        radial_parameter2=fourth.radial_parameter2,
        parameter3=fourth.parameter3,
        radial2_parameter2=fourth.radial2_parameter2,
        radial_parameter3=fourth.radial_parameter3,
        parameter4=fourth.parameter4,
    )


def test_axial_bridge_consumes_real_compact_repaired_u_fifth_jet():
    hierarchy = _hierarchy()
    X = 0.5 * 2.17**2  # inside the first active U-repair bump
    eta = -0.23
    order = 2

    repaired = hierarchy.axial_fifth_mixed_jet(1, X, eta)
    assert repaired != _base_u_fifth_mixed_jet(X, eta)

    actual = hierarchy_owned_axial_preceding_diffusion_third_parameter_jet(
        hierarchy,
        order,
        X,
        eta,
    )
    lower = preceding_diffusion_second_parameter_jet(
        H,
        -0.5 - H,
        order,
        X,
        eta,
        _profile_fourth_from_axial(repaired),
    )
    assert actual.value == lower.value
    assert actual.parameter == lower.parameter
    assert actual.parameter2 == lower.parameter2

    epsilon = 2.0e-5
    plus = preceding_diffusion_second_parameter_jet(
        H,
        -0.5 - H,
        order,
        X,
        eta + epsilon,
        _profile_fourth_from_axial(
            hierarchy.axial_fifth_mixed_jet(1, X, eta + epsilon)
        ),
    ).parameter2
    minus = preceding_diffusion_second_parameter_jet(
        H,
        -0.5 - H,
        order,
        X,
        eta - epsilon,
        _profile_fourth_from_axial(
            hierarchy.axial_fifth_mixed_jet(1, X, eta - epsilon)
        ),
    ).parameter2
    oracle = (plus - minus) / (2.0 * epsilon)
    assert actual.parameter3 == pytest.approx(oracle, rel=2.0e-5, abs=2.0e-7)


def test_axial_bridge_fails_closed_when_strict_lower_order_is_missing():
    hierarchy = Section5LowerHistoryFifthMixedHierarchy(
        H,
        C,
        (_leading_source(),),
    )
    with pytest.raises(ValueError, match="missing one or more strict-lower"):
        hierarchy_owned_axial_preceding_diffusion_third_parameter_jet(
            hierarchy,
            2,
            0.0,
            0.1,
        )


def test_axial_bridge_requires_fifth_mixed_hierarchy():
    with pytest.raises(TypeError, match="Section5LowerHistoryFifthMixedHierarchy"):
        hierarchy_owned_axial_preceding_diffusion_third_parameter_jet(
            object(),
            1,
            0.0,
            0.0,
        )
