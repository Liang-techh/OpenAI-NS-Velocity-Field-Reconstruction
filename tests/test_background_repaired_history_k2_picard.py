import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_inner_solver import picard_map_eq_5_7
from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_preceding_diffusion_second_parameter_jet import (
    ProfileFourthMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_fourth_mixed_jets import (
    AxialFifthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history_forcing_second_parameter import (
    hierarchy_owned_positive_axis_forcing_second_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_k1_picard_parameter import (
    hierarchy_owned_k1_picard_parameter_jet_eq_5_7,
)
from openai_ns_reconstruction.background_repaired_history_k2_picard import (
    hierarchy_owned_k2_picard_value_eq_5_7,
)
from openai_ns_reconstruction.background_repaired_history_phi_fourth_mixed import (
    Section5LowerHistoryPhiFourthMixedHierarchy,
    Section5PhiFourthMixedCoefficientJetSource,
)
from openai_ns_reconstruction.profiles import toy_gaussian_profile


C = 2.75
H = 0.005
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
    value = 1
    for j in range(derivative):
        value *= power - j
    return value if derivative <= power else 0


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
                [_poly_derivative(row, eta, derivative) for row in _MOMENT_POLYS]
            )
            for derivative in range(6)
        ]
    )


def _base_phi_fourth_mixed_jet(X, eta):
    e = math.exp(-float(X))
    coefficients = [1.0, 0.04, 0.05, 0.01, 0.002]
    value, parameter, parameter2, parameter3, parameter4 = [
        e * _poly_derivative(coefficients, eta, derivative)
        for derivative in range(5)
    ]
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
    return Section5PhiFourthMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _base_phi_fourth_mixed_jet(
            X, eta
        ).third().second(),
        axial_third_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth().third(),
        provenance="analytic strong leading test fixture; Issue #1 remains upstream",
        axial_fourth_mixed_jet_provider=lambda X, eta: _base_u_fifth_mixed_jet(
            X, eta
        ).fourth(),
        phi_third_mixed_jet_provider=lambda X, eta: _base_phi_fourth_mixed_jet(
            X, eta
        ).third(),
        axial_fifth_mixed_jet_provider=_base_u_fifth_mixed_jet,
        phi_fourth_mixed_jet_provider=_base_phi_fourth_mixed_jet,
    )


def _hierarchy():
    repair = Lemma52MomentRepair.from_intervals(
        0.2,
        ((2.0, 2.35), (3.0, 3.35)),
        ((4.0, 4.35), (5.0, 5.35), (6.0, 6.35)),
        quadrature_points=64,
    )
    profile = Lemma52RepairedProfileAdapter(
        repair,
        toy_gaussian_profile(),
        _moment_jets,
        lambda _eta: np.array([1.5, 0.0, 0.0, 0.0, 0.0, 0.0]),
    )
    repaired = Section5PhiFourthMixedCoefficientJetSource.from_lemma52_repair_full_strong_mixed(
        1,
        C=C,
        profile=profile,
        base_phi_fourth_mixed_jet=_base_phi_fourth_mixed_jet,
        base_U_fifth_mixed_jet=_base_u_fifth_mixed_jet,
    )
    return Section5LowerHistoryPhiFourthMixedHierarchy(
        H,
        C,
        (_leading_source(), repaired),
        quadrature_points=32,
    )


@pytest.mark.parametrize("xi", [0.0, math.sqrt(0.5 * 2.17**2)])
def test_k2_picard_matches_direct_second_application_on_compact_repaired_history(xi):
    hierarchy = _hierarchy()
    order = 2
    eta = 0.23
    quadrature_points = 8
    fields = hierarchy.positive_axis_fields(order)

    def previous(s, eta_value):
        return hierarchy_owned_k1_picard_parameter_jet_eq_5_7(
            hierarchy,
            order,
            s,
            eta_value,
            quadrature_points=quadrature_points,
        ).value

    def previous_parameter(s, eta_value):
        return hierarchy_owned_k1_picard_parameter_jet_eq_5_7(
            hierarchy,
            order,
            s,
            eta_value,
            quadrature_points=quadrature_points,
        ).parameter

    def forcing(s, eta_value):
        return hierarchy_owned_positive_axis_forcing_second_parameter_jet(
            hierarchy,
            order,
            s,
            eta_value,
        ).value

    expected = picard_map_eq_5_7(
        order,
        xi,
        eta,
        previous,
        previous_parameter,
        fields.A0,
        fields.A1,
        forcing,
        quadrature_points=quadrature_points,
    )
    actual = hierarchy_owned_k2_picard_value_eq_5_7(
        hierarchy,
        order,
        xi,
        eta,
        quadrature_points=quadrature_points,
    )

    assert actual == pytest.approx(expected, rel=2.0e-12, abs=2.0e-12)
    if xi == 0.0:
        assert np.array_equal(actual, np.zeros(6))
    else:
        previous_value = previous(xi, eta)
        assert np.linalg.norm(actual - previous_value) > 1.0e-10


def test_k2_picard_requires_strong_phi_fourth_hierarchy():
    with pytest.raises(TypeError, match="Section5LowerHistoryPhiFourthMixedHierarchy"):
        hierarchy_owned_k2_picard_value_eq_5_7(object(), 1, 0.5, 0.1)
