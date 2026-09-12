import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_inner_solver import (
    first_picard_term_eq_5_7,
    singular_inverse_eq_5_7,
)
from openai_ns_reconstruction.background_lower_history_source import (
    positive_axis_point_data_from_lower_history,
)
from openai_ns_reconstruction.background_moment_repair import Lemma52MomentRepair
from openai_ns_reconstruction.background_moment_repair_profile import (
    Lemma52RepairedProfileAdapter,
)
from openai_ns_reconstruction.background_positive_axis import positive_axis_forcing
from openai_ns_reconstruction.background_preceding_diffusion_parameter_jet import (
    ProfileThirdMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_third_mixed_jets import (
    AxialFourthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history_first_picard_parameter import (
    hierarchy_owned_first_picard_parameter_jet_eq_5_7,
)
from openai_ns_reconstruction.background_repaired_history_forcing_parameter import (
    hierarchy_owned_positive_axis_forcing_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_k1_picard import (
    hierarchy_owned_k1_picard_value_eq_5_7,
)
from openai_ns_reconstruction.background_repaired_history_phi_third_mixed import (
    Section5LowerHistoryPhiThirdMixedHierarchy,
    Section5PhiThirdMixedCoefficientJetSource,
)
from openai_ns_reconstruction.background_repaired_history_source_parameter import (
    hierarchy_owned_lower_source_parameter_jet,
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


def _leading_source(*, with_fourth=True):
    kwargs = {}
    if with_fourth:
        kwargs["axial_fourth_mixed_jet_provider"] = _base_u_fourth_mixed_jet
    return Section5PhiThirdMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _base_phi_third_mixed_jet(
            X, eta
        ).second(),
        axial_third_mixed_jet_provider=lambda X, eta: _base_u_fourth_mixed_jet(
            X, eta
        ).third(),
        provenance="analytic leading-order test fixture; Issue #1 remains upstream",
        phi_third_mixed_jet_provider=_base_phi_third_mixed_jet,
        **kwargs,
    )


def _hierarchy():
    repaired = Section5PhiThirdMixedCoefficientJetSource.from_lemma52_repair_full_mixed(
        1,
        C=C,
        profile=_functional_repair(),
        base_phi_third_mixed_jet=_base_phi_third_mixed_jet,
        base_U_fourth_mixed_jet=_base_u_fourth_mixed_jet,
    )
    return Section5LowerHistoryPhiThirdMixedHierarchy(
        H,
        C,
        (_leading_source(), repaired),
        quadrature_points=64,
    )


def _landed_source(hierarchy, order, X, eta):
    phi = [hierarchy.phi_second_jet(j, X, eta) for j in range(order)]
    axial = [hierarchy.axial_second_jet(j, X, eta) for j in range(order)]
    beta = [hierarchy.beta_second_jet(j, X, eta) for j in range(order)]
    return positive_axis_point_data_from_lower_history(
        hierarchy.h,
        order,
        X,
        eta,
        phi,
        axial,
        beta,
    ).source


def _landed_forcing(hierarchy, order, xi, eta):
    return positive_axis_forcing(
        hierarchy.h,
        hierarchy.C,
        xi,
        eta,
        _landed_source(hierarchy, order, xi * xi, eta),
    )


@pytest.mark.parametrize("X", [0.0, 0.5 * 2.17**2])
def test_full_lower_source_eta_jet_matches_landed_value_finite_difference(X):
    hierarchy = _hierarchy()
    order = 2
    eta = 0.23
    eps = 2.0e-6

    actual = hierarchy_owned_lower_source_parameter_jet(
        hierarchy,
        order,
        X,
        eta,
    )
    center = _landed_source(hierarchy, order, X, eta)
    plus = _landed_source(hierarchy, order, X, eta + eps)
    minus = _landed_source(hierarchy, order, X, eta - eps)

    assert actual.value == center
    for field in ("angular", "axial", "pressure_product", "omega_quotient"):
        finite_difference = (
            getattr(plus, field) - getattr(minus, field)
        ) / (2.0 * eps)
        assert getattr(actual.parameter, field) == pytest.approx(
            finite_difference,
            rel=2.0e-5,
            abs=2.0e-7,
        )


@pytest.mark.parametrize("xi", [0.0, math.sqrt(0.5) * 2.17])
def test_hierarchy_owned_forcing_eta_jet_matches_landed_value_finite_difference(xi):
    hierarchy = _hierarchy()
    order = 2
    eta = 0.23
    eps = 2.0e-6

    actual = hierarchy_owned_positive_axis_forcing_parameter_jet(
        hierarchy,
        order,
        xi,
        eta,
    )

    center = _landed_forcing(hierarchy, order, xi, eta)
    plus = _landed_forcing(hierarchy, order, xi, eta + eps)
    minus = _landed_forcing(hierarchy, order, xi, eta - eps)
    finite_difference = (plus - minus) / (2.0 * eps)

    np.testing.assert_array_equal(actual.value, center)
    np.testing.assert_allclose(
        actual.parameter,
        finite_difference,
        rtol=2.0e-5,
        atol=2.0e-7,
    )


@pytest.mark.parametrize("xi", [0.0, math.sqrt(0.5) * 2.17])
def test_hierarchy_owned_first_picard_parameter_jet_matches_value_finite_difference(xi):
    hierarchy = _hierarchy()
    order = 2
    eta = 0.23
    eps = 2.0e-6
    quadrature_points = 32

    actual = hierarchy_owned_first_picard_parameter_jet_eq_5_7(
        hierarchy,
        order,
        xi,
        eta,
        quadrature_points=quadrature_points,
    )

    def landed_first_picard(e):
        return first_picard_term_eq_5_7(
            order,
            xi,
            e,
            lambda s, eta_value: _landed_forcing(
                hierarchy,
                order,
                s,
                eta_value,
            ),
            quadrature_points=quadrature_points,
        )

    center = landed_first_picard(eta)
    plus = landed_first_picard(eta + eps)
    minus = landed_first_picard(eta - eps)
    finite_difference = (plus - minus) / (2.0 * eps)

    np.testing.assert_allclose(actual.value, center, rtol=0.0, atol=1.0e-13)
    np.testing.assert_allclose(
        actual.parameter,
        finite_difference,
        rtol=3.0e-5,
        atol=3.0e-7,
    )


def test_hierarchy_owned_k1_picard_matches_independent_rhs_composition():
    hierarchy = _hierarchy()
    order = 2
    xi = math.sqrt(0.5) * 2.17
    eta = 0.23
    quadrature_points = 6

    actual = hierarchy_owned_k1_picard_value_eq_5_7(
        hierarchy,
        order,
        xi,
        eta,
        quadrature_points=quadrature_points,
    )
    fields = hierarchy.positive_axis_fields(order)

    def rhs(s, eta_value):
        first = hierarchy_owned_first_picard_parameter_jet_eq_5_7(
            hierarchy,
            order,
            s,
            eta_value,
            quadrature_points=quadrature_points,
        )
        forcing = hierarchy_owned_positive_axis_forcing_parameter_jet(
            hierarchy,
            order,
            s,
            eta_value,
        ).value
        return (
            fields.A0(s, eta_value) @ first.value
            + fields.A1(s, eta_value) @ first.parameter
            + forcing
        )

    expected = singular_inverse_eq_5_7(
        xi,
        eta,
        rhs,
        quadrature_points=quadrature_points,
    )
    first = hierarchy_owned_first_picard_parameter_jet_eq_5_7(
        hierarchy,
        order,
        xi,
        eta,
        quadrature_points=quadrature_points,
    ).value

    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=2.0e-13)
    assert not np.allclose(actual, first, rtol=1.0e-10, atol=1.0e-12)


def test_hierarchy_owned_k1_picard_axis_is_zero_after_strong_preflight():
    actual = hierarchy_owned_k1_picard_value_eq_5_7(
        _hierarchy(),
        2,
        0.0,
        0.23,
        quadrature_points=6,
    )
    np.testing.assert_array_equal(actual, np.zeros(6))


def test_order_one_fails_closed_without_leading_fourth_mixed_u():
    hierarchy = Section5LowerHistoryPhiThirdMixedHierarchy(
        H,
        C,
        (_leading_source(with_fourth=False),),
        quadrature_points=32,
    )

    with pytest.raises(ValueError, match="fourth-mixed U"):
        hierarchy_owned_lower_source_parameter_jet(
            hierarchy,
            1,
            0.5,
            0.1,
        )

    with pytest.raises(ValueError, match="fourth-mixed U"):
        hierarchy_owned_positive_axis_forcing_parameter_jet(
            hierarchy,
            1,
            math.sqrt(0.5),
            0.1,
        )

    # G(0)=0 must not bypass the hierarchy truth boundary at the axis.
    with pytest.raises(ValueError, match="fourth-mixed U"):
        hierarchy_owned_first_picard_parameter_jet_eq_5_7(
            hierarchy,
            1,
            0.0,
            0.1,
            quadrature_points=16,
        )

    # The nontrivial Picard bridge must preserve the same fail-closed boundary.
    with pytest.raises(ValueError, match="fourth-mixed U"):
        hierarchy_owned_k1_picard_value_eq_5_7(
            hierarchy,
            1,
            0.0,
            0.1,
            quadrature_points=6,
        )


def test_source_eta_jet_requires_strong_phi_history():
    with pytest.raises(TypeError, match="Section5LowerHistoryPhiThirdMixedHierarchy"):
        hierarchy_owned_lower_source_parameter_jet(object(), 1, 0.5, 0.1)


def test_forcing_eta_jet_requires_strong_phi_history():
    with pytest.raises(TypeError, match="Section5LowerHistoryPhiThirdMixedHierarchy"):
        hierarchy_owned_positive_axis_forcing_parameter_jet(object(), 1, 0.5, 0.1)


def test_first_picard_eta_jet_requires_strong_phi_history():
    with pytest.raises(TypeError, match="Section5LowerHistoryPhiThirdMixedHierarchy"):
        hierarchy_owned_first_picard_parameter_jet_eq_5_7(object(), 1, 0.5, 0.1)


def test_k1_picard_requires_strong_phi_history():
    with pytest.raises(TypeError, match="Section5LowerHistoryPhiThirdMixedHierarchy"):
        hierarchy_owned_k1_picard_value_eq_5_7(object(), 1, 0.5, 0.1)
