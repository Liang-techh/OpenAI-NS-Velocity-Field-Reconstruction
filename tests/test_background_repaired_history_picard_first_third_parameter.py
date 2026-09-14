import numpy as np
import pytest

from openai_ns_reconstruction.background_inner_solver import (
    first_picard_term_eq_5_7,
    singular_inverse_eq_5_7,
)
from openai_ns_reconstruction.background_moment_repair_phi_fifth_mixed_jets import (
    ProfileFifthMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_fifth_mixed_jets import (
    AxialSixthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history_forcing_third_parameter import (
    hierarchy_owned_positive_axis_forcing_third_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_picard_first_third_parameter import (
    hierarchy_owned_first_picard_third_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_sixth_mixed import (
    Section5LowerHistorySixthMixedHierarchy,
    Section5SixthMixedCoefficientJetSource,
)


H = 0.005
C = 2.75


def _phi(X, eta):
    del eta
    value = 0.31 * (1.0 + 0.17 * X + 0.03 * X * X)
    radial = 0.31 * (0.17 + 0.06 * X)
    radial2 = 0.31 * 0.06
    return ProfileFifthMixedJet(
        value=value,
        radial=radial,
        radial2=radial2,
        parameter=0.0,
        radial_parameter=0.0,
        parameter2=0.0,
        radial2_parameter=0.0,
        radial_parameter2=0.0,
        parameter3=0.0,
        radial2_parameter2=0.0,
        radial_parameter3=0.0,
        parameter4=0.0,
        radial2_parameter3=0.0,
        radial_parameter4=0.0,
        parameter5=0.0,
    )


def _axial(X, eta):
    del eta
    value = 0.22 * (1.0 + 0.11 * X + 0.025 * X * X)
    radial = 0.22 * (0.11 + 0.05 * X)
    radial2 = 0.22 * 0.05
    return AxialSixthMixedJet(
        value=value,
        radial=radial,
        radial2=radial2,
        parameter=0.0,
        radial_parameter=0.0,
        parameter2=0.0,
        radial2_parameter=0.0,
        radial_parameter2=0.0,
        parameter3=0.0,
        radial2_parameter2=0.0,
        radial_parameter3=0.0,
        parameter4=0.0,
        radial2_parameter3=0.0,
        radial_parameter4=0.0,
        parameter5=0.0,
        radial2_parameter4=0.0,
        radial_parameter5=0.0,
        parameter6=0.0,
    )


def _source0():
    return Section5SixthMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _phi(X, eta).fourth().third().second(),
        axial_third_mixed_jet_provider=lambda X, eta: _axial(X, eta).fifth().fourth().third(),
        provenance="analytic Picard regression fixture; not paper coefficient data",
        normalization_C=None,
        axial_fourth_mixed_jet_provider=lambda X, eta: _axial(X, eta).fifth().fourth(),
        phi_third_mixed_jet_provider=lambda X, eta: _phi(X, eta).fourth().third(),
        axial_fifth_mixed_jet_provider=lambda X, eta: _axial(X, eta).fifth(),
        phi_fourth_mixed_jet_provider=lambda X, eta: _phi(X, eta).fourth(),
        phi_fifth_mixed_jet_provider=_phi,
        axial_sixth_mixed_jet_provider=_axial,
    )


def _hierarchy():
    return Section5LowerHistorySixthMixedHierarchy(
        H,
        C,
        (_source0(),),
        quadrature_points=6,
    )


def test_first_picard_third_eta_value_delegates_to_landed_Gf():
    hierarchy = _hierarchy()
    xi = 0.24
    eta = 0.13
    actual = hierarchy_owned_first_picard_third_parameter_jet(
        hierarchy, 1, xi, eta
    )
    expected = first_picard_term_eq_5_7(
        1,
        xi,
        eta,
        lambda s, e: hierarchy_owned_positive_axis_forcing_third_parameter_jet(
            hierarchy, 1, s, e
        ).value,
        quadrature_points=hierarchy.quadrature_points,
    )
    np.testing.assert_array_equal(actual.value, expected)


def test_first_picard_third_eta_commutes_with_singular_inverse():
    hierarchy = _hierarchy()
    xi = 0.21
    eta = -0.17
    actual = hierarchy_owned_first_picard_third_parameter_jet(
        hierarchy, 1, xi, eta
    )
    expected = singular_inverse_eq_5_7(
        xi,
        eta,
        lambda s, e: hierarchy_owned_positive_axis_forcing_third_parameter_jet(
            hierarchy, 1, s, e
        ).parameter3,
        quadrature_points=hierarchy.quadrature_points,
    )
    np.testing.assert_array_equal(actual.parameter3, expected)


def test_first_picard_third_eta_axis_value_and_all_derivatives_are_exactly_zero():
    hierarchy = _hierarchy()
    actual = hierarchy_owned_first_picard_third_parameter_jet(
        hierarchy, 1, 0.0, 0.19
    )
    for row in (actual.value, actual.parameter, actual.parameter2, actual.parameter3):
        np.testing.assert_array_equal(row, np.zeros(6))


def test_first_picard_third_eta_fails_closed_on_wrong_hierarchy_type():
    with pytest.raises(TypeError, match="Section5LowerHistorySixthMixedHierarchy"):
        hierarchy_owned_first_picard_third_parameter_jet(object(), 1, 0.0, 0.0)
