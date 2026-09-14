import numpy as np
import pytest

from openai_ns_reconstruction.background_inner_solver import singular_inverse_eq_5_7
from openai_ns_reconstruction.background_moment_repair_phi_fifth_mixed_jets import (
    ProfileFifthMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_fifth_mixed_jets import (
    AxialSixthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history_picard_fourth import (
    hierarchy_owned_fourth_picard_value,
)
from openai_ns_reconstruction.background_repaired_history_picard_third_first_parameter import (
    hierarchy_owned_third_picard_first_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_positive_axis_matrix_jets import (
    hierarchy_owned_positive_axis_matrix_second_parameter_jet,
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
        provenance="analytic fourth-Picard regression fixture; not paper coefficient data",
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
        quadrature_points=3,
    )


def _rhs(hierarchy, order, s, eta):
    third = hierarchy_owned_third_picard_first_parameter_jet(
        hierarchy,
        order,
        s,
        eta,
    )
    matrices = hierarchy_owned_positive_axis_matrix_second_parameter_jet(
        hierarchy,
        order,
        s,
        eta,
    )
    return matrices.A0 @ third.value + matrices.A1 @ third.parameter


def test_fourth_picard_value_is_exact_K_of_third_picard_term():
    hierarchy = _hierarchy()
    order = 1
    xi = 0.13
    eta = 0.09
    actual = hierarchy_owned_fourth_picard_value(
        hierarchy,
        order,
        xi,
        eta,
    )
    expected = singular_inverse_eq_5_7(
        xi,
        eta,
        lambda s, e: _rhs(hierarchy, order, s, e),
        quadrature_points=hierarchy.quadrature_points,
    )
    np.testing.assert_array_equal(actual.value, expected)


def test_fourth_picard_value_closes_current_eta_derivative_triangle():
    hierarchy = _hierarchy()
    third = hierarchy_owned_third_picard_first_parameter_jet(
        hierarchy,
        1,
        0.12,
        -0.08,
    )
    actual = hierarchy_owned_fourth_picard_value(
        hierarchy,
        1,
        0.12,
        -0.08,
    )
    assert third.value.shape == (6,)
    assert third.parameter.shape == (6,)
    assert actual.value.shape == (6,)
    assert np.all(np.isfinite(actual.value))
    assert not actual.value.flags.writeable


def test_fourth_picard_value_axis_is_exactly_zero_and_read_only():
    actual = hierarchy_owned_fourth_picard_value(
        _hierarchy(),
        1,
        0.0,
        0.18,
    )
    np.testing.assert_array_equal(actual.value, np.zeros(6))
    assert not actual.value.flags.writeable


def test_fourth_picard_value_fails_closed_on_wrong_hierarchy_type():
    with pytest.raises(TypeError, match="Section5LowerHistorySixthMixedHierarchy"):
        hierarchy_owned_fourth_picard_value(object(), 1, 0.0, 0.0)
