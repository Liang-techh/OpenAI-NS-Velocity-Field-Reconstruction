import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair_phi_fifth_mixed_jets import (
    ProfileFifthMixedJet,
)
from openai_ns_reconstruction.background_positive_axis import (
    PositiveAxisBaseJet,
    RadialParameterJet,
    positive_axis_A0,
    positive_axis_A1,
)
from openai_ns_reconstruction.background_regular_flux_fifth_mixed_jets import (
    AxialSixthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history_positive_axis_base import (
    hierarchy_owned_positive_axis_base_jet,
    hierarchy_owned_positive_axis_matrices,
)
from openai_ns_reconstruction.background_repaired_history_sixth_mixed import (
    Section5LowerHistorySixthMixedHierarchy,
    Section5SixthMixedCoefficientJetSource,
)


H = 0.005
C = 2.75


def _phi(X, eta):
    value = 0.31 * (1.0 + 0.17 * X + 0.03 * X * X) * (1.0 + 0.04 * eta)
    radial = 0.31 * (0.17 + 0.06 * X) * (1.0 + 0.04 * eta)
    radial2 = 0.31 * 0.06 * (1.0 + 0.04 * eta)
    return ProfileFifthMixedJet(
        value=value,
        radial=radial,
        radial2=radial2,
        parameter=0.04 * 0.31 * (1.0 + 0.17 * X + 0.03 * X * X),
        radial_parameter=0.04 * 0.31 * (0.17 + 0.06 * X),
        parameter2=0.0,
        radial2_parameter=0.04 * 0.31 * 0.06,
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
    value = 0.22 * (1.0 + 0.11 * X + 0.025 * X * X) * (1.0 - 0.03 * eta)
    radial = 0.22 * (0.11 + 0.05 * X) * (1.0 - 0.03 * eta)
    radial2 = 0.22 * 0.05 * (1.0 - 0.03 * eta)
    return AxialSixthMixedJet(
        value=value,
        radial=radial,
        radial2=radial2,
        parameter=-0.03 * 0.22 * (1.0 + 0.11 * X + 0.025 * X * X),
        radial_parameter=-0.03 * 0.22 * (0.11 + 0.05 * X),
        parameter2=0.0,
        radial2_parameter=-0.03 * 0.22 * 0.05,
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
        provenance="analytic PositiveAxis-base regression fixture; not paper coefficient data",
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


def test_positive_axis_base_is_exact_projection_of_order_zero_hierarchy():
    hierarchy = _hierarchy()
    X = 0.19
    eta = -0.23
    actual = hierarchy_owned_positive_axis_base_jet(hierarchy, X, eta)
    phi = hierarchy.phi_fifth_mixed_jet(0, X, eta)
    axial = hierarchy.axial_sixth_mixed_jet(0, X, eta)
    beta = hierarchy.beta_fifth_mixed_jet(0, X, eta)
    expected = PositiveAxisBaseJet(
        phi=RadialParameterJet(phi.value, phi.radial, phi.radial2, phi.parameter),
        axial=RadialParameterJet(axial.value, axial.radial, axial.radial2, axial.parameter),
        beta=beta.value,
    )
    assert actual == expected


def test_positive_axis_matrices_delegate_to_landed_exact_formulas():
    hierarchy = _hierarchy()
    order = 1
    xi = 0.37
    eta = 0.18
    actual = hierarchy_owned_positive_axis_matrices(hierarchy, order, xi, eta)
    base = hierarchy_owned_positive_axis_base_jet(hierarchy, xi * xi, eta)
    expected_A0 = positive_axis_A0(H, order, C, xi, eta, base)
    expected_A1 = positive_axis_A1(H, xi, eta, base)
    np.testing.assert_array_equal(actual.A0, expected_A0)
    np.testing.assert_array_equal(actual.A1, expected_A1)
    assert not actual.A0.flags.writeable
    assert not actual.A1.flags.writeable


def test_positive_axis_base_uses_hierarchy_eq_5_2_beta_not_independent_table():
    hierarchy = _hierarchy()
    X = 0.07
    eta = 0.29
    actual = hierarchy_owned_positive_axis_base_jet(hierarchy, X, eta)
    assert actual.beta == hierarchy.beta_fifth_mixed_jet(0, X, eta).value


def test_positive_axis_base_fails_closed_on_wrong_hierarchy_type():
    with pytest.raises(TypeError, match="Section5LowerHistorySixthMixedHierarchy"):
        hierarchy_owned_positive_axis_base_jet(object(), 0.0, 0.0)
    with pytest.raises(TypeError, match="Section5LowerHistorySixthMixedHierarchy"):
        hierarchy_owned_positive_axis_matrices(object(), 1, 0.0, 0.0)
