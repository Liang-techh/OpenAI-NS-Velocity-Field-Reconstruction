import numpy as np
import pytest

from openai_ns_reconstruction.background_moment_repair_phi_fifth_mixed_jets import (
    ProfileFifthMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_fifth_mixed_jets import (
    AxialSixthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history_picard_enclosure import (
    HierarchyBoundPicardMajorantInputs,
    hierarchy_owned_conditional_picard_value_enclosure,
    hierarchy_owned_picard_prefix4_value,
)
from openai_ns_reconstruction.background_repaired_history_picard_first_third_parameter import (
    hierarchy_owned_first_picard_third_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_picard_second_second_parameter import (
    hierarchy_owned_second_picard_second_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_picard_third_first_parameter import (
    hierarchy_owned_third_picard_first_parameter_jet,
)
from openai_ns_reconstruction.background_repaired_history_picard_fourth import (
    hierarchy_owned_fourth_picard_value,
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
        provenance="analytic Picard-enclosure regression fixture; not paper coefficient data",
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
        quadrature_points=1,
    )


def _majorant(hierarchy):
    return HierarchyBoundPicardMajorantInputs(
        hierarchy=hierarchy,
        order=1,
        C_n=0.1,
        radial_extent=0.1,
        rho=0.8,
        rho_prime=0.2,
        provenance=(
            "test-only analytic majorant inputs; exercises Eq. (5.8) binding, "
            "not a paper coefficient bound"
        ),
    )


def test_prefix4_is_exact_sum_of_landed_hierarchy_owned_terms():
    hierarchy = _hierarchy()
    xi = 0.035
    eta = 0.04
    actual = hierarchy_owned_picard_prefix4_value(hierarchy, 1, xi, eta)

    first = hierarchy_owned_first_picard_third_parameter_jet(hierarchy, 1, xi, eta)
    second = hierarchy_owned_second_picard_second_parameter_jet(hierarchy, 1, xi, eta)
    third = hierarchy_owned_third_picard_first_parameter_jet(hierarchy, 1, xi, eta)
    fourth = hierarchy_owned_fourth_picard_value(hierarchy, 1, xi, eta)
    expected = first.value + second.value + third.value + fourth.value

    np.testing.assert_array_equal(actual, expected)
    assert not actual.flags.writeable


def test_conditional_enclosure_uses_same_hierarchy_prefix_and_eq_5_8_tail():
    hierarchy = _hierarchy()
    majorant = _majorant(hierarchy)
    actual = hierarchy_owned_conditional_picard_value_enclosure(
        hierarchy,
        1,
        0.03,
        -0.05,
        majorant,
    )
    expected_prefix = hierarchy_owned_picard_prefix4_value(hierarchy, 1, 0.03, -0.05)

    np.testing.assert_array_equal(actual.prefix_value, expected_prefix)
    assert actual.start_order == 4
    assert actual.tail_upper_bound == majorant.tail_certificate.tail_upper_bound
    assert np.isfinite(actual.tail_upper_bound)
    assert not actual.prefix_value.flags.writeable
    assert actual.certifies_tolerance(actual.tail_upper_bound * 1.01)


def test_conditional_enclosure_rejects_cross_wired_hierarchy_and_order():
    hierarchy = _hierarchy()
    other = _hierarchy()
    majorant = _majorant(hierarchy)

    with pytest.raises(ValueError, match="different hierarchy"):
        hierarchy_owned_conditional_picard_value_enclosure(
            other, 1, 0.02, 0.0, majorant
        )
    with pytest.raises(ValueError, match="different recursive order"):
        hierarchy_owned_conditional_picard_value_enclosure(
            hierarchy, 2, 0.02, 0.0, majorant
        )


def test_conditional_enclosure_fails_closed_outside_majorant_domain():
    hierarchy = _hierarchy()
    majorant = _majorant(hierarchy)

    with pytest.raises(ValueError, match="radial extent"):
        hierarchy_owned_conditional_picard_value_enclosure(
            hierarchy, 1, 0.11, 0.0, majorant
        )
    with pytest.raises(ValueError, match="smaller parameter strip"):
        hierarchy_owned_conditional_picard_value_enclosure(
            hierarchy, 1, 0.02, 0.21, majorant
        )


def test_picard_majorant_and_enclosure_fail_closed_on_wrong_types():
    with pytest.raises(TypeError, match="Section5LowerHistorySixthMixedHierarchy"):
        HierarchyBoundPicardMajorantInputs(
            hierarchy=object(),
            order=1,
            C_n=0.1,
            radial_extent=0.1,
            rho=0.8,
            rho_prime=0.2,
            provenance="wrong-type regression",
        )
    with pytest.raises(TypeError, match="HierarchyBoundPicardMajorantInputs"):
        hierarchy_owned_conditional_picard_value_enclosure(
            _hierarchy(), 1, 0.0, 0.0, object()
        )
