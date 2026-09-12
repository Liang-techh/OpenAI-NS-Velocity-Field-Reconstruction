import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_preceding_diffusion_parameter_jet import (
    ProfileThirdMixedJet,
)
from openai_ns_reconstruction.background_regular_flux_third_mixed_jets import (
    AxialFourthMixedJet,
)
from openai_ns_reconstruction.background_repaired_history_matrix_parameter import (
    hierarchy_owned_positive_axis_matrix_parameter_jets,
)
from openai_ns_reconstruction.background_repaired_history_phi_third_mixed import (
    Section5LowerHistoryPhiThirdMixedHierarchy,
    Section5PhiThirdMixedCoefficientJetSource,
)


C = 2.75
H = 0.005


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


def _hierarchy():
    leading = Section5PhiThirdMixedCoefficientJetSource(
        order=0,
        phi_second_jet_provider=lambda X, eta: _base_phi_third_mixed_jet(
            X, eta
        ).second(),
        axial_third_mixed_jet_provider=lambda X, eta: _base_u_fourth_mixed_jet(
            X, eta
        ).third(),
        axial_fourth_mixed_jet_provider=_base_u_fourth_mixed_jet,
        phi_third_mixed_jet_provider=_base_phi_third_mixed_jet,
        provenance="analytic leading-order matrix-jet test fixture; Issue #1 remains upstream",
    )
    return Section5LowerHistoryPhiThirdMixedHierarchy(
        H,
        C,
        (leading,),
        quadrature_points=64,
    )


@pytest.mark.parametrize("xi", [0.0, math.sqrt(0.5) * 2.17])
def test_hierarchy_owned_matrix_eta_jets_match_landed_fields_finite_difference(xi):
    hierarchy = _hierarchy()
    order = 1
    eta = 0.23
    eps = 2.0e-6

    actual = hierarchy_owned_positive_axis_matrix_parameter_jets(
        hierarchy,
        order,
        xi,
        eta,
    )
    fields = hierarchy.positive_axis_fields(order)

    for jet, field in ((actual.A0, fields.A0), (actual.A1, fields.A1)):
        center = field(xi, eta)
        plus = field(xi, eta + eps)
        minus = field(xi, eta - eps)
        finite_difference = (plus - minus) / (2.0 * eps)

        np.testing.assert_array_equal(jet.value, center)
        np.testing.assert_allclose(
            jet.parameter,
            finite_difference,
            rtol=3.0e-5,
            atol=3.0e-7,
        )


def test_matrix_eta_jets_reject_unowned_recursive_order():
    with pytest.raises(ValueError, match="strict lower coefficients"):
        hierarchy_owned_positive_axis_matrix_parameter_jets(
            _hierarchy(),
            2,
            0.5,
            0.1,
        )


def test_matrix_eta_jets_require_strong_history_type():
    with pytest.raises(TypeError, match="Section5LowerHistoryPhiThirdMixedHierarchy"):
        hierarchy_owned_positive_axis_matrix_parameter_jets(
            object(),
            1,
            0.5,
            0.1,
        )
