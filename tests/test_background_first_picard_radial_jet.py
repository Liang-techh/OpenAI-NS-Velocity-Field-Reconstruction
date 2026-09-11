import math

import numpy as np
import pytest

from openai_ns_reconstruction.background_first_picard_radial_jet import (
    first_picard_radial_jet_from_hierarchy,
)
from openai_ns_reconstruction.background_lower_history_source import ProfileSecondJet
from openai_ns_reconstruction.background_regular_flux_second_jets import AxialThirdMixedJet
from openai_ns_reconstruction.background_repaired_history import (
    Section5CoefficientJetSource,
    Section5LowerHistoryJetHierarchy,
)


H = 0.005
C = 2.4


def _phi_second_jet(X, eta):
    X = float(X)
    eta = float(eta)
    e = math.exp(-X)
    p = 1.0 + 0.1 * eta + 0.05 * eta * eta
    dp = 0.1 + 0.1 * eta
    d2p = 0.1
    value = C * e * p
    parameter = C * e * dp
    return ProfileSecondJet(
        value=value,
        radial=-value,
        radial2=value,
        parameter=parameter,
        radial_parameter=-parameter,
        parameter2=C * e * d2p,
    )


def _u_third_mixed_jet(X, eta):
    X = float(X)
    eta = float(eta)
    e = math.exp(-X)
    q = eta + 0.1 * eta * eta
    dq = 1.0 + 0.2 * eta
    d2q = 0.2
    value = e * q
    parameter = e * dq
    parameter2 = e * d2q
    return AxialThirdMixedJet(
        value=value,
        radial=-value,
        radial2=value,
        parameter=parameter,
        radial_parameter=-parameter,
        parameter2=parameter2,
        radial2_parameter=parameter,
        radial_parameter2=-parameter2,
        parameter3=0.0,
    )


def _hierarchy():
    leading = Section5CoefficientJetSource(
        order=0,
        phi_second_jet_provider=_phi_second_jet,
        axial_third_mixed_jet_provider=_u_third_mixed_jet,
        provenance="analytic leading fixture for radial-jet regression; Issue #1 unresolved",
    )
    return Section5LowerHistoryJetHierarchy(
        H,
        C,
        (leading,),
        quadrature_points=96,
    )


def test_radial_jet_matches_independent_centered_difference_of_integral_term():
    hierarchy = _hierarchy()
    xi = 0.91
    eta = 0.17
    jet = first_picard_radial_jet_from_hierarchy(
        hierarchy,
        1,
        xi,
        eta,
        quadrature_points=128,
    )

    step = 2.0e-5
    plus = hierarchy.first_picard_term(
        1, xi + step, eta, quadrature_points=128
    )
    minus = hierarchy.first_picard_term(
        1, xi - step, eta, quadrature_points=128
    )
    finite_difference = (plus - minus) / (2.0 * step)

    assert np.allclose(jet.d_xi, finite_difference, rtol=3e-7, atol=3e-9)
    assert np.allclose(
        jet.value,
        hierarchy.first_picard_term(1, xi, eta, quadrature_points=128),
        rtol=0.0,
        atol=0.0,
    )


def test_radial_jet_closes_positive_radius_singular_ode_identity():
    hierarchy = _hierarchy()
    xi = 1.13
    eta = -0.21
    jet = first_picard_radial_jet_from_hierarchy(hierarchy, 1, xi, eta)

    weights = np.array([0.0, 0.0, 2.0, 0.0, 3.0, 1.0])
    independent_residual = (
        jet.d_xi + weights * jet.value / xi - hierarchy.positive_axis_fields(1).forcing(xi, eta)
    )
    assert np.max(np.abs(independent_residual)) < 2e-13
    assert jet.max_abs_residual < 2e-13
    assert not jet.value.flags.writeable
    assert not jet.d_xi.flags.writeable


def test_radial_jet_fails_closed_at_axis_and_without_complete_lower_history():
    hierarchy = _hierarchy()
    with pytest.raises(ValueError):
        first_picard_radial_jet_from_hierarchy(hierarchy, 1, 0.0, 0.0)
    with pytest.raises(ValueError):
        first_picard_radial_jet_from_hierarchy(hierarchy, 1, -0.1, 0.0)
    with pytest.raises(ValueError):
        first_picard_radial_jet_from_hierarchy(hierarchy, 2, 0.8, 0.0)
    with pytest.raises(TypeError):
        first_picard_radial_jet_from_hierarchy(object(), 1, 0.8, 0.0)
