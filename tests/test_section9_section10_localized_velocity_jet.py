from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.local_field import LocalField
from openai_ns_reconstruction.section9_eq921_prefix_jet import (
    Section9FinitePrefixJetCertificate,
    required_spacetime_multiindices,
)
from openai_ns_reconstruction.section9_section10_localized_velocity_jet import (
    localize_section9_prefix_jet,
)
from openai_ns_reconstruction.spatial_localization import (
    section10_localized_field,
    section10_spatial_cutoff,
    section10_spatial_cutoff_gradient,
)


def _A(x: float, y: float, z: float, t: float) -> np.ndarray:
    del t
    return np.array([0.0, 0.0, x * y + 0.3 * x + 0.2 * z])


def _B(x: float, y: float, z: float, t: float) -> float:
    del x, y, z, t
    return 0.3


def _prefix_at(x: float, y: float, z: float, *, order: int = 1, **overrides):
    required = required_spacetime_multiindices(order)
    zero_vector = np.zeros(3)
    A = {alpha: np.array(zero_vector, copy=True) for alpha in required}
    B = {alpha: 0.0 for alpha in required}
    p = {alpha: 0.0 for alpha in required}

    A[(0, 0, 0, 0)] = _A(x, y, z, 0.6)
    B[(0, 0, 0, 0)] = 0.3
    p[(0, 0, 0, 0)] = 1.25
    if order >= 1:
        A[(0, 1, 0, 0)] = np.array([0.0, 0.0, y + 0.3])
        A[(0, 0, 1, 0)] = np.array([0.0, 0.0, x])
        A[(0, 0, 0, 1)] = np.array([0.0, 0.0, 0.2])

    kwargs = dict(
        q=Fraction(1, 2),
        prefix_order=1,
        stages=(1,),
        derivative_order=order,
        A_derivatives=A,
        B_derivatives=B,
        p_derivatives=p,
        provider_id="manufactured-prefix",
        provider_revision="test-r1",
        provider_provenance="analytic polynomial fixture",
    )
    kwargs.update(overrides)
    return Section9FinitePrefixJetCertificate(**kwargs)


def test_prefix_jet_localization_matches_independent_numeric_curl_of_cA():
    x, y, z, t = 0.19, 0.06, 0.02, 0.6
    prefix = _prefix_at(x, y, z)
    certificate = localize_section9_prefix_jet(prefix, x=x, y=y, z=z, t=t)

    # Independent hand curl for A=(0,0,xy+0.3x+0.2z).
    np.testing.assert_allclose(certificate.curl_A, np.array([x, -(y + 0.3), 0.0]))
    assert certificate.cutoff_value == section10_spatial_cutoff(x, y, z, t)
    np.testing.assert_allclose(
        certificate.cutoff_gradient,
        section10_spatial_cutoff_gradient(x, y, z, t),
    )

    # This path finite-differences curl(cA) as a whole because analytic_curl is
    # deliberately absent.  It therefore does not reuse the production
    # product-rule evaluation under test.
    numerical = section10_localized_field(LocalField(_A, _B)).velocity(
        x, y, z, t, eps=2.0e-6
    )
    np.testing.assert_allclose(
        certificate.localized_velocity,
        numerical,
        rtol=4.0e-5,
        atol=4.0e-6,
    )
    assert certificate.formal_localized_velocity_ready
    assert not certificate.paper_exact_velocity_available
    assert not certificate.residual_artifact_ready


def test_fixed_section10_support_forces_pointwise_zero_outside_cylinder():
    x, y, z, t = 0.30, 0.0, 0.0, 0.6
    certificate = localize_section9_prefix_jet(
        _prefix_at(x, y, z), x=x, y=y, z=z, t=t
    )
    assert certificate.cutoff_value == 0.0
    np.testing.assert_array_equal(certificate.cutoff_gradient, np.zeros(3))
    np.testing.assert_array_equal(certificate.localized_velocity, np.zeros(3))


def test_requires_first_spatial_derivatives_and_fail_closed_prefix_truth():
    with pytest.raises(ValueError, match="at least one"):
        localize_section9_prefix_jet(_prefix_at(0.1, 0.1, 0.0, order=0), x=0.1, y=0.1, z=0.0, t=0.6)

    untrusted = _prefix_at(0.1, 0.1, 0.0, paper_exact_velocity_available=True)
    with pytest.raises(ValueError, match="truth gate"):
        localize_section9_prefix_jet(untrusted, x=0.1, y=0.1, z=0.0, t=0.6)


def test_nonzero_azimuthal_scalar_on_axis_fails_closed():
    with pytest.raises(ValueError, match="vanish on the axis"):
        localize_section9_prefix_jet(_prefix_at(0.0, 0.0, 0.0), x=0.0, y=0.0, z=0.0, t=0.6)
