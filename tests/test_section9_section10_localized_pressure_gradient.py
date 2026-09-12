from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.section9_eq921_prefix_jet import (
    Section9FinitePrefixJetCertificate,
    required_spacetime_multiindices,
)
from openai_ns_reconstruction.section9_section10_localized_pressure_gradient import (
    localize_section9_prefix_pressure_gradient,
)
from openai_ns_reconstruction.spatial_localization import (
    section10_spatial_cutoff,
)


def _linear_pressure(x: float, y: float, z: float, t: float) -> float:
    return 5.0 + 2.0 * x - 3.0 * y + 4.0 * z + 0.5 * t


def _prefix_at(
    x: float,
    y: float,
    z: float,
    t: float,
    *,
    derivative_order: int = 1,
) -> Section9FinitePrefixJetCertificate:
    indices = required_spacetime_multiindices(derivative_order)
    A = {alpha: np.zeros(3) for alpha in indices}
    B = {alpha: 0.0 for alpha in indices}
    p = {alpha: 0.0 for alpha in indices}
    p[(0, 0, 0, 0)] = _linear_pressure(x, y, z, t)
    if derivative_order >= 1:
        p[(1, 0, 0, 0)] = 0.5
        p[(0, 1, 0, 0)] = 2.0
        p[(0, 0, 1, 0)] = -3.0
        p[(0, 0, 0, 1)] = 4.0
    return Section9FinitePrefixJetCertificate(
        q=Fraction(1, 2),
        prefix_order=1,
        stages=(1,),
        derivative_order=derivative_order,
        A_derivatives=A,
        B_derivatives=B,
        p_derivatives=p,
        provider_id="manufactured-pressure-jet",
        provider_revision="test-v1",
        provider_provenance="tests/test_section9_section10_localized_pressure_gradient.py",
    )


def test_transition_collar_gradient_matches_independent_finite_difference() -> None:
    x, y, z, t = 0.21, 0.02, 0.18, 0.6
    prefix = _prefix_at(x, y, z, t)
    certificate = localize_section9_prefix_pressure_gradient(
        prefix, x=x, y=y, z=z, t=t
    )

    def cut_pressure(xp: float, yp: float, zp: float) -> float:
        return (
            section10_spatial_cutoff(xp, yp, zp, t)
            * _linear_pressure(xp, yp, zp, t)
        )

    eps = 1.0e-6
    numeric = np.empty(3)
    numeric[0] = (
        cut_pressure(x + eps, y, z) - cut_pressure(x - eps, y, z)
    ) / (2.0 * eps)
    numeric[1] = (
        cut_pressure(x, y + eps, z) - cut_pressure(x, y - eps, z)
    ) / (2.0 * eps)
    numeric[2] = (
        cut_pressure(x, y, z + eps) - cut_pressure(x, y, z - eps)
    ) / (2.0 * eps)

    np.testing.assert_allclose(
        certificate.localized_pressure_gradient,
        numeric,
        rtol=2.0e-5,
        atol=2.0e-7,
    )
    assert certificate.formal_localized_pressure_gradient_ready
    assert certificate.paper_exact_velocity_available is False
    assert certificate.residual_artifact_ready is False


def test_plateau_preserves_pressure_gradient_exactly() -> None:
    x, y, z, t = 0.05, 0.02, 0.03, 0.4
    certificate = localize_section9_prefix_pressure_gradient(
        _prefix_at(x, y, z, t), x=x, y=y, z=z, t=t
    )
    assert certificate.cutoff_value == 1.0
    np.testing.assert_array_equal(certificate.cutoff_gradient, np.zeros(3))
    np.testing.assert_array_equal(
        certificate.localized_pressure_gradient,
        np.array([2.0, -3.0, 4.0]),
    )
    assert certificate.localized_pressure == _linear_pressure(x, y, z, t)


def test_official_support_exterior_has_zero_cut_pressure_and_gradient() -> None:
    x, y, z, t = 0.30, 0.0, 0.0, 0.7
    certificate = localize_section9_prefix_pressure_gradient(
        _prefix_at(x, y, z, t), x=x, y=y, z=z, t=t
    )
    assert certificate.cutoff_value == 0.0
    np.testing.assert_array_equal(certificate.cutoff_gradient, np.zeros(3))
    assert certificate.localized_pressure == 0.0
    np.testing.assert_array_equal(
        certificate.localized_pressure_gradient,
        np.zeros(3),
    )


def test_derivative_order_zero_fails_closed() -> None:
    prefix = _prefix_at(0.1, 0.0, 0.0, 0.5, derivative_order=0)
    assert prefix.formal_prefix_jet_ready
    with pytest.raises(ValueError, match="derivative_order must be at least one"):
        localize_section9_prefix_pressure_gradient(
            prefix, x=0.1, y=0.0, z=0.0, t=0.5
        )
