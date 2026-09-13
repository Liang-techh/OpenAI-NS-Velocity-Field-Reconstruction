import numpy as np
import pytest
from fractions import Fraction

from openai_ns_reconstruction.local_field import LocalField, LocalizedField
from openai_ns_reconstruction.section9_eq921_prefix_jet import (
    Section9FinitePrefixJetCertificate,
    required_spacetime_multiindices,
)
from openai_ns_reconstruction.section9_section10_convective_acceleration import (
    localize_section9_prefix_convective_acceleration,
)
from openai_ns_reconstruction.spatial_localization import section10_spatial_cutoff


def _A(x: float, y: float, z: float, t: float) -> np.ndarray:
    del t
    return np.array([0.0, 0.0, x * y + 0.3 * x + 0.2 * z])


def _B(x: float, y: float, z: float, t: float) -> float:
    del x, y, z, t
    return 0.0


def _prefix_at(
    x: float,
    y: float,
    z: float,
    t: float,
    *,
    order: int = 2,
    **overrides,
) -> Section9FinitePrefixJetCertificate:
    del t
    required = required_spacetime_multiindices(order)
    A = {alpha: np.zeros(3) for alpha in required}
    B = {alpha: 0.0 for alpha in required}
    p = {alpha: 0.0 for alpha in required}

    A[(0, 0, 0, 0)] = _A(x, y, z, 0.0)
    if order >= 1:
        A[(0, 1, 0, 0)] = np.array([0.0, 0.0, y + 0.3])
        A[(0, 0, 1, 0)] = np.array([0.0, 0.0, x])
        A[(0, 0, 0, 1)] = np.array([0.0, 0.0, 0.2])
    if order >= 2:
        A[(0, 1, 1, 0)] = np.array([0.0, 0.0, 1.0])

    kwargs = dict(
        q=Fraction(1, 2),
        prefix_order=1,
        stages=(1,),
        derivative_order=order,
        A_derivatives=A,
        B_derivatives=B,
        p_derivatives=p,
        provider_id="manufactured-convection-prefix",
        provider_revision="test-r1",
        provider_provenance="analytic polynomial fixture",
    )
    kwargs.update(overrides)
    return Section9FinitePrefixJetCertificate(**kwargs)


def test_convective_acceleration_matches_independent_directional_difference():
    # Both radial and axial factors lie in the executable fixed-cutoff
    # transition collar, so the oracle exercises the full localization rather
    # than only the plateau identity.
    x, y, z, t = 0.20, 0.02, 0.16, 0.6
    certificate = localize_section9_prefix_convective_acceleration(
        _prefix_at(x, y, z, t), x=x, y=y, z=z, t=t
    )

    # Independent oracle: finite-difference the complete LocalizedField
    # velocity in the frozen direction u(x0).  No analytic cutoff gradient or
    # analytic spatial Jacobian is supplied here, so the poloidal velocity is
    # obtained by numerically curling cA as a whole.
    field = LocalizedField(
        LocalField(_A, _B),
        section10_spatial_cutoff,
        cutoff_gradient=None,
    )
    inner_eps = 1.5e-6
    u0 = field.velocity(x, y, z, t, eps=inner_eps)
    h = 3.0e-5

    plus = field.velocity(
        x + h * u0[0],
        y + h * u0[1],
        z + h * u0[2],
        t,
        eps=inner_eps,
    )
    minus = field.velocity(
        x - h * u0[0],
        y - h * u0[1],
        z - h * u0[2],
        t,
        eps=inner_eps,
    )
    independent = (plus - minus) / (2.0 * h)

    np.testing.assert_allclose(
        certificate.convective_acceleration,
        independent,
        rtol=4.0e-3,
        atol=4.0e-3,
    )
    assert certificate.formal_convective_acceleration_ready
    assert not certificate.axis_convective_acceleration_covered
    assert not certificate.residual_artifact_ready
    assert not certificate.paper_exact_velocity_available


def test_convective_acceleration_is_exactly_zero_outside_fixed_support():
    x, y, z, t = 0.30, 0.02, 0.01, 0.6
    certificate = localize_section9_prefix_convective_acceleration(
        _prefix_at(x, y, z, t), x=x, y=y, z=z, t=t
    )
    np.testing.assert_array_equal(
        certificate.spatial.base.localized_velocity,
        np.zeros(3),
    )
    np.testing.assert_array_equal(
        certificate.spatial.localized_velocity_spatial_jacobian,
        np.zeros((3, 3)),
    )
    np.testing.assert_array_equal(
        certificate.convective_acceleration,
        np.zeros(3),
    )


def test_convective_acceleration_inherits_second_order_and_truth_gates():
    x, y, z, t = 0.20, 0.02, 0.16, 0.6
    with pytest.raises(ValueError, match="at least two"):
        localize_section9_prefix_convective_acceleration(
            _prefix_at(x, y, z, t, order=1), x=x, y=y, z=z, t=t
        )

    with pytest.raises(ValueError, match="truth gate"):
        localize_section9_prefix_convective_acceleration(
            _prefix_at(x, y, z, t, paper_exact_velocity_available=True),
            x=x,
            y=y,
            z=z,
            t=t,
        )
