import math
from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.local_field import LocalField
from openai_ns_reconstruction.section9_eq921_prefix_jet import (
    Section9FinitePrefixJetCertificate,
    required_spacetime_multiindices,
)
from openai_ns_reconstruction.section9_section10_formal_ns_operator import (
    assemble_section9_section10_formal_ns_operator,
)
from openai_ns_reconstruction.spatial_localization import (
    section10_localized_field,
    section10_spatial_cutoff,
)
from openai_ns_reconstruction.verify import (
    gradient_scalar_numeric,
    jacobian_numeric,
    laplacian_vector_numeric,
    time_derivative_numeric,
)


def _monomial_cube_derivative(value: float, order: int) -> float:
    if order == 0:
        return value**3
    if order == 1:
        return 3.0 * value**2
    if order == 2:
        return 6.0 * value
    if order == 3:
        return 6.0
    return 0.0


def _F_spatial_derivative(
    dx: int, dy: int, dz: int, x: float, y: float, z: float
) -> float:
    active = sum(value > 0 for value in (dx, dy, dz))
    if active == 0:
        return x**3 + y**3 + z**3
    if active > 1:
        return 0.0
    if dx:
        return _monomial_cube_derivative(x, dx)
    if dy:
        return _monomial_cube_derivative(y, dy)
    return _monomial_cube_derivative(z, dz)


def _G_spatial_derivative(
    dx: int, dy: int, dz: int, x: float, y: float, z: float
) -> float:
    degree = dx + dy + dz
    if degree == 0:
        return 0.4 * x - 0.2 * y + 0.1 * z + 0.3
    if degree != 1:
        return 0.0
    if dx == 1:
        return 0.4
    if dy == 1:
        return -0.2
    if dz == 1:
        return 0.1
    return 0.0


def _A(x: float, y: float, z: float, t: float) -> np.ndarray:
    return np.array([0.0, 0.0, (1.0 + t) * (x**3 + y**3 + z**3)])


def _B(x: float, y: float, z: float, t: float) -> float:
    del x, y, z, t
    return 0.0


def _p(x: float, y: float, z: float, t: float) -> float:
    return (1.0 + 0.1 * t) * (0.4 * x - 0.2 * y + 0.1 * z + 0.3)


def _analytic_curl(x: float, y: float, z: float, t: float) -> np.ndarray:
    del z
    return (1.0 + t) * np.array([3.0 * y * y, -3.0 * x * x, 0.0])


def _prefix_at(
    x: float,
    y: float,
    z: float,
    t: float,
    *,
    order: int = 3,
    **overrides,
) -> Section9FinitePrefixJetCertificate:
    required = required_spacetime_multiindices(order)
    A = {}
    B = {}
    p = {}
    for dt, dx, dy, dz in required:
        spatial_A = _F_spatial_derivative(dx, dy, dz, x, y, z)
        if dt == 0:
            A_factor = 1.0 + t
        elif dt == 1:
            A_factor = 1.0
        else:
            A_factor = 0.0
        A[(dt, dx, dy, dz)] = np.array(
            [0.0, 0.0, A_factor * spatial_A]
        )
        B[(dt, dx, dy, dz)] = 0.0

        spatial_p = _G_spatial_derivative(dx, dy, dz, x, y, z)
        if dt == 0:
            p_factor = 1.0 + 0.1 * t
        elif dt == 1:
            p_factor = 0.1
        else:
            p_factor = 0.0
        p[(dt, dx, dy, dz)] = p_factor * spatial_p

    kwargs = dict(
        q=Fraction(1, 2),
        prefix_order=1,
        stages=(1,),
        derivative_order=order,
        A_derivatives=A,
        B_derivatives=B,
        p_derivatives=p,
        provider_id="manufactured-formal-ns-operator-prefix",
        provider_revision="test-r1",
        provider_provenance="analytic time-dependent polynomial fixture",
    )
    kwargs.update(overrides)
    return Section9FinitePrefixJetCertificate(**kwargs)


def test_formal_ns_operator_matches_independent_final_field_differences():
    # Both radial and axial cutoff factors are in their transition collars.
    # The manufactured field has nonzero u_t, convection, Delta u and grad(cp).
    x, y, z, t = 0.20, 0.02, 0.16, 0.6
    nu = 0.7
    certificate = assemble_section9_section10_formal_ns_operator(
        _prefix_at(x, y, z, t),
        x=x,
        y=y,
        z=z,
        t=t,
        viscosity=nu,
    )

    field = section10_localized_field(
        LocalField(_A, _B, analytic_curl=_analytic_curl)
    )
    velocity = field.velocity
    localized_pressure = lambda xx, yy, zz, tt: (
        section10_spatial_cutoff(xx, yy, zz, tt) * _p(xx, yy, zz, tt)
    )

    u0 = velocity(x, y, z, t)
    independent = (
        time_derivative_numeric(velocity, x, y, z, t, eps=2.0e-6)
        + jacobian_numeric(velocity, x, y, z, t, eps=2.0e-6) @ u0
        - nu * laplacian_vector_numeric(velocity, x, y, z, t, eps=2.0e-5)
        + gradient_scalar_numeric(localized_pressure, x, y, z, t, eps=2.0e-6)
    )

    np.testing.assert_allclose(
        certificate.formal_operator_sum,
        independent,
        rtol=2.0e-5,
        atol=4.0e-3,
    )
    np.testing.assert_allclose(
        certificate.viscous_term,
        -nu * certificate.laplacian.localized_velocity_laplacian,
        rtol=0.0,
        atol=0.0,
    )
    assert certificate.formal_ns_operator_ready
    assert certificate.component_source_identity_verified
    assert certificate.fixed_section10_cutoff_coherent
    assert not certificate.residual_artifact_ready
    assert not certificate.forcing_artifact_ready
    assert not certificate.paper_exact_velocity_available


def test_formal_ns_operator_is_exactly_zero_outside_fixed_spatial_support():
    x, y, z, t = 0.30, 0.02, 0.01, 0.6
    certificate = assemble_section9_section10_formal_ns_operator(
        _prefix_at(x, y, z, t),
        x=x,
        y=y,
        z=z,
        t=t,
        viscosity=0.7,
    )
    np.testing.assert_array_equal(certificate.formal_operator_sum, np.zeros(3))
    np.testing.assert_array_equal(
        certificate.time_derivative.localized_velocity_time_derivative,
        np.zeros(3),
    )
    np.testing.assert_array_equal(
        certificate.convection.convective_acceleration,
        np.zeros(3),
    )
    np.testing.assert_array_equal(
        certificate.laplacian.localized_velocity_laplacian,
        np.zeros(3),
    )
    np.testing.assert_array_equal(
        certificate.pressure.localized_pressure_gradient,
        np.zeros(3),
    )


def test_formal_ns_operator_keeps_truth_and_derivative_order_gates_closed():
    x, y, z, t = 0.20, 0.02, 0.16, 0.6
    with pytest.raises(ValueError, match="at least three"):
        assemble_section9_section10_formal_ns_operator(
            _prefix_at(x, y, z, t, order=2),
            x=x,
            y=y,
            z=z,
            t=t,
        )

    with pytest.raises(ValueError, match="nonnegative"):
        assemble_section9_section10_formal_ns_operator(
            _prefix_at(x, y, z, t),
            x=x,
            y=y,
            z=z,
            t=t,
            viscosity=-0.1,
        )

    with pytest.raises(ValueError, match="truth gate"):
        assemble_section9_section10_formal_ns_operator(
            _prefix_at(x, y, z, t, paper_exact_velocity_available=True),
            x=x,
            y=y,
            z=z,
            t=t,
        )
