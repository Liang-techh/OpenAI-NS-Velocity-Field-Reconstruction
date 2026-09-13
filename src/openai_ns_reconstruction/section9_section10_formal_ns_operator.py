"""Compose the landed finite-prefix Section 10 Navier--Stokes operator terms.

This module is deliberately a *formal operator adapter*, not a residual or
forcing constructor.  It takes one admitted finite Section 9 Eq. (9.21) prefix
jet and, at one off-axis spacetime point, combines the already-landed analytic
Section 10 pieces

    u_t,
    (u . grad) u,
    Delta u,
    grad(c p)

with the repository sign convention

    u_t + (u . grad) u - nu Delta u + grad(c p).

All four pieces are derived atomically from the same prefix certificate, point,
and fixed Section 10 spatial cutoff path.  The adapter checks their source and
point identity before summing them; callers cannot inject precomputed terms or
an alternate cutoff.

The result is intentionally named ``formal_operator_sum``.  The input is still
a provider-supplied finite prefix and the axis-regularity path is still absent.
This file therefore does not establish the actual Section 7/8 correction
sequence, the infinite/locally-finite Eq. (9.21) field, the all-order smooth
extension through t=1, a genuine Navier--Stokes residual/forcing artifact,
finite-energy closure, blow-up closure, or paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .section9_eq921_prefix_jet import Section9FinitePrefixJetCertificate
from .section9_section10_convective_acceleration import (
    Section9Section10ConvectiveAccelerationCertificate,
    localize_section9_prefix_convective_acceleration,
)
from .section9_section10_localized_pressure_gradient import (
    Section9Section10LocalizedPressureGradientCertificate,
    localize_section9_prefix_pressure_gradient,
)
from .section9_section10_localized_velocity_laplacian import (
    Section9Section10LocalizedVelocityLaplacianCertificate,
    localize_section9_prefix_laplacian,
)
from .section9_section10_localized_velocity_time_derivative import (
    Section9Section10LocalizedVelocityTimeDerivativeCertificate,
    localize_section9_prefix_time_derivative,
)


def _nonnegative_finite_scalar(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a nonnegative finite real scalar")
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"{name} must be a nonnegative finite real scalar")
    return out


def _readonly_vec3(value: object, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (3,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite three-vector")
    frozen = np.array(out, dtype=float, copy=True)
    frozen.setflags(write=False)
    return frozen


def _base_identity(base: object) -> tuple[object, ...]:
    """Return the identity fields shared by localized velocity certificates."""
    return (
        base.point,
        base.q,
        base.prefix_order,
        base.derivative_order,
        base.provider_id,
        base.provider_revision,
        base.provider_provenance,
    )


def _pressure_identity(
    pressure: Section9Section10LocalizedPressureGradientCertificate,
) -> tuple[object, ...]:
    return (
        pressure.point,
        pressure.q,
        pressure.prefix_order,
        pressure.derivative_order,
        pressure.provider_id,
        pressure.provider_revision,
        pressure.provider_provenance,
    )


@dataclass(frozen=True)
class Section9Section10FormalNSOperatorCertificate:
    """One-point finite-prefix formal NS operator sum, explicitly non-residual."""

    viscosity: float
    time_derivative: Section9Section10LocalizedVelocityTimeDerivativeCertificate
    convection: Section9Section10ConvectiveAccelerationCertificate
    laplacian: Section9Section10LocalizedVelocityLaplacianCertificate
    pressure: Section9Section10LocalizedPressureGradientCertificate
    viscous_term: np.ndarray
    formal_operator_sum: np.ndarray
    component_source_identity_verified: bool = True
    fixed_section10_cutoff_coherent: bool = True
    formal_ns_operator_combination_verified: bool = True
    actual_section9_sequence_verified: bool = False
    eq_9_21_infinite_sum_constructed: bool = False
    section9_field_smooth_extension_through_t1_constructed: bool = False
    residual_artifact_ready: bool = False
    forcing_artifact_ready: bool = False
    endpoint_residual_closure_verified: bool = False
    finite_energy_closure_verified: bool = False
    blow_up_closure_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_ns_operator_ready(self) -> bool:
        return (
            self.time_derivative.formal_localized_time_derivative_ready
            and self.convection.formal_convective_acceleration_ready
            and self.laplacian.formal_localized_laplacian_ready
            and self.pressure.formal_localized_pressure_gradient_ready
            and self.laplacian.spatial.base.derivative_order >= 3
            and self.component_source_identity_verified
            and self.fixed_section10_cutoff_coherent
            and self.formal_ns_operator_combination_verified
            and not self.actual_section9_sequence_verified
            and not self.eq_9_21_infinite_sum_constructed
            and not self.section9_field_smooth_extension_through_t1_constructed
            and not self.residual_artifact_ready
            and not self.forcing_artifact_ready
            and not self.endpoint_residual_closure_verified
            and not self.finite_energy_closure_verified
            and not self.blow_up_closure_verified
            and not self.paper_exact_velocity_available
        )


def assemble_section9_section10_formal_ns_operator(
    prefix: Section9FinitePrefixJetCertificate,
    *,
    x: float,
    y: float,
    z: float,
    t: float,
    viscosity: float = 1.0,
) -> Section9Section10FormalNSOperatorCertificate:
    """Compose the analytic finite-prefix NS operator terms at one point.

    This function deliberately stops one semantic level before a residual.  A
    genuine residual artifact may only be produced after the actual Section 7/8
    correction sequence, locally-finite Eq. (9.21) construction, endpoint
    extension, and source/provenance gates have been established.
    """
    if not isinstance(prefix, Section9FinitePrefixJetCertificate):
        raise TypeError("prefix must be a Section9FinitePrefixJetCertificate")
    if not prefix.formal_prefix_jet_ready:
        raise ValueError("prefix must pass the finite-prefix jet truth gate")
    if prefix.derivative_order < 3:
        raise ValueError(
            "prefix derivative_order must be at least three to compose the NS operator"
        )
    nu = _nonnegative_finite_scalar(viscosity, "viscosity")

    time_derivative = localize_section9_prefix_time_derivative(
        prefix, x=x, y=y, z=z, t=t
    )
    convection = localize_section9_prefix_convective_acceleration(
        prefix, x=x, y=y, z=z, t=t
    )
    laplacian = localize_section9_prefix_laplacian(
        prefix, x=x, y=y, z=z, t=t
    )
    pressure = localize_section9_prefix_pressure_gradient(
        prefix, x=x, y=y, z=z, t=t
    )

    expected_identity = _base_identity(time_derivative.base)
    component_identities = (
        _base_identity(convection.spatial.base),
        _base_identity(laplacian.spatial.base),
        _pressure_identity(pressure),
    )
    if any(identity != expected_identity for identity in component_identities):
        raise ValueError(
            "localized NS operator components do not share one point/source identity"
        )

    base = time_derivative.base
    cutoff_pairs = (
        (convection.spatial.base.cutoff_value, convection.spatial.base.cutoff_gradient),
        (laplacian.spatial.base.cutoff_value, laplacian.spatial.base.cutoff_gradient),
        (pressure.cutoff_value, pressure.cutoff_gradient),
    )
    for cutoff_value, cutoff_gradient in cutoff_pairs:
        if cutoff_value != base.cutoff_value or not np.array_equal(
            cutoff_gradient, base.cutoff_gradient
        ):
            raise ValueError(
                "localized NS operator components do not share the fixed Section 10 cutoff"
            )

    viscous_term = _readonly_vec3(
        -nu * laplacian.localized_velocity_laplacian,
        "formal viscous term",
    )
    formal_operator_sum = _readonly_vec3(
        time_derivative.localized_velocity_time_derivative
        + convection.convective_acceleration
        + viscous_term
        + pressure.localized_pressure_gradient,
        "formal NS operator sum",
    )

    certificate = Section9Section10FormalNSOperatorCertificate(
        viscosity=nu,
        time_derivative=time_derivative,
        convection=convection,
        laplacian=laplacian,
        pressure=pressure,
        viscous_term=viscous_term,
        formal_operator_sum=formal_operator_sum,
    )
    if not certificate.formal_ns_operator_ready:
        raise ValueError("formal NS operator certificate failed its truth gate")
    return certificate
