"""Analytic convective acceleration for the finite-prefix Section 9 -> Section 10 bridge.

This module is one bounded Navier--Stokes operator adapter.  It combines the
already-landed fixed-Section-10 localized velocity with its analytic Cartesian
spatial Jacobian to evaluate

    (u_cut . grad) u_cut = (grad u_cut) u_cut.

The input remains an admitted, provider-supplied *finite* Eq. (9.21) prefix.
Production code performs no finite differencing and does not assemble a full
Navier--Stokes residual or forcing.  The symmetry axis remains fail-closed
through the underlying spatial-Jacobian adapter pending the actual paper
field's axis-regularity construction.

Accordingly this file does not establish the actual Section 7/8 correction
sequence, the infinite/locally-finite Eq. (9.21) field, all-order smooth
extension through t=1, a genuine residual/forcing closure, finite-energy
control, blow-up closure, or paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .section9_eq921_prefix_jet import Section9FinitePrefixJetCertificate
from .section9_section10_localized_velocity_spatial_jacobian import (
    Section9Section10LocalizedVelocitySpatialJacobianCertificate,
    localize_section9_prefix_spatial_jacobian,
)


def _readonly_vec3(value: object, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (3,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite three-vector")
    frozen = np.array(out, dtype=float, copy=True)
    frozen.setflags(write=False)
    return frozen


@dataclass(frozen=True)
class Section9Section10ConvectiveAccelerationCertificate:
    """One-point off-axis analytic ``(u_cut . grad) u_cut`` certificate."""

    spatial: Section9Section10LocalizedVelocitySpatialJacobianCertificate
    convective_acceleration: np.ndarray
    analytic_convection_from_velocity_and_jacobian_verified: bool = True
    axis_convective_acceleration_covered: bool = False
    actual_section9_sequence_verified: bool = False
    eq_9_21_infinite_sum_constructed: bool = False
    section9_field_smooth_extension_through_t1_constructed: bool = False
    residual_artifact_ready: bool = False
    endpoint_residual_closure_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_convective_acceleration_ready(self) -> bool:
        return (
            self.spatial.formal_localized_spatial_jacobian_ready
            and self.analytic_convection_from_velocity_and_jacobian_verified
            and not self.axis_convective_acceleration_covered
            and not self.actual_section9_sequence_verified
            and not self.eq_9_21_infinite_sum_constructed
            and not self.section9_field_smooth_extension_through_t1_constructed
            and not self.residual_artifact_ready
            and not self.endpoint_residual_closure_verified
            and not self.paper_exact_velocity_available
        )


def localize_section9_prefix_convective_acceleration(
    prefix: Section9FinitePrefixJetCertificate,
    *,
    x: float,
    y: float,
    z: float,
    t: float,
) -> Section9Section10ConvectiveAccelerationCertificate:
    """Evaluate the analytic off-axis convection term for one finite prefix.

    The spatial Jacobian stores ``partial_j u`` in column ``j``.  Therefore
    matrix-vector multiplication by the localized velocity gives exactly

        ((grad u) u)_i = sum_j u_j partial_j u_i.

    No caller-supplied cutoff, cutoff derivative, velocity, or Jacobian is
    accepted: all are obtained atomically from the existing fixed Section 10
    finite-prefix bridge.  At least second-order prefix jets are consequently
    required by the delegated spatial-Jacobian truth gate.
    """
    if not isinstance(prefix, Section9FinitePrefixJetCertificate):
        raise TypeError("prefix must be a Section9FinitePrefixJetCertificate")

    spatial = localize_section9_prefix_spatial_jacobian(
        prefix, x=x, y=y, z=z, t=t
    )
    convection = _readonly_vec3(
        spatial.localized_velocity_spatial_jacobian
        @ spatial.base.localized_velocity,
        "localized convective acceleration",
    )

    certificate = Section9Section10ConvectiveAccelerationCertificate(
        spatial=spatial,
        convective_acceleration=convection,
    )
    if not certificate.formal_convective_acceleration_ready:
        raise ValueError("convective-acceleration certificate failed its truth gate")
    return certificate
