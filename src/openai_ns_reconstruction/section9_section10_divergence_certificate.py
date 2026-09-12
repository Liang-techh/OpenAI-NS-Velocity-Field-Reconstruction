"""Paper-specific divergence check for the finite Section 9 -> Section 10 bridge.

Section 10 uses

    u_cut = curl(c A) + c B e_theta,

where the fixed spatial cutoff ``c`` is axisymmetric.  Away from the symmetry
axis, ``div(c B e_theta)`` reduces to ``c e_theta . grad(B)`` because
``e_theta . grad(c)=0`` and ``div(e_theta)=0``.  Thus the direct swirl remains
divergence free exactly when the supplied Section 9 scalar jet satisfies the
axisymmetric first-derivative relation.

This module evaluates that relation against the already-landed analytic spatial
Jacobian.  It is deliberately a finite-prefix, one-point certificate boundary:
it does not manufacture missing Section 7/8 corrections, an infinite Eq. (9.21)
field, a t=1 extension, a Navier--Stokes residual, forcing, finite-energy proof,
or paper-exact velocity.  Floating comparisons are explicit and therefore do
not promote a numerical pass to a theorem-level claim.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .section9_eq921_prefix_jet import Section9FinitePrefixJetCertificate
from .section9_section10_localized_velocity_spatial_jacobian import (
    Section9Section10LocalizedVelocitySpatialJacobianCertificate,
    localize_section9_prefix_spatial_jacobian,
)


@dataclass(frozen=True)
class Section9Section10DivergenceCertificate:
    """One-point analytic divergence ledger with fail-closed paper truth flags."""

    spatial: Section9Section10LocalizedVelocitySpatialJacobianCertificate
    tolerance: float
    e_theta: np.ndarray
    cutoff_azimuthal_derivative: float
    B_azimuthal_derivative: float
    poloidal_divergence: float
    localized_swirl_divergence: float
    expected_localized_swirl_divergence: float
    localized_divergence: float
    analytic_trace_consistency_error: float
    support_exterior_exact_zero: bool
    fixed_cutoff_axisymmetry_at_point_verified: bool
    prefix_swirl_axisymmetry_at_point_verified: bool
    analytic_trace_consistency_verified: bool
    divergence_free_at_point_verified: bool
    axis_covered: bool = False
    actual_section9_sequence_verified: bool = False
    eq_9_21_infinite_sum_constructed: bool = False
    section9_field_smooth_extension_through_t1_constructed: bool = False
    residual_artifact_ready: bool = False
    endpoint_residual_closure_verified: bool = False
    finite_energy_verified: bool = False
    blowup_path_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_divergence_evaluation_ready(self) -> bool:
        """Whether the finite-prefix divergence evaluation itself is admissible.

        A non-axisymmetric supplied ``B`` jet may legitimately make
        ``divergence_free_at_point_verified`` false; that is a useful failed
        certificate rather than a malformed evaluation.  The formal gate thus
        checks the analytic bookkeeping and truth ceiling, not the PDE claim.
        """
        return (
            self.spatial.formal_localized_spatial_jacobian_ready
            and self.fixed_cutoff_axisymmetry_at_point_verified
            and self.analytic_trace_consistency_verified
            and not self.axis_covered
            and not self.actual_section9_sequence_verified
            and not self.eq_9_21_infinite_sum_constructed
            and not self.section9_field_smooth_extension_through_t1_constructed
            and not self.residual_artifact_ready
            and not self.endpoint_residual_closure_verified
            and not self.finite_energy_verified
            and not self.blowup_path_verified
            and not self.paper_exact_velocity_available
        )


def _finite_nonnegative(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite nonnegative real")
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"{name} must be a finite nonnegative real")
    return out


def certify_section9_section10_divergence(
    prefix: Section9FinitePrefixJetCertificate,
    *,
    x: float,
    y: float,
    z: float,
    t: float,
    tolerance: float = 1.0e-11,
) -> Section9Section10DivergenceCertificate:
    """Evaluate the off-axis finite-prefix divergence identity at one point.

    The production path uses only analytic jets already present in ``prefix``
    and the fixed Section 10 cutoff derivative path.  It performs no finite
    differencing.  ``tolerance`` controls only floating identity comparisons;
    it is not an analytic-error or convergence bound.
    """
    tol = _finite_nonnegative(tolerance, "tolerance")
    spatial = localize_section9_prefix_spatial_jacobian(
        prefix, x=x, y=y, z=z, t=t
    )
    x0, y0, _z0, _t0 = spatial.base.point
    r = math.hypot(x0, y0)
    if r == 0.0:
        # The spatial-Jacobian constructor already rejects the axis, but keep
        # the truth boundary explicit here as well.
        raise ValueError(
            "divergence certificate requires an off-axis point; axis regularity "
            "must come from the actual paper field"
        )

    e_theta = np.array([-y0 / r, x0 / r, 0.0], dtype=float)
    e_theta.setflags(write=False)

    cutoff_azimuthal = float(e_theta @ spatial.base.cutoff_gradient)
    B_azimuthal = float(e_theta @ spatial.B_spatial_gradient)

    swirl = spatial.base.local_velocity - spatial.base.curl_A
    localized_swirl_jacobian = (
        np.outer(swirl, spatial.base.cutoff_gradient)
        + spatial.base.cutoff_value * spatial.swirl_spatial_jacobian
    )
    poloidal_jacobian = (
        spatial.localized_velocity_spatial_jacobian - localized_swirl_jacobian
    )

    poloidal_divergence = float(np.trace(poloidal_jacobian))
    localized_swirl_divergence = float(np.trace(localized_swirl_jacobian))
    expected_swirl_divergence = float(
        spatial.base.B_value * cutoff_azimuthal
        + spatial.base.cutoff_value * B_azimuthal
    )
    localized_divergence = float(
        np.trace(spatial.localized_velocity_spatial_jacobian)
    )
    consistency_error = abs(
        localized_divergence
        - (poloidal_divergence + expected_swirl_divergence)
    )

    exterior_zero = bool(
        spatial.base.cutoff_value == 0.0
        and np.array_equal(spatial.base.cutoff_gradient, np.zeros(3))
        and np.array_equal(spatial.cutoff_hessian, np.zeros((3, 3)))
        and np.array_equal(
            spatial.localized_velocity_spatial_jacobian, np.zeros((3, 3))
        )
    )
    cutoff_axisymmetric = abs(cutoff_azimuthal) <= tol
    B_axisymmetric = abs(B_azimuthal) <= tol
    trace_consistent = (
        abs(poloidal_divergence) <= tol
        and abs(localized_swirl_divergence - expected_swirl_divergence) <= tol
        and consistency_error <= tol
    )
    divergence_free = exterior_zero or (
        cutoff_axisymmetric
        and B_axisymmetric
        and trace_consistent
        and abs(localized_divergence) <= tol
    )

    certificate = Section9Section10DivergenceCertificate(
        spatial=spatial,
        tolerance=tol,
        e_theta=e_theta,
        cutoff_azimuthal_derivative=cutoff_azimuthal,
        B_azimuthal_derivative=B_azimuthal,
        poloidal_divergence=poloidal_divergence,
        localized_swirl_divergence=localized_swirl_divergence,
        expected_localized_swirl_divergence=expected_swirl_divergence,
        localized_divergence=localized_divergence,
        analytic_trace_consistency_error=consistency_error,
        support_exterior_exact_zero=exterior_zero,
        fixed_cutoff_axisymmetry_at_point_verified=cutoff_axisymmetric,
        prefix_swirl_axisymmetry_at_point_verified=B_axisymmetric,
        analytic_trace_consistency_verified=trace_consistent,
        divergence_free_at_point_verified=divergence_free,
    )
    if not certificate.formal_divergence_evaluation_ready:
        raise ValueError("finite-prefix divergence certificate failed its truth gate")
    return certificate
