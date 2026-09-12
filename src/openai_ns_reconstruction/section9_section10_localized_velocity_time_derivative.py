"""Analytic time derivative for the finite-prefix Section 9 -> Section 10 bridge.

This module advances the already-landed one-point localized-velocity adapter by
one derivative needed by the Navier--Stokes residual: ``d_t u_cut``.  For the
fixed Section 10 spatial cutoff ``c(x,y,z)`` the cutoff is time independent, so

    d_t u_cut
      = c curl(d_t A) + grad(c) x d_t A + c (d_t B) e_theta.

The input remains a finite, provider-supplied ``Section9FinitePrefixJetCertificate``.
At least total derivative order two is required because ``curl(d_t A)`` uses the
mixed derivatives ``D_t D_x A``, ``D_t D_y A``, and ``D_t D_z A``.  Production
code performs no finite differencing of the Section 9 fields or of the cutoff.

This is not a residual or forcing constructor.  It does not prove that the
finite prefix is the manuscript correction sequence, construct the infinite
Eq. (9.21) sum, extend it smoothly through ``t=1``, or promote the velocity to
paper-exact status.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .local_field import azimuthal_vector
from .section9_eq921_prefix_jet import Section9FinitePrefixJetCertificate
from .section9_section10_localized_velocity_jet import (
    Section9Section10LocalizedVelocityJetCertificate,
    localize_section9_prefix_jet,
)

_DT = (1, 0, 0, 0)
_DTDX = (1, 1, 0, 0)
_DTDY = (1, 0, 1, 0)
_DTDZ = (1, 0, 0, 1)


def _readonly_vec3(value: object, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (3,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite three-vector")
    frozen = np.array(out, dtype=float, copy=True)
    frozen.setflags(write=False)
    return frozen


def _finite_scalar(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real scalar")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


def _curl_time_derivative_from_prefix_jet(
    prefix: Section9FinitePrefixJetCertificate,
) -> np.ndarray:
    """Return ``d_t curl(A)=curl(d_t A)`` from mixed analytic derivatives."""
    dAt_dx = _readonly_vec3(prefix.A_derivatives[_DTDX], "D_t D_x A")
    dAt_dy = _readonly_vec3(prefix.A_derivatives[_DTDY], "D_t D_y A")
    dAt_dz = _readonly_vec3(prefix.A_derivatives[_DTDZ], "D_t D_z A")
    curl_t = np.array(
        [
            dAt_dy[2] - dAt_dz[1],
            dAt_dz[0] - dAt_dx[2],
            dAt_dx[1] - dAt_dy[0],
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(curl_t)):
        raise ArithmeticError("analytic time derivative of curl(A) is not finite")
    curl_t.setflags(write=False)
    return curl_t


@dataclass(frozen=True)
class Section9Section10LocalizedVelocityTimeDerivativeCertificate:
    """One-point analytic ``d_t u_cut`` with the current fail-closed boundary."""

    base: Section9Section10LocalizedVelocityJetCertificate
    A_time_derivative: np.ndarray
    B_time_derivative: float
    curl_A_time_derivative: np.ndarray
    local_velocity_time_derivative: np.ndarray
    localized_velocity_time_derivative: np.ndarray
    analytic_time_derivative_from_prefix_jet_verified: bool = True
    section10_spatial_cutoff_time_independent_verified: bool = True
    actual_section9_sequence_verified: bool = False
    eq_9_21_infinite_sum_constructed: bool = False
    section9_field_smooth_extension_through_t1_constructed: bool = False
    residual_artifact_ready: bool = False
    endpoint_residual_closure_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_localized_time_derivative_ready(self) -> bool:
        return (
            self.base.formal_localized_velocity_ready
            and self.base.derivative_order >= 2
            and self.analytic_time_derivative_from_prefix_jet_verified
            and self.section10_spatial_cutoff_time_independent_verified
            and not self.actual_section9_sequence_verified
            and not self.eq_9_21_infinite_sum_constructed
            and not self.section9_field_smooth_extension_through_t1_constructed
            and not self.residual_artifact_ready
            and not self.endpoint_residual_closure_verified
            and not self.paper_exact_velocity_available
        )


def localize_section9_prefix_time_derivative(
    prefix: Section9FinitePrefixJetCertificate,
    *,
    x: float,
    y: float,
    z: float,
    t: float,
) -> Section9Section10LocalizedVelocityTimeDerivativeCertificate:
    """Evaluate the analytic time derivative of the fixed-cutoff localized field.

    The same fixed Section 10 cutoff/gradient pair used by
    ``localize_section9_prefix_jet`` is reused through the returned base
    certificate.  No caller-supplied temporal cutoff is accepted here: the
    official Section 10 spatial cutoff is time independent, and the separately
    landed paper ``timeSwitch`` is plateaued at ``t=1`` rather than being an
    arbitrary window for this local algebra.
    """
    if not isinstance(prefix, Section9FinitePrefixJetCertificate):
        raise TypeError("prefix must be a Section9FinitePrefixJetCertificate")
    if not prefix.formal_prefix_jet_ready:
        raise ValueError("prefix must pass the finite-prefix jet truth gate")
    if prefix.derivative_order < 2:
        raise ValueError(
            "prefix derivative_order must be at least two to form d_t curl(A)"
        )

    base = localize_section9_prefix_jet(prefix, x=x, y=y, z=z, t=t)
    A_t = _readonly_vec3(prefix.A_derivatives[_DT], "D_t A")
    B_t = _finite_scalar(prefix.B_derivatives[_DT], "D_t B")
    curl_A_t = _curl_time_derivative_from_prefix_jet(prefix)

    x0, y0, _z0, _t0 = base.point
    swirl_t = _readonly_vec3(azimuthal_vector(B_t, x0, y0), "D_t(B e_theta)")
    local_t = np.array(curl_A_t + swirl_t, dtype=float)
    if not np.all(np.isfinite(local_t)):
        raise ArithmeticError("local velocity time derivative is not finite")
    local_t.setflags(write=False)

    localized_t = np.array(
        base.cutoff_value * curl_A_t
        + np.cross(base.cutoff_gradient, A_t)
        + base.cutoff_value * swirl_t,
        dtype=float,
    )
    if not np.all(np.isfinite(localized_t)):
        raise ArithmeticError("localized velocity time derivative is not finite")
    localized_t.setflags(write=False)

    certificate = Section9Section10LocalizedVelocityTimeDerivativeCertificate(
        base=base,
        A_time_derivative=A_t,
        B_time_derivative=B_t,
        curl_A_time_derivative=curl_A_t,
        local_velocity_time_derivative=local_t,
        localized_velocity_time_derivative=localized_t,
    )
    if not certificate.formal_localized_time_derivative_ready:
        raise ValueError("localized time-derivative certificate failed its truth gate")
    return certificate
