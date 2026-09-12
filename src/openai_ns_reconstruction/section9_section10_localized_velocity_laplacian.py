"""Analytic off-axis Laplacian for the finite-prefix Section 9 -> Section 10 bridge.

This module is a bounded derivative adapter only.  It differentiates the
already-landed finite-prefix localization

    u_cut = c (curl(A) + B e_theta) + grad(c) x A

through one spatial Laplacian, using the fixed Section 10 cutoff geometry and
its analytic gradient/Hessian/Delta/grad(Delta) path.  A complete finite
Section 9 spacetime jet through total order three is required because
``Delta curl(A) = curl(Delta A)`` uses third spatial derivatives of ``A``.

Production code uses no finite differences.  The symmetry axis remains
fail-closed because Cartesian derivatives of ``B e_theta`` at ``r=0`` require
the actual paper field's axis-regularity structure.  The input is still a
provider-supplied finite Eq. (9.21) prefix, so this file does not construct the
actual Section 7/8 correction sequence, the infinite/local-finite Eq. (9.21)
field, the smooth extension through ``t=1``, a genuine Navier--Stokes residual
or forcing, finite-energy closure, blow-up closure, or paper-exact velocity.
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
from .section10_spatial_cutoff_viscous_adapter import (
    section10_spatial_cutoff_laplacian,
    section10_spatial_cutoff_laplacian_gradient,
)

_DXX = (0, 2, 0, 0)
_DYY = (0, 0, 2, 0)
_DZZ = (0, 0, 0, 2)
_DXXX = (0, 3, 0, 0)
_DXXY = (0, 2, 1, 0)
_DXXZ = (0, 2, 0, 1)
_DXYY = (0, 1, 2, 0)
_DXZZ = (0, 1, 0, 2)
_DYYY = (0, 0, 3, 0)
_DYYZ = (0, 0, 2, 1)
_DYZZ = (0, 0, 1, 2)
_DZZZ = (0, 0, 0, 3)


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


def _A_laplacian(prefix: Section9FinitePrefixJetCertificate) -> np.ndarray:
    return _readonly_vec3(
        np.asarray(prefix.A_derivatives[_DXX], dtype=float)
        + np.asarray(prefix.A_derivatives[_DYY], dtype=float)
        + np.asarray(prefix.A_derivatives[_DZZ], dtype=float),
        "Delta A",
    )


def _curl_A_laplacian(prefix: Section9FinitePrefixJetCertificate) -> np.ndarray:
    """Return ``Delta curl(A)=curl(Delta A)`` from third spatial A jets."""
    d_dx_lap_A = (
        np.asarray(prefix.A_derivatives[_DXXX], dtype=float)
        + np.asarray(prefix.A_derivatives[_DXYY], dtype=float)
        + np.asarray(prefix.A_derivatives[_DXZZ], dtype=float)
    )
    d_dy_lap_A = (
        np.asarray(prefix.A_derivatives[_DXXY], dtype=float)
        + np.asarray(prefix.A_derivatives[_DYYY], dtype=float)
        + np.asarray(prefix.A_derivatives[_DYZZ], dtype=float)
    )
    d_dz_lap_A = (
        np.asarray(prefix.A_derivatives[_DXXZ], dtype=float)
        + np.asarray(prefix.A_derivatives[_DYYZ], dtype=float)
        + np.asarray(prefix.A_derivatives[_DZZZ], dtype=float)
    )
    return _readonly_vec3(
        np.array(
            [
                d_dy_lap_A[2] - d_dz_lap_A[1],
                d_dz_lap_A[0] - d_dx_lap_A[2],
                d_dx_lap_A[1] - d_dy_lap_A[0],
            ],
            dtype=float,
        ),
        "Delta curl(A)",
    )


def _B_laplacian(prefix: Section9FinitePrefixJetCertificate) -> float:
    return _finite_scalar(
        prefix.B_derivatives[_DXX]
        + prefix.B_derivatives[_DYY]
        + prefix.B_derivatives[_DZZ],
        "Delta B",
    )


def _swirl_laplacian(
    *,
    x: float,
    y: float,
    B_value: float,
    B_gradient: np.ndarray,
    B_laplacian: float,
) -> np.ndarray:
    """Return ``Delta(B e_theta)`` in Cartesian coordinates for ``r>0``."""
    r = math.hypot(x, y)
    if r == 0.0:
        raise ValueError(
            "Laplacian of B e_theta requires an off-axis point; "
            "axis regularity must be supplied separately"
        )
    r3 = r * r * r
    e_theta = np.array([-y / r, x / r, 0.0], dtype=float)
    basis_jacobian = np.array(
        [
            [x * y / r3, -(x * x) / r3, 0.0],
            [(y * y) / r3, -x * y / r3, 0.0],
            [0.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    basis_laplacian = -e_theta / (r * r)
    value = (
        B_laplacian * e_theta
        + 2.0
        * (
            B_gradient[0] * basis_jacobian[:, 0]
            + B_gradient[1] * basis_jacobian[:, 1]
        )
        + B_value * basis_laplacian
    )
    return _readonly_vec3(value, "Delta(B e_theta)")


@dataclass(frozen=True)
class Section9Section10LocalizedVelocityLaplacianCertificate:
    """One-point off-axis analytic ``Delta u_cut`` with fail-closed truth flags."""

    spatial: Section9Section10LocalizedVelocitySpatialJacobianCertificate
    cutoff_laplacian: float
    cutoff_laplacian_gradient: np.ndarray
    A_laplacian: np.ndarray
    curl_A_laplacian: np.ndarray
    B_laplacian: float
    swirl_laplacian: np.ndarray
    local_velocity_laplacian: np.ndarray
    localized_velocity_laplacian: np.ndarray
    analytic_laplacian_from_prefix_jet_verified: bool = True
    section10_fixed_cutoff_third_order_applied: bool = True
    azimuthal_basis_laplacian_off_axis_verified: bool = True
    axis_laplacian_covered: bool = False
    actual_section9_sequence_verified: bool = False
    eq_9_21_infinite_sum_constructed: bool = False
    section9_field_smooth_extension_through_t1_constructed: bool = False
    residual_artifact_ready: bool = False
    endpoint_residual_closure_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_localized_laplacian_ready(self) -> bool:
        return (
            self.spatial.formal_localized_spatial_jacobian_ready
            and self.spatial.base.derivative_order >= 3
            and self.analytic_laplacian_from_prefix_jet_verified
            and self.section10_fixed_cutoff_third_order_applied
            and self.azimuthal_basis_laplacian_off_axis_verified
            and not self.axis_laplacian_covered
            and not self.actual_section9_sequence_verified
            and not self.eq_9_21_infinite_sum_constructed
            and not self.section9_field_smooth_extension_through_t1_constructed
            and not self.residual_artifact_ready
            and not self.endpoint_residual_closure_verified
            and not self.paper_exact_velocity_available
        )


def localize_section9_prefix_laplacian(
    prefix: Section9FinitePrefixJetCertificate,
    *,
    x: float,
    y: float,
    z: float,
    t: float,
) -> Section9Section10LocalizedVelocityLaplacianCertificate:
    """Evaluate the analytic off-axis Laplacian of the localized finite prefix.

    With ``V=curl(A)+B e_theta``, ``g=grad(c)``, ``H=Hess(c)``, and
    ``L=Delta c``, the identity used is

    ``Delta u_cut = L V + 2 grad(V) g + c Delta V
                     + grad(L) x A
                     + 2 sum_j H[:,j] x partial_j A
                     + g x Delta A``.

    The fixed cutoff derivatives are not caller-overridable.  The result is a
    reusable formal-structure derivative certificate, not a genuine residual or
    forcing artifact.
    """
    if not isinstance(prefix, Section9FinitePrefixJetCertificate):
        raise TypeError("prefix must be a Section9FinitePrefixJetCertificate")
    if not prefix.formal_prefix_jet_ready:
        raise ValueError("prefix must pass the finite-prefix jet truth gate")
    if prefix.derivative_order < 3:
        raise ValueError(
            "prefix derivative_order must be at least three to form Delta curl(A)"
        )

    spatial = localize_section9_prefix_spatial_jacobian(
        prefix, x=x, y=y, z=z, t=t
    )
    base = spatial.base
    x0, y0, z0, t0 = base.point

    A_laplacian = _A_laplacian(prefix)
    curl_A_laplacian = _curl_A_laplacian(prefix)
    B_laplacian = _B_laplacian(prefix)
    swirl_laplacian = _swirl_laplacian(
        x=x0,
        y=y0,
        B_value=base.B_value,
        B_gradient=spatial.B_spatial_gradient,
        B_laplacian=B_laplacian,
    )
    local_laplacian = _readonly_vec3(
        curl_A_laplacian + swirl_laplacian,
        "Delta local velocity",
    )

    cutoff_laplacian = _finite_scalar(
        section10_spatial_cutoff_laplacian(x0, y0, z0, t0),
        "Section 10 cutoff Laplacian",
    )
    cutoff_laplacian_gradient = _readonly_vec3(
        section10_spatial_cutoff_laplacian_gradient(x0, y0, z0, t0),
        "Section 10 grad(Delta c)",
    )

    mixed_cross = np.zeros(3, dtype=float)
    for axis in range(3):
        mixed_cross += np.cross(
            spatial.cutoff_hessian[:, axis],
            spatial.A_spatial_jacobian[:, axis],
        )

    localized_laplacian = _readonly_vec3(
        cutoff_laplacian * base.local_velocity
        + 2.0
        * (spatial.local_velocity_spatial_jacobian @ base.cutoff_gradient)
        + base.cutoff_value * local_laplacian
        + np.cross(cutoff_laplacian_gradient, base.A_value)
        + 2.0 * mixed_cross
        + np.cross(base.cutoff_gradient, A_laplacian),
        "Delta localized velocity",
    )

    certificate = Section9Section10LocalizedVelocityLaplacianCertificate(
        spatial=spatial,
        cutoff_laplacian=cutoff_laplacian,
        cutoff_laplacian_gradient=cutoff_laplacian_gradient,
        A_laplacian=A_laplacian,
        curl_A_laplacian=curl_A_laplacian,
        B_laplacian=B_laplacian,
        swirl_laplacian=swirl_laplacian,
        local_velocity_laplacian=local_laplacian,
        localized_velocity_laplacian=localized_laplacian,
    )
    if not certificate.formal_localized_laplacian_ready:
        raise ValueError("localized Laplacian certificate failed its truth gate")
    return certificate
