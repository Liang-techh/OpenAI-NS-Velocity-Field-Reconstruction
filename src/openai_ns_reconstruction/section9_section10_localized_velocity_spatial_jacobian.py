"""Analytic spatial Jacobian for the finite-prefix Section 9 -> Section 10 bridge.

Paper/formal location
---------------------
Section 10 localizes the Section 9 field by Eq. (10.4),

    u_cut = curl(c A) + c B e_theta
          = c curl(A) + grad(c) x A + c B e_theta,

with the fixed spatial cutoff geometry pinned by the official
``NavierStokes/SpatialLocalization.lean`` source.  This module differentiates
that already-landed algebra once in Cartesian space.  It consumes the existing
finite ``Section9FinitePrefixJetCertificate`` and the exact analytic gradient
and Hessian of this repository's geometry-faithful executable cutoff
representative.

At least total derivative order two is required because ``grad curl(A)`` uses
second spatial derivatives of ``A``.  Away from the symmetry axis the direct
swirl derivative is computed analytically from ``B`` and ``grad B`` together
with the Cartesian derivatives of ``e_theta``.  The axis is deliberately
fail-closed here: pointwise derivatives of ``B e_theta`` at ``r=0`` require the
separate smooth axis-regularity structure of the actual paper field and cannot
be inferred from a scalar ``B`` jet alone.

Production code performs no finite differencing.  This is only a reusable
finite-prefix derivative adapter toward the convection term; it is not a
Navier--Stokes residual/forcing constructor and does not establish the actual
Section 7/8 correction sequence, infinite Eq. (9.21), the smooth extension
through ``t=1``, bounded energy, blow-up closure, or paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .section9_eq921_prefix_jet import Section9FinitePrefixJetCertificate
from .section9_section10_localized_velocity_jet import (
    Section9Section10LocalizedVelocityJetCertificate,
    localize_section9_prefix_jet,
)
from .spatial_localization import section10_spatial_cutoff_hessian

_DX = (0, 1, 0, 0)
_DY = (0, 0, 1, 0)
_DZ = (0, 0, 0, 1)
_DXX = (0, 2, 0, 0)
_DXY = (0, 1, 1, 0)
_DXZ = (0, 1, 0, 1)
_DYY = (0, 0, 2, 0)
_DYZ = (0, 0, 1, 1)
_DZZ = (0, 0, 0, 2)


def _readonly_vec3(value: object, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (3,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite three-vector")
    frozen = np.array(out, dtype=float, copy=True)
    frozen.setflags(write=False)
    return frozen


def _readonly_matrix3(value: object, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (3, 3) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite 3x3 matrix")
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


def _A_spatial_jacobian(prefix: Section9FinitePrefixJetCertificate) -> np.ndarray:
    """Return the matrix whose j-th column is ``partial_j A``."""
    jacobian = np.column_stack(
        [
            _readonly_vec3(prefix.A_derivatives[_DX], "D_x A"),
            _readonly_vec3(prefix.A_derivatives[_DY], "D_y A"),
            _readonly_vec3(prefix.A_derivatives[_DZ], "D_z A"),
        ]
    )
    return _readonly_matrix3(jacobian, "spatial Jacobian of A")


def _curl_A_spatial_jacobian(
    prefix: Section9FinitePrefixJetCertificate,
) -> np.ndarray:
    """Return ``grad(curl A)`` from the second Cartesian derivatives of ``A``."""
    A_xx = _readonly_vec3(prefix.A_derivatives[_DXX], "D_xx A")
    A_xy = _readonly_vec3(prefix.A_derivatives[_DXY], "D_xy A")
    A_xz = _readonly_vec3(prefix.A_derivatives[_DXZ], "D_xz A")
    A_yy = _readonly_vec3(prefix.A_derivatives[_DYY], "D_yy A")
    A_yz = _readonly_vec3(prefix.A_derivatives[_DYZ], "D_yz A")
    A_zz = _readonly_vec3(prefix.A_derivatives[_DZZ], "D_zz A")

    d_dx = np.array(
        [
            A_xy[2] - A_xz[1],
            A_xz[0] - A_xx[2],
            A_xx[1] - A_xy[0],
        ],
        dtype=float,
    )
    d_dy = np.array(
        [
            A_yy[2] - A_yz[1],
            A_yz[0] - A_xy[2],
            A_xy[1] - A_yy[0],
        ],
        dtype=float,
    )
    d_dz = np.array(
        [
            A_yz[2] - A_zz[1],
            A_zz[0] - A_xz[2],
            A_xz[1] - A_yz[0],
        ],
        dtype=float,
    )
    return _readonly_matrix3(
        np.column_stack([d_dx, d_dy, d_dz]),
        "spatial Jacobian of curl(A)",
    )


def _azimuthal_basis_and_spatial_jacobian(
    x: float, y: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``e_theta`` and its Cartesian spatial Jacobian for ``r>0``."""
    r = math.hypot(x, y)
    if r == 0.0:
        raise ValueError(
            "spatial derivative of B e_theta requires an off-axis point; "
            "axis regularity must be supplied separately"
        )
    r3 = r * r * r
    e_theta = _readonly_vec3(np.array([-y / r, x / r, 0.0]), "e_theta")
    basis_jacobian = np.array(
        [
            [x * y / r3, -(x * x) / r3, 0.0],
            [(y * y) / r3, -x * y / r3, 0.0],
            [0.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    return e_theta, _readonly_matrix3(
        basis_jacobian, "spatial Jacobian of e_theta"
    )


@dataclass(frozen=True)
class Section9Section10LocalizedVelocitySpatialJacobianCertificate:
    """One-point off-axis analytic ``grad u_cut`` with fail-closed truth flags."""

    base: Section9Section10LocalizedVelocityJetCertificate
    cutoff_hessian: np.ndarray
    A_spatial_jacobian: np.ndarray
    B_spatial_gradient: np.ndarray
    curl_A_spatial_jacobian: np.ndarray
    swirl_spatial_jacobian: np.ndarray
    local_velocity_spatial_jacobian: np.ndarray
    localized_velocity_spatial_jacobian: np.ndarray
    analytic_spatial_jacobian_from_prefix_jet_verified: bool = True
    section10_fixed_cutoff_hessian_applied: bool = True
    azimuthal_basis_derivative_off_axis_verified: bool = True
    axis_spatial_derivative_covered: bool = False
    actual_section9_sequence_verified: bool = False
    eq_9_21_infinite_sum_constructed: bool = False
    section9_field_smooth_extension_through_t1_constructed: bool = False
    residual_artifact_ready: bool = False
    endpoint_residual_closure_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_localized_spatial_jacobian_ready(self) -> bool:
        return (
            self.base.formal_localized_velocity_ready
            and self.base.derivative_order >= 2
            and self.analytic_spatial_jacobian_from_prefix_jet_verified
            and self.section10_fixed_cutoff_hessian_applied
            and self.azimuthal_basis_derivative_off_axis_verified
            and not self.axis_spatial_derivative_covered
            and not self.actual_section9_sequence_verified
            and not self.eq_9_21_infinite_sum_constructed
            and not self.section9_field_smooth_extension_through_t1_constructed
            and not self.residual_artifact_ready
            and not self.endpoint_residual_closure_verified
            and not self.paper_exact_velocity_available
        )


def localize_section9_prefix_spatial_jacobian(
    prefix: Section9FinitePrefixJetCertificate,
    *,
    x: float,
    y: float,
    z: float,
    t: float,
) -> Section9Section10LocalizedVelocitySpatialJacobianCertificate:
    """Evaluate the analytic off-axis spatial Jacobian of the localized velocity.

    With ``C=curl(A)``, ``S=B e_theta``, ``g=grad(c)``, and ``H=Hess(c)``,
    differentiating Eq. (10.4) in coordinate ``x_j`` gives

        partial_j u_cut
          = g_j (C+S) + c partial_j(C+S)
            + H[:,j] x A + g x partial_j A.

    The fixed cutoff, its gradient, and its Hessian are not caller-overridable.
    The input prefix remains finite/provider-supplied and the returned object is
    therefore only a formal-structure derivative certificate.
    """
    if not isinstance(prefix, Section9FinitePrefixJetCertificate):
        raise TypeError("prefix must be a Section9FinitePrefixJetCertificate")
    if not prefix.formal_prefix_jet_ready:
        raise ValueError("prefix must pass the finite-prefix jet truth gate")
    if prefix.derivative_order < 2:
        raise ValueError(
            "prefix derivative_order must be at least two to form grad curl(A)"
        )

    base = localize_section9_prefix_jet(prefix, x=x, y=y, z=z, t=t)
    x0, y0, z0, t0 = base.point
    A_jacobian = _A_spatial_jacobian(prefix)
    B_gradient = _readonly_vec3(
        np.array(
            [
                _finite_scalar(prefix.B_derivatives[_DX], "D_x B"),
                _finite_scalar(prefix.B_derivatives[_DY], "D_y B"),
                _finite_scalar(prefix.B_derivatives[_DZ], "D_z B"),
            ],
            dtype=float,
        ),
        "spatial gradient of B",
    )
    curl_jacobian = _curl_A_spatial_jacobian(prefix)
    cutoff_hessian = _readonly_matrix3(
        section10_spatial_cutoff_hessian(x0, y0, z0, t0),
        "Section 10 cutoff Hessian",
    )

    e_theta, basis_jacobian = _azimuthal_basis_and_spatial_jacobian(x0, y0)
    swirl_jacobian = _readonly_matrix3(
        np.outer(e_theta, B_gradient) + base.B_value * basis_jacobian,
        "spatial Jacobian of B e_theta",
    )
    local_jacobian = _readonly_matrix3(
        curl_jacobian + swirl_jacobian,
        "spatial Jacobian of local velocity",
    )

    localized_jacobian = np.zeros((3, 3), dtype=float)
    for axis in range(3):
        localized_jacobian[:, axis] = (
            base.cutoff_gradient[axis] * base.local_velocity
            + base.cutoff_value * local_jacobian[:, axis]
            + np.cross(cutoff_hessian[:, axis], base.A_value)
            + np.cross(base.cutoff_gradient, A_jacobian[:, axis])
        )
    localized_jacobian = _readonly_matrix3(
        localized_jacobian,
        "spatial Jacobian of localized velocity",
    )

    certificate = Section9Section10LocalizedVelocitySpatialJacobianCertificate(
        base=base,
        cutoff_hessian=cutoff_hessian,
        A_spatial_jacobian=A_jacobian,
        B_spatial_gradient=B_gradient,
        curl_A_spatial_jacobian=curl_jacobian,
        swirl_spatial_jacobian=swirl_jacobian,
        local_velocity_spatial_jacobian=local_jacobian,
        localized_velocity_spatial_jacobian=localized_jacobian,
    )
    if not certificate.formal_localized_spatial_jacobian_ready:
        raise ValueError("localized spatial-Jacobian certificate failed its truth gate")
    return certificate
