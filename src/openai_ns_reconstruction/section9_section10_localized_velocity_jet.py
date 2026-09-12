"""Finite-order analytic bridge from a Section 9 prefix jet to Section 10 velocity.

This module deliberately stops before any Navier--Stokes residual or forcing is
constructed.  It consumes an already-admitted ``Section9FinitePrefixJetCertificate``
with at least first spatial derivatives and evaluates, at one caller-supplied
Cartesian point, the algebra

    u_loc = curl(A) + B e_theta,
    u_cut = c curl(A) + grad(c) x A + c B e_theta,

using the repository's fixed Section 10 support/plateau representative and its
matching analytic gradient.  Thus the Section 9 derivative plumbing is connected
to the already-landed Section 10 product-rule path without finite-differencing
``A`` or the cutoff in production code.

The point itself is not carried by ``Section9FinitePrefixJetCertificate`` and
therefore remains caller-supplied metadata.  This adapter does not prove that
the supplied finite prefix is the actual manuscript correction sequence, does
not construct the infinite Eq. (9.21) field, does not cover ``t=1``, and does
not produce a residual, forcing, or paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .local_field import azimuthal_vector
from .section9_eq921_prefix_jet import Section9FinitePrefixJetCertificate
from .spatial_localization import (
    section10_spatial_cutoff,
    section10_spatial_cutoff_gradient,
)

_ZERO = (0, 0, 0, 0)
_DX = (0, 1, 0, 0)
_DY = (0, 0, 1, 0)
_DZ = (0, 0, 0, 1)


def _finite_coordinate(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real coordinate")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


def _readonly_vec3(value: object, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (3,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite three-vector")
    frozen = np.array(out, dtype=float, copy=True)
    frozen.setflags(write=False)
    return frozen


def _curl_from_prefix_jet(prefix: Section9FinitePrefixJetCertificate) -> np.ndarray:
    """Read ``curl(A)`` from the first Cartesian derivatives in ``prefix``."""
    dA_dx = _readonly_vec3(prefix.A_derivatives[_DX], "D_x A")
    dA_dy = _readonly_vec3(prefix.A_derivatives[_DY], "D_y A")
    dA_dz = _readonly_vec3(prefix.A_derivatives[_DZ], "D_z A")
    curl = np.array(
        [
            dA_dy[2] - dA_dz[1],
            dA_dz[0] - dA_dx[2],
            dA_dx[1] - dA_dy[0],
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(curl)):
        raise ArithmeticError("analytic curl from Section 9 prefix jet is not finite")
    curl.setflags(write=False)
    return curl


@dataclass(frozen=True)
class Section9Section10LocalizedVelocityJetCertificate:
    """One-point finite-prefix velocity bridge with explicit truth boundaries."""

    point: tuple[float, float, float, float]
    q: object
    prefix_order: int
    derivative_order: int
    provider_id: str
    provider_revision: str
    provider_provenance: str
    A_value: np.ndarray
    B_value: float
    curl_A: np.ndarray
    local_velocity: np.ndarray
    cutoff_value: float
    cutoff_gradient: np.ndarray
    localized_velocity: np.ndarray
    analytic_curl_from_prefix_jet_verified: bool = True
    section10_fixed_cutoff_product_rule_applied: bool = True
    point_binding_to_actual_correction_field_verified: bool = False
    actual_section9_sequence_verified: bool = False
    eq_9_21_infinite_sum_constructed: bool = False
    section9_field_smooth_extension_through_t1_constructed: bool = False
    residual_artifact_ready: bool = False
    endpoint_residual_closure_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_localized_velocity_ready(self) -> bool:
        return (
            self.derivative_order >= 1
            and self.analytic_curl_from_prefix_jet_verified
            and self.section10_fixed_cutoff_product_rule_applied
            and not self.point_binding_to_actual_correction_field_verified
            and not self.actual_section9_sequence_verified
            and not self.eq_9_21_infinite_sum_constructed
            and not self.section9_field_smooth_extension_through_t1_constructed
            and not self.residual_artifact_ready
            and not self.endpoint_residual_closure_verified
            and not self.paper_exact_velocity_available
        )


def localize_section9_prefix_jet(
    prefix: Section9FinitePrefixJetCertificate,
    *,
    x: float,
    y: float,
    z: float,
    t: float,
) -> Section9Section10LocalizedVelocityJetCertificate:
    """Apply the fixed Section 10 velocity localization to one Section 9 prefix jet.

    ``prefix`` must already pass its own finite-prefix truth gate and must carry
    first spatial derivatives.  The function intentionally accepts no generic
    cutoff or caller-supplied cutoff gradient: the fixed Section 10 pair is used
    atomically so the support/plateau geometry cannot drift from the analytic
    product-rule path.
    """
    if not isinstance(prefix, Section9FinitePrefixJetCertificate):
        raise TypeError("prefix must be a Section9FinitePrefixJetCertificate")
    if not prefix.formal_prefix_jet_ready:
        raise ValueError("prefix must pass the finite-prefix jet truth gate")
    if prefix.derivative_order < 1:
        raise ValueError("prefix derivative_order must be at least one to form curl(A)")

    x = _finite_coordinate(x, "x")
    y = _finite_coordinate(y, "y")
    z = _finite_coordinate(z, "z")
    t = _finite_coordinate(t, "t")

    A_value = _readonly_vec3(prefix.A_derivatives[_ZERO], "A")
    B_value = float(prefix.B_derivatives[_ZERO])
    if not math.isfinite(B_value):
        raise ValueError("B must be finite")

    curl_A = _curl_from_prefix_jet(prefix)
    swirl = _readonly_vec3(azimuthal_vector(B_value, x, y), "B e_theta")
    local_velocity = np.array(curl_A + swirl, dtype=float, copy=True)
    local_velocity.setflags(write=False)

    cutoff_value = float(section10_spatial_cutoff(x, y, z, t))
    cutoff_gradient = _readonly_vec3(
        section10_spatial_cutoff_gradient(x, y, z, t),
        "Section 10 cutoff gradient",
    )
    localized_velocity = np.array(
        cutoff_value * curl_A
        + np.cross(cutoff_gradient, A_value)
        + cutoff_value * swirl,
        dtype=float,
    )
    if not np.all(np.isfinite(localized_velocity)):
        raise ArithmeticError("localized velocity is not finite")
    localized_velocity.setflags(write=False)

    return Section9Section10LocalizedVelocityJetCertificate(
        point=(x, y, z, t),
        q=prefix.q,
        prefix_order=prefix.prefix_order,
        derivative_order=prefix.derivative_order,
        provider_id=prefix.provider_id,
        provider_revision=prefix.provider_revision,
        provider_provenance=prefix.provider_provenance,
        A_value=A_value,
        B_value=B_value,
        curl_A=curl_A,
        local_velocity=local_velocity,
        cutoff_value=cutoff_value,
        cutoff_gradient=cutoff_gradient,
        localized_velocity=localized_velocity,
    )
