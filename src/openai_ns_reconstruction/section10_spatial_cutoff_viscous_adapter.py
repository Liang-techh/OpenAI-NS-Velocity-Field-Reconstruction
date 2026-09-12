"""Analytic viscous-operator derivatives for the fixed Section 10 cutoff.

This module stays on the already-admitted Section 10 spatial-localization path.
The geometry is the official ``chi(16 r^2) chi(4 z)`` support/plateau geometry,
while transition-collar values remain those of the repository's explicit
C-infinity representative rather than a pointwise claim about Mathlib's
noncomputable ``ContDiffBump``.

The main purpose is to expose ``grad(Delta c)`` analytically.  A future
Laplacian of

    u_cut = c curl(A) + grad(c) x A + c B e_theta

needs this third-order cutoff derivative through ``Delta(grad(c) x A)``.
Nothing here constructs the missing Section 7/8 correction field, the infinite
Eq. (9.21) sum, a t=1 extension, or a Navier--Stokes residual/forcing.
"""
from __future__ import annotations

import math
import numpy as np

from .cutoffs import standard_cutoff_third_derivative
from .spatial_localization import (
    radial_square,
    section10_spatial_cutoff_hessian,
    symmetric_smooth_bump,
    symmetric_smooth_bump_derivative,
    symmetric_smooth_bump_second_derivative,
)


def _finite_coordinates(
    x: float, y: float, z: float, t: float
) -> tuple[float, float, float, float]:
    values = tuple(float(v) for v in (x, y, z, t))
    if not all(math.isfinite(v) for v in values):
        raise ValueError("spatial-localization coordinates must be finite")
    return values


def symmetric_smooth_bump_third_derivative(s: float) -> float:
    """Third derivative of the executable even bump representative."""
    s = float(s)
    if not math.isfinite(s):
        raise ValueError("spatial-localization coordinates must be finite")
    if s == 0.0:
        return 0.0
    sign = 1.0 if s > 0.0 else -1.0
    return sign * standard_cutoff_third_derivative(abs(s))


def section10_spatial_cutoff_laplacian(
    x: float, y: float, z: float, t: float = 0.0
) -> float:
    """Analytic ``Delta c`` from the already-landed cutoff Hessian."""
    x, y, z, t = _finite_coordinates(x, y, z, t)
    value = float(np.trace(section10_spatial_cutoff_hessian(x, y, z, t)))
    if not math.isfinite(value):
        raise ArithmeticError("spatial cutoff Laplacian is not finite")
    return value


def section10_spatial_cutoff_laplacian_gradient(
    x: float, y: float, z: float, t: float = 0.0
) -> np.ndarray:
    """Exact Cartesian ``grad(Delta c)`` for the executable fixed cutoff.

    For ``c=b(16 r^2)b(4z)`` and ``r^2=x^2+y^2``, differentiating the trace of
    the existing analytic Hessian gives

    ``d_x Delta c = x[(4096 b_r'' + 32768 r^2 b_r''')b_z
                       + 512 b_r' b_z'']``

    with the analogous ``y`` expression, and

    ``d_z Delta c = (256 b_r' + 4096 r^2 b_r'')b_z'
                    + 64 b_r b_z'''``.

    The formula uses no finite differences in production.  It is exact only
    for the repository's geometry-faithful C-infinity representative.
    """
    x, y, z, _t = _finite_coordinates(x, y, z, t)
    r2 = radial_square(x, y)
    r_arg = 16.0 * r2
    z_arg = 4.0 * z

    radial = symmetric_smooth_bump(r_arg)
    axial = symmetric_smooth_bump(z_arg)
    radial_d1 = symmetric_smooth_bump_derivative(r_arg)
    axial_d1 = symmetric_smooth_bump_derivative(z_arg)
    radial_d2 = symmetric_smooth_bump_second_derivative(r_arg)
    axial_d2 = symmetric_smooth_bump_second_derivative(z_arg)
    radial_d3 = symmetric_smooth_bump_third_derivative(r_arg)
    axial_d3 = symmetric_smooth_bump_third_derivative(z_arg)

    radial_common = (
        (4096.0 * radial_d2 + 32768.0 * r2 * radial_d3) * axial
        + 512.0 * radial_d1 * axial_d2
    )
    gradient = np.array(
        [
            x * radial_common,
            y * radial_common,
            (256.0 * radial_d1 + 4096.0 * r2 * radial_d2) * axial_d1
            + 64.0 * radial * axial_d3,
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(gradient)):
        raise ArithmeticError("spatial cutoff Laplacian gradient is not finite")
    return gradient
