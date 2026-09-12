"""Section 10 spatial-localization geometry, separated from paper-exact values.

Provenance
----------
OpenAI ``NavierStokes/SpatialLocalization.lean`` at the repository-pinned
formalization commit defines

    radialSquare(x) = x_0^2 + x_1^2
    cutoffProfile(r2, z) = chi(16 r2) chi(4 z)

where ``chi`` is Mathlib's C-infinity ``ContDiffBump`` with inner radius 1/2
and outer radius 1. Consequently the cutoff is one on
``r^2 < 1/32, |z| < 1/8`` and supported in
``r^2 <= 1/16, |z| <= 1/4``.

The exact numerical transition values of Mathlib's ``ContDiffBump`` are not
reimplemented here. ``section10_spatial_cutoff`` uses this repository's
explicit C-infinity representative with the *same plateau and support
geometry*. ``section10_spatial_cutoff_gradient`` and
``section10_spatial_cutoff_hessian`` are exact analytic derivatives of that
representative, so product-rule paths do not need to finite-difference the
cutoff. None of these functions is promoted to a paper-exact scalar
realization in the transition collar.
"""
from __future__ import annotations

import math
import numpy as np

from .cutoffs import (
    smooth_cutoff,
    standard_cutoff_derivative,
    standard_cutoff_second_derivative,
)

SUPPORT_RADIUS_SQUARED = 1.0 / 16.0
SUPPORT_HALF_HEIGHT = 1.0 / 4.0
PLATEAU_RADIUS_SQUARED = 1.0 / 32.0
PLATEAU_HALF_HEIGHT = 1.0 / 8.0
PERIOD_HALF_WIDTH = 1.0 / 2.0


def _finite(*values: float) -> tuple[float, ...]:
    out = tuple(float(v) for v in values)
    if not all(math.isfinite(v) for v in out):
        raise ValueError("spatial-localization coordinates must be finite")
    return out


def radial_square(x: float, y: float) -> float:
    """Squared distance to the symmetry axis; mirrors the Lean definition."""
    x, y = _finite(x, y)
    return x * x + y * y


def symmetric_smooth_bump(s: float) -> float:
    """C-infinity even bump with plateau |s|<=1/2 and support |s|<1.

    This has the same support/plateau specification as the Mathlib bump used by
    the official formalization, but not necessarily the same collar values.
    """
    (s,) = _finite(s)
    return smooth_cutoff(abs(s))


def symmetric_smooth_bump_derivative(s: float) -> float:
    """Exact derivative of this repository's even representative.

    ``smooth_cutoff`` is flat on a neighborhood of zero, so the derivative at
    the absolute-value cusp is exactly zero. In the transition collar the sign
    factor implements d/ds smooth_cutoff(|s|).
    """
    (s,) = _finite(s)
    if s == 0.0:
        return 0.0
    sign = 1.0 if s > 0.0 else -1.0
    return sign * standard_cutoff_derivative(abs(s))


def symmetric_smooth_bump_second_derivative(s: float) -> float:
    """Exact second derivative of this repository's even representative."""
    (s,) = _finite(s)
    return standard_cutoff_second_derivative(abs(s))


def support_cylinder_contains(x: float, y: float, z: float) -> bool:
    """Closed support cylinder from ``SpatialLocalization.supportCylinder``."""
    x, y, z = _finite(x, y, z)
    return radial_square(x, y) <= SUPPORT_RADIUS_SQUARED and abs(z) <= SUPPORT_HALF_HEIGHT


def plateau_contains(x: float, y: float, z: float) -> bool:
    """Open unit-cutoff cylinder from ``SpatialLocalization.plateau``."""
    x, y, z = _finite(x, y, z)
    return radial_square(x, y) < PLATEAU_RADIUS_SQUARED and abs(z) < PLATEAU_HALF_HEIGHT


def section10_spatial_cutoff(x: float, y: float, z: float, t: float = 0.0) -> float:
    """Geometry-faithful representative of the official fixed spatial cutoff.

    ``t`` is accepted so the function can be passed directly as a repository
    ``ScalarField``. The Section 10 spatial cutoff itself is time-independent.
    """
    x, y, z, _t = _finite(x, y, z, t)
    r2 = radial_square(x, y)
    return symmetric_smooth_bump(16.0 * r2) * symmetric_smooth_bump(4.0 * z)


def section10_spatial_cutoff_gradient(
    x: float, y: float, z: float, t: float = 0.0
) -> np.ndarray:
    """Analytic Cartesian gradient of ``section10_spatial_cutoff``.

    For ``c(x,y,z)=b(16(x^2+y^2)) b(4z)`` with even ``b``, the chain rule gives
    ``(32 x b'_r b_z, 32 y b'_r b_z, 4 b_r b'_z)``. This is the derivative
    needed by the Section 10 product rule. It is exact for the executable bump
    representative, not a claim that its transition values equal Mathlib's
    noncomputable ``ContDiffBump`` pointwise.
    """
    x, y, z, _t = _finite(x, y, z, t)
    r_arg = 16.0 * radial_square(x, y)
    z_arg = 4.0 * z
    radial = symmetric_smooth_bump(r_arg)
    axial = symmetric_smooth_bump(z_arg)
    radial_derivative = symmetric_smooth_bump_derivative(r_arg)
    axial_derivative = symmetric_smooth_bump_derivative(z_arg)
    gradient = np.array(
        [
            32.0 * x * radial_derivative * axial,
            32.0 * y * radial_derivative * axial,
            4.0 * radial * axial_derivative,
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(gradient)):
        raise ArithmeticError("spatial cutoff gradient is not finite")
    return gradient


def section10_spatial_cutoff_hessian(
    x: float, y: float, z: float, t: float = 0.0
) -> np.ndarray:
    """Analytic Cartesian Hessian of the fixed Section 10 cutoff representative.

    This is a derivative adapter for later analytic ``grad(u)``/``Delta u``
    assembly. It preserves the already-fixed support/plateau geometry and does
    not promote the executable transition collar to the paper's noncomputable
    Mathlib ``ContDiffBump`` values.
    """
    x, y, z, _t = _finite(x, y, z, t)
    r_arg = 16.0 * radial_square(x, y)
    z_arg = 4.0 * z
    radial = symmetric_smooth_bump(r_arg)
    axial = symmetric_smooth_bump(z_arg)
    radial_d1 = symmetric_smooth_bump_derivative(r_arg)
    axial_d1 = symmetric_smooth_bump_derivative(z_arg)
    radial_d2 = symmetric_smooth_bump_second_derivative(r_arg)
    axial_d2 = symmetric_smooth_bump_second_derivative(z_arg)

    h_xx = (32.0 * radial_d1 + 1024.0 * x * x * radial_d2) * axial
    h_yy = (32.0 * radial_d1 + 1024.0 * y * y * radial_d2) * axial
    h_xy = 1024.0 * x * y * radial_d2 * axial
    h_xz = 128.0 * x * radial_d1 * axial_d1
    h_yz = 128.0 * y * radial_d1 * axial_d1
    h_zz = 16.0 * radial * axial_d2
    hessian = np.array(
        [
            [h_xx, h_xy, h_xz],
            [h_xy, h_yy, h_yz],
            [h_xz, h_yz, h_zz],
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(hessian)):
        raise ArithmeticError("spatial cutoff Hessian is not finite")
    return hessian


def section10_localized_field(local: "LocalField") -> "LocalizedField":
    """Bind a local field to the fixed Section 10 cutoff and exact gradient.

    This adapter intentionally makes the fixed support/plateau representative
    and its matching analytic derivative an atomic choice.  In particular it
    prevents a caller from accidentally pairing the Section 10 cutoff with a
    generic or mismatched gradient while using ``LocalizedField``'s analytic
    product-rule path.

    The adapter does *not* certify ``local`` as the paper's completed profile,
    and the executable transition collar remains only geometry-faithful to the
    formal ``ContDiffBump``.  It therefore carries no paper-exact promotion.
    """
    # Import lazily to keep spatial geometry independent of the local-field
    # module at import time while still enforcing the expected constructor.
    from .local_field import LocalField, LocalizedField

    if not isinstance(local, LocalField):
        raise TypeError("local must be a LocalField")
    return LocalizedField(
        local=local,
        cutoff=section10_spatial_cutoff,
        cutoff_gradient=section10_spatial_cutoff_gradient,
    )


def section10_axisymmetric_cutoff(r: float, z: float, t: float = 0.0) -> float:
    """Axisymmetric adapter for ``LocalizedField.from_axisymmetric``."""
    r, z, t = _finite(r, z, t)
    if r < 0:
        raise ValueError("cylindrical radius must be nonnegative")
    return section10_spatial_cutoff(r, 0.0, z, t)


def support_is_strictly_inside_period_cube(x: float, y: float, z: float) -> bool:
    """Executable counterpart of the formalization's strict cube separation.

    Returns ``False`` outside the support cylinder. Every supported point has
    |x_i| <= 1/4 < 1/2, so translated unit-period copies cannot overlap in the
    central inner cube used by the formalization.
    """
    x, y, z = _finite(x, y, z)
    if not support_cylinder_contains(x, y, z):
        return False
    return max(abs(x), abs(y), abs(z)) < PERIOD_HALF_WIDTH
