"""Section 10 spatial-localization geometry, separated from paper-exact values.

Provenance
----------
OpenAI ``NavierStokes/SpatialLocalization.lean`` at the repository-pinned
formalization commit defines

    radialSquare(x) = x_0^2 + x_1^2
    cutoffProfile(r2, z) = chi(16 r2) chi(4 z)

where ``chi`` is Mathlib's C-infinity ``ContDiffBump`` with inner radius 1/2
and outer radius 1.  Consequently the cutoff is one on
``r^2 < 1/32, |z| < 1/8`` and supported in
``r^2 <= 1/16, |z| <= 1/4``.

The exact numerical transition values of Mathlib's ``ContDiffBump`` are not
reimplemented here.  ``section10_spatial_cutoff`` uses this repository's
explicit C-infinity representative with the *same plateau and support
geometry*.  It must therefore be treated as ``formal-geometry`` / diagnostic
in the transition collar, not as a paper-exact scalar cutoff.
"""
from __future__ import annotations

import math

from .cutoffs import smooth_cutoff

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
    ``ScalarField``.  The Section 10 spatial cutoff itself is time-independent.
    """
    x, y, z, _t = _finite(x, y, z, t)
    r2 = radial_square(x, y)
    return symmetric_smooth_bump(16.0 * r2) * symmetric_smooth_bump(4.0 * z)


def section10_axisymmetric_cutoff(r: float, z: float, t: float = 0.0) -> float:
    """Axisymmetric adapter for ``LocalizedField.from_axisymmetric``."""
    r, z, t = _finite(r, z, t)
    if r < 0:
        raise ValueError("cylindrical radius must be nonnegative")
    return section10_spatial_cutoff(r, 0.0, z, t)


def support_is_strictly_inside_period_cube(x: float, y: float, z: float) -> bool:
    """Executable counterpart of the formalization's strict cube separation.

    Returns ``False`` outside the support cylinder.  Every supported point has
    |x_i| <= 1/4 < 1/2, so translated unit-period copies cannot overlap in the
    central inner cube used by the formalization.
    """
    x, y, z = _finite(x, y, z)
    if not support_cylinder_contains(x, y, z):
        return False
    return max(abs(x), abs(y), abs(z)) < PERIOD_HALF_WIDTH
