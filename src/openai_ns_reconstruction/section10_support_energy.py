"""Section 10 compact-support and fixed-time kinetic-energy certificates.

This module records consequences of the *fixed* Section 10 support cylinder.
It deliberately does not certify the unresolved local velocity field, a global
in-time L-infinity bound, or bounded energy along the blow-up limit.

For the official support geometry

    r^2 <= 1/16,  |z| <= 1/4,

the enclosing cylinder has radius 1/4, height 1/2, and volume pi/32.  Hence,
at any fixed time at which an independently certified velocity bound
``|u| <= M`` holds on the localized field,

    (1/2) * integral |u|^2 dx <= pi * M^2 / 64.

The implication is exact at the geometry/algebra level.  Supplying ``M`` is an
upstream proof obligation; sampled maxima are not promoted to certificates.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math

from .coordinates import _finite
from .spatial_localization import (
    SUPPORT_HALF_HEIGHT,
    SUPPORT_RADIUS_SQUARED,
    support_cylinder_contains,
)

# Exact rational geometry pinned by SpatialLocalization.supportCylinder.
SUPPORT_RADIUS_SQUARED_EXACT = Fraction(1, 16)
SUPPORT_HALF_HEIGHT_EXACT = Fraction(1, 4)
SUPPORT_RADIUS_EXACT = Fraction(1, 4)
SUPPORT_HEIGHT_EXACT = Fraction(1, 2)
SUPPORT_VOLUME_OVER_PI_EXACT = Fraction(1, 32)
FIXED_TIME_ENERGY_COEFFICIENT_OVER_PI_EXACT = Fraction(1, 64)


def section10_support_volume() -> float:
    """Volume of the enclosing support cylinder, exactly ``pi/32``."""
    return math.pi * float(SUPPORT_VOLUME_OVER_PI_EXACT)


def fixed_time_energy_bound_from_speed_sup(speed_sup_bound: float) -> float:
    """Return the support-volume bound ``E(t) <= pi*M^2/64``.

    ``speed_sup_bound`` must be an independently certified upper bound for
    ``|u(t,x)|`` on the actual localized field at the time under study.  This
    function only propagates that bound through the fixed support geometry; it
    does not infer a supremum from samples.

    A finite bound here proves finite kinetic energy at that fixed time.  It
    does *not* prove a uniform energy bound as ``t -> 1-`` unless the supplied
    upstream bound is itself strong enough to do so.
    """
    speed_sup_bound = _finite(speed_sup_bound, "speed_sup_bound")
    if speed_sup_bound < 0:
        raise ValueError("speed_sup_bound must be nonnegative")
    return (
        math.pi
        * float(FIXED_TIME_ENERGY_COEFFICIENT_OVER_PI_EXACT)
        * speed_sup_bound
        * speed_sup_bound
    )


@dataclass(frozen=True)
class Section10SupportEnergyCertificate:
    """Geometry-level certificate parameterized by a certified speed bound.

    The object is intentionally narrow: it certifies the exact enclosing
    support geometry and the resulting fixed-time energy implication.  It does
    not certify smoothness, the value of the speed bound, the blow-up path, or
    a uniform-in-time energy estimate.
    """

    speed_sup_bound: float

    def __post_init__(self) -> None:
        value = _finite(self.speed_sup_bound, "speed_sup_bound")
        if value < 0:
            raise ValueError("speed_sup_bound must be nonnegative")
        object.__setattr__(self, "speed_sup_bound", value)

    @property
    def support_radius(self) -> float:
        return float(SUPPORT_RADIUS_EXACT)

    @property
    def support_half_height(self) -> float:
        return float(SUPPORT_HALF_HEIGHT_EXACT)

    @property
    def support_volume(self) -> float:
        return section10_support_volume()

    @property
    def fixed_time_energy_upper_bound(self) -> float:
        return fixed_time_energy_bound_from_speed_sup(self.speed_sup_bound)

    def enclosing_support_contains(self, x: float, y: float, z: float) -> bool:
        """Membership in the official closed support cylinder."""
        return support_cylinder_contains(x, y, z)


def geometry_matches_spatial_localization_constants() -> bool:
    """Guard against silent drift between the certificate and cutoff module."""
    return (
        SUPPORT_RADIUS_SQUARED == float(SUPPORT_RADIUS_SQUARED_EXACT)
        and SUPPORT_HALF_HEIGHT == float(SUPPORT_HALF_HEIGHT_EXACT)
    )
