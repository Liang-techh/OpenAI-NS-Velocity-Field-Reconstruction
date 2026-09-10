"""Velocity evaluators for the leading OpenAI blow-up field.

Implements Eqs. (4.3), (4.5), and (4.7) of the paper.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from .coordinates import similarity_coordinates
from .profiles import LeadingProfile


@dataclass(frozen=True)
class CylindricalVelocity:
    u_r: float
    u_theta: float
    u_z: float
    q: float
    eta: float
    X: float

    @property
    def vector(self) -> np.ndarray:
        return np.array([self.u_r, self.u_theta, self.u_z], dtype=float)


def leading_velocity_cylindrical(
    r: float,
    z: float,
    t: float,
    profile: LeadingProfile,
    *,
    h: float = 0.005,
    quadrature_points: int = 801,
) -> CylindricalVelocity:
    """Evaluate the leading axisymmetric field u^(0).

    Paper equations:
        u_theta^(0) = q^(-A) E(X,eta)
        u_z^(0)     = q^(-A) U(X,eta)
        r u_r^(0)   = V0(X,eta)
    where V0 is fixed by incompressibility, Eq. (4.7).
    """
    s = similarity_coordinates(r, z, t, h)
    E = float(profile.E(s.X, s.eta))
    U = float(profile.U(s.X, s.eta))
    V0 = profile.V0(s.X, s.eta, h, n=quadrature_points)

    u_theta = s.q ** (-s.A) * E
    u_z = s.q ** (-s.A) * U
    if abs(r) < 1e-15:
        # Smooth axisymmetric Cartesian field has u_r=u_theta=0 on the axis.
        u_r = 0.0
        u_theta = 0.0
    else:
        u_r = V0 / float(r)
    return CylindricalVelocity(u_r, u_theta, u_z, s.q, s.eta, s.X)


def leading_velocity_cartesian(
    x: float,
    y: float,
    z: float,
    t: float,
    profile: LeadingProfile,
    *,
    h: float = 0.005,
    quadrature_points: int = 801,
) -> np.ndarray:
    """Evaluate the leading velocity in Cartesian coordinates.

    This is the cylindrical-to-Cartesian form of Eq. (4.5).  It is equivalent to

        u1 = (v0/(2q)) x - q^(-A-1/2) F y
        u2 = (v0/(2q)) y + q^(-A-1/2) F x
        u3 = q^(-A) U,

    with E=sqrt(2X)F and V0=X v0.
    """
    x = float(x)
    y = float(y)
    r = math.hypot(x, y)
    cyl = leading_velocity_cylindrical(
        r, z, t, profile, h=h, quadrature_points=quadrature_points
    )
    if r == 0.0:
        return np.array([0.0, 0.0, cyl.u_z], dtype=float)
    c = x / r
    s = y / r
    u_x = cyl.u_r * c - cyl.u_theta * s
    u_y = cyl.u_r * s + cyl.u_theta * c
    return np.array([u_x, u_y, cyl.u_z], dtype=float)


def blowup_probe(
    tau: float,
    X_in: float,
    profile: LeadingProfile,
    *,
    h: float = 0.005,
) -> CylindricalVelocity:
    """Evaluate the leading field on the paper's singular path z=0, r=sqrt(2 X_in tau)."""
    if tau <= 0.0:
        raise ValueError("tau must be positive")
    r = math.sqrt(2.0 * X_in * tau)
    return leading_velocity_cylindrical(r, 0.0, 1.0 - tau, profile, h=h)
