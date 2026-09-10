"""Leading velocity/pressure, Eqs. (4.3)-(4.7), with explicit tau interfaces."""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from .coordinates import (SimilarityPoint, _finite, similarity_coordinates,
                          similarity_coordinates_from_tau)
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


def _evaluate(r: float, s: SimilarityPoint, profile: LeadingProfile,
              quadrature_points: int) -> CylindricalVelocity:
    U = _finite(profile.U(s.X, s.eta), "U")
    v0 = profile.radial_flux_factor(s.X, s.eta, s.h, n=quadrature_points, d=s.d, L=s.L)
    E = (math.sqrt(2*s.X)*profile.smooth_swirl_factor(s.X, s.eta)
         if profile.F is not None else _finite(profile.E(s.X, s.eta), "E"))
    if r == 0 and E != 0:
        raise ValueError("a smooth axisymmetric profile must have E(0,eta)=0")
    # No arbitrary 1e-15 radius cutoff: nonzero radii retain their velocity.
    ur = (r / math.sqrt(s.q)) * v0 / (2 * math.sqrt(s.q))
    ut, uz = s.q**(-s.A)*E, s.q**(-s.A)*U
    if not all(math.isfinite(v) for v in (ur, ut, uz)):
        raise OverflowError("velocity is outside floating-point range")
    return CylindricalVelocity(ur, ut, uz, s.q, s.eta, s.X)


def leading_velocity_cylindrical(r: float, z: float, t: float, profile: LeadingProfile, *,
                                 h: float = 0.005, quadrature_points: int = 32) -> CylindricalVelocity:
    return _evaluate(float(r), similarity_coordinates(r,z,t,h), profile, quadrature_points)


def leading_velocity_cylindrical_from_tau(r: float, z: float, tau: float,
                                          profile: LeadingProfile, *, h: float = 0.005,
                                          quadrature_points: int = 32) -> CylindricalVelocity:
    return _evaluate(float(r), similarity_coordinates_from_tau(r,z,tau,h), profile, quadrature_points)


def _cartesian(x: float, y: float, v: CylindricalVelocity) -> np.ndarray:
    r = math.hypot(x, y)
    if r == 0:
        return np.array([0., 0., v.u_z])
    return np.array([v.u_r*(x/r)-v.u_theta*(y/r),
                     v.u_r*(y/r)+v.u_theta*(x/r), v.u_z])


def leading_velocity_cartesian(x: float, y: float, z: float, t: float,
                               profile: LeadingProfile, *, h: float = 0.005,
                               quadrature_points: int = 32) -> np.ndarray:
    """Cartesian form equivalent to Eq. (4.5), with the regular axis value."""
    x, y = _finite(x,"x"), _finite(y,"y")
    return _cartesian(x,y,leading_velocity_cylindrical(math.hypot(x,y),z,t,profile,
                       h=h,quadrature_points=quadrature_points))


def leading_velocity_cartesian_from_tau(x: float, y: float, z: float, tau: float,
                                        profile: LeadingProfile, *, h: float = 0.005,
                                        quadrature_points: int = 32) -> np.ndarray:
    x, y = _finite(x,"x"), _finite(y,"y")
    return _cartesian(x,y,leading_velocity_cylindrical_from_tau(math.hypot(x,y),z,tau,
                        profile,h=h,quadrature_points=quadrature_points))


def leading_pressure(x: float, y: float, z: float, t: float,
                     profile: LeadingProfile, *, h: float = 0.005) -> float:
    if profile.Pi is None:
        raise ValueError("pressure profile Pi is required; do not silently use zero pressure")
    s = similarity_coordinates(math.hypot(x,y),z,t,h)
    return _finite(s.q**(-2*s.A)*profile.Pi(s.X,s.eta), "pressure")


def blowup_probe(tau: float, X_in: float, profile: LeadingProfile, *,
                 h: float = 0.005) -> CylindricalVelocity:
    tau, X_in = _finite(tau,"tau"), _finite(X_in,"X_in")
    if tau <= 0 or X_in < 0:
        raise ValueError("tau must be positive and X_in nonnegative")
    r = math.sqrt(2*X_in) * math.sqrt(tau)
    return leading_velocity_cylindrical_from_tau(r,0.,tau,profile,h=h)
