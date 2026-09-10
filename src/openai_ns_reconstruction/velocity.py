"""Leading kinematics, Eqs. (4.3)-(4.7); profile validity is a separate task."""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np
from .coordinates import similarity_coordinates, similarity_coordinates_from_tau, SimilarityPoint
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
    E, U = float(profile.E(s.X, s.eta)), float(profile.U(s.X, s.eta))
    if not math.isfinite(E) or not math.isfinite(U):
        raise ArithmeticError("profile returned nonfinite velocity")
    if r == 0 and E != 0:
        raise ValueError("a regular axisymmetric profile must satisfy E(0,eta)=0")
    flux = profile.V0(s.X, s.eta, s.h, n=quadrature_points, d=s.d)
    velocity = CylindricalVelocity(0.0 if r == 0 else flux / r,
                                  s.q**(-s.A) * E, s.q**(-s.A) * U,
                                  s.q, s.eta, s.X)
    if not np.all(np.isfinite(velocity.vector)):
        raise ArithmeticError("velocity lies outside binary64 range")
    return velocity


def leading_velocity_cylindrical(r: float, z: float, t: float, profile: LeadingProfile,
                                 *, h: float = 0.005,
                                 quadrature_points: int = 32) -> CylindricalVelocity:
    return _evaluate(r, similarity_coordinates(r, z, t, h), profile, quadrature_points)


def leading_velocity_from_tau(r: float, z: float, tau: float, profile: LeadingProfile,
                               *, h: float = 0.005,
                               quadrature_points: int = 32) -> CylindricalVelocity:
    """Leading cylindrical field without forming the lossy subtraction 1-tau."""
    return _evaluate(r, similarity_coordinates_from_tau(r, z, tau, h),
                     profile, quadrature_points)


def leading_velocity_cartesian(x: float, y: float, z: float, t: float,
                               profile: LeadingProfile, *, h: float = 0.005,
                               quadrature_points: int = 32) -> np.ndarray:
    r = math.hypot(x, y)
    v = leading_velocity_cylindrical(r, z, t, profile, h=h,
                                     quadrature_points=quadrature_points)
    if r == 0:
        return np.array([0.0, 0.0, v.u_z])
    return np.array([(v.u_r * x - v.u_theta * y) / r,
                     (v.u_r * y + v.u_theta * x) / r, v.u_z])


def leading_pressure_cartesian(x: float, y: float, z: float, t: float,
                               profile: LeadingProfile, *, h: float = 0.005) -> float:
    if profile.Pi is None:
        raise ValueError("pressure profile Pi is missing; zero is not an implicit default")
    s = similarity_coordinates(math.hypot(x, y), z, t, h)
    value = s.q**(-2 * s.A) * float(profile.Pi(s.X, s.eta))
    if not math.isfinite(value):
        raise ArithmeticError("pressure is not finite")
    return value


def blowup_probe(tau: float, X_in: float, profile: LeadingProfile, *,
                  h: float = 0.005) -> CylindricalVelocity:
    if not math.isfinite(tau) or tau <= 0 or not math.isfinite(X_in) or X_in < 0:
        raise ValueError("tau must be positive and X_in nonnegative, both finite")
    r = math.sqrt(tau) * math.sqrt(2 * X_in)
    return leading_velocity_from_tau(r, 0.0, tau, profile, h=h)
