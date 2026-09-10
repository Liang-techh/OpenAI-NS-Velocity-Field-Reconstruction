"""All-order axisymmetric background architecture from Section 5.

The paper writes the nth coefficient at lambda_n = 2 n h as

    u_theta,n = q^(-A+lambda_n) E_n,
    u_z,n     = q^(-A+lambda_n) U_n,
    r u_r,n   = q^(lambda_n) V_n,
    p_n       = q^(-2A+lambda_n) Pi_n.                 (5.1)

The actual smooth background is assembled with q-dependent cutoffs applied to vector
potentials before curl.  This module encodes that representation without inventing the
paper's recursively constructed coefficient profiles.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence
import math
import numpy as np

from .coordinates import similarity_coordinates
from .profiles import LeadingProfile

ScalarCutoff = Callable[[float], float]


def standard_cutoff(s: float) -> float:
    """Simple C1 compact cutoff for experiments, not the paper's fixed C-infinity cutoff."""
    s = float(s)
    if s <= 0.5:
        return 1.0
    if s >= 1.0:
        return 0.0
    # smoothstep on [1/2,1], enough for numerical experiments only
    u = 2.0 * (s - 0.5)
    return 1.0 - (3.0 * u * u - 2.0 * u * u * u)


@dataclass(frozen=True)
class BackgroundCoefficient:
    n: int
    profile: LeadingProfile
    cutoff_scale: float = 1.0

    def lambda_n(self, h: float) -> float:
        return 2.0 * self.n * h


def coefficient_radial_flux(
    X: float,
    eta: float,
    coefficient: BackgroundCoefficient,
    *,
    h: float,
    quadrature_points: int = 801,
) -> float:
    """Return the Section 5 radial-flux coefficient ``V_n(X, eta)``.

    This is the second identity in Eq. (5.2):

        V_n / X = [2 eta U_n
                   - 2 eta (D + lambda_n) A_X(U_n)
                   - d partial_eta A_X(U_n)] / L.

    It is the exact kinematic consequence of incompressibility for the nth formal
    coefficient.  It does *not* construct the recursive positive-order profiles
    ``U_n`` themselves; callers must provide those through ``coefficient.profile``.
    For n=0 the formula reduces to the leading-flow identity (4.7).
    """
    X = float(X)
    eta = float(eta)
    if X < 0.0:
        raise ValueError("X must be nonnegative")
    if X == 0.0:
        return 0.0
    if not (0.0 < h < 0.5):
        raise ValueError("h must satisfy 0 < h < 1/2")

    D = 0.5 - h
    d = 1.0 - eta * eta
    L = 1.0 - 2.0 * h * eta * eta
    lam = coefficient.lambda_n(h)
    U = float(coefficient.profile.U(X, eta))
    AU = coefficient.profile.radial_average_U(X, eta, n=quadrature_points)
    dAU = coefficient.profile.radial_average_dU_deta(X, eta, n=quadrature_points)
    return (X / L) * (2.0 * eta * U - 2.0 * eta * (D + lam) * AU - d * dAU)


def coefficient_streamfunction(
    r: float,
    z: float,
    t: float,
    coefficient: BackgroundCoefficient,
    *,
    h: float,
    quadrature_points: int = 801,
) -> float:
    """Physical Stokes streamfunction S_n in Eq. (5.27).

    F_n = X A_X(U_n),
    S_n = q^(1-A+lambda_n) F_n.
    """
    s = similarity_coordinates(r, z, t, h)
    AU = coefficient.profile.radial_average_U(s.X, s.eta, n=quadrature_points)
    F_n = s.X * AU
    return s.q ** (1.0 - s.A + coefficient.lambda_n(h)) * F_n


def coefficient_vector_potential_cartesian(
    x: float,
    y: float,
    z: float,
    t: float,
    coefficient: BackgroundCoefficient,
    *,
    h: float,
    quadrature_points: int = 801,
) -> np.ndarray:
    """A_n=(S_n/r)e_theta, written smoothly away from the axis, Eq. (5.27)."""
    r2 = float(x) ** 2 + float(y) ** 2
    if r2 == 0.0:
        return np.zeros(3)
    r = math.sqrt(r2)
    S_n = coefficient_streamfunction(
        r, z, t, coefficient, h=h, quadrature_points=quadrature_points
    )
    factor = S_n / r2
    return np.array([-factor * y, factor * x, 0.0], dtype=float)


def background_potential_cartesian(
    x: float,
    y: float,
    z: float,
    t: float,
    coefficients: Sequence[BackgroundCoefficient],
    *,
    h: float = 0.005,
    cutoff: ScalarCutoff = standard_cutoff,
) -> np.ndarray:
    """Cutoff-summed axisymmetric vector potential used before taking curl.

    This mirrors the structural requirement in Sections 5 and 10: cut off the potential,
    not the poloidal velocity, so divergence-freeness survives after taking curl.
    """
    r = math.hypot(x, y)
    s = similarity_coordinates(r, z, t, h)
    total = np.zeros(3)
    for c in coefficients:
        weight = 1.0 if c.n == 0 else cutoff(c.cutoff_scale * s.q)
        total += weight * coefficient_vector_potential_cartesian(x, y, z, t, c, h=h)
    return total


def background_swirl_cartesian(
    x: float,
    y: float,
    z: float,
    t: float,
    coefficients: Sequence[BackgroundCoefficient],
    *,
    h: float = 0.005,
    cutoff: ScalarCutoff = standard_cutoff,
) -> np.ndarray:
    """Cutoff-summed direct azimuthal part B e_theta of the background."""
    r = math.hypot(x, y)
    if r == 0.0:
        return np.zeros(3)
    s = similarity_coordinates(r, z, t, h)
    B = 0.0
    for c in coefficients:
        weight = 1.0 if c.n == 0 else cutoff(c.cutoff_scale * s.q)
        lam = c.lambda_n(h)
        B += weight * s.q ** (-s.A + lam) * float(c.profile.E(s.X, s.eta))
    e_theta = np.array([-y / r, x / r, 0.0])
    return B * e_theta
