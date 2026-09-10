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
