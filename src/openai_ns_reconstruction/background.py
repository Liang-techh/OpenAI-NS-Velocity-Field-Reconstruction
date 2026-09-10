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

from .coordinates import similarity_coordinates, _finite, _validate_h
from .cutoffs import standard_cutoff, standard_cutoff_derivative
from .profiles import LeadingProfile

ScalarCutoff = Callable[[float], float]


@dataclass(frozen=True)
class BackgroundCoefficient:
    n: int
    profile: LeadingProfile
    cutoff_scale: float = 1.0

    def __post_init__(self) -> None:
        if isinstance(self.n, bool) or not isinstance(self.n, int) or self.n < 0:
            raise ValueError("coefficient order n must be a nonnegative integer")
        if _finite(self.cutoff_scale, "cutoff_scale") <= 0:
            raise ValueError("cutoff_scale must be positive")

    def lambda_n(self, h: float) -> float:
        return 2.0 * self.n * _validate_h(h)


def coefficient_streamfunction(
    r: float,
    z: float,
    t: float,
    coefficient: BackgroundCoefficient,
    *,
    h: float,
    quadrature_points: int = 32,
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
    quadrature_points: int = 32,
) -> np.ndarray:
    """A_n=(S_n/r)e_theta, in axis-regular Cartesian form, Eq. (5.27)."""
    x, y = _finite(x,"x"), _finite(y,"y")
    s = similarity_coordinates(math.hypot(x,y), z, t, h)
    AU = coefficient.profile.radial_average_U(s.X,s.eta,n=quadrature_points)
    # S/r^2 = q^(-A+lambda_n) A_X(U_n)/2, regular even on the axis.
    factor = 0.5 * s.q**(-s.A+coefficient.lambda_n(h)) * AU
    return factor * np.array([-y,x,0.0])


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
        weight = 1.0 if c.n == 0 else _finite(cutoff(c.cutoff_scale * s.q),"cutoff")
        if weight == 0:
            continue
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
    s = similarity_coordinates(r, z, t, h)
    B = 0.0
    for c in coefficients:
        weight = 1.0 if c.n == 0 else _finite(cutoff(c.cutoff_scale * s.q),"cutoff")
        if weight == 0:
            continue
        E = _finite(c.profile.E(s.X,s.eta),"E")
        if r == 0 and E != 0:
            raise ValueError("E must vanish on the axis")
        lam = c.lambda_n(h)
        B += weight * s.q ** (-s.A + lam) * E
    B = _finite(B,"background swirl")
    if r == 0:
        return np.zeros(3)
    e_theta = np.array([-y / r, x / r, 0.0])
    return B * e_theta


def background_velocity_cartesian(
    x: float, y: float, z: float, t: float,
    coefficients: Sequence[BackgroundCoefficient], *, h: float = 0.005,
    cutoff: ScalarCutoff = standard_cutoff,
    cutoff_derivative: ScalarCutoff | None = None,
    quadrature_points: int = 32,
) -> np.ndarray:
    """Analytic curl of the cutoff-summed potential, plus direct swirl.

    The derivatives of chi(a_n*q) are essential for incompressibility in the
    transition zone. This instantiates the algebra, not the coefficient solver.
    Custom cutoffs must supply their derivative explicitly.
    """
    if cutoff_derivative is None:
        if cutoff is not standard_cutoff:
            raise ValueError("custom cutoff requires cutoff_derivative")
        cutoff_derivative = standard_cutoff_derivative
    x, y = _finite(x,"x"), _finite(y,"y")
    r = math.hypot(x,y)
    s = similarity_coordinates(r,z,t,h)
    total = np.zeros(3)
    for c in coefficients:
        w = 1.0 if c.n == 0 else _finite(cutoff(c.cutoff_scale*s.q),"cutoff")
        wq = 0.0 if c.n == 0 else c.cutoff_scale * _finite(
            cutoff_derivative(c.cutoff_scale*s.q),"cutoff derivative")
        if w == 0.0 and wq == 0.0:
            continue
        lam, p = c.lambda_n(h), c.profile
        AU = p.radial_average_U(s.X,s.eta,n=quadrature_points)
        v = p.radial_flux_factor(s.X,s.eta,h,lam=lam,n=quadrature_points,d=s.d,L=s.L)
        radial = 0.5*w*s.q**(lam-1)*v - wq*s.eta*s.q**lam*AU/s.L
        uz = w*s.q**(-s.A+lam)*_finite(p.U(s.X,s.eta),"U")
        swirl = np.zeros(3)
        if r != 0:
            E = _finite(p.E(s.X,s.eta),"E")
            swirl = w*s.q**(-s.A+lam)*E*np.array([-y/r,x/r,0.0])
        elif _finite(p.E(0,s.eta),"E") != 0:
            raise ValueError("E must vanish on the axis")
        total += np.array([radial*x,radial*y,uz]) + swirl
    if not np.all(np.isfinite(total)):
        raise OverflowError("background velocity is outside floating-point range")
    return total


def background_pressure(
    x: float, y: float, z: float, t: float,
    coefficients: Sequence[BackgroundCoefficient], *, h: float = 0.005,
    cutoff: ScalarCutoff = standard_cutoff,
) -> float:
    """Finite cutoff pressure sum, requiring every active Pi_n."""
    s = similarity_coordinates(math.hypot(x,y),z,t,h)
    total = 0.0
    for c in coefficients:
        w = 1.0 if c.n == 0 else _finite(cutoff(c.cutoff_scale*s.q),"cutoff")
        if w == 0:
            continue
        if c.profile.Pi is None:
            raise ValueError(f"pressure profile Pi_{c.n} is missing")
        total += w*s.q**(-2*s.A+c.lambda_n(h))*c.profile.Pi(s.X,s.eta)
    return _finite(total,"background pressure")
