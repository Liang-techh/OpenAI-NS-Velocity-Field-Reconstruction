"""Section 5 background assembly; recursive coefficients are still caller inputs.

lambda_n=2nh. Cutoffs apply to potentials BEFORE curl. The supplied C-infinity
cutoff is an experimental choice, not the paper's recursively certified schedule.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Sequence
import math
import numpy as np
from .coordinates import similarity_coordinates, validate_h
from .profiles import LeadingProfile
from .cutoffs import smooth_cutoff

ScalarCutoff = Callable[[float], float]
standard_cutoff = smooth_cutoff  # backwards-compatible name; formerly only C1


@dataclass(frozen=True)
class BackgroundCoefficient:
    n: int
    profile: LeadingProfile
    cutoff_scale: float = 1.0

    def __post_init__(self):
        if isinstance(self.n, bool) or not isinstance(self.n, int) or self.n < 0:
            raise ValueError("coefficient index n must be a nonnegative integer")
        if not math.isfinite(self.cutoff_scale) or self.cutoff_scale <= 0:
            raise ValueError("cutoff_scale must be finite and positive")

    def lambda_n(self, h: float) -> float:
        return 2 * self.n * validate_h(h)


def coefficient_streamfunction(r: float, z: float, t: float,
                                coefficient: BackgroundCoefficient, *, h: float,
                                quadrature_points: int = 32) -> float:
    """S_n=q^(1-A+lambda_n) X A_X(U_n), Eq. (5.27)."""
    s = similarity_coordinates(r, z, t, h)
    AU = coefficient.profile.radial_average_U(s.X, s.eta, n=quadrature_points)
    return s.q**(1 - s.A + coefficient.lambda_n(h)) * s.X * AU


def coefficient_vector_potential_cartesian(x: float, y: float, z: float, t: float,
                                            coefficient: BackgroundCoefficient, *,
                                            h: float, quadrature_points: int = 32) -> np.ndarray:
    # S_n/r^2 = 1/2 q^(-A+lambda_n) A_X(U_n), regular even at r=0.
    s = similarity_coordinates(math.hypot(x, y), z, t, h)
    AU = coefficient.profile.radial_average_U(s.X, s.eta, n=quadrature_points)
    factor = 0.5 * s.q**(-s.A + coefficient.lambda_n(h)) * AU
    return factor * np.array([-y, x, 0.0])


def background_potential_cartesian(x: float, y: float, z: float, t: float,
                                    coefficients: Sequence[BackgroundCoefficient], *,
                                    h: float = 0.005,
                                    cutoff: ScalarCutoff = standard_cutoff) -> np.ndarray:
    s = similarity_coordinates(math.hypot(x, y), z, t, h)
    total = np.zeros(3)
    for c in coefficients:
        weight = 1.0 if c.n == 0 else float(cutoff(c.cutoff_scale * s.q))
        if weight:
            total += weight * coefficient_vector_potential_cartesian(x, y, z, t, c, h=h)
    return total


def background_swirl_cartesian(x: float, y: float, z: float, t: float,
                                coefficients: Sequence[BackgroundCoefficient], *,
                                h: float = 0.005,
                                cutoff: ScalarCutoff = standard_cutoff) -> np.ndarray:
    r = math.hypot(x, y)
    s = similarity_coordinates(r, z, t, h)
    B = 0.0
    for c in coefficients:
        weight = 1.0 if c.n == 0 else float(cutoff(c.cutoff_scale * s.q))
        if weight:
            E = float(c.profile.E(s.X, s.eta))
            if r == 0 and E != 0:
                raise ValueError("regular swirl must vanish on the axis")
            B += weight * s.q**(-s.A + c.lambda_n(h)) * E
    return np.zeros(3) if r == 0 else B * np.array([-y / r, x / r, 0.0])


def coefficient_radial_flux(
    X: float,
    eta: float,
    coefficient: BackgroundCoefficient,
    *,
    h: float,
    quadrature_points: int = 32,
) -> float:
    """Return the Section 5 radial-flux coefficient ``V_n(X, eta)``.

    This is the second identity in Eq. (5.2):

        V_n / X = [2 eta U_n
                   - 2 eta (D + lambda_n) A_X(U_n)
                   - d partial_eta A_X(U_n)] / L.

    It is the exact kinematic consequence of incompressibility for the nth formal
    coefficient. It does not construct the recursive positive-order profiles.
    The concurrent coefficient formula is preserved, using the improved average
    backend and its default quadrature order. For n=0 this reduces to Eq. (4.7).
    """
    from .profiles import _point

    X, eta = _point(X, eta)
    h = validate_h(h)
    if X == 0.0:
        return 0.0
    D = 0.5 - h
    d = 1.0 - eta * eta
    L = 1.0 - 2.0 * h * eta * eta
    lam = coefficient.lambda_n(h)
    U = float(coefficient.profile.U(X, eta))
    AU = coefficient.profile.radial_average_U(X, eta, n=quadrature_points)
    dAU = coefficient.profile.radial_average_dU_deta(X, eta, n=quadrature_points)
    value = (X / L) * (2.0 * eta * U - 2.0 * eta * (D + lam) * AU - d * dAU)
    if not math.isfinite(value):
        raise ArithmeticError("coefficient radial flux is not finite")
    return value
