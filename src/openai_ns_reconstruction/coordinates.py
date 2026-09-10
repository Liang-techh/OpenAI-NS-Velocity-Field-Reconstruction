"""Similarity coordinates used in OpenAI, Eq. (4.1).

For tau = 1-t, A = 1/2+h, D = 1/2-h, the paper introduces

    z   = q**D * eta
    tau = q * (1-eta**2)
    X   = r**2/(2*q).

Eliminating eta gives the scalar equation

    q - z**2 * q**(2*h) = tau,

whose derivative is L = 1 - 2*h*eta**2 > 0 for h < 1/2.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class SimilarityPoint:
    tau: float
    q: float
    eta: float
    X: float
    h: float
    A: float
    D: float
    d: float
    L: float


def _f(q: float, z: float, tau: float, h: float) -> float:
    return q - z * z * q ** (2.0 * h) - tau


def solve_q(z: float, t: float, h: float, *, rtol: float = 1e-13, max_iter: int = 200) -> float:
    """Solve the paper's implicit concentration scale q(z,t).

    Valid for t < 1 and 0 < h < 1/2.  The root is unique because
    d/dq [q-z^2 q^(2h)] = 1-2h eta^2 > 0 on the physical branch.
    A monotone bisection is used so the implementation stays robust close
    to the singular time.
    """
    if not (0.0 < h < 0.5):
        raise ValueError("h must satisfy 0 < h < 1/2")
    tau = 1.0 - float(t)
    if tau <= 0.0:
        raise ValueError("the similarity chart is defined for t < 1")

    z = float(z)
    if z == 0.0:
        return tau

    D = 0.5 - h
    q_axis = abs(z) ** (1.0 / D)
    lo = max(tau, q_axis)
    # The physical root lies strictly above both tau and the eta=+-1 endpoint.
    # f(lo) can be <= 0, so grow a safe upper bracket.
    hi = max(2.0 * lo, lo + 1.0)
    while _f(hi, z, tau, h) <= 0.0:
        hi *= 2.0

    flo = _f(lo, z, tau, h)
    if flo > 0.0:
        # Roundoff can make the theoretical lower bracket slightly positive.
        lo *= max(0.5, 1.0 - 1e-12)

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        fm = _f(mid, z, tau, h)
        if fm <= 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= rtol * max(1.0, abs(mid)):
            return 0.5 * (lo + hi)
    raise RuntimeError("solve_q did not converge")


def similarity_coordinates(r: float, z: float, t: float, h: float) -> SimilarityPoint:
    """Return (tau,q,eta,X,A,D,d,L) for the paper's similarity chart."""
    r = abs(float(r))
    q = solve_q(z, t, h)
    tau = 1.0 - float(t)
    A = 0.5 + h
    D = 0.5 - h
    eta = float(z) / (q ** D)
    X = (r * r) / (2.0 * q)
    d = 1.0 - eta * eta
    L = 1.0 - 2.0 * h * eta * eta

    # Numerical guard: exact mathematics has |eta| < 1 for tau > 0.
    if abs(eta) > 1.0 + 1e-10:
        raise ArithmeticError("computed eta left the physical chart")
    return SimilarityPoint(tau, q, eta, X, h, A, D, d, L)


def coordinate_identity_error(r: float, z: float, t: float, h: float) -> tuple[float, float]:
    """Absolute errors in z=q^D eta and tau=q(1-eta^2)."""
    s = similarity_coordinates(r, z, t, h)
    e_z = abs(float(z) - s.q ** s.D * s.eta)
    e_tau = abs(s.tau - s.q * (1.0 - s.eta * s.eta))
    return e_z, e_tau
