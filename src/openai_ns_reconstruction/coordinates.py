"""Similarity geometry, Eq. (4.1); binary64 numerics are not proof certificates.

Use the ``*_from_tau`` entry points near the singularity: forming t=1-tau
loses tau altogether below machine resolution. d=tau/q is stored separately
from eta, since 1-eta**2 can also suffer cancellation near the chart endpoints.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import sys


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


def validate_h(h: float) -> float:
    h = float(h)
    if not math.isfinite(h) or not 0 < h < 0.5:
        raise ValueError("chart requires finite 0 < h < 1/2; full paper uses h < 1/100")
    return h


def solve_q_from_tau(z: float, tau: float, h: float, *, rtol: float = 1e-13,
                     max_iter: int = 200) -> float:
    """Solve q-z^2 q^(2h)=tau using a dimensionless, bracketed bisection.

    Scale by max(tau, |z|^(1/D)); this avoids a unit-sized bracket and an
    absolute 1e-13 stopping threshold for roots many orders smaller than 1.
    expm1 evaluates the equation without subtracting almost equal terms.
    """
    z, tau, h = float(z), float(tau), validate_h(h)
    if not math.isfinite(z) or not math.isfinite(tau) or tau <= 0:
        raise ValueError("z must be finite and tau must be finite and positive")
    if not math.isfinite(rtol) or not 4 * sys.float_info.epsilon <= rtol < 1:
        raise ValueError("rtol must be at least 4 machine epsilons and smaller than 1")
    if isinstance(max_iter, bool) or not isinstance(max_iter, int) or max_iter < 1:
        raise ValueError("max_iter must be a positive integer")
    if z == 0:
        return tau
    p = 2 * h
    log_endpoint = 2 * math.log(abs(z)) / (1 - p)
    log_scale = max(math.log(tau), log_endpoint)
    try:
        scale = math.exp(log_scale)
    except OverflowError as exc:
        raise ArithmeticError("q lies outside binary64 range") from exc
    if not math.isfinite(scale) or scale <= 0:
        raise ArithmeticError("q scale is not representable")
    a = tau / scale
    log_b = min(0.0, 2 * math.log(abs(z)) - (1 - p) * math.log(scale))

    def residual(y: float) -> float:
        return y * (-math.expm1(log_b + (p - 1) * math.log(y))) - a

    lo, hi = 1.0, 2.0
    if residual(lo) >= 0:  # endpoint and root coincide to floating precision
        return scale
    for _ in range(max_iter):
        if residual(hi) > 0:
            break
        hi *= 2
        if not math.isfinite(hi):
            raise ArithmeticError("unable to bracket q in binary64")
    else:
        raise RuntimeError("solve_q could not bracket the root")
    for _ in range(max_iter):
        mid = lo + (hi - lo) / 2
        if residual(mid) <= 0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= rtol * mid:
            q = scale * (lo + (hi - lo) / 2)
            if not math.isfinite(q) or q <= 0:
                raise ArithmeticError("q is not representable")
            return q
    raise RuntimeError("solve_q did not converge")


def solve_q(z: float, t: float, h: float, *, rtol: float = 1e-13,
            max_iter: int = 200) -> float:
    """Compatibility entry point for finite t<1."""
    return solve_q_from_tau(z, 1 - float(t), h, rtol=rtol, max_iter=max_iter)


def similarity_coordinates_from_tau(r: float, z: float, tau: float,
                                    h: float) -> SimilarityPoint:
    r, z, tau = float(r), float(z), float(tau)
    h = validate_h(h)
    if not math.isfinite(r) or r < 0:
        raise ValueError("cylindrical radius must be finite and nonnegative")
    q = solve_q_from_tau(z, tau, h)
    A, D = 0.5 + h, 0.5 - h
    eta = z / q**D
    if abs(eta) > 1 + 1e-12:
        raise ArithmeticError("computed eta left the physical chart")
    eta = min(1.0, max(-1.0, eta))  # endpoint rounding only
    X = 0.5 * (r / math.sqrt(q))**2
    if not math.isfinite(X):
        raise ArithmeticError("X lies outside binary64 range")
    return SimilarityPoint(tau, q, eta, X, h, A, D, tau / q,
                           1 - 2 * h * eta**2)


def similarity_coordinates(r: float, z: float, t: float, h: float) -> SimilarityPoint:
    return similarity_coordinates_from_tau(r, z, 1 - float(t), h)


def coordinate_identity_error(r: float, z: float, t: float,
                              h: float) -> tuple[float, float]:
    """Absolute errors; the second includes endpoint cancellation in 1-eta^2."""
    s = similarity_coordinates(r, z, t, h)
    return abs(z - s.q**s.D * s.eta), abs(s.tau - s.q * (1 - s.eta**2))
