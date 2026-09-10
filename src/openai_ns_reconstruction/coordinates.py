"""Similarity chart, Eqs. (4.1)-(4.2); floating-point evaluation, not a proof.

Use the ``*_from_tau`` functions near t=1: forming t=1-tau loses tau below
machine precision. We solve on the physical branch in log-scaled coordinates
and keep d=tau/q, rather than subtracting two almost equal numbers.
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


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _validate_h(h: float) -> float:
    h = _finite(h, "h")
    if not 0.0 < h < 0.5:
        raise ValueError("h must satisfy 0 < h < 1/2 (geometry domain)")
    return h


# Public compatibility name used by the newer Stage-1/2 modules on main.
def validate_h(h: float) -> float:
    return _validate_h(h)


def _tau(t: float) -> float:
    t = _finite(t, "t")
    if t >= 1.0:
        raise ValueError("t must be less than 1; use the tau API near t=1")
    return 1.0 - t


def solve_q_from_tau(z: float, tau: float, h: float, *,
                     rtol: float = 1e-13, max_iter: int = 200) -> float:
    """Solve q-z^2*q^(2h)=tau, with a relative (not unit-scale) tolerance."""
    z, tau, h = _finite(z, "z"), _finite(tau, "tau"), _validate_h(h)
    rtol = _finite(rtol, "rtol")
    if tau <= 0:
        raise ValueError("tau must be positive")
    if not 4 * sys.float_info.epsilon <= rtol < 1:
        raise ValueError("rtol must be between 4*machine_epsilon and 1")
    if isinstance(max_iter, bool) or not isinstance(max_iter, int) or max_iter < 1:
        raise ValueError("max_iter must be a positive integer")
    if z == 0:
        return tau
    k = 1.0 - 2.0 * h
    log_tau = math.log(tau)
    log_axis = 2.0 * math.log(abs(z)) / k
    log_scale = max(log_tau, log_axis)
    if log_scale > math.log(sys.float_info.max):
        raise OverflowError("q is outside floating-point range")
    scale = tau if log_tau >= log_axis else math.exp(log_axis)
    log_gamma = k * (log_axis - log_scale)
    log_a = log_tau - log_scale

    def residual(y: float) -> float:
        return -math.expm1(log_gamma - k * y) - math.exp(log_a - y)

    if residual(0.0) == 0.0:
        return scale
    lo, hi = 0.0, math.log(2.0)
    for _ in range(2048):
        if residual(hi) > 0:
            break
        hi *= 2.0
        if log_scale + hi > math.log(sys.float_info.max):
            hi = math.log(sys.float_info.max) - log_scale
            if hi <= 0 or residual(hi) <= 0:
                raise OverflowError("q is outside floating-point range")
            break
    else:
        raise RuntimeError("could not bracket q")
    for _ in range(max_iter):
        mid = lo + (hi - lo) * 0.5
        if residual(mid) <= 0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= math.log1p(rtol):
            y = (lo + hi) * 0.5
            q = scale * math.exp(y)
            if not math.isfinite(q) or q <= 0:
                raise OverflowError("q is outside floating-point range")
            return q
    raise RuntimeError("solve_q did not converge; increase max_iter")


def solve_q(z: float, t: float, h: float, *,
            rtol: float = 1e-13, max_iter: int = 200) -> float:
    return solve_q_from_tau(z, _tau(t), h, rtol=rtol, max_iter=max_iter)


def similarity_coordinates_from_tau(r: float, z: float, tau: float,
                                    h: float) -> SimilarityPoint:
    r = _finite(r, "r")
    if r < 0:
        raise ValueError("cylindrical radius r must be nonnegative")
    q = solve_q_from_tau(z, tau, h)
    h, tau, z = float(h), float(tau), float(z)
    A, D = 0.5 + h, 0.5 - h
    eta = z / q**D
    if abs(eta) > 1.0 + 1e-10:
        raise ArithmeticError("computed eta left the physical chart")
    d = tau / q
    L = (1.0 - 2.0 * h) + 2.0 * h * d
    X = 0.5 * (r / math.sqrt(q))**2
    if not math.isfinite(X):
        raise OverflowError("X is outside floating-point range")
    return SimilarityPoint(tau, q, eta, X, h, A, D, d, L)


def similarity_coordinates(r: float, z: float, t: float, h: float) -> SimilarityPoint:
    return similarity_coordinates_from_tau(r, z, _tau(t), h)


def coordinate_identity_error(r: float, z: float, t: float, h: float) -> tuple[float, float]:
    s = similarity_coordinates(r, z, t, h)
    return abs(float(z) - s.q**s.D * s.eta), abs(s.tau - s.q * (1 - s.eta**2))


def coordinate_derivatives(s: SimilarityPoint) -> dict[str, float]:
    """Analytic q, eta and X derivatives from Lemma 4.1, at fixed r,z,t."""
    qL, qDL = s.q * s.L, s.q**s.D * s.L
    return {"q_t": -1.0 / s.L, "q_z": 2 * s.eta * s.q / qDL,
            "eta_t": s.D * s.eta / qL, "eta_z": s.d / qDL,
            "X_t": s.X / qL, "X_z": -2 * s.eta * s.X / qDL,
            "X_r": math.sqrt(2 * s.X) / math.sqrt(s.q)}


def weighted_profile_derivatives(s: SimilarityPoint, b: float, f: float,
                                 f_X: float, f_eta: float) -> dict[str, float]:
    """Derivatives of q**b*f(X,eta), given analytic profile partials, Eq. (4.2)."""
    b, f = _finite(b, "b"), _finite(f, "f")
    f_X, f_eta = _finite(f_X, "f_X"), _finite(f_eta, "f_eta")
    T = (-b * f + s.D * s.eta * f_eta + s.X * f_X) / s.L
    Z = (2 * b * s.eta * f + s.d * f_eta - 2 * s.eta * s.X * f_X) / s.L
    return {"t": s.q**(b - 1) * T, "z": s.q**(b - s.D) * Z,
            "r": s.q**(b - 0.5) * math.sqrt(2 * s.X) * f_X}
