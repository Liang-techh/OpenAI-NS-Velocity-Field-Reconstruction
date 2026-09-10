"""Finite-difference diagnostics, NOT a proof of smooth forcing or NS breakdown.

R(u,p)=u_t+(u.grad)u-nu*Delta(u)+grad(p). Computing f=R and then
checking R-f with the same stencil is a tautology, not independent evidence.
The manufactured-solution tests instead supply an independently derived force.
"""
from __future__ import annotations
from typing import Callable
import math
import numpy as np

VectorField = Callable[[float, float, float, float], np.ndarray]
ScalarField = Callable[[float, float, float, float], float]
TimeBounds = tuple[float | None, float | None] | None


def _vector(value) -> np.ndarray:
    a = np.asarray(value, dtype=float)
    if a.shape != (3,) or not np.all(np.isfinite(a)):
        raise ValueError("field must return a finite length-3 vector")
    return a


def _step(eps: float) -> float:
    eps = float(eps)
    if not math.isfinite(eps) or eps <= 0:
        raise ValueError("finite-difference step must be finite and positive")
    return eps


def _point(x, y, z, t):
    values = tuple(float(v) for v in (x, y, z, t))
    if not all(math.isfinite(v) for v in values):
        raise ValueError("evaluation point must be finite")
    return values


def _shift(xyz, axis, eps):
    plus, minus = list(xyz), list(xyz)
    plus[axis] += eps
    minus[axis] -= eps
    if plus[axis] == xyz[axis] or minus[axis] == xyz[axis]:
        raise ValueError("finite-difference step is below coordinate resolution")
    return plus, minus


def jacobian_numeric(u: VectorField, x: float, y: float, z: float, t: float,
                     *, eps: float = 1e-5) -> np.ndarray:
    x, y, z, t = _point(x, y, z, t)
    eps = _step(eps)
    J = np.empty((3, 3))
    for axis in range(3):
        p, m = _shift((x, y, z), axis, eps)
        J[:, axis] = (_vector(u(*p, t)) - _vector(u(*m, t))) / (2 * eps)
    return J


def divergence_numeric(u: VectorField, x: float, y: float, z: float, t: float,
                       *, eps: float = 1e-5) -> float:
    return float(np.trace(jacobian_numeric(u, x, y, z, t, eps=eps)))


def laplacian_vector_numeric(u: VectorField, x: float, y: float, z: float, t: float,
                             *, eps: float = 1e-4) -> np.ndarray:
    x, y, z, t = _point(x, y, z, t)
    eps = _step(eps)
    center, out = _vector(u(x, y, z, t)), np.zeros(3)
    for axis in range(3):
        p, m = _shift((x, y, z), axis, eps)
        out += (_vector(u(*p, t)) - 2 * center + _vector(u(*m, t))) / eps**2
    return out


def time_derivative_numeric(u: VectorField, x: float, y: float, z: float, t: float,
                            *, eps: float = 1e-5,
                            time_bounds: TimeBounds = None) -> np.ndarray:
    """Second-order stencil; bounds are [lower,upper), None means unbounded.

    Uses a forward stencil at the included lower endpoint and a shortened
    central stencil inside the domain. Never evaluates at the excluded upper
    endpoint when bounds are supplied. Too-small floating steps fail loudly.
    """
    x, y, z, t = _point(x, y, z, t)
    dt = _step(eps)
    forward = False
    if time_bounds is not None:
        lower, upper = time_bounds
        if any(b is not None and not math.isfinite(b) for b in time_bounds):
            raise ValueError("time bounds must be finite or None")
        if lower is not None and upper is not None and lower >= upper:
            raise ValueError("lower time bound must be smaller than upper")
        if (lower is not None and t < lower) or (upper is not None and t >= upper):
            raise ValueError("time lies outside [lower, upper)")
        if upper is not None:
            dt = min(dt, (upper - t) / 4)
        if lower is not None:
            forward = t == lower
            if not forward:
                dt = min(dt, (t - lower) / 4)
    if t + dt == t or t - dt == t:
        raise ValueError("time stencil is below floating resolution; use tau/analytic data")
    if forward:
        return (-3 * _vector(u(x, y, z, t)) + 4 * _vector(u(x, y, z, t + dt))
                - _vector(u(x, y, z, t + 2 * dt))) / (2 * dt)
    return (_vector(u(x, y, z, t + dt)) - _vector(u(x, y, z, t - dt))) / (2 * dt)


def gradient_scalar_numeric(p: ScalarField, x: float, y: float, z: float, t: float,
                            *, eps: float = 1e-5) -> np.ndarray:
    x, y, z, t = _point(x, y, z, t)
    eps = _step(eps)
    g = np.empty(3)
    for axis in range(3):
        xp, xm = _shift((x, y, z), axis, eps)
        vp, vm = float(p(*xp, t)), float(p(*xm, t))
        if not math.isfinite(vp) or not math.isfinite(vm):
            raise ValueError("pressure must be finite")
        g[axis] = (vp - vm) / (2 * eps)
    return g


def navier_stokes_residual_numeric(u: VectorField, p: ScalarField,
                                   x: float, y: float, z: float, t: float, *,
                                   viscosity: float = 1.0, eps_space: float = 1e-4,
                                   eps_time: float = 1e-5,
                                   time_bounds: TimeBounds = (0.0, 1.0)) -> np.ndarray:
    """Required force R(u,p); default time domain matches this repository.

    Set time_bounds=None for manufactured fields defined on all real times.
    Step-refinement checks are required before interpreting a numerical value.
    """
    if not math.isfinite(viscosity) or viscosity <= 0:
        raise ValueError("Navier-Stokes viscosity must be finite and positive")
    point = _point(x, y, z, t)
    uv = _vector(u(*point))
    ut = time_derivative_numeric(u, *point, eps=eps_time, time_bounds=time_bounds)
    J = jacobian_numeric(u, *point, eps=eps_space)
    lap = laplacian_vector_numeric(u, *point, eps=eps_space)
    gp = gradient_scalar_numeric(p, *point, eps=eps_space)
    return ut + J @ uv - viscosity * lap + gp


def reconstruct_forcing_numeric(u: VectorField, p: ScalarField,
                                 x: float, y: float, z: float, t: float, *,
                                 viscosity: float = 1.0, eps_space: float = 1e-4,
                                 eps_time: float = 1e-5,
                                 time_bounds: TimeBounds = (0.0, 1.0)) -> np.ndarray:
    return navier_stokes_residual_numeric(
        u, p, x, y, z, t, viscosity=viscosity, eps_space=eps_space,
        eps_time=eps_time, time_bounds=time_bounds)


def forced_ns_closure_error_numeric(u: VectorField, p: ScalarField,
                                    forcing: VectorField, x: float, y: float,
                                    z: float, t: float, *, viscosity: float = 1.0,
                                    eps_space: float = 1e-4, eps_time: float = 1e-5,
                                    time_bounds: TimeBounds = (0.0, 1.0)) -> np.ndarray:
    return navier_stokes_residual_numeric(
        u, p, x, y, z, t, viscosity=viscosity, eps_space=eps_space,
        eps_time=eps_time, time_bounds=time_bounds) - _vector(forcing(x, y, z, t))


def loglog_slope(xs: np.ndarray, ys: np.ndarray) -> float:
    xs, ys = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    if xs.ndim != 1 or xs.shape != ys.shape or xs.size < 2:
        raise ValueError("xs and ys must be equal-length 1D arrays with at least 2 samples")
    if (not np.all(np.isfinite(xs)) or not np.all(np.isfinite(ys))
            or np.any(xs <= 0) or np.any(ys == 0) or np.unique(xs).size < 2):
        raise ValueError("need finite xs>0, nonzero ys, and distinct xs")
    return float(np.polyfit(np.log(xs), np.log(np.abs(ys)), 1)[0])
