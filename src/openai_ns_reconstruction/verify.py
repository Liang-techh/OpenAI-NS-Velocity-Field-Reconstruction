"""Numerical diagnostics for reconstructed fields.

These checks are diagnostics, not substitutes for the Lean proof. They are intended to catch
implementation mistakes while the executable reconstruction is being built.

Forced Navier--Stokes sign convention used here:

    u_t + (u . grad) u - nu Delta u + grad p = f.

Accordingly, ``navier_stokes_residual_numeric`` is the force required by a supplied ``(u,p)``
pair, and ``forced_ns_closure_error_numeric`` measures residual minus a candidate force.
"""

from __future__ import annotations

from typing import Callable
import numpy as np

VectorField = Callable[[float, float, float, float], np.ndarray]
ScalarField = Callable[[float, float, float, float], float]


def divergence_numeric(u: VectorField, x: float, y: float, z: float, t: float, *, eps: float = 1e-5) -> float:
    """Centered-difference divergence."""
    out = 0.0
    for axis in range(3):
        p = [x, y, z]
        m = [x, y, z]
        p[axis] += eps
        m[axis] -= eps
        up = np.asarray(u(*p, t), dtype=float)
        um = np.asarray(u(*m, t), dtype=float)
        out += (up[axis] - um[axis]) / (2.0 * eps)
    return float(out)


def laplacian_vector_numeric(u: VectorField, x: float, y: float, z: float, t: float, *, eps: float = 1e-4) -> np.ndarray:
    center = np.asarray(u(x, y, z, t), dtype=float)
    out = np.zeros(3)
    for axis in range(3):
        p = [x, y, z]
        m = [x, y, z]
        p[axis] += eps
        m[axis] -= eps
        out += (np.asarray(u(*p, t), float) - 2.0 * center + np.asarray(u(*m, t), float)) / (eps * eps)
    return out


def time_derivative_numeric(u: VectorField, x: float, y: float, z: float, t: float, *, eps: float = 1e-5) -> np.ndarray:
    return (
        np.asarray(u(x, y, z, t + eps), dtype=float)
        - np.asarray(u(x, y, z, t - eps), dtype=float)
    ) / (2.0 * eps)


def jacobian_numeric(u: VectorField, x: float, y: float, z: float, t: float, *, eps: float = 1e-5) -> np.ndarray:
    J = np.zeros((3, 3))
    for axis in range(3):
        p = [x, y, z]
        m = [x, y, z]
        p[axis] += eps
        m[axis] -= eps
        J[:, axis] = (np.asarray(u(*p, t), float) - np.asarray(u(*m, t), float)) / (2.0 * eps)
    return J


def gradient_scalar_numeric(p: ScalarField, x: float, y: float, z: float, t: float, *, eps: float = 1e-5) -> np.ndarray:
    g = np.zeros(3)
    for axis in range(3):
        xp = [x, y, z]
        xm = [x, y, z]
        xp[axis] += eps
        xm[axis] -= eps
        g[axis] = (float(p(*xp, t)) - float(p(*xm, t))) / (2.0 * eps)
    return g


def navier_stokes_residual_numeric(
    u: VectorField,
    p: ScalarField,
    x: float,
    y: float,
    z: float,
    t: float,
    *,
    viscosity: float = 1.0,
    eps_space: float = 1e-4,
    eps_time: float = 1e-5,
) -> np.ndarray:
    """Return ``u_t + (u.grad)u - nu Delta u + grad p``.

    Under the repository's forced-NS convention this quantity is exactly the force ``f``
    required by the supplied velocity and pressure. The routine is numerical and therefore
    intended as a reconstruction diagnostic rather than a proof object.
    """
    uv = np.asarray(u(x, y, z, t), dtype=float)
    ut = time_derivative_numeric(u, x, y, z, t, eps=eps_time)
    J = jacobian_numeric(u, x, y, z, t, eps=eps_space)
    adv = J @ uv
    lap = laplacian_vector_numeric(u, x, y, z, t, eps=eps_space)
    gp = gradient_scalar_numeric(p, x, y, z, t, eps=eps_space)
    return ut + adv - float(viscosity) * lap + gp


def reconstruct_forcing_numeric(
    u: VectorField,
    p: ScalarField,
    x: float,
    y: float,
    z: float,
    t: float,
    *,
    viscosity: float = 1.0,
    eps_space: float = 1e-4,
    eps_time: float = 1e-5,
) -> np.ndarray:
    """Numerically reconstruct the force required by a candidate ``(u,p)`` field.

    This is a named Stage-8 entry point for Section 10 work. It deliberately delegates to the
    residual implementation so the force convention cannot silently diverge from verification.
    """
    return navier_stokes_residual_numeric(
        u,
        p,
        x,
        y,
        z,
        t,
        viscosity=viscosity,
        eps_space=eps_space,
        eps_time=eps_time,
    )


def forced_ns_closure_error_numeric(
    u: VectorField,
    p: ScalarField,
    forcing: VectorField,
    x: float,
    y: float,
    z: float,
    t: float,
    *,
    viscosity: float = 1.0,
    eps_space: float = 1e-4,
    eps_time: float = 1e-5,
) -> np.ndarray:
    """Return the defect in the forced equation, ``R(u,p) - f``.

    Once the paper-exact localized field and smooth force are instantiated, this provides an
    independent numerical closure check without changing either field.
    """
    residual = navier_stokes_residual_numeric(
        u,
        p,
        x,
        y,
        z,
        t,
        viscosity=viscosity,
        eps_space=eps_space,
        eps_time=eps_time,
    )
    f = np.asarray(forcing(x, y, z, t), dtype=float)
    if f.shape != (3,):
        raise ValueError("forcing must return a length-3 vector")
    return residual - f


def loglog_slope(xs: np.ndarray, ys: np.ndarray) -> float:
    """Least-squares slope of log|y| against log x."""
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if np.any(xs <= 0) or np.any(np.abs(ys) <= 0):
        raise ValueError("xs and |ys| must be positive")
    return float(np.polyfit(np.log(xs), np.log(np.abs(ys)), 1)[0])
