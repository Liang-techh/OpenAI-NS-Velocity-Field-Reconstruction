"""Executable Section 10 time-localization adapter.

Provenance
----------
The pinned OpenAI formalization defines in ``SmoothCutoffs.lean``

    timeSwitch(t) = 1 - cutoff((4/3) t),

where ``cutoff`` is Mathlib's noncomputable C-infinity ``ContDiffBump`` with
inner radius 1/2 and outer radius 1. ``TimeLocalization.lean`` then sets
``v = chi u`` and ``q = chi p`` and proves

    R(v,q) = chi R(u,p) + chi' u + (chi^2-chi) (u.grad)u.

This module uses the repository's explicit C-infinity representative of the
same cutoff support/plateau geometry. Therefore ``time_switch`` is exactly zero
for |t| <= 3/8 and exactly one for |t| >= 3/4 (in particular for t >= 3/4),
but its transition-collar values are not claimed equal to Mathlib's
noncomputable bump pointwise.

This is only the time-switch part of the construction. It does not construct
the unresolved incoming singular fields and does not provide the separate
smooth global force extension through t=1.
"""
from __future__ import annotations

from typing import Callable
import math
import numpy as np

from .coordinates import _finite
from .cutoffs import smooth_cutoff, standard_cutoff_derivative
from .verify import (
    finite_vector,
    jacobian_numeric,
    navier_stokes_residual_numeric,
)

VectorField = Callable[[float, float, float, float], np.ndarray]
ScalarField = Callable[[float, float, float, float], float]

TIME_SCALE = 4.0 / 3.0
EARLY_HALF_WIDTH = 3.0 / 8.0
LATE_START = 3.0 / 4.0


def _bump_even(s: float) -> float:
    s = _finite(s, "time cutoff argument")
    return smooth_cutoff(abs(s))


def _bump_even_derivative(s: float) -> float:
    s = _finite(s, "time cutoff argument")
    if s == 0.0:
        return 0.0
    sign = 1.0 if s > 0.0 else -1.0
    return sign * standard_cutoff_derivative(abs(s))


def time_switch(t: float) -> float:
    """Executable representative of the pinned formalization's ``timeSwitch``."""
    t = _finite(t, "t")
    value = 1.0 - _bump_even(TIME_SCALE * t)
    if not math.isfinite(value):
        raise ArithmeticError("time switch is not finite")
    return value


def time_switch_derivative(t: float) -> float:
    """Analytic derivative of this repository's explicit switch representative."""
    t = _finite(t, "t")
    value = -TIME_SCALE * _bump_even_derivative(TIME_SCALE * t)
    if not math.isfinite(value):
        raise ArithmeticError("time-switch derivative is not finite")
    return value


def activated_velocity(u: VectorField) -> VectorField:
    """Return the time-localized velocity ``chi(t) u``."""
    def field(x: float, y: float, z: float, t: float) -> np.ndarray:
        chi = time_switch(t)
        if chi == 0.0:
            return np.zeros(3)
        return chi * finite_vector(u(x, y, z, t))
    return field


def activated_pressure(p: ScalarField) -> ScalarField:
    """Return the time-localized pressure ``chi(t) p``."""
    def field(x: float, y: float, z: float, t: float) -> float:
        chi = time_switch(t)
        if chi == 0.0:
            return 0.0
        value = chi * _finite(p(x, y, z, t), "pressure")
        if not math.isfinite(value):
            raise ArithmeticError("activated pressure is not finite")
        return value
    return field


def activated_residual_formula_numeric(
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
    """Evaluate the exact time-activation residual identity's right-hand side.

    Spatial and original-field time derivatives are numerical diagnostics, but
    ``chi`` and ``chi'`` are evaluated analytically. This function is kept
    separate from evaluating the residual of the activated fields so tests can
    cross-check the formal identity without defining the expected answer as the
    same residual computation.
    """
    x, y, z, t = (_finite(v, name) for v, name in zip(
        (x, y, z, t), ("x", "y", "z", "t")
    ))
    chi = time_switch(t)
    dchi = time_switch_derivative(t)
    uv = finite_vector(u(x, y, z, t))
    original_residual = navier_stokes_residual_numeric(
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
    advection = jacobian_numeric(u, x, y, z, t, eps=eps_space) @ uv
    rhs = chi * original_residual + dchi * uv + (chi * chi - chi) * advection
    return finite_vector(rhs)
