"""Closed-past preparation for the Section 10 endpoint construction.

The pinned OpenAI formalization's ``PastExtension.lean`` does not assume any
regularity of the incoming fields for negative time.  Instead it first applies
the concrete ``TimeLocalization.timeSwitch`` and then replaces all negative-
time values by zero.  The zero germ around ``t=0`` makes that replacement
smooth in the theorem, and ``CandidateFromLimits`` subsequently takes the
actual Navier--Stokes residual of these closed-past fields.

This module makes the algebraic branch construction executable.  It deliberately
does not claim a smoothness proof: the repository's time switch uses an explicit
C-infinity representative with the same support/plateau geometry as Mathlib's
noncomputable ``ContDiffBump``, and ``past_residual_numeric`` is a finite-
difference diagnostic rather than an all-order residual-derivative object.
"""
from __future__ import annotations

from collections.abc import Callable

import numpy as np

from .coordinates import _finite
from .time_localization import activated_pressure, activated_velocity
from .verify import finite_vector, navier_stokes_residual_numeric

VectorField = Callable[[float, float, float, float], np.ndarray]
ScalarField = Callable[[float, float, float, float], float]


def zero_before_vector(field: VectorField) -> VectorField:
    """Preserve ``t>=0`` and return exactly zero for ``t<0``.

    This is the executable coordinate-order analogue of
    ``PastExtension.zeroBefore`` for velocity-valued fields.  The negative
    branch is short-circuited, so no unavailable negative-time input is sampled.
    """
    if not callable(field):
        raise ValueError("field must be callable")

    def extended(x: float, y: float, z: float, t: float) -> np.ndarray:
        t = _finite(t, "t")
        if t < 0.0:
            return np.zeros(3)
        return finite_vector(field(x, y, z, t))

    return extended


def zero_before_scalar(field: ScalarField) -> ScalarField:
    """Scalar counterpart of :func:`zero_before_vector`."""
    if not callable(field):
        raise ValueError("field must be callable")

    def extended(x: float, y: float, z: float, t: float) -> float:
        t = _finite(t, "t")
        if t < 0.0:
            return 0.0
        return _finite(field(x, y, z, t), "scalar field")

    return extended


def past_velocity(u: VectorField) -> VectorField:
    """Pinned ``pastVelocity`` branch: ``zeroBefore (activatedVelocity u)``."""
    return zero_before_vector(activated_velocity(u))


def past_pressure(p: ScalarField) -> ScalarField:
    """Pinned ``pastPressure`` branch: ``zeroBefore (activatedPressure p)``."""
    return zero_before_scalar(activated_pressure(p))


def past_residual_numeric(
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
    endpoint: float = 1.0,
) -> np.ndarray:
    """Numerically evaluate the actual residual of the closed-past fields.

    ``PastExtension.pastResidual`` is an exact theorem-level definition.  This
    routine evaluates the same operator with the repository's independent
    finite-difference diagnostics.  It is therefore useful for regression and
    for later residual-jet plumbing, but it is not evidence of all-order
    smoothness or of ``CandidateFromLimits``' locally uniform endpoint limits.

    For ``t<=0`` the residual is returned exactly as zero without sampling the
    unresolved input fields.  This is justified structurally by the negative
    zero branch together with the time-switch zero germ around ``t=0``.  For
    ``0<t<endpoint`` the diagnostic respects the open singular endpoint when
    choosing its temporal stencil.
    """
    x = _finite(x, "x")
    y = _finite(y, "y")
    z = _finite(z, "z")
    t = _finite(t, "t")
    endpoint = _finite(endpoint, "endpoint")
    if endpoint <= 0.0:
        raise ValueError("endpoint must be positive")
    if t <= 0.0:
        return np.zeros(3)
    if t >= endpoint:
        raise ValueError("past residual is defined only for t<endpoint")

    return navier_stokes_residual_numeric(
        past_velocity(u),
        past_pressure(p),
        x,
        y,
        z,
        t,
        viscosity=viscosity,
        eps_space=eps_space,
        eps_time=eps_time,
        time_domain=(None, endpoint),
    )
