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

The endpoint-transfer certificate below records one exact structural
consequence needed by the t=1 Borel glue: any residual derivative majorant
whose whole pre-endpoint validity window starts in the official unit plateau
can be reused for the activated residual without adding time-cutoff correction
terms.  This is a fail-closed interval/identity transfer only; it does not prove
the source majorant or manufacture endpoint jets.

This module still does not construct the unresolved incoming singular fields
or provide the separate smooth global force extension through t=1.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import math
import numpy as np

from .coordinates import _finite
from .cutoffs import smooth_cutoff, standard_cutoff_derivative
from .endpoint_limit_majorant import EndpointPowerLawMajorant
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
SECTION10_ENDPOINT = 1.0


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


@dataclass(frozen=True)
class EndpointLocalizationTransferCertificate:
    """Conditional transfer of an endpoint majorant through the official switch.

    ``certified`` means only that the supplied majorant's whole validity window
    ``[valid_from, 1)`` lies in the exact unit plateau ``t>=3/4`` of the
    executable Section 10 time switch.  On that interval ``chi=1``, ``chi'=0``
    and ``chi^2-chi=0`` exactly, so the activated residual equals the original
    residual as a function.  Hence, whenever the relevant derivatives exist,
    all spacetime derivative bounds transfer unchanged on the open endpoint
    regime.

    The source ``EndpointPowerLawMajorant`` remains a caller-supplied analytic
    hypothesis.  This certificate does not establish that hypothesis, locally
    uniform convergence, endpoint jet existence, or smooth Borel gluing.
    """

    derivative_degree: int
    spatial_window: int
    valid_from: float
    endpoint: float

    @property
    def certified(self) -> bool:
        values = (self.valid_from, self.endpoint)
        return (
            self.derivative_degree >= 0
            and self.spatial_window >= 0
            and all(math.isfinite(value) for value in values)
            and self.valid_from >= LATE_START
            and self.valid_from < self.endpoint
            and self.endpoint == SECTION10_ENDPOINT
        )

    def correction_coefficients(self, t: float) -> tuple[float, float, float]:
        """Return ``(chi, chi', chi^2-chi)`` inside the certified interval.

        The result is required to be exactly ``(1,0,0)``.  Any future change to
        the executable cutoff that violates the official late plateau therefore
        fails closed instead of silently transferring an endpoint majorant
        across a transition collar.
        """
        if not self.certified:
            raise ValueError("endpoint localization transfer is not certified")
        t = _finite(t, "t")
        if t < self.valid_from or t >= self.endpoint:
            raise ValueError("time must lie in the certified pre-endpoint interval")
        chi = time_switch(t)
        dchi = time_switch_derivative(t)
        nonlinear = chi * chi - chi
        if chi != 1.0 or dchi != 0.0 or nonlinear != 0.0:
            raise ArithmeticError("time localization is not the identity on the claimed endpoint window")
        return chi, dchi, nonlinear


def section10_endpoint_localization_transfer(
    majorant: EndpointPowerLawMajorant,
) -> EndpointLocalizationTransferCertificate:
    """Fail closed unless an endpoint majorant lives in the official late plateau.

    This prevents an arbitrary transition-window majorant from being treated as
    input to the t=1 endpoint glue.  No numerical sampling is used to enlarge
    the admissible interval: the gate is the pinned exact threshold ``3/4`` and
    the pinned Section 10 endpoint ``1``.
    """
    if not isinstance(majorant, EndpointPowerLawMajorant):
        raise ValueError("majorant must be an EndpointPowerLawMajorant")
    if majorant.endpoint != SECTION10_ENDPOINT:
        raise ValueError("Section 10 endpoint localization requires endpoint=1")
    if majorant.valid_from < LATE_START:
        raise ValueError("endpoint majorant validity must start in the t>=3/4 unit plateau")
    certificate = EndpointLocalizationTransferCertificate(
        derivative_degree=majorant.derivative_degree,
        spatial_window=majorant.spatial_window,
        valid_from=majorant.valid_from,
        endpoint=majorant.endpoint,
    )
    if not certificate.certified:
        raise ArithmeticError("endpoint localization transfer invariant failed")
    certificate.correction_coefficients(certificate.valid_from)
    return certificate


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
