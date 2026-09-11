"""Finite-interval numerical/Volterra adapter for Section 7 primary amplitudes.

The pinned ``PrimaryODE.solution`` is a finite-interval Volterra solution of
the time-dependent two-mode linear ODE constructed from the moving frame.
This module supplies a parameterized numerical path and an independent
integral-equation diagnostic around the already-landed pointwise algebra.

It is deliberately **formal-structure** infrastructure.  The datum/forcing
callbacks remain upstream obligations, SciPy integration is numerical rather
than a proof of the Lean existence theorem, and no returned trajectory is
paper-exact unless separate provenance certifies all of its inputs.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable, Iterable

import numpy as np
from scipy.integrate import quad_vec, solve_ivp

from .primary_amplitude_ode import PrimaryModalDatum, modal_rhs


DatumPath = Callable[[float], PrimaryModalDatum]
ForcingPath = Callable[[float], Iterable[float]]


def _finite(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


def _state2(value: Iterable[float], name: str) -> np.ndarray:
    out = np.asarray(tuple(value), dtype=float)
    if out.shape != (2,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite two-vector")
    return out


def _forcing3(value: Iterable[float]) -> tuple[float, float, float]:
    out = np.asarray(tuple(value), dtype=float)
    if out.shape != (3,) or not np.all(np.isfinite(out)):
        raise ValueError("forcing callback must return three finite scalars")
    return float(out[0]), float(out[1]), float(out[2])


def _harmonic(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError("harmonic must be an integer")
    return value


@dataclass(frozen=True)
class PrimaryAmplitudeInterval:
    """Closed interval and harmonic for a pinned ``PrimaryODE.solution`` slice."""

    start: float
    end: float
    harmonic: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "start", _finite(self.start, "start"))
        object.__setattr__(self, "end", _finite(self.end, "end"))
        object.__setattr__(self, "harmonic", _harmonic(self.harmonic))
        if self.end < self.start:
            raise ValueError("end must satisfy end >= start")

    def require_time(self, time: float) -> float:
        t = _finite(time, "time")
        if t < self.start or t > self.end:
            raise ValueError("time must lie in the closed amplitude interval")
        return t


@dataclass(frozen=True)
class VolterraDiagnostic:
    """Numerical defect for ``z(t)=z0+integral(Az+g)``.

    ``quadrature_error`` is SciPy's error estimate; it is not an interval or a
    rigorous theorem certificate.
    """

    time: float
    defect: np.ndarray
    defect_norm: float
    quadrature_error: float

    def __post_init__(self) -> None:
        defect = _state2(self.defect, "defect")
        object.__setattr__(self, "defect", defect)
        object.__setattr__(self, "defect_norm", _finite(self.defect_norm, "defect_norm"))
        object.__setattr__(
            self, "quadrature_error", _finite(self.quadrature_error, "quadrature_error")
        )
        if self.defect_norm < 0.0 or self.quadrature_error < 0.0:
            raise ValueError("diagnostic norms/errors must be nonnegative")


class NumericalPrimaryAmplitudePath:
    """Dense numerical trajectory returned by ``solve_primary_amplitude_ivp``."""

    def __init__(
        self,
        interval: PrimaryAmplitudeInterval,
        initial_state: Iterable[float],
        dense_solution: Callable[[float], np.ndarray],
    ) -> None:
        self.interval = interval
        self.initial_state = _state2(initial_state, "initial_state").copy()
        self._dense_solution = dense_solution

    def __call__(self, time: float) -> np.ndarray:
        t = self.interval.require_time(time)
        value = np.asarray(self._dense_solution(t), dtype=float).reshape(-1)
        if value.shape != (2,) or not np.all(np.isfinite(value)):
            raise RuntimeError("numerical amplitude path returned a non-finite two-vector")
        return value.copy()


def primary_modal_rhs_at(
    *,
    time: float,
    state: Iterable[float],
    interval: PrimaryAmplitudeInterval,
    datum_at: DatumPath,
    forcing_at: ForcingPath,
) -> np.ndarray:
    """Evaluate the official modal ODE using time-dependent supplied inputs."""
    t = interval.require_time(time)
    data = datum_at(t)
    if not isinstance(data, PrimaryModalDatum):
        raise TypeError("datum_at must return PrimaryModalDatum")
    f_radial, f_K, f_N = _forcing3(forcing_at(t))
    return modal_rhs(
        data,
        interval.harmonic,
        state,
        f_radial=f_radial,
        f_K=f_K,
        f_N=f_N,
    )


def solve_primary_amplitude_ivp(
    *,
    interval: PrimaryAmplitudeInterval,
    initial_state: Iterable[float],
    datum_at: DatumPath,
    forcing_at: ForcingPath,
    rtol: float = 1e-10,
    atol: float = 1e-12,
) -> NumericalPrimaryAmplitudePath:
    """Numerically solve the finite-interval modal IVP with dense output.

    This is execution plumbing for the theorem-shaped ODE, not the
    noncomputable Lean existence proof and not a paper-exact pulse certificate.
    """
    z0 = _state2(initial_state, "initial_state")
    rtol = _finite(rtol, "rtol")
    atol = _finite(atol, "atol")
    if rtol <= 0.0 or atol <= 0.0:
        raise ValueError("rtol and atol must be positive")

    # Fail closed at the initial point before entering SciPy.
    primary_modal_rhs_at(
        time=interval.start,
        state=z0,
        interval=interval,
        datum_at=datum_at,
        forcing_at=forcing_at,
    )

    if interval.start == interval.end:
        return NumericalPrimaryAmplitudePath(interval, z0, lambda _t: z0.copy())

    def rhs(time: float, state: np.ndarray) -> np.ndarray:
        return primary_modal_rhs_at(
            time=time,
            state=state,
            interval=interval,
            datum_at=datum_at,
            forcing_at=forcing_at,
        )

    result = solve_ivp(
        rhs,
        (interval.start, interval.end),
        z0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=True,
    )
    if not result.success or result.sol is None:
        raise RuntimeError(f"primary amplitude IVP solve failed: {result.message}")
    end_state = np.asarray(result.y[:, -1], dtype=float)
    if end_state.shape != (2,) or not np.all(np.isfinite(end_state)):
        raise RuntimeError("primary amplitude IVP solve produced non-finite state")
    return NumericalPrimaryAmplitudePath(interval, z0, result.sol)


def volterra_diagnostic(
    *,
    interval: PrimaryAmplitudeInterval,
    initial_state: Iterable[float],
    state_at: Callable[[float], Iterable[float]],
    datum_at: DatumPath,
    forcing_at: ForcingPath,
    time: float,
    epsabs: float = 1e-10,
    epsrel: float = 1e-10,
) -> VolterraDiagnostic:
    """Check the finite-interval Volterra equation by independent quadrature.

    The state trajectory is treated as a black box.  The checker integrates
    the ODE RHS independently and compares against ``state_at(time)-z0``.
    """
    z0 = _state2(initial_state, "initial_state")
    t = interval.require_time(time)
    epsabs = _finite(epsabs, "epsabs")
    epsrel = _finite(epsrel, "epsrel")
    if epsabs <= 0.0 or epsrel <= 0.0:
        raise ValueError("epsabs and epsrel must be positive")

    start_value = _state2(state_at(interval.start), "state_at(start)")
    if not np.array_equal(start_value, z0):
        raise ValueError("state_at(start) must equal initial_state exactly")

    target = _state2(state_at(t), "state_at(time)")
    if t == interval.start:
        defect = target - z0
        return VolterraDiagnostic(
            time=t,
            defect=defect,
            defect_norm=float(np.linalg.norm(defect, ord=2)),
            quadrature_error=0.0,
        )

    def integrand(s: float) -> np.ndarray:
        state = _state2(state_at(float(s)), "state_at(quadrature point)")
        return primary_modal_rhs_at(
            time=float(s),
            state=state,
            interval=interval,
            datum_at=datum_at,
            forcing_at=forcing_at,
        )

    integral, error = quad_vec(
        integrand,
        interval.start,
        t,
        epsabs=epsabs,
        epsrel=epsrel,
    )
    integral = _state2(integral, "quadrature integral")
    defect = target - z0 - integral
    return VolterraDiagnostic(
        time=t,
        defect=defect,
        defect_norm=float(np.linalg.norm(defect, ord=2)),
        quadrature_error=float(error),
    )
