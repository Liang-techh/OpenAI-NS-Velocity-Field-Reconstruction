"""Theorem-shaped primary amplitude ODE algebra for Section 7.

This module transcribes the finite-dimensional coefficient and forcing maps in
the pinned official Lean files ``NavierStokes/MovingFrameODE.lean``,
``NavierStokes/GrowingMode.lean`` and ``NavierStokes/PrimaryODE.lean`` at
``openai/NavierStokesAndEuler@f9e8bc5``.

The physical tangent coordinates satisfy

``x' = (a-d)x + b y + forceX`` and ``y' = c x - d y + forceY``.

With ``x=p+q`` and ``y=h(p-q)``, where ``h'/h=rate``, the official moving
basis rewrites this as the two-mode system

``z' = A_modal z + force_modal``.

Only this algebraic adapter is implemented here.  It does not construct the
paper's Proposition 5.5 base field, certify the phase/frame hypotheses, solve
the finite-interval Volterra equation, construct the pulse amplitudes, or
promote caller-supplied coefficients to a paper-exact oscillatory wave.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np


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


def _harmonic(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError("harmonic must be an integer")
    return value


@dataclass(frozen=True)
class PrimaryModalDatum:
    """Scalar projections needed by the pinned ``PrimaryODE.FrameData`` map.

    ``g_K`` and ``g_N`` are the shear projections onto the moving tangent
    frame, and ``N_theta`` is the theta component of the second frame vector.
    The actual frame/base-field construction remains an upstream obligation.
    """

    rho: float
    rho_dot: float
    g_K: float
    F: float
    N_theta: float
    g_N: float
    rotation: float
    eigenvalue: float
    eigenvector: float
    eigen_rate: float
    viscosity: float

    def __post_init__(self) -> None:
        names = (
            "rho",
            "rho_dot",
            "g_K",
            "F",
            "N_theta",
            "g_N",
            "rotation",
            "eigenvalue",
            "eigenvector",
            "eigen_rate",
            "viscosity",
        )
        for name in names:
            object.__setattr__(self, name, _finite(getattr(self, name), name))
        if self.eigenvector == 0.0:
            raise ValueError("eigenvector must be nonzero for the moving eigenbasis")


@dataclass(frozen=True)
class MovingCoefficients:
    """The ``a,b,c`` coefficients before scalar viscous damping."""

    a: float
    b: float
    c: float


@dataclass(frozen=True)
class ModalErrors:
    """The four moving-eigenbasis error entries used by ``modalOperator``."""

    e11: float
    e12: float
    e21: float
    e22: float


def moving_coefficients(data: PrimaryModalDatum) -> MovingCoefficients:
    """Pinned ``coeff11/coeff12/coeff21`` after taking frame projections."""
    denom = 1.0 + data.rho * data.rho
    a = data.rho * (data.g_K - data.rho_dot) / denom
    b = (2.0 * data.F * data.N_theta - data.rho * data.rotation) / denom
    c = -(2.0 * data.F * data.N_theta + data.g_N) + data.rho * data.rotation
    return MovingCoefficients(a=a, b=b, c=c)


def modal_errors(data: PrimaryModalDatum) -> ModalErrors:
    """Pinned ``modal11``--``modal22`` after subtracting the reference mode.

    The official ``PrimaryODE`` first forms
    ``B=b-lambda/h`` and ``C=c-lambda*h`` and then applies the moving-basis
    formulas, including the logarithmic eigenvector derivative ``rate``.
    """
    coeff = moving_coefficients(data)
    h = data.eigenvector
    B = coeff.b - data.eigenvalue / h
    C = coeff.c - data.eigenvalue * h
    rate = data.eigen_rate
    e11 = (coeff.a + h * B + C / h - rate) / 2.0
    e12 = (coeff.a - h * B + C / h + rate) / 2.0
    e21 = (coeff.a + h * B - C / h + rate) / 2.0
    e22 = (coeff.a - h * B - C / h - rate) / 2.0
    return ModalErrors(e11=e11, e12=e12, e21=e21, e22=e22)


def harmonic_damping(data: PrimaryModalDatum, harmonic: int) -> float:
    """Pinned ``damping j = j^2 * viscosity``."""
    j = _harmonic(harmonic)
    return float(j * j) * data.viscosity


def modal_operator(data: PrimaryModalDatum, harmonic: int) -> np.ndarray:
    """Return the exact 2x2 matrix shape of ``GrowingMode.modalOperator``."""
    d = harmonic_damping(data, harmonic)
    e = modal_errors(data)
    lam = data.eigenvalue
    return np.array(
        [
            [lam - d + e.e11, e.e12],
            [e.e21, -lam - d + e.e22],
        ],
        dtype=float,
    )


def physical_forcing_components(
    *, rho: float, f_radial: float, f_K: float, f_N: float
) -> np.ndarray:
    """Pinned ``forceX`` and ``forceY`` from projected physical forcing data.

    ``f_K`` and ``f_N`` denote the independently supplied frame projections
    ``<K,tail(f)>`` and ``<N,tail(f)>``.  This function does not infer them from
    a surrogate ambient field.
    """
    rho = _finite(rho, "rho")
    f_radial = _finite(f_radial, "f_radial")
    f_K = _finite(f_K, "f_K")
    f_N = _finite(f_N, "f_N")
    force_x = -(f_radial - rho * f_K) / (1.0 + rho * rho)
    force_y = -f_N
    return np.array([force_x, force_y], dtype=float)


def modal_forcing(
    data: PrimaryModalDatum, *, f_radial: float, f_K: float, f_N: float
) -> np.ndarray:
    """Pinned ``PrimaryODE.FrameData.forcing`` in the moving eigenbasis."""
    force_x, force_y = physical_forcing_components(
        rho=data.rho, f_radial=f_radial, f_K=f_K, f_N=f_N
    )
    h = data.eigenvector
    return np.array(
        [(force_x + force_y / h) / 2.0, (force_x - force_y / h) / 2.0],
        dtype=float,
    )


def modal_rhs(
    data: PrimaryModalDatum,
    harmonic: int,
    state: Iterable[float],
    *,
    f_radial: float = 0.0,
    f_K: float = 0.0,
    f_N: float = 0.0,
) -> np.ndarray:
    """Evaluate the pinned two-mode ODE right-hand side at one point."""
    z = _state2(state, "state")
    return modal_operator(data, harmonic) @ z + modal_forcing(
        data, f_radial=f_radial, f_K=f_K, f_N=f_N
    )


def modal_to_physical(data: PrimaryModalDatum, state: Iterable[float]) -> np.ndarray:
    """Return ``(x,y)=(p+q, h(p-q))`` from the official moving eigenbasis."""
    p, q = _state2(state, "state")
    return np.array([p + q, data.eigenvector * (p - q)], dtype=float)


def physical_to_modal(data: PrimaryModalDatum, physical_state: Iterable[float]) -> np.ndarray:
    """Inverse map ``p=(x+y/h)/2``, ``q=(x-y/h)/2``."""
    x, y = _state2(physical_state, "physical_state")
    h = data.eigenvector
    return np.array([(x + y / h) / 2.0, (x - y / h) / 2.0], dtype=float)


def physical_rhs(
    data: PrimaryModalDatum,
    harmonic: int,
    physical_state: Iterable[float],
    *,
    f_radial: float = 0.0,
    f_K: float = 0.0,
    f_N: float = 0.0,
) -> np.ndarray:
    """Evaluate the independent ``rhsX/rhsY`` moving-frame equations."""
    x, y = _state2(physical_state, "physical_state")
    coeff = moving_coefficients(data)
    damping = harmonic_damping(data, harmonic)
    force_x, force_y = physical_forcing_components(
        rho=data.rho, f_radial=f_radial, f_K=f_K, f_N=f_N
    )
    return np.array(
        [
            (coeff.a - damping) * x + coeff.b * y + force_x,
            coeff.c * x - damping * y + force_y,
        ],
        dtype=float,
    )


def modal_derivative_to_physical(
    data: PrimaryModalDatum,
    state: Iterable[float],
    state_derivative: Iterable[float],
) -> np.ndarray:
    """Differentiate ``x=p+q, y=h(p-q)`` using ``h'=rate*h``.

    This is useful for independently checking that a modal RHS really
    reconstructs the pinned physical moving-frame equations.
    """
    p, q = _state2(state, "state")
    p_dot, q_dot = _state2(state_derivative, "state_derivative")
    h = data.eigenvector
    h_dot = data.eigen_rate * h
    return np.array(
        [p_dot + q_dot, h_dot * (p - q) + h * (p_dot - q_dot)],
        dtype=float,
    )
