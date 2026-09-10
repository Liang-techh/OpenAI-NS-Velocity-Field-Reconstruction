"""Local and localized velocity representations from Eqs. (9.21) and (10.4).

Paper structure:

    A     = A0 + sum_j chi(a_j q) A_j
    B eθ  = B0 eθ + sum_j chi(a_j q) B_j
    u_loc = curl(A) + B eθ                                  (9.21)

and after spatial/time localization

    u = curl(c A) + c B eθ.                                 (10.4)

This file provides an executable composition layer.  Individual correction potentials A_j
and direct azimuthal increments B_j are supplied as callables; constructing the exact paper
sequence is a separate stage rather than being silently approximated.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence
import math
import numpy as np

VectorField = Callable[[float, float, float, float], np.ndarray]
ScalarField = Callable[[float, float, float, float], float]


def curl_numeric(A: VectorField, x: float, y: float, z: float, t: float, *, eps: float = 1e-5) -> np.ndarray:
    """Centered-difference curl of a Cartesian vector potential."""
    def d(axis: int) -> np.ndarray:
        p = [x, y, z]
        m = [x, y, z]
        p[axis] += eps
        m[axis] -= eps
        return (np.asarray(A(*p, t), float) - np.asarray(A(*m, t), float)) / (2.0 * eps)

    d_dx = d(0)
    d_dy = d(1)
    d_dz = d(2)
    return np.array([
        d_dy[2] - d_dz[1],
        d_dz[0] - d_dx[2],
        d_dx[1] - d_dy[0],
    ])


def azimuthal_vector(B: float, x: float, y: float) -> np.ndarray:
    r = math.hypot(x, y)
    if r == 0.0:
        return np.zeros(3)
    return float(B) * np.array([-y / r, x / r, 0.0])


@dataclass(frozen=True)
class LocalField:
    """Executable form of u_loc = curl(A) + B e_theta."""

    vector_potential: VectorField
    azimuthal_scalar: ScalarField

    def velocity(self, x: float, y: float, z: float, t: float, *, eps: float = 1e-5) -> np.ndarray:
        return curl_numeric(self.vector_potential, x, y, z, t, eps=eps) + azimuthal_vector(
            self.azimuthal_scalar(x, y, z, t), x, y
        )


@dataclass(frozen=True)
class LocalizedField:
    """Executable form of Eq. (10.4): u = curl(c A) + c B e_theta."""

    local: LocalField
    cutoff: ScalarField

    def localized_potential(self, x: float, y: float, z: float, t: float) -> np.ndarray:
        return float(self.cutoff(x, y, z, t)) * np.asarray(
            self.local.vector_potential(x, y, z, t), dtype=float
        )

    def velocity(self, x: float, y: float, z: float, t: float, *, eps: float = 1e-5) -> np.ndarray:
        c = float(self.cutoff(x, y, z, t))
        return curl_numeric(self.localized_potential, x, y, z, t, eps=eps) + c * azimuthal_vector(
            self.local.azimuthal_scalar(x, y, z, t), x, y
        )


def sum_vector_fields(fields: Sequence[VectorField]) -> VectorField:
    def total(x: float, y: float, z: float, t: float) -> np.ndarray:
        out = np.zeros(3)
        for f in fields:
            out += np.asarray(f(x, y, z, t), dtype=float)
        return out
    return total


def sum_scalar_fields(fields: Sequence[ScalarField]) -> ScalarField:
    def total(x: float, y: float, z: float, t: float) -> float:
        return float(sum(float(f(x, y, z, t)) for f in fields))
    return total
