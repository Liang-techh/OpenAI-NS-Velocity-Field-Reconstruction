"""Composition layers for (9.21), (10.4), not the missing correction solver.

curl(A) is divergence-free analytically. The extra B e_theta is divergence-free
ONLY if B is axisymmetric; localization requires c*B to stay axisymmetric too.
Arbitrary Cartesian callables do NOT establish either assumption. Factory
methods below express the scalar axisymmetry explicitly. Numerical curl is only
a finite-difference approximation to this structural representation.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Sequence
import math
import numpy as np
from .verify import jacobian_numeric, _vector

VectorField = Callable[[float, float, float, float], np.ndarray]
ScalarField = Callable[[float, float, float, float], float]
AxisymmetricScalar = Callable[[float, float, float], float]


def curl_numeric(A: VectorField, x: float, y: float, z: float, t: float,
                 *, eps: float = 1e-5) -> np.ndarray:
    J = jacobian_numeric(A, x, y, z, t, eps=eps)
    return np.array([J[2, 1] - J[1, 2], J[0, 2] - J[2, 0], J[1, 0] - J[0, 1]])


def azimuthal_vector(B: float, x: float, y: float) -> np.ndarray:
    B, r = float(B), math.hypot(x, y)
    if not math.isfinite(B) or not math.isfinite(r):
        raise ValueError("azimuthal data must be finite")
    if r == 0:
        if B != 0:
            raise ValueError("regular azimuthal scalar must vanish on the axis")
        return np.zeros(3)
    return B * np.array([-y / r, x / r, 0.0])


@dataclass(frozen=True)
class LocalField:
    vector_potential: VectorField
    azimuthal_scalar: ScalarField

    @classmethod
    def from_axisymmetric(cls, vector_potential: VectorField,
                          B: AxisymmetricScalar) -> LocalField:
        return cls(vector_potential, lambda x, y, z, t: B(math.hypot(x, y), z, t))

    def velocity(self, x: float, y: float, z: float, t: float,
                 *, eps: float = 1e-5) -> np.ndarray:
        return curl_numeric(self.vector_potential, x, y, z, t, eps=eps) + azimuthal_vector(
            self.azimuthal_scalar(x, y, z, t), x, y)


@dataclass(frozen=True)
class LocalizedField:
    local: LocalField
    cutoff: ScalarField

    @classmethod
    def from_axisymmetric(cls, local: LocalField,
                          cutoff: AxisymmetricScalar) -> LocalizedField:
        return cls(local, lambda x, y, z, t: cutoff(math.hypot(x, y), z, t))

    def localized_potential(self, x: float, y: float, z: float, t: float) -> np.ndarray:
        c = float(self.cutoff(x, y, z, t))
        if not math.isfinite(c):
            raise ValueError("cutoff must be finite")
        if c == 0:  # do not call an undefined/singular local field outside support
            return np.zeros(3)
        return c * _vector(self.local.vector_potential(x, y, z, t))

    def velocity(self, x: float, y: float, z: float, t: float,
                 *, eps: float = 1e-5) -> np.ndarray:
        c = float(self.cutoff(x, y, z, t))
        if not math.isfinite(c):
            raise ValueError("cutoff must be finite")
        out = curl_numeric(self.localized_potential, x, y, z, t, eps=eps)
        if c != 0:
            out += c * azimuthal_vector(self.local.azimuthal_scalar(x, y, z, t), x, y)
        return out


def sum_vector_fields(fields: Sequence[VectorField]) -> VectorField:
    fields = tuple(fields)
    def total(x, y, z, t):
        out = np.zeros(3)
        for field in fields:
            out += _vector(field(x, y, z, t))
        return out
    return total


def sum_scalar_fields(fields: Sequence[ScalarField]) -> ScalarField:
    fields = tuple(fields)
    def total(x, y, z, t):
        return float(sum(float(field(x, y, z, t)) for field in fields))
    return total
