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
from .coordinates import _finite
from .verify import finite_vector, jacobian_numeric

VectorField = Callable[[float, float, float, float], np.ndarray]
ScalarField = Callable[[float, float, float, float], float]


def curl_numeric(A: VectorField, x: float, y: float, z: float, t: float, *, eps: float = 1e-5) -> np.ndarray:
    """Diagnostic curl; prefer an analytic callback for nested residuals."""
    J = jacobian_numeric(A,x,y,z,t,eps=eps)
    return np.array([J[2,1]-J[1,2],J[0,2]-J[2,0],J[1,0]-J[0,1]])


def azimuthal_vector(B: float, x: float, y: float) -> np.ndarray:
    B, x, y = _finite(B,"B"), _finite(x,"x"), _finite(y,"y")
    r = math.hypot(x, y)
    if r == 0.0:
        if B != 0:
            raise ValueError("smooth azimuthal velocity must vanish on the axis")
        return np.zeros(3)
    return float(B) * np.array([-y / r, x / r, 0.0])


@dataclass(frozen=True)
class LocalField:
    """u_loc=curl(A)+B*e_theta, with an optional analytic curl.

    Divergence-freeness requires B to be axisymmetric and axis-regular.
    Arbitrary Cartesian callbacks are NOT automatically incompressible.
    """

    vector_potential: VectorField
    azimuthal_scalar: ScalarField
    analytic_curl: VectorField | None = None

    def velocity(self, x: float, y: float, z: float, t: float, *, eps: float = 1e-5) -> np.ndarray:
        poloidal = (finite_vector(self.analytic_curl(x,y,z,t)) if self.analytic_curl is not None
                    else curl_numeric(self.vector_potential,x,y,z,t,eps=eps))
        return poloidal + azimuthal_vector(
            self.azimuthal_scalar(x, y, z, t), x, y
        )


@dataclass(frozen=True)
class LocalizedField:
    """Eq. (10.4), retaining grad(c) cross A.

    For the direct swirl to remain divergence-free, c*B must be
    independent of theta. Prefer an axisymmetric cutoff, not an arbitrary one.
    """

    local: LocalField
    cutoff: ScalarField
    cutoff_gradient: VectorField | None = None

    def localized_potential(self, x: float, y: float, z: float, t: float) -> np.ndarray:
        return _finite(self.cutoff(x,y,z,t),"cutoff") * finite_vector(self.local.vector_potential(x,y,z,t))

    def velocity(self, x: float, y: float, z: float, t: float, *, eps: float = 1e-5) -> np.ndarray:
        c = _finite(self.cutoff(x,y,z,t),"cutoff")
        if self.local.analytic_curl is not None and self.cutoff_gradient is not None:
            poloidal = c*finite_vector(self.local.analytic_curl(x,y,z,t)) + np.cross(
                finite_vector(self.cutoff_gradient(x,y,z,t)),
                finite_vector(self.local.vector_potential(x,y,z,t)))
        else:
            poloidal = curl_numeric(self.localized_potential,x,y,z,t,eps=eps)
        return poloidal + c * azimuthal_vector(
            self.local.azimuthal_scalar(x, y, z, t), x, y
        )


def sum_vector_fields(fields: Sequence[VectorField]) -> VectorField:
    def total(x: float, y: float, z: float, t: float) -> np.ndarray:
        out = np.zeros(3)
        for f in fields:
            out += finite_vector(f(x,y,z,t))
        return out
    return total


def sum_scalar_fields(fields: Sequence[ScalarField]) -> ScalarField:
    def total(x: float, y: float, z: float, t: float) -> float:
        return float(sum(_finite(f(x,y,z,t),"scalar field") for f in fields))
    return total
