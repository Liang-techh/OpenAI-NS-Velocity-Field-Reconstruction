"""Leading-profile interface, Eqs. (4.3), (4.6), (4.7).

No arbitrary E,U pair establishes the paper's profile construction. paper_exact
is caller metadata, NOT a certification gate. Numerical quadrature is diagnostic.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Callable
import math
import numpy as np
from .coordinates import validate_h

ScalarFn = Callable[[float, float], float]


@lru_cache(maxsize=16)
def _gauss_rule(n: int) -> tuple[np.ndarray, np.ndarray]:
    if isinstance(n, bool) or not isinstance(n, int) or not 2 <= n <= 2048:
        raise ValueError("quadrature order must be an integer in [2, 2048]")
    x, w = np.polynomial.legendre.leggauss(n)
    x, w = (x + 1) / 2, w / 2
    x.setflags(write=False)
    w.setflags(write=False)
    return x, w


def _point(X: float, eta: float) -> tuple[float, float]:
    X, eta = float(X), float(eta)
    if not math.isfinite(X) or X < 0 or not math.isfinite(eta) or abs(eta) > 1:
        raise ValueError("profile domain is finite X>=0, |eta|<=1")
    return X, eta


@dataclass(frozen=True)
class LeadingProfile:
    E: ScalarFn
    U: ScalarFn
    dU_deta: ScalarFn
    Pi: ScalarFn | None = None
    name: str = "unnamed-profile"
    paper_exact: bool = False
    # Optional exact average/regular-axis functions avoid all quadrature work.
    average_U: ScalarFn | None = None
    average_dU_deta: ScalarFn | None = None
    F: ScalarFn | None = None  # E=sqrt(2X) F

    def _average(self, fn: ScalarFn, exact: ScalarFn | None,
                 X: float, eta: float, n: int) -> float:
        X, eta = _point(X, eta)
        if X == 0:
            value = float(fn(0, eta))
        elif exact is not None:
            value = float(exact(X, eta))
        else:
            nodes, weights = _gauss_rule(n)
            values = np.asarray([fn(float(X * s), eta) for s in nodes], dtype=float)
            value = float(weights @ values)
        if not math.isfinite(value):
            raise ArithmeticError("profile average is not finite")
        return value

    def radial_average_U(self, X: float, eta: float, *, n: int = 32) -> float:
        """A_X(U)=integral_0^1 U(Xs,eta) ds, with its exact axis limit."""
        return self._average(self.U, self.average_U, X, eta, n)

    def radial_average_dU_deta(self, X: float, eta: float, *, n: int = 32) -> float:
        return self._average(self.dU_deta, self.average_dU_deta, X, eta, n)

    def V0(self, X: float, eta: float, h: float, *, n: int = 32,
           d: float | None = None) -> float:
        X, eta = _point(X, eta)
        h = validate_h(h)
        if X == 0:
            return 0.0
        d = 1 - eta**2 if d is None else float(d)
        if not math.isfinite(d) or not 0 <= d <= 1 + 1e-12:
            raise ValueError("d must be in [0,1]")
        AU = self.radial_average_U(X, eta, n=n)
        dAU = self.radial_average_dU_deta(X, eta, n=n)
        value = X / (1 - 2 * h * eta**2) * (
            2 * eta * float(self.U(X, eta)) - 2 * (0.5 - h) * eta * AU - d * dAU)
        if not math.isfinite(value):
            raise ArithmeticError("radial flux is not finite")
        return value

    def pressure_radial_derivative(self, X: float, eta: float) -> float:
        X, eta = _point(X, eta)
        if X == 0:
            if self.F is None:
                raise ValueError("axis pressure derivative requires regular profile F")
            return float(self.F(0, eta))**2
        return float(self.E(X, eta))**2 / (2 * X)
