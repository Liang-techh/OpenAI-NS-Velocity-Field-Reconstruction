"""Leading profile interface; Eqs. (4.3)-(4.7).

F=E/sqrt(2X) is the smooth quantity on the axis. The exact primitives
can be supplied; otherwise cached Gaussian quadrature is used on [0,1].
No concrete profile in this module instantiates Theorem 4.6.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable
import math
import numpy as np
from .coordinates import _finite, _validate_h
from .quadrature import unit_rule

ScalarFn = Callable[[float, float], float]


@dataclass(frozen=True)
class LeadingProfile:
    E: ScalarFn
    U: ScalarFn
    dU_deta: ScalarFn
    Pi: ScalarFn | None = None
    name: str = "unnamed-profile"
    paper_exact: bool = False
    F: ScalarFn | None = None
    average_U: ScalarFn | None = None
    average_dU_deta: ScalarFn | None = None
    provenance: str | None = None

    def __post_init__(self) -> None:
        for name in ("E", "U", "dU_deta"):
            if not callable(getattr(self, name)):
                raise TypeError(f"{name} must be callable")
        for name in ("Pi", "F", "average_U", "average_dU_deta"):
            value = getattr(self, name)
            if value is not None and not callable(value):
                raise TypeError(f"{name} must be callable or None")
        if self.paper_exact and not (self.provenance and self.provenance.strip()):
            raise ValueError("paper_exact requires explicit provenance; a flag is not a proof")

    @staticmethod
    def _point(X: float, eta: float) -> tuple[float, float]:
        X, eta = _finite(X, "X"), _finite(eta, "eta")
        if X < 0 or abs(eta) > 1 + 1e-10:
            raise ValueError("profile domain requires X>=0 and |eta|<=1")
        return X, eta

    def _average(self, fn: ScalarFn, exact: ScalarFn | None,
                 X: float, eta: float, n: int) -> float:
        X, eta = self._point(X, eta)
        nodes, weights = unit_rule(n)
        if X == 0:
            return _finite(fn(0.0, eta), "axis profile value")
        if exact is not None:
            return _finite(exact(X, eta), "radial average")
        values = np.array([fn(float(X * node), eta) for node in nodes], dtype=float)
        if values.shape != nodes.shape or not np.all(np.isfinite(values)):
            raise ValueError("profile must return finite scalars")
        return float(weights @ values)

    def radial_average_U(self, X: float, eta: float, *, n: int = 32) -> float:
        return self._average(self.U, self.average_U, X, eta, n)

    def radial_average_dU_deta(self, X: float, eta: float, *, n: int = 32) -> float:
        return self._average(self.dU_deta, self.average_dU_deta, X, eta, n)

    def radial_flux_factor(self, X: float, eta: float, h: float, *,
                           lam: float = 0.0, n: int = 32,
                           d: float | None = None, L: float | None = None) -> float:
        """V_n/X from Eqs. (4.7), (5.2), (5.27)."""
        X, eta = self._point(X, eta)
        h, lam = _validate_h(h), _finite(lam, "lam")
        d = 1 - eta**2 if d is None else _finite(d, "d")
        L = 1 - 2*h*eta**2 if L is None else _finite(L, "L")
        if d < -1e-10 or L <= 0:
            raise ValueError("invalid physical-chart d or L")
        U = _finite(self.U(X, eta), "U")
        AU = self.radial_average_U(X, eta, n=n)
        dAU = self.radial_average_dU_deta(X, eta, n=n)
        return (2*eta*U - 2*(0.5-h+lam)*eta*AU - d*dAU) / L

    def V0(self, X: float, eta: float, h: float, *, n: int = 32,
           d: float | None = None, L: float | None = None) -> float:
        X, eta = self._point(X, eta)
        return X * self.radial_flux_factor(X, eta, h, n=n, d=d, L=L)

    def smooth_swirl_factor(self, X: float, eta: float) -> float:
        X, eta = self._point(X, eta)
        if self.F is not None:
            return _finite(self.F(X, eta), "F")
        if X == 0:
            raise ValueError("supply smooth F=E/sqrt(2X) to evaluate an axis limit")
        return _finite(self.E(X, eta), "E") / math.sqrt(2*X)

    def pressure_radial_derivative(self, X: float, eta: float) -> float:
        return self.smooth_swirl_factor(X, eta)**2


def _point(X: float, eta: float) -> tuple[float, float]:
    """Compatibility helper retained for Stage-2 code written against main."""
    return LeadingProfile._point(X, eta)


def toy_gaussian_profile() -> LeadingProfile:
    """Explicit analytic TEST fixture, not the paper's constructed leading profile."""
    def avg(X: float) -> float:
        return 1.0 if X == 0 else -math.expm1(-X) / X
    def F(X: float, eta: float) -> float:
        return math.exp(-X) * (1 + 0.1*eta**2)
    return LeadingProfile(
        E=lambda X,e: math.sqrt(2*X)*F(X,e),
        U=lambda X,e: e*math.exp(-X), dU_deta=lambda X,e: math.exp(-X),
        Pi=lambda X,e: -0.5*math.exp(-2*X)*(1+0.1*e**2)**2,
        name="toy-gaussian-not-openai", F=F,
        average_U=lambda X,e: e*avg(X), average_dU_deta=lambda X,e: avg(X),
        provenance="Repository analytic diagnostic fixture, not Theorem 4.6.",
    )
