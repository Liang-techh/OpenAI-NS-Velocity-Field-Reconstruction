"""Profile interfaces for the leading OpenAI Navier--Stokes vortex.

The paper proves existence of specially constructed smooth profiles E(X,eta), U(X,eta),
and Pi(X,eta).  They are not a single elementary closed-form tuple.  This module therefore
keeps the paper-exact kinematic formulas separate from whatever concrete profile constructor
is plugged in.

Equation references: (4.3), (4.6), (4.7).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import numpy as np

ScalarFn = Callable[[float, float], float]


@dataclass(frozen=True)
class LeadingProfile:
    """A concrete leading profile with the derivatives needed by Eq. (4.7).

    Parameters
    ----------
    E:
        Azimuthal similarity profile E(X, eta).
    U:
        Axial similarity profile U(X, eta).
    dU_deta:
        Partial derivative of U with respect to eta.
    Pi:
        Optional pressure profile. It is not needed to evaluate velocity, but is
        carried so the same object can later feed the residual verifier.
    name:
        Provenance/status label.
    paper_exact:
        True only for a profile built from the paper's complete profile construction.
    """

    E: ScalarFn
    U: ScalarFn
    dU_deta: ScalarFn
    Pi: ScalarFn | None = None
    name: str = "unnamed-profile"
    paper_exact: bool = False

    def radial_average_U(self, X: float, eta: float, *, n: int = 801) -> float:
        """A_X(U) = X^-1 integral_0^X U(x,eta) dx, Eq. (4.6)."""
        X = float(X)
        if X < 0:
            raise ValueError("X must be nonnegative")
        if X == 0.0:
            return float(self.U(0.0, eta))
        n = max(3, int(n))
        if n % 2 == 0:
            n += 1
        xs = np.linspace(0.0, X, n)
        vals = np.array([self.U(float(x), eta) for x in xs], dtype=float)
        return float(np.trapz(vals, xs) / X)

    def radial_average_dU_deta(self, X: float, eta: float, *, n: int = 801) -> float:
        """d_eta A_X(U), evaluated by differentiating under the radial integral."""
        X = float(X)
        if X < 0:
            raise ValueError("X must be nonnegative")
        if X == 0.0:
            return float(self.dU_deta(0.0, eta))
        n = max(3, int(n))
        if n % 2 == 0:
            n += 1
        xs = np.linspace(0.0, X, n)
        vals = np.array([self.dU_deta(float(x), eta) for x in xs], dtype=float)
        return float(np.trapz(vals, xs) / X)

    def V0(self, X: float, eta: float, h: float, *, n: int = 801) -> float:
        """Radial flux V0 = r u_r from the exact incompressibility identity (4.7)."""
        X = float(X)
        eta = float(eta)
        if X == 0.0:
            return 0.0
        D = 0.5 - h
        d = 1.0 - eta * eta
        L = 1.0 - 2.0 * h * eta * eta
        U = float(self.U(X, eta))
        AU = self.radial_average_U(X, eta, n=n)
        dAU = self.radial_average_dU_deta(X, eta, n=n)
        return (X / L) * (2.0 * eta * U - 2.0 * D * eta * AU - d * dAU)

    def pressure_radial_derivative(self, X: float, eta: float) -> float:
        """Pi_X = E^2/(2X), Eq. (4.7), with the smooth-axis limit left to profile data."""
        X = float(X)
        if X <= 0.0:
            raise ValueError("use the smooth profile-specific axis limit at X=0")
        e = float(self.E(X, eta))
        return e * e / (2.0 * X)
