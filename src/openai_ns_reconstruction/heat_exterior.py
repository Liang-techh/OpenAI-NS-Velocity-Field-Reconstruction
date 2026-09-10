"""Executable exterior swirl from Appendix A.6, Eqs. (A.32)-(A.38).

This is a parameterized paper formula evaluated by adaptive quadrature, NOT
the complete counterexample: the axis profile, moment matching, stress cone,
recursive coefficients and oscillatory corrections are not supplied here.
In particular this exterior-only field is not regular at r=0.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from scipy.integrate import quad
from .coordinates import validate_h


def _quad_checked(fn, a, b, *, rtol=2e-11, atol=1e-12):
    result = quad(fn, a, b, epsabs=atol, epsrel=rtol, limit=250, full_output=1)
    value, error = result[:2]
    if len(result) > 3 or not math.isfinite(value) or not math.isfinite(error):
        raise ArithmeticError("heat quadrature failed to converge: " + str(result[3:]))
    if error > 10 * max(atol, rtol * abs(value)):
        raise ArithmeticError("heat quadrature error estimate exceeds tolerance")
    return float(value)


def heat_factor(Z: float, h: float = 0.005, *, derivative: int = 0,
                rtol: float = 2e-11, atol: float = 1e-12) -> float:
    """H^(m)(Z) by (A.34); exact endpoint formula (A.35) at Z=0.

    Quadrature error estimates are numerical diagnostics, not interval bounds.
    Orders 0..12 are supported; higher-order certification is not claimed.
    """
    Z, h = float(Z), validate_h(h)
    if not math.isfinite(Z) or Z < 0:
        raise ValueError("heat factor is defined only for finite Z>=0")
    if (isinstance(derivative, bool) or not isinstance(derivative, int)
            or not 0 <= derivative <= 12):
        raise ValueError("derivative order must be an integer in [0,12]")
    if not (math.isfinite(rtol) and 1e-13 <= rtol < 1
            and math.isfinite(atol) and atol > 0):
        raise ValueError("invalid quadrature tolerance")
    m = derivative
    rising = math.prod(h + k for k in range(m))
    if Z == 0:
        return (-1)**m * rising * math.prod(1 + h + k for k in range(m))
    coefficient = (-1)**m * rising / math.gamma(1 + h)
    log_Z = math.log(Z)
    power = h + m

    def integrand(v):
        if v == 0:
            return 0.0
        lv = math.log(v)
        w = log_Z + lv
        # log(1+Zv), without overflowing the product Zv.
        log_den = w + math.log1p(math.exp(-w)) if w > 0 else math.log1p(math.exp(w))
        return math.exp(-v + power * (lv - log_den))

    integral = _quad_checked(integrand, 0.0, np.inf, rtol=rtol,
                             atol=atol / max(1.0, abs(coefficient)))
    return coefficient * integral


@dataclass(frozen=True)
class HeatExterior:
    h: float = 0.005
    c_inf: float = 1.0

    def __post_init__(self):
        validate_h(self.h)
        if not math.isfinite(self.c_inf) or self.c_inf <= 0:
            raise ValueError("c_inf must be finite and positive")

    @property
    def A(self):
        return 0.5 + self.h

    def profile_E(self, X: float, eta: float) -> float:
        """E_heat=c_inf X^(-A) H(2(1-eta^2)/X), including eta=+-1."""
        if not math.isfinite(X) or X <= 0 or not math.isfinite(eta) or abs(eta) > 1:
            raise ValueError("exterior profile requires finite X>0 and |eta|<=1")
        return self.c_inf * X**(-self.A) * heat_factor(2 * (1 - eta**2) / X, self.h)

    def _scales(self, r, tau):
        r, tau = float(r), float(tau)
        if not math.isfinite(r) or r <= 0 or not math.isfinite(tau) or tau < 0:
            raise ValueError("exterior requires finite r>0 and tau>=0")
        s = r * r / 2
        if not math.isfinite(s) or s <= 0:
            raise ArithmeticError("exterior radial scale is not representable")
        return s, 2 * (tau / s)

    def swirl_from_tau(self, r: float, tau: float) -> float:
        s, Z = self._scales(r, tau)
        return self.c_inf * s**(-self.A) * heat_factor(Z, self.h)

    def swirl(self, r: float, t: float) -> float:
        return self.swirl_from_tau(r, 1 - float(t))

    def pressure_from_tau(self, r: float, tau: float) -> float:
        """p=-integral_r^infinity K(rho,t)^2/rho d(rho), zero at infinity.

        w=(r/rho)^2 maps the integral to [0,1]. This is the exterior
        centrifugal pressure only, not the global pressure construction.
        """
        s, Z = self._scales(r, tau)
        integral = _quad_checked(
            lambda w: w**(2 * self.A - 1) * heat_factor(Z * w, self.h)**2,
            0.0, 1.0, rtol=1e-9, atol=1e-11)
        return -0.5 * self.c_inf**2 * s**(-2 * self.A) * integral

    def velocity(self, x: float, y: float, z: float, t: float) -> np.ndarray:
        r = math.hypot(x, y)
        K = self.swirl(r, t)
        return K * np.array([-y / r, x / r, 0.0])

    def pressure(self, x: float, y: float, z: float, t: float) -> float:
        return self.pressure_from_tau(math.hypot(x, y), 1 - float(t))

    def ode_defect(self, Z: float) -> float:
        """Numerical check of (A.37), using independently integrated H,H',H''."""
        H, H1, H2 = (heat_factor(Z, self.h, derivative=m) for m in range(3))
        a = 1 + self.h
        return Z**2 * H2 + (1 + 2 * a * Z) * H1 + a * (a - 1) * H
