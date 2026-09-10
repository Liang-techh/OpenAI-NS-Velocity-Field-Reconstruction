"""Natural-axis data and profile rescaling from the official Lean construction.

This module implements the *explicit* part of the Stage-1 construction that is
already materialized in OpenAI's Lean files ``NaturalAxisData.lean``,
``NaturalAxisCoefficients.lean`` and ``NaturalProfile.lean``.

It does **not** solve the upstream coefficient-space fixed-point problem.  Thus
objects assembled here are executable realizations of the exact rescaling
formulas, but are not by themselves the paper's completed Theorem 4.6 profile.

Official Lean source (pinned by ``references/provenance_manifest.json``):

* ``NaturalAxisData``: D, A, d, L, U, H, W, chi.
* ``NaturalAxisCoefficients``: realGradient, realPhase, realAmplitude.
* ``NaturalProfile``: rescalePoint, angularProfile, affineProfile.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp
from typing import Callable

import numpy as np

Scalar1D = Callable[[float], float]
Scalar2D = Callable[[float, float], float]


def D(h: float) -> float:
    return 0.5 - float(h)


def A(h: float) -> float:
    return 0.5 + float(h)


def d(eta: float) -> float:
    eta = float(eta)
    return 1.0 - eta * eta


def L(h: float, eta: float) -> float:
    eta = float(eta)
    return 1.0 - 2.0 * float(h) * eta * eta


def axis_U(j: float, eta: float) -> float:
    return 4.0 * float(eta) + float(j)


def H(h: float, j: float, eta: float) -> float:
    eta = float(eta)
    return D(h) * eta + d(eta) * axis_U(j, eta)


def W(h: float, j: float, eta: float) -> float:
    eta = float(eta)
    return 1.0 - 4.0 * d(eta) - 2.0 * D(h) * eta * axis_U(j, eta)


def chi(h: float, j: float, sigma: float, eta: float) -> float:
    hh = H(h, j, eta)
    s = float(sigma)
    return hh * hh / (hh * hh + s * s)


def real_gradient(h: float, j: float, sigma: float, eta: float) -> float:
    """Lean ``realGradient`` = ``-L H / (H^2 + sigma^2)``."""

    hh = H(h, j, eta)
    s = float(sigma)
    return -L(h, eta) * hh / (hh * hh + s * s)


def real_phase(
    h: float,
    j: float,
    sigma: float,
    eta: float,
    *,
    samples: int = 4001,
) -> float:
    """Numerically evaluate Lean's exact integral definition of ``realPhase``.

    Lean defines ``realPhase(x) = x * integral_0^1 realGradient(t*x) dt``.
    This is equivalently ``integral_0^x realGradient(s) ds``.  The defining
    formula is paper/formalization exact; this routine evaluates its integral
    by trapezoidal quadrature, so the returned floating-point number is a
    numerical approximation rather than an exact certificate.
    """

    eta = float(eta)
    if eta == 0.0:
        return 0.0
    n = max(3, int(samples))
    if n % 2 == 0:
        n += 1
    xs = np.linspace(0.0, eta, n)
    vals = np.array([real_gradient(h, j, sigma, float(x)) for x in xs], dtype=float)
    return float(np.trapezoid(vals, xs))


def real_amplitude(
    h: float,
    j: float,
    sigma: float,
    Lambda: float,
    C: float,
    eta: float,
    *,
    samples: int = 4001,
) -> float:
    """Lean ``realAmplitude = exp(Lambda * realPhase) / C``."""

    C = float(C)
    if C == 0.0:
        raise ValueError("C must be nonzero")
    return exp(float(Lambda) * real_phase(h, j, sigma, eta, samples=samples)) / C


@dataclass(frozen=True)
class NaturalAxisParameters:
    """Numerical parameters constrained to the manuscript's proved small range."""

    h: float
    j: float
    sigma: float
    Lambda: float
    C: float
    phase_samples: int = 4001

    def __post_init__(self) -> None:
        if not (0.0 < self.h <= 1.0 / 1000.0):
            raise ValueError("h must satisfy 0 < h <= 1/1000")
        if not (0.0 < self.j <= 1.0 / 1000.0):
            raise ValueError("j must satisfy 0 < j <= 1/1000")
        if self.sigma <= 0.0:
            raise ValueError("sigma must be positive")
        if self.Lambda <= 0.0:
            raise ValueError("Lambda must be positive")
        if self.C <= 0.0:
            raise ValueError("C must be positive")

    def amplitude(self, eta: float) -> float:
        return real_amplitude(
            self.h,
            self.j,
            self.sigma,
            self.Lambda,
            self.C,
            eta,
            samples=self.phase_samples,
        )


@dataclass(frozen=True)
class NaturalProfileAssembly:
    """Exact Lean rescaling of supplied scaled fixed-point fields.

    The supplied ``phi``, ``u``, ``average`` and ``pressure`` are the scaled
    fields appearing in ``NaturalProfile.ProfileFamily``.  Supplying arbitrary
    functions here does not make a Theorem 4.6 witness; the missing fixed-point
    construction remains tracked as a Stage-1 blocker.
    """

    parameters: NaturalAxisParameters
    phi: Scalar2D
    u: Scalar2D
    du_deta: Scalar2D
    average: Scalar2D
    pressure: Scalar2D
    axis_pressure: Scalar1D

    def rescale_point(self, X: float, eta: float) -> tuple[float, float]:
        return self.parameters.Lambda * float(X), float(eta)

    def E(self, X: float, eta: float) -> float:
        """Angular profile: ``a(eta) * phi(Lambda*X, eta)``."""

        Y, eta = self.rescale_point(X, eta)
        return self.parameters.amplitude(eta) * float(self.phi(Y, eta))

    def U(self, X: float, eta: float) -> float:
        """Axial profile: ``(4 eta + j) + Lambda^-1 u(Lambda X,eta)``."""

        Y, eta = self.rescale_point(X, eta)
        return axis_U(self.parameters.j, eta) + float(self.u(Y, eta)) / self.parameters.Lambda

    def dU_deta(self, X: float, eta: float) -> float:
        """Eta derivative of the natural axial profile."""

        Y, eta = self.rescale_point(X, eta)
        return 4.0 + float(self.du_deta(Y, eta)) / self.parameters.Lambda

    def radial_average(self, X: float, eta: float) -> float:
        """Natural radial-average field reconstructed by the same affine rule."""

        Y, eta = self.rescale_point(X, eta)
        return axis_U(self.parameters.j, eta) + float(self.average(Y, eta)) / self.parameters.Lambda

    def Pi(self, X: float, eta: float) -> float:
        """Pressure profile: ``P0(eta) + Lambda^-1 P(Lambda X,eta)``."""

        Y, eta = self.rescale_point(X, eta)
        return float(self.axis_pressure(eta)) + float(self.pressure(Y, eta)) / self.parameters.Lambda

    def to_leading_profile(self):
        """Adapt to ``LeadingProfile`` without falsely claiming paper exactness."""

        from .profiles import LeadingProfile

        return LeadingProfile(
            E=self.E,
            U=self.U,
            dU_deta=self.dU_deta,
            Pi=self.Pi,
            name="natural-profile-rescaling (fixed-point inputs unresolved)",
            paper_exact=False,
        )
