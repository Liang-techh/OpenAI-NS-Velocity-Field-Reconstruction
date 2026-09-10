"""Exact pressure row of the Section 5 coefficient recursion.

Paper provenance
----------------
For the formal expansion (5.1), OpenAI's paper writes the order-n pressure
identity, for n >= 1,

    d_X Pi_n = C^-2 sum_{i+j=n} phi_i phi_j - Omega_{n-1}/(2X).   (5.5)

Immediately after (5.6) the paper separates this into the term linear in the
current unknown plus already-known lower-order data,

    d_X Pi_n = 2 C^-2 phi_0 phi_n
             + C^-2 sum_{i=1}^{n-1} phi_i phi_{n-i}
             - (Omega_{n-1}/X)/2.

Lemma 5.1 uses Pi_n(0, eta)=0 on the fixed inner interval.  The paper also
proves after (5.6) that Omega_k is divisible by X.  This module therefore asks
for the *regular quotient* Omega_{n-1}/X rather than numerically dividing by X
at the axis.

This is one exact row of the coupled recursion, not the missing solver for
(5.2)--(5.6): phi_n and the regular Omega quotient are inputs until the angular,
axial and radial rows are constructed.  Numerical quadrature used to integrate
Pi_n is an executable approximation, not a paper-exact certificate.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence
import math

from .quadrature import integrate

ScalarProfile = Callable[[float, float], float]


def _point(X: float, eta: float) -> tuple[float, float]:
    X, eta = float(X), float(eta)
    if not math.isfinite(X) or X < 0:
        raise ValueError("X must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1:
        raise ValueError("eta must be finite with |eta| <= 1")
    return X, eta


def _positive_order(n: int) -> int:
    if isinstance(n, bool) or not isinstance(n, int) or n < 1:
        raise ValueError("recurrence order n must be a positive integer")
    return n


def _finite_value(fn: ScalarProfile, X: float, eta: float, name: str) -> float:
    value = float(fn(X, eta))
    if not math.isfinite(value):
        raise ValueError(f"{name} must return finite values")
    return value


@dataclass(frozen=True)
class PressureRowTerms:
    """The three terms in the post-(5.6) linearized form of Eq. (5.5)."""

    current_linear: float
    lower_order_convolution: float
    radial_source: float

    @property
    def total(self) -> float:
        return self.current_linear + self.lower_order_convolution + self.radial_source


def pressure_row_terms_eq_5_5(
    n: int,
    X: float,
    eta: float,
    phi: Sequence[ScalarProfile],
    omega_prev_over_x: ScalarProfile,
    *,
    C: float,
) -> PressureRowTerms:
    """Evaluate the exact source split of the order-n pressure row (5.5).

    ``phi`` must contain ``phi_0, ..., phi_n``.  ``omega_prev_over_x`` is the
    regular extension of ``Omega_{n-1}/X`` whose existence is established in
    the paragraph following Eq. (5.6).  Passing the quotient directly keeps the
    axis X=0 regular and prevents a numerical 0/0 surrogate.
    """

    n = _positive_order(n)
    X, eta = _point(X, eta)
    C = float(C)
    if not math.isfinite(C) or C <= 0:
        raise ValueError("C must be finite and positive")
    if len(phi) < n + 1:
        raise ValueError("phi must contain coefficients phi_0 through phi_n")
    if not callable(omega_prev_over_x):
        raise TypeError("omega_prev_over_x must be callable")

    values = [_finite_value(phi[i], X, eta, f"phi_{i}") for i in range(n + 1)]
    inv_c2 = 1.0 / (C * C)
    current = 2.0 * inv_c2 * values[0] * values[n]
    lower = inv_c2 * sum(values[i] * values[n - i] for i in range(1, n))
    omega_over_x = _finite_value(omega_prev_over_x, X, eta, "Omega_{n-1}/X")
    return PressureRowTerms(current, lower, -0.5 * omega_over_x)


def pressure_derivative_eq_5_5(
    n: int,
    X: float,
    eta: float,
    phi: Sequence[ScalarProfile],
    omega_prev_over_x: ScalarProfile,
    *,
    C: float,
) -> float:
    """Return ``partial_X Pi_n`` from Eq. (5.5)."""

    return pressure_row_terms_eq_5_5(
        n, X, eta, phi, omega_prev_over_x, C=C
    ).total


def solve_pressure_eq_5_5(
    n: int,
    X: float,
    eta: float,
    phi: Sequence[ScalarProfile],
    omega_prev_over_x: ScalarProfile,
    *,
    C: float,
    quadrature_points: int = 32,
) -> float:
    """Integrate Eq. (5.5) with the Lemma 5.1 datum ``Pi_n(0,eta)=0``.

    This is a genuine positive-order pressure solve once ``phi_n`` and the
    already-determined regular source ``Omega_{n-1}/X`` are supplied.  It does
    not solve the coupled angular/axial equations (5.3)--(5.4) that determine
    ``phi_n`` and ``U_n``.
    """

    n = _positive_order(n)
    X, eta = _point(X, eta)
    return integrate(
        lambda x: pressure_derivative_eq_5_5(
            n, x, eta, phi, omega_prev_over_x, C=C
        ),
        0.0,
        X,
        n=quadrature_points,
    )
