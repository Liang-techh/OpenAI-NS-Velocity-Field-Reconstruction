"""Eq. (5.15) forward reconstruction after the Lemma 5.2 repair.

Section 5.2 first repairs the compact moments by modifying ``U_n`` and ``E_n``.
The remaining positive-order profiles are then *reconstructed* rather than
independently fitted:

    F_n(X,eta)  = integral_0^X U_n(x,eta) dx,
    V_n(X,eta)  = -integral_0^X Z_{-A,n} U_n(x,eta) dx,
    Pi_n(X,eta) = integral_0^X [ C^-2 sum_{i+j=n} phi_i phi_j
                                 - Omega_{n-1}/(2x) ] dx.          (5.15)

Using Eq. (5.2), the ``V_n`` integral can be evaluated without differentiating
or dividing by ``X`` at the axis.  If ``dF_deta = integral_0^X d_eta U_n dx``
then

    V_n = [2 eta X U_n - 2 eta (D+lambda_n) F_n - d dF_deta] / L.

This module implements that exact algebraic bridge for *supplied repaired
profiles*.  It does not manufacture the repaired ``U_n``/``phi_n`` or the
lower-order ``Omega`` source, and finite Gauss-Legendre quadrature is numerical
verification infrastructure rather than a paper-exact certificate.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Callable, Sequence
import math

from .background_recurrence import solve_pressure_eq_5_5
from .coordinates import validate_h
from .quadrature import integrate

ScalarProfile = Callable[[float, float], float]


def _positive_order(n: int) -> int:
    if isinstance(n, bool) or not isinstance(n, Integral) or n < 1:
        raise ValueError("coefficient order n must be a positive integer")
    return int(n)


def _point(X: float, eta: float) -> tuple[float, float]:
    X, eta = float(X), float(eta)
    if not math.isfinite(X) or X < 0.0:
        raise ValueError("X must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return X, eta


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


@dataclass(frozen=True)
class Eq515Reconstruction:
    """The three reconstructed profiles in Eq. (5.15) at one ``(X,eta)``.

    ``dF_deta`` is carried as an audit quantity because Eq. (5.2) uses it in
    the regular axis formula for ``V_n``.
    """

    F: float
    V: float
    Pi: float
    dF_deta: float

    def __post_init__(self) -> None:
        for value, name in (
            (self.F, "F_n"),
            (self.V, "V_n"),
            (self.Pi, "Pi_n"),
            (self.dF_deta, "d_eta F_n"),
        ):
            _finite(value, name)


def reconstruct_eq_5_15(
    n: int,
    X: float,
    eta: float,
    U_n: ScalarProfile,
    dU_n_deta: ScalarProfile,
    phi: Sequence[ScalarProfile],
    omega_prev_over_x: ScalarProfile,
    *,
    h: float,
    C: float,
    quadrature_points: int = 32,
) -> Eq515Reconstruction:
    """Reconstruct ``F_n,V_n,Pi_n`` exactly in the structure of Eq. (5.15).

    ``U_n`` and ``phi[n]`` are expected to be the already-extended/repaired
    profiles from Lemma 5.2.  ``phi`` must contain ``phi_0,...,phi_n`` and
    ``omega_prev_over_x`` is the smooth quotient supplied by Eq. (5.6).

    The pressure path delegates to the already landed Eq. (5.5) implementation;
    the radial-flux path uses Eq. (5.2), algebraically equivalent to the second
    integral in Eq. (5.15), so it remains regular at ``X=0``.
    """

    n = _positive_order(n)
    X, eta = _point(X, eta)
    h = validate_h(h)
    C = _finite(C, "C")
    if C <= 0.0:
        raise ValueError("C must be positive")
    if not callable(U_n) or not callable(dU_n_deta):
        raise TypeError("U_n and dU_n_deta must be callable")
    if not callable(omega_prev_over_x):
        raise TypeError("omega_prev_over_x must be callable")
    if len(phi) < n + 1 or not all(callable(phi[j]) for j in range(n + 1)):
        raise ValueError("phi must contain callable profiles phi_0 through phi_n")

    F = integrate(
        lambda x: _finite(U_n(x, eta), "U_n"),
        0.0,
        X,
        n=quadrature_points,
    )
    dF_deta = integrate(
        lambda x: _finite(dU_n_deta(x, eta), "d_eta U_n"),
        0.0,
        X,
        n=quadrature_points,
    )

    D = 0.5 - h
    lam = 2.0 * n * h
    d = 1.0 - eta * eta
    L = 1.0 - 2.0 * h * eta * eta
    if L <= 0.0:
        raise ValueError("similarity factor L must be positive")
    endpoint_u = _finite(U_n(X, eta), "U_n")
    V = (
        2.0 * eta * X * endpoint_u
        - 2.0 * eta * (D + lam) * F
        - d * dF_deta
    ) / L

    Pi = solve_pressure_eq_5_5(
        n,
        X,
        eta,
        phi,
        omega_prev_over_x,
        C=C,
        quadrature_points=quadrature_points,
    )

    return Eq515Reconstruction(
        F=_finite(F, "F_n"),
        V=_finite(V, "V_n"),
        Pi=_finite(Pi, "Pi_n"),
        dF_deta=_finite(dF_deta, "d_eta F_n"),
    )
