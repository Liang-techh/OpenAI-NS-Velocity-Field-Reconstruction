"""Exact paper-derived rows of the Section 5 coefficient recursion.

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
proves after (5.6) that every Omega_k is divisible by X because Eq. (5.2) gives
V_j = X v_j with smooth v_j.  This module therefore represents the radial
coefficient by the regular jet of v_j and evaluates Omega_k/X directly, never
forming a numerical 0/0 at the axis.

The Eq. (5.6) implementation is still a pointwise recurrence row, not the
missing coupled solver for (5.2)--(5.6): the smooth v_j/U_j jets must come from
already constructed coefficient profiles.  Numerical quadrature used to
integrate Pi_n is an executable approximation, not a paper-exact certificate.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence
import math

from .coordinates import validate_h
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


def _nonnegative_order(k: int) -> int:
    if isinstance(k, bool) or not isinstance(k, int) or k < 0:
        raise ValueError("coefficient order k must be a nonnegative integer")
    return k


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _finite_value(fn: ScalarProfile, X: float, eta: float, name: str) -> float:
    return _finite(fn(X, eta), name)


@dataclass(frozen=True)
class RegularFluxJet:
    """Second jet of ``v`` in the paper's regular factorization ``V = X v``.

    All entries are evaluated at one fixed ``(X, eta)``.  The second derivatives
    are exactly what is needed to evaluate the shifted axial-viscosity operator
    ``Z^[2]`` in Eq. (5.6) without differentiating a numerically divided ``V/X``.
    """

    value: float
    dX: float
    dEta: float
    dXX: float
    dXdEta: float
    dEtaEta: float

    def checked(self, name: str = "v") -> "RegularFluxJet":
        values = (
            _finite(self.value, name),
            _finite(self.dX, f"dX {name}"),
            _finite(self.dEta, f"dEta {name}"),
            _finite(self.dXX, f"dXX {name}"),
            _finite(self.dXdEta, f"dXdEta {name}"),
            _finite(self.dEtaEta, f"dEtaEta {name}"),
        )
        return RegularFluxJet(*values)


@dataclass(frozen=True)
class OmegaRowTerms:
    """The five regular contributions to ``Omega_k/X`` from Eq. (5.6)."""

    scaling: float
    radial_transport: float
    axial_transport: float
    radial_viscosity: float
    shifted_axial_viscosity: float

    @property
    def total(self) -> float:
        return (
            self.scaling
            + self.radial_transport
            + self.axial_transport
            + self.radial_viscosity
            + self.shifted_axial_viscosity
        )


def _geometry(eta: float, h: float) -> tuple[float, float, float]:
    h = validate_h(h)
    D = 0.5 - h
    d = 1.0 - eta * eta
    L = 1.0 - 2.0 * h * eta * eta
    if L <= 0.0:
        raise ValueError("similarity factor L must be positive")
    return D, d, L


def regular_T_on_Xv_over_X(
    jet: RegularFluxJet,
    X: float,
    eta: float,
    *,
    b: float,
    h: float,
) -> float:
    """Evaluate ``T_b(X v)/X`` using Eq. (4.2), regularly at ``X=0``."""

    X, eta = _point(X, eta)
    jet = jet.checked()
    b = _finite(b, "b")
    D, _, L = _geometry(eta, h)
    return ((1.0 - b) * jet.value + D * eta * jet.dEta + X * jet.dX) / L


def regular_Z_on_Xv_over_X(
    jet: RegularFluxJet,
    X: float,
    eta: float,
    *,
    b: float,
    h: float,
) -> float:
    """Evaluate ``Z_b(X v)/X`` from Eq. (4.2), including the axis value."""

    X, eta = _point(X, eta)
    jet = jet.checked()
    b = _finite(b, "b")
    _, d, L = _geometry(eta, h)
    return (
        2.0 * eta * (b - 1.0) * jet.value
        + d * jet.dEta
        - 2.0 * eta * X * jet.dX
    ) / L


def regular_Z2_on_Xv_over_X(
    jet: RegularFluxJet,
    X: float,
    eta: float,
    *,
    b: float,
    h: float,
) -> float:
    """Evaluate ``Z_{b-D} Z_b (X v) / X`` with no axis division.

    Section 5 abbreviates ``Z^[2]_{a,n}=Z_{a+lambda_n-D} Z_{a+lambda_n}``.
    Here ``b`` is the inner operator index ``a+lambda_n``.  Writing
    ``Z_b(Xv)=Xw`` keeps both applications regular.  The formulas for ``w_X``
    and ``w_eta`` below are analytic derivatives of Eq. (4.2), so only the
    second jet of ``v`` is required.
    """

    X, eta = _point(X, eta)
    jet = jet.checked()
    b = _finite(b, "b")
    h = validate_h(h)
    D, d, L = _geometry(eta, h)

    numerator = (
        2.0 * eta * (b - 1.0) * jet.value
        + d * jet.dEta
        - 2.0 * eta * X * jet.dX
    )
    numerator_X = (
        2.0 * eta * (b - 2.0) * jet.dX
        + d * jet.dXdEta
        - 2.0 * eta * X * jet.dXX
    )
    numerator_eta = (
        2.0 * (b - 1.0) * jet.value
        + 2.0 * eta * (b - 2.0) * jet.dEta
        + d * jet.dEtaEta
        - 2.0 * X * jet.dX
        - 2.0 * eta * X * jet.dXdEta
    )

    w = numerator / L
    w_X = numerator_X / L
    # L_eta = -4 h eta.
    w_eta = numerator_eta / L + 4.0 * h * eta * numerator / (L * L)
    return (
        2.0 * eta * (b - D - 1.0) * w
        + d * w_eta
        - 2.0 * eta * X * w_X
    ) / L


def omega_row_terms_eq_5_6(
    k: int,
    X: float,
    eta: float,
    radial_flux_jets: Sequence[RegularFluxJet],
    U_values: Sequence[float],
    *,
    h: float,
) -> OmegaRowTerms:
    """Evaluate the exact regular quotient ``Omega_k/X`` from Eq. (5.6).

    Eq. (5.2) gives ``V_j = X v_j``.  Dividing Eq. (5.6) by ``X`` then gives
    regular terms

    ``T_{0,k}V_k/X``,
    ``sum v_i (v_j/2 + X d_X v_j)``,
    ``sum U_i Z_{0,j}V_j/X``,
    ``-2(2 d_X v_k + X d_XX v_k)``, and
    ``-Z^[2]_{0,k-1}V_{k-1}/X``.

    The last term is zero for ``k=0`` by the paper's negative-index convention.
    The inputs therefore consist only of already-known pointwise coefficient
    jets.  This function does not construct those profiles or solve (5.3)-(5.4).
    """

    k = _nonnegative_order(k)
    X, eta = _point(X, eta)
    h = validate_h(h)
    if len(radial_flux_jets) < k + 1:
        raise ValueError("radial_flux_jets must contain v_0 through v_k")
    if len(U_values) < k + 1:
        raise ValueError("U_values must contain U_0 through U_k")

    jets = [radial_flux_jets[i].checked(f"v_{i}") for i in range(k + 1)]
    U = [_finite(U_values[i], f"U_{i}") for i in range(k + 1)]
    lam = lambda n: 2.0 * n * h

    scaling = regular_T_on_Xv_over_X(jets[k], X, eta, b=lam(k), h=h)
    radial_transport = 0.0
    axial_transport = 0.0
    for i in range(k + 1):
        j = k - i
        radial_transport += jets[i].value * (
            0.5 * jets[j].value + X * jets[j].dX
        )
        axial_transport += U[i] * regular_Z_on_Xv_over_X(
            jets[j], X, eta, b=lam(j), h=h
        )

    radial_viscosity = -2.0 * (2.0 * jets[k].dX + X * jets[k].dXX)
    shifted_axial_viscosity = 0.0
    if k >= 1:
        shifted_axial_viscosity = -regular_Z2_on_Xv_over_X(
            jets[k - 1], X, eta, b=lam(k - 1), h=h
        )

    return OmegaRowTerms(
        scaling=scaling,
        radial_transport=radial_transport,
        axial_transport=axial_transport,
        radial_viscosity=radial_viscosity,
        shifted_axial_viscosity=shifted_axial_viscosity,
    )


def omega_over_x_eq_5_6(
    k: int,
    X: float,
    eta: float,
    radial_flux_jets: Sequence[RegularFluxJet],
    U_values: Sequence[float],
    *,
    h: float,
) -> float:
    """Return the regular source ``Omega_k/X`` required by Eq. (5.5)."""

    return omega_row_terms_eq_5_6(
        k, X, eta, radial_flux_jets, U_values, h=h
    ).total


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
