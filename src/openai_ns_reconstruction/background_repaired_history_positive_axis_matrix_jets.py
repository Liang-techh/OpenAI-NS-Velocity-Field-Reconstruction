"""Hierarchy-owned eta jets of the Section 5 PositiveAxis matrices.

For positive coefficient order ``n``, the displayed Eq. (5.7) matrices ``A0``
and ``A1`` depend only on the order-zero leading/base profile.  The landed
hierarchy bridge already owns their value rows.  This module differentiates
those same displayed formulas analytically through ``partial_eta^2`` while
taking every profile row from coefficient order zero of
``Section5LowerHistorySixthMixedHierarchy``.

The implementation uses a small exact derivative algebra for value/first/second
eta derivatives.  Quotients are differentiated by the identity ``b*q=a``;
there is no finite difference, fitted matrix derivative, sampled ``V/X``
division, generic cutoff, or caller-maintained matrix/source derivative table.

This remains Stage-2 ``formal-structure`` infrastructure.  The underlying
order-zero strong profile is still an Issue-1/upstream input, and these finite
matrix jets do not establish a Picard fixed point, an all-order coefficient
hierarchy, convergence, or paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import math

import numpy as np

from .background_repaired_history_positive_axis_base import (
    hierarchy_owned_positive_axis_matrices,
)
from .background_repaired_history_sixth_mixed import (
    Section5LowerHistorySixthMixedHierarchy,
)


@dataclass(frozen=True)
class _EtaJet2:
    """Scalar value and first two exact eta derivatives."""

    value: float
    d1: float = 0.0
    d2: float = 0.0

    def __post_init__(self) -> None:
        for name in ("value", "d1", "d2"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)

    @staticmethod
    def coerce(value: "_EtaJet2 | float") -> "_EtaJet2":
        if isinstance(value, _EtaJet2):
            return value
        return _EtaJet2(float(value))

    def __add__(self, other: "_EtaJet2 | float") -> "_EtaJet2":
        other = self.coerce(other)
        return _EtaJet2(
            self.value + other.value,
            self.d1 + other.d1,
            self.d2 + other.d2,
        )

    __radd__ = __add__

    def __neg__(self) -> "_EtaJet2":
        return _EtaJet2(-self.value, -self.d1, -self.d2)

    def __sub__(self, other: "_EtaJet2 | float") -> "_EtaJet2":
        return self + (-self.coerce(other))

    def __rsub__(self, other: "_EtaJet2 | float") -> "_EtaJet2":
        return self.coerce(other) - self

    def __mul__(self, other: "_EtaJet2 | float") -> "_EtaJet2":
        other = self.coerce(other)
        return _EtaJet2(
            self.value * other.value,
            self.d1 * other.value + self.value * other.d1,
            self.d2 * other.value + 2.0 * self.d1 * other.d1 + self.value * other.d2,
        )

    __rmul__ = __mul__

    def __truediv__(self, other: "_EtaJet2 | float") -> "_EtaJet2":
        other = self.coerce(other)
        if other.value == 0.0:
            raise ZeroDivisionError("eta-jet denominator has zero value")
        q0 = self.value / other.value
        q1 = (self.d1 - other.d1 * q0) / other.value
        q2 = (
            self.d2
            - 2.0 * other.d1 * q1
            - other.d2 * q0
        ) / other.value
        return _EtaJet2(q0, q1, q2)

    def __rtruediv__(self, other: "_EtaJet2 | float") -> "_EtaJet2":
        return self.coerce(other) / self


def _positive_order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("order must be a positive integer")
    return int(value)


def _point(xi: float, eta: float) -> tuple[float, float]:
    xi = float(xi)
    eta = float(eta)
    if not math.isfinite(xi) or xi < 0.0:
        raise ValueError("xi must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return xi, eta


def _matrix6(value: np.ndarray, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (6, 6) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite matrix of shape (6, 6)")
    out = out.copy()
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class HierarchyOwnedPositiveAxisMatrixSecondParameterJet:
    """Value and first two eta derivatives of the exact Eq. (5.7) matrices."""

    A0: np.ndarray
    A0_parameter: np.ndarray
    A0_parameter2: np.ndarray
    A1: np.ndarray
    A1_parameter: np.ndarray
    A1_parameter2: np.ndarray

    def __post_init__(self) -> None:
        for name in (
            "A0",
            "A0_parameter",
            "A0_parameter2",
            "A1",
            "A1_parameter",
            "A1_parameter2",
        ):
            object.__setattr__(self, name, _matrix6(getattr(self, name), name))


def hierarchy_owned_positive_axis_matrix_second_parameter_jet(
    hierarchy: Section5LowerHistorySixthMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
) -> HierarchyOwnedPositiveAxisMatrixSecondParameterJet:
    """Build ``A0/A1`` through ``partial_eta^2`` from one strong hierarchy.

    The value matrices delegate to the landed hierarchy-owned value bridge.
    Derivative rows replay the displayed Eq. (5.7) formulas using only
    order-zero ``phi/U`` mixed jets and the order-zero analytic Eq. (5.2)
    ``beta=V/X`` jet from the same hierarchy.
    """

    if not isinstance(hierarchy, Section5LowerHistorySixthMixedHierarchy):
        raise TypeError("hierarchy must be a Section5LowerHistorySixthMixedHierarchy")
    order = _positive_order(order)
    xi, eta = _point(xi, eta)
    X = xi * xi

    values = hierarchy_owned_positive_axis_matrices(hierarchy, order, xi, eta)
    phi = hierarchy.phi_fifth_mixed_jet(0, X, eta)
    axial = hierarchy.axial_sixth_mixed_jet(0, X, eta)
    beta = hierarchy.beta_fifth_mixed_jet(0, X, eta)

    h = float(hierarchy.h)
    C = float(hierarchy.C)
    lam = 2.0 * order * h
    a = 0.5 + h
    d = 0.5 - h
    angular_power = -a - 0.5
    axial_power = -a
    inv_c2 = 1.0 / (C * C)

    eta_j = _EtaJet2(eta, 1.0, 0.0)
    edge = 1.0 - eta_j * eta_j
    ell = 1.0 - 2.0 * h * eta_j * eta_j
    if ell.value <= 0.0:
        raise ValueError("ell(h,eta) must be positive")

    phi_value = _EtaJet2(phi.value, phi.parameter, phi.parameter2)
    phi_radial = _EtaJet2(
        phi.radial,
        phi.radial_parameter,
        phi.radial_parameter2,
    )
    phi_parameter = _EtaJet2(
        phi.parameter,
        phi.parameter2,
        phi.parameter3,
    )

    axial_value = _EtaJet2(axial.value, axial.parameter, axial.parameter2)
    axial_radial = _EtaJet2(
        axial.radial,
        axial.radial_parameter,
        axial.radial_parameter2,
    )
    axial_parameter = _EtaJet2(
        axial.parameter,
        axial.parameter2,
        axial.parameter3,
    )
    beta_value = _EtaJet2(beta.value, beta.parameter, beta.parameter2)

    M = 1.0 - 2.0 * eta_j * axial_value
    R = M / ell + beta_value
    Gphi = X * phi_radial + phi_value
    Gu = X * axial_radial

    def transported_value(
        power: float,
        value: _EtaJet2,
        parameter: _EtaJet2,
        radial: _EtaJet2,
    ) -> _EtaJet2:
        return (
            2.0 * eta_j * power * value
            + edge * parameter
            - 2.0 * eta_j * X * radial
        ) / ell

    A0_parameter = np.zeros((6, 6), dtype=float)
    A0_parameter2 = np.zeros((6, 6), dtype=float)
    A1_parameter = np.zeros((6, 6), dtype=float)
    A1_parameter2 = np.zeros((6, 6), dtype=float)

    def assign(
        first: np.ndarray,
        second: np.ndarray,
        row: int,
        column: int,
        jet: _EtaJet2,
    ) -> None:
        first[row, column] = jet.d1
        second[row, column] = jet.d2

    assign(
        A0_parameter,
        A0_parameter2,
        3,
        0,
        4.0 * xi * inv_c2 * phi_value,
    )
    assign(
        A0_parameter,
        A0_parameter2,
        4,
        0,
        2.0 * (beta_value - (angular_power + lam) * M / ell),
    )
    assign(
        A0_parameter,
        A0_parameter2,
        4,
        1,
        2.0
        * (
            transported_value(
                angular_power,
                phi_value,
                phi_parameter,
                phi_radial,
            )
            + 2.0 * eta_j * (a - lam) * Gphi / ell
        ),
    )
    assign(
        A0_parameter,
        A0_parameter2,
        4,
        2,
        -4.0 * eta_j * (d + lam) * Gphi / ell,
    )
    assign(
        A0_parameter,
        A0_parameter2,
        4,
        4,
        xi * R,
    )
    assign(
        A0_parameter,
        A0_parameter2,
        5,
        0,
        -8.0 * eta_j * X * inv_c2 * phi_value / ell,
    )
    assign(
        A0_parameter,
        A0_parameter2,
        5,
        1,
        2.0
        * (
            -(axial_power + lam) * M / ell
            + transported_value(
                axial_power,
                axial_value,
                axial_parameter,
                axial_radial,
            )
            + 2.0 * eta_j * (a - lam) * Gu / ell
        ),
    )
    assign(
        A0_parameter,
        A0_parameter2,
        5,
        2,
        -4.0 * eta_j * (d + lam) * Gu / ell,
    )
    assign(
        A0_parameter,
        A0_parameter2,
        5,
        3,
        4.0 * eta_j * (-2.0 * a + lam) / ell,
    )
    assign(
        A0_parameter,
        A0_parameter2,
        5,
        5,
        xi * R,
    )

    H = d * eta_j + edge * axial_value
    assign(
        A1_parameter,
        A1_parameter2,
        4,
        0,
        2.0 * H / ell,
    )
    assign(
        A1_parameter,
        A1_parameter2,
        4,
        1,
        -2.0 * edge * Gphi / ell,
    )
    assign(
        A1_parameter,
        A1_parameter2,
        4,
        2,
        -2.0 * edge * Gphi / ell,
    )
    assign(
        A1_parameter,
        A1_parameter2,
        5,
        1,
        2.0 * (H - edge * Gu) / ell,
    )
    assign(
        A1_parameter,
        A1_parameter2,
        5,
        2,
        -2.0 * edge * Gu / ell,
    )
    assign(
        A1_parameter,
        A1_parameter2,
        5,
        3,
        2.0 * edge / ell,
    )

    return HierarchyOwnedPositiveAxisMatrixSecondParameterJet(
        A0=values.A0,
        A0_parameter=A0_parameter,
        A0_parameter2=A0_parameter2,
        A1=values.A1,
        A1_parameter=A1_parameter,
        A1_parameter2=A1_parameter2,
    )
