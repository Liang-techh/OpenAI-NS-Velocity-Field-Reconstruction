"""Hierarchy-owned analytic second-eta jets of the Eq. (5.7) matrices.

The landed matrix-parameter bridge owns ``(A0, partial_eta A0)`` and
``(A1, partial_eta A1)`` from the strong Section-5 lower-history hierarchy.
Advancing the differentiated Picard recursion one more level requires the
second eta derivatives of those same matrices.

This module computes only that new derivative layer.  Values and first eta
derivatives are delegated unchanged to the landed matrix-parameter bridge, so
there is no duplicate lower-order implementation.  The second derivatives are
obtained by exact second-order scalar jet arithmetic applied to the displayed
PositiveAxis formulas, using only hierarchy-owned third-mixed ``phi_0/U_0``
data and the Eq. (5.2)-derived second jet of ``beta_0=V_0/X``.  No production
finite difference, sampled matrix table, generic cutoff, fitted coefficient, or
caller-supplied matrix derivative is admitted.

This remains Stage-2 ``formal-structure`` infrastructure.  Genuine Issue-#1
leading profile data remain upstream, and this matrix jet alone does not prove
Picard convergence or materialize a paper-exact positive-order coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .background_repaired_history_matrix_parameter import (
    hierarchy_owned_positive_axis_matrix_parameter_jets,
)
from .background_repaired_history_phi_third_mixed import (
    Section5LowerHistoryPhiThirdMixedHierarchy,
)


def _matrix6(value: np.ndarray, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (6, 6) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite matrix of shape (6, 6)")
    out = out.copy()
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class PositiveAxisMatrixSecondParameterJet:
    """Value and first two analytic eta derivatives of one Eq. (5.7) matrix."""

    value: np.ndarray
    parameter: np.ndarray
    parameter2: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", _matrix6(self.value, "value"))
        object.__setattr__(self, "parameter", _matrix6(self.parameter, "parameter"))
        object.__setattr__(self, "parameter2", _matrix6(self.parameter2, "parameter2"))


@dataclass(frozen=True)
class PositiveAxisMatrixSecondParameterJets:
    """Second eta jets for hierarchy-owned ``A0`` and ``A1``."""

    A0: PositiveAxisMatrixSecondParameterJet
    A1: PositiveAxisMatrixSecondParameterJet

    def __post_init__(self) -> None:
        if not isinstance(self.A0, PositiveAxisMatrixSecondParameterJet) or not isinstance(
            self.A1, PositiveAxisMatrixSecondParameterJet
        ):
            raise TypeError(
                "A0 and A1 must be PositiveAxisMatrixSecondParameterJet instances"
            )


@dataclass(frozen=True)
class _Eta2:
    """Scalar value/first/second derivative jet in eta."""

    value: float
    parameter: float = 0.0
    parameter2: float = 0.0

    def __post_init__(self) -> None:
        for name in ("value", "parameter", "parameter2"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)

    @staticmethod
    def coerce(value: "_Eta2 | float") -> "_Eta2":
        if isinstance(value, _Eta2):
            return value
        return _Eta2(float(value))

    def __add__(self, other: "_Eta2 | float") -> "_Eta2":
        other = self.coerce(other)
        return _Eta2(
            self.value + other.value,
            self.parameter + other.parameter,
            self.parameter2 + other.parameter2,
        )

    __radd__ = __add__

    def __neg__(self) -> "_Eta2":
        return _Eta2(-self.value, -self.parameter, -self.parameter2)

    def __sub__(self, other: "_Eta2 | float") -> "_Eta2":
        return self + (-self.coerce(other))

    def __rsub__(self, other: "_Eta2 | float") -> "_Eta2":
        return self.coerce(other) - self

    def __mul__(self, other: "_Eta2 | float") -> "_Eta2":
        other = self.coerce(other)
        return _Eta2(
            self.value * other.value,
            self.parameter * other.value + self.value * other.parameter,
            self.parameter2 * other.value
            + 2.0 * self.parameter * other.parameter
            + self.value * other.parameter2,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "_Eta2":
        if self.value == 0.0:
            raise ZeroDivisionError("cannot invert a zero eta jet")
        inverse = 1.0 / self.value
        return _Eta2(
            inverse,
            -self.parameter * inverse * inverse,
            2.0 * self.parameter * self.parameter * inverse**3
            - self.parameter2 * inverse * inverse,
        )

    def __truediv__(self, other: "_Eta2 | float") -> "_Eta2":
        return self * self.coerce(other).reciprocal()

    def __rtruediv__(self, other: "_Eta2 | float") -> "_Eta2":
        return self.coerce(other) * self.reciprocal()


def _entry_matrix(entries: list[list[_Eta2]]) -> np.ndarray:
    out = np.array(
        [[entry.parameter2 for entry in row] for row in entries],
        dtype=float,
    )
    if out.shape != (6, 6) or not np.all(np.isfinite(out)):
        raise OverflowError("positive-axis matrix second eta jet is outside binary64 range")
    return out


def hierarchy_owned_positive_axis_matrix_second_parameter_jets(
    hierarchy: Section5LowerHistoryPhiThirdMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
) -> PositiveAxisMatrixSecondParameterJets:
    """Return hierarchy-owned ``(A,d_eta A,d_eta^2 A)`` for ``A0`` and ``A1``.

    The existing hierarchy-owned matrix-parameter implementation remains
    authoritative for value and first derivative.  Only the second derivative
    is new here.  Its scalar inputs are the order-zero third-mixed ``phi/U``
    jets and the order-zero Eq. (5.2)-derived ``beta`` second jet, which are
    sufficient because the displayed matrices contain at most one explicit
    ``partial_eta`` of a base profile.
    """

    if not isinstance(hierarchy, Section5LowerHistoryPhiThirdMixedHierarchy):
        raise TypeError(
            "hierarchy must be a Section5LowerHistoryPhiThirdMixedHierarchy"
        )

    lower = hierarchy_owned_positive_axis_matrix_parameter_jets(
        hierarchy,
        order,
        xi,
        eta,
    )

    xi = float(xi)
    eta = float(eta)
    X = xi * xi
    h = hierarchy.h
    C = hierarchy.C

    phi = hierarchy.phi_third_mixed_jet(0, X, eta)
    axial = hierarchy.axial_third_mixed_jet(0, X, eta)
    beta = hierarchy.beta_second_jet(0, X, eta)

    E = _Eta2(eta, 1.0, 0.0)
    edge = 1.0 - E * E
    ell = 1.0 - 2.0 * h * E * E
    if ell.value <= 0.0:
        raise ValueError("ell(h,eta) must be positive")

    phi_value = _Eta2(phi.value, phi.parameter, phi.parameter2)
    phi_radial = _Eta2(
        phi.radial,
        phi.radial_parameter,
        phi.radial_parameter2,
    )
    phi_parameter = _Eta2(
        phi.parameter,
        phi.parameter2,
        phi.parameter3,
    )
    axial_value = _Eta2(axial.value, axial.parameter, axial.parameter2)
    axial_radial = _Eta2(
        axial.radial,
        axial.radial_parameter,
        axial.radial_parameter2,
    )
    axial_parameter = _Eta2(
        axial.parameter,
        axial.parameter2,
        axial.parameter3,
    )
    beta_value = _Eta2(beta.value, beta.parameter, beta.parameter2)

    lam = 2.0 * int(order) * h
    a = 0.5 + h
    d = 0.5 - h
    angular_power = -a - 0.5
    axial_power = -a
    inverse_c2 = 1.0 / (C * C)

    M = 1.0 - 2.0 * E * axial_value
    R = M / ell + beta_value
    Gphi = X * phi_radial + phi_value
    Gu = X * axial_radial

    def axial_eta_value(
        power: float,
        value: _Eta2,
        parameter: _Eta2,
        radial: _Eta2,
    ) -> _Eta2:
        return (
            2.0 * E * power * value
            + edge * parameter
            - 2.0 * E * X * radial
        ) / ell

    zero = _Eta2(0.0)
    A0 = [[zero for _ in range(6)] for _ in range(6)]
    A1 = [[zero for _ in range(6)] for _ in range(6)]

    A0[0][4] = _Eta2(1.0)
    A0[1][5] = _Eta2(1.0)
    A0[2][5] = _Eta2(-1.0)
    A0[3][0] = 4.0 * xi * inverse_c2 * phi_value
    A0[4][0] = 2.0 * (beta_value - (angular_power + lam) * M / ell)
    A0[4][1] = 2.0 * (
        axial_eta_value(angular_power, phi_value, phi_parameter, phi_radial)
        + 2.0 * E * (a - lam) * Gphi / ell
    )
    A0[4][2] = -4.0 * E * (d + lam) * Gphi / ell
    A0[4][4] = xi * R
    A0[5][0] = -8.0 * E * X * inverse_c2 * phi_value / ell
    A0[5][1] = 2.0 * (
        -(axial_power + lam) * M / ell
        + axial_eta_value(axial_power, axial_value, axial_parameter, axial_radial)
        + 2.0 * E * (a - lam) * Gu / ell
    )
    A0[5][2] = -4.0 * E * (d + lam) * Gu / ell
    A0[5][3] = 4.0 * E * (-2.0 * a + lam) / ell
    A0[5][5] = xi * R

    H = d * E + edge * axial_value
    A1[4][0] = 2.0 * H / ell
    A1[4][1] = -2.0 * edge * Gphi / ell
    A1[4][2] = A1[4][1]
    A1[5][1] = 2.0 * (H - edge * Gu) / ell
    A1[5][2] = -2.0 * edge * Gu / ell
    A1[5][3] = 2.0 * edge / ell

    return PositiveAxisMatrixSecondParameterJets(
        A0=PositiveAxisMatrixSecondParameterJet(
            lower.A0.value,
            lower.A0.parameter,
            _entry_matrix(A0),
        ),
        A1=PositiveAxisMatrixSecondParameterJet(
            lower.A1.value,
            lower.A1.parameter,
            _entry_matrix(A1),
        ),
    )
