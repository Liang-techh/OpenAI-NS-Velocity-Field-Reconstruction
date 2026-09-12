"""Hierarchy-owned analytic ``eta`` jets of the Eq. (5.7) matrices.

PR #184 landed the first nontrivial hierarchy-owned Picard value

    W_n^(1) = G(A0 W_n^(0) + A1 partial_eta W_n^(0) + f_n).

Differentiating that iterate requires, among other terms, analytic
``partial_eta A0`` and ``partial_eta A1``.  The displayed PositiveAxis matrices
only depend on the order-zero base profiles through the second ``(X, eta)``
jets of ``phi_0``, ``U_0`` and ``beta_0=V_0/X``.  Those jets are already owned
by ``Section5LowerHistoryPhiThirdMixedHierarchy``.

This module differentiates the landed matrix formulas analytically at fixed
``xi`` (hence fixed ``X=xi^2``).  No finite difference, sampled matrix table,
generic cutoff, fitted coefficient, or caller-supplied matrix derivative is
accepted on the production path.

The result is still Stage-2 ``formal-structure`` infrastructure.  In
particular it does not by itself provide ``partial_eta W_n^(1)``: that derivative
also contains ``A1 partial_eta^2 W_n^(0)``, so a hierarchy-owned second eta jet
of the first Picard iterate/forcing is still required.  The genuine Issue-#1
leading profile remains upstream and ``paper_exact_velocity_available=false``.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import math

import numpy as np

from .background_lower_history_source import ProfileSecondJet
from .background_positive_axis import (
    PositiveAxisBaseJet,
    positive_axis_A0,
    positive_axis_A1,
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
class PositiveAxisMatrixParameterJet:
    """Value and analytic ``eta`` derivative of one Eq. (5.7) matrix field."""

    value: np.ndarray
    parameter: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", _matrix6(self.value, "value"))
        object.__setattr__(self, "parameter", _matrix6(self.parameter, "parameter"))


@dataclass(frozen=True)
class PositiveAxisMatrixParameterJets:
    """Analytic parameter jets for the hierarchy-owned ``A0`` and ``A1``."""

    A0: PositiveAxisMatrixParameterJet
    A1: PositiveAxisMatrixParameterJet

    def __post_init__(self) -> None:
        if not isinstance(self.A0, PositiveAxisMatrixParameterJet) or not isinstance(
            self.A1, PositiveAxisMatrixParameterJet
        ):
            raise TypeError("A0 and A1 must be PositiveAxisMatrixParameterJet instances")


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


def _axial_value_parameter(
    power: float,
    X: float,
    eta: float,
    jet: ProfileSecondJet,
    ell: float,
    inverse_ell_parameter: float,
) -> float:
    """Differentiate the PositiveAxis ``axialValue`` quotient in ``eta``."""

    edge = 1.0 - eta * eta
    numerator = (
        2.0 * eta * power * jet.value
        + edge * jet.parameter
        - 2.0 * eta * X * jet.radial
    )
    numerator_parameter = (
        2.0 * power * jet.value
        + 2.0 * eta * power * jet.parameter
        - 2.0 * eta * jet.parameter
        + edge * jet.parameter2
        - 2.0 * X * jet.radial
        - 2.0 * eta * X * jet.radial_parameter
    )
    return (
        numerator_parameter / ell
        + numerator * inverse_ell_parameter
    )


def hierarchy_owned_positive_axis_matrix_parameter_jets(
    hierarchy: Section5LowerHistoryPhiThirdMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
) -> PositiveAxisMatrixParameterJets:
    """Return hierarchy-owned ``(A0,d_eta A0)`` and ``(A1,d_eta A1)``.

    The matrix values are evaluated by the already-landed PositiveAxis formulas.
    Their derivatives are derived only from the hierarchy-owned order-zero
    second jets.  ``beta_0`` and ``d_eta beta_0`` therefore continue to come
    from the analytic Eq. (5.2) regular-flux adapter rather than an independent
    derivative table.
    """

    if not isinstance(hierarchy, Section5LowerHistoryPhiThirdMixedHierarchy):
        raise TypeError(
            "hierarchy must be a Section5LowerHistoryPhiThirdMixedHierarchy"
        )
    order = _positive_order(order)
    xi, eta = _point(xi, eta)

    # Public hierarchy validation: recursive order n is legal only when every
    # strict-lower coefficient 0,...,n-1 is owned.  Constructing the field
    # callbacks performs that check without evaluating unrelated source terms.
    hierarchy.positive_axis_fields(order)

    X = xi * xi
    phi = hierarchy.phi_second_jet(0, X, eta)
    axial = hierarchy.axial_second_jet(0, X, eta)
    beta = hierarchy.beta_second_jet(0, X, eta)

    base = PositiveAxisBaseJet(
        phi=phi.first_jet(),
        axial=axial.first_jet(),
        beta=beta.value,
    )
    value_A0 = positive_axis_A0(
        hierarchy.h,
        order,
        hierarchy.C,
        xi,
        eta,
        base,
    )
    value_A1 = positive_axis_A1(
        hierarchy.h,
        xi,
        eta,
        base,
    )

    h = hierarchy.h
    C = hierarchy.C
    lam = 2.0 * order * h
    a = 0.5 + h
    d = 0.5 - h
    edge = 1.0 - eta * eta
    edge_parameter = -2.0 * eta
    ell = 1.0 - 2.0 * h * eta * eta
    if ell <= 0.0:
        raise ValueError("ell(h,eta) must be positive")
    inverse_ell_parameter = 4.0 * h * eta / (ell * ell)
    inverse_c2 = 1.0 / (C * C)

    angular_power = -a - 0.5
    axial_power = -a

    M = 1.0 - 2.0 * eta * axial.value
    M_parameter = -2.0 * axial.value - 2.0 * eta * axial.parameter
    M_over_ell_parameter = (
        M_parameter / ell + M * inverse_ell_parameter
    )

    R_parameter = M_over_ell_parameter + beta.parameter

    Gphi = X * phi.radial + phi.value
    Gphi_parameter = X * phi.radial_parameter + phi.parameter
    Gu = X * axial.radial
    Gu_parameter = X * axial.radial_parameter

    phi_axial_value_parameter = _axial_value_parameter(
        angular_power,
        X,
        eta,
        phi,
        ell,
        inverse_ell_parameter,
    )
    axial_axial_value_parameter = _axial_value_parameter(
        axial_power,
        X,
        eta,
        axial,
        ell,
        inverse_ell_parameter,
    )

    def eta_times_over_ell_parameter(value: float, parameter: float) -> float:
        return (
            (value + eta * parameter) / ell
            + eta * value * inverse_ell_parameter
        )

    dA0 = np.zeros((6, 6), dtype=float)
    dA0[3, 0] = 4.0 * xi * inverse_c2 * phi.parameter
    dA0[4, 0] = 2.0 * (
        beta.parameter - (angular_power + lam) * M_over_ell_parameter
    )
    dA0[4, 1] = 2.0 * (
        phi_axial_value_parameter
        + 2.0
        * (a - lam)
        * eta_times_over_ell_parameter(Gphi, Gphi_parameter)
    )
    dA0[4, 2] = -4.0 * (d + lam) * eta_times_over_ell_parameter(
        Gphi,
        Gphi_parameter,
    )
    dA0[4, 4] = xi * R_parameter
    dA0[5, 0] = -8.0 * X * inverse_c2 * eta_times_over_ell_parameter(
        phi.value,
        phi.parameter,
    )
    dA0[5, 1] = 2.0 * (
        -(axial_power + lam) * M_over_ell_parameter
        + axial_axial_value_parameter
        + 2.0
        * (a - lam)
        * eta_times_over_ell_parameter(Gu, Gu_parameter)
    )
    dA0[5, 2] = -4.0 * (d + lam) * eta_times_over_ell_parameter(
        Gu,
        Gu_parameter,
    )
    dA0[5, 3] = 4.0 * (-2.0 * a + lam) * (
        1.0 / ell + eta * inverse_ell_parameter
    )
    dA0[5, 5] = xi * R_parameter

    H = d * eta + edge * axial.value
    H_parameter = d + edge_parameter * axial.value + edge * axial.parameter

    dA1 = np.zeros((6, 6), dtype=float)
    dA1[4, 0] = 2.0 * (
        H_parameter / ell + H * inverse_ell_parameter
    )
    edge_Gphi_parameter = edge_parameter * Gphi + edge * Gphi_parameter
    dA1[4, 1] = -2.0 * (
        edge_Gphi_parameter / ell + edge * Gphi * inverse_ell_parameter
    )
    dA1[4, 2] = dA1[4, 1]

    H_minus_edge_Gu = H - edge * Gu
    H_minus_edge_Gu_parameter = (
        H_parameter - edge_parameter * Gu - edge * Gu_parameter
    )
    dA1[5, 1] = 2.0 * (
        H_minus_edge_Gu_parameter / ell
        + H_minus_edge_Gu * inverse_ell_parameter
    )
    edge_Gu_parameter = edge_parameter * Gu + edge * Gu_parameter
    dA1[5, 2] = -2.0 * (
        edge_Gu_parameter / ell + edge * Gu * inverse_ell_parameter
    )
    dA1[5, 3] = 2.0 * (
        edge_parameter / ell + edge * inverse_ell_parameter
    )

    if not np.all(np.isfinite(dA0)) or not np.all(np.isfinite(dA1)):
        raise OverflowError("positive-axis matrix parameter jet is outside binary64 range")

    return PositiveAxisMatrixParameterJets(
        A0=PositiveAxisMatrixParameterJet(value_A0, dA0),
        A1=PositiveAxisMatrixParameterJet(value_A1, dA1),
    )
