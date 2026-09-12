"""Hierarchy-owned analytic eta jet for the Section 5 strict-lower source.

The landed ``actualLowerSource`` bridge already determines the four source
scalars from strict-lower coefficient second jets.  Separate landed adapters
now also provide every stronger derivative needed to differentiate that source:

* repaired ``phi`` third-mixed jets for angular preceding diffusion;
* repaired ``U`` third/fourth-mixed jets for axial preceding diffusion and
  Eq. (5.6);
* analytic Eq. (5.2) regular-flux jets; and
* the hierarchy-owned analytic ``partial_eta(Omega/X)`` bridge.

This module combines those existing pieces into one fail-closed analytic
``eta`` jet of ``actualLowerSource``.  It accepts only the strong repaired
history hierarchy, so the derivative cannot be paired with unrelated
caller-supplied beta/Omega tables.  No finite difference is used in production.

Issue #1 still owns the genuine leading mixed profile.  Therefore this is
``formal-structure`` solver infrastructure, not a paper-exact source and not yet
a complete ``partial_eta f_n`` or ``k=1`` Eq. (5.7) Picard application.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import math

from .background_lower_history_source import (
    PositiveAxisPointData,
    ProfileSecondJet,
    positive_axis_point_data_from_lower_history,
)
from .background_positive_axis import PositiveAxisSourceJet
from .background_preceding_diffusion_parameter_jet import (
    ProfileThirdMixedJet,
    preceding_diffusion_parameter_jet,
)
from .background_regular_flux_second_jets import AxialThirdMixedJet
from .background_repaired_history_omega_parameter import (
    hierarchy_owned_omega_parameter_jet,
)
from .background_repaired_history_phi_third_mixed import (
    Section5LowerHistoryPhiThirdMixedHierarchy,
)


@dataclass(frozen=True)
class LowerSourceParameterJet:
    """Value and analytic ``eta`` derivative of ``actualLowerSource``.

    ``parameter`` reuses ``PositiveAxisSourceJet`` as a named four-vector: each
    field stores the eta derivative of the correspondingly named source field.
    """

    value: PositiveAxisSourceJet
    parameter: PositiveAxisSourceJet

    def __post_init__(self) -> None:
        if not isinstance(self.value, PositiveAxisSourceJet):
            raise TypeError("value must be a PositiveAxisSourceJet")
        if not isinstance(self.parameter, PositiveAxisSourceJet):
            raise TypeError("parameter must be a PositiveAxisSourceJet")


def _positive_order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("order must be a positive integer")
    return int(value)


def _point(X: float, eta: float) -> tuple[float, float]:
    X = float(X)
    eta = float(eta)
    if not math.isfinite(X) or X < 0.0:
        raise ValueError("X must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return X, eta


def _axial_value_parameter(
    h: float,
    power: float,
    X: float,
    eta: float,
    jet: ProfileSecondJet,
) -> tuple[float, float]:
    """Return ``(Z_power F, partial_eta Z_power F)`` analytically."""

    d = 1.0 - eta * eta
    ell = 1.0 - 2.0 * h * eta * eta
    numerator = (
        2.0 * eta * power * jet.value
        + d * jet.parameter
        - 2.0 * eta * X * jet.radial
    )
    numerator_eta = (
        2.0 * power * jet.value
        + 2.0 * eta * (power - 1.0) * jet.parameter
        + d * jet.parameter2
        - 2.0 * X * jet.radial
        - 2.0 * eta * X * jet.radial_parameter
    )
    value = numerator / ell
    parameter = numerator_eta / ell + 4.0 * h * eta * numerator / (ell * ell)
    return value, parameter


def _profile_third_from_axial(jet: AxialThirdMixedJet) -> ProfileThirdMixedJet:
    if not isinstance(jet, AxialThirdMixedJet):
        raise TypeError("jet must be an AxialThirdMixedJet")
    return ProfileThirdMixedJet(
        value=jet.value,
        radial=jet.radial,
        radial2=jet.radial2,
        parameter=jet.parameter,
        radial_parameter=jet.radial_parameter,
        parameter2=jet.parameter2,
        radial2_parameter=jet.radial2_parameter,
        radial_parameter2=jet.radial_parameter2,
        parameter3=jet.parameter3,
    )


def hierarchy_owned_lower_source_parameter_jet(
    hierarchy: Section5LowerHistoryPhiThirdMixedHierarchy,
    order: int,
    X: float,
    eta: float,
) -> LowerSourceParameterJet:
    """Return hierarchy-owned ``(actualLowerSource, partial_eta source)``.

    ``order`` is the positive coefficient order being solved, so only
    coefficients ``0,...,order-1`` are consumed.  The ordinary source value is
    delegated to the already-landed ``actualLowerSource`` implementation.  Its
    eta derivative is then evaluated analytically from the same hierarchy-owned
    coefficient data.

    The strong hierarchy requirement is deliberate.  In particular, order one
    requires genuine strong order-zero phi and U data; until Issue #1 provides
    those data, production use fails closed instead of fabricating a leading
    derivative record.
    """

    if not isinstance(hierarchy, Section5LowerHistoryPhiThirdMixedHierarchy):
        raise TypeError(
            "hierarchy must be a Section5LowerHistoryPhiThirdMixedHierarchy"
        )
    order = _positive_order(order)
    X, eta = _point(X, eta)
    if len(hierarchy.sources) < order:
        raise ValueError(
            "requested source order is missing one or more strict-lower coefficients"
        )

    phi = [hierarchy.phi_second_jet(j, X, eta) for j in range(order)]
    axial = [hierarchy.axial_second_jet(j, X, eta) for j in range(order)]
    beta = [hierarchy.beta_second_jet(j, X, eta) for j in range(order)]

    point: PositiveAxisPointData = positive_axis_point_data_from_lower_history(
        hierarchy.h,
        order,
        X,
        eta,
        phi,
        axial,
        beta,
    )

    angular_power = -1.0 - hierarchy.h
    axial_power = -0.5 - hierarchy.h
    angular_parameter = 0.0
    axial_parameter = 0.0
    pressure_parameter = 0.0

    for i in range(1, order):
        j = order - i
        lam_j = 2.0 * j * hierarchy.h

        angular_gradient = X * phi[j].radial + phi[j].value
        angular_gradient_parameter = (
            X * phi[j].radial_parameter + phi[j].parameter
        )
        angular_parameter += (
            beta[i].parameter * angular_gradient
            + beta[i].value * angular_gradient_parameter
        )
        angular_z, angular_z_parameter = _axial_value_parameter(
            hierarchy.h,
            angular_power + lam_j,
            X,
            eta,
            phi[j],
        )
        angular_parameter += (
            axial[i].parameter * angular_z
            + axial[i].value * angular_z_parameter
        )

        axial_parameter += (
            beta[i].parameter * X * axial[j].radial
            + beta[i].value * X * axial[j].radial_parameter
        )
        axial_z, axial_z_parameter = _axial_value_parameter(
            hierarchy.h,
            axial_power + lam_j,
            X,
            eta,
            axial[j],
        )
        axial_parameter += (
            axial[i].parameter * axial_z
            + axial[i].value * axial_z_parameter
        )

        pressure_parameter += (
            phi[i].parameter * phi[j].value
            + phi[i].value * phi[j].parameter
        )

    angular_preceding = preceding_diffusion_parameter_jet(
        hierarchy.h,
        angular_power,
        order,
        X,
        eta,
        hierarchy.phi_third_mixed_jet(order - 1, X, eta),
    )
    angular_parameter -= angular_preceding.parameter

    axial_preceding = preceding_diffusion_parameter_jet(
        hierarchy.h,
        axial_power,
        order,
        X,
        eta,
        _profile_third_from_axial(
            hierarchy.axial_third_mixed_jet(order - 1, X, eta)
        ),
    )
    axial_parameter -= axial_preceding.parameter

    omega = hierarchy_owned_omega_parameter_jet(
        hierarchy,
        order - 1,
        X,
        eta,
    )

    parameter = PositiveAxisSourceJet(
        angular=angular_parameter,
        axial=axial_parameter,
        pressure_product=pressure_parameter,
        omega_quotient=omega.parameter,
    )
    return LowerSourceParameterJet(value=point.source, parameter=parameter)
