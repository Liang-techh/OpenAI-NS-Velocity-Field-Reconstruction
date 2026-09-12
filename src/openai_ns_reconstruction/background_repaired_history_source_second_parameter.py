"""Hierarchy-owned analytic second eta jet for the Section 5 strict-lower source.

The landed first-parameter bridge already owns
``(actualLowerSource, partial_eta actualLowerSource)`` from the strong repaired
Section-5 hierarchy.  Subsequent increments landed exactly the stronger data
needed for one more eta derivative:

* repaired fourth-mixed ``phi`` jets;
* repaired fifth-mixed ``U`` jets and their fourth-mixed projections;
* fourth-mixed ``beta=V/X`` jets derived only through analytic Eq. (5.2);
* analytic second-eta preceding-diffusion jets; and
* the hierarchy-owned analytic second-eta ``Omega/X`` bridge.

This module composes those existing pieces.  The already-landed first-parameter
bridge remains authoritative for the value and first derivative; only the new
second derivative is assembled here.  No finite difference, caller-maintained
beta/Omega derivative table, generic cutoff, or fitted coefficient is admitted
in production.

Genuine order-zero strong data remain an Issue-#1 input.  This is therefore
Stage-2 ``formal-structure`` solver infrastructure, not a paper-exact recursive
coefficient or a convergence claim.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import math

from .background_positive_axis import PositiveAxisSourceJet
from .background_preceding_diffusion_second_parameter_jet import (
    ProfileFourthMixedJet,
    preceding_diffusion_second_parameter_jet,
)
from .background_regular_flux_third_mixed_jets import AxialFourthMixedJet
from .background_repaired_history_omega_second_parameter import (
    hierarchy_owned_omega_second_parameter_jet,
)
from .background_repaired_history_phi_fourth_mixed import (
    Section5LowerHistoryPhiFourthMixedHierarchy,
)
from .background_repaired_history_source_parameter import (
    hierarchy_owned_lower_source_parameter_jet,
)


@dataclass(frozen=True)
class LowerSourceSecondParameterJet:
    """Value and first two analytic eta derivatives of ``actualLowerSource``."""

    value: PositiveAxisSourceJet
    parameter: PositiveAxisSourceJet
    parameter2: PositiveAxisSourceJet

    def __post_init__(self) -> None:
        for name in ("value", "parameter", "parameter2"):
            if not isinstance(getattr(self, name), PositiveAxisSourceJet):
                raise TypeError(f"{name} must be a PositiveAxisSourceJet")


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


def _profile_fourth_from_axial(jet: AxialFourthMixedJet) -> ProfileFourthMixedJet:
    if not isinstance(jet, AxialFourthMixedJet):
        raise TypeError("jet must be an AxialFourthMixedJet")
    return ProfileFourthMixedJet(
        value=jet.value,
        radial=jet.radial,
        radial2=jet.radial2,
        parameter=jet.parameter,
        radial_parameter=jet.radial_parameter,
        parameter2=jet.parameter2,
        radial2_parameter=jet.radial2_parameter,
        radial_parameter2=jet.radial_parameter2,
        parameter3=jet.parameter3,
        radial2_parameter2=jet.radial2_parameter2,
        radial_parameter3=jet.radial_parameter3,
        parameter4=jet.parameter4,
    )


def _axial_value_parameter2(
    h: float,
    power: float,
    X: float,
    eta: float,
    jet: ProfileFourthMixedJet,
) -> tuple[float, float, float]:
    """Return ``(Z_power F, partial_eta Z_power F, partial_eta^2 Z_power F)``."""

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
    numerator_eta2 = (
        (4.0 * power - 2.0) * jet.parameter
        + 2.0 * eta * (power - 2.0) * jet.parameter2
        + d * jet.parameter3
        - 4.0 * X * jet.radial_parameter
        - 2.0 * eta * X * jet.radial_parameter2
    )
    value = numerator / ell
    parameter = numerator_eta / ell + 4.0 * h * eta * numerator / (ell * ell)
    parameter2 = (
        numerator_eta2 / ell
        + 8.0 * h * eta * numerator_eta / (ell * ell)
        + 4.0 * h * numerator / (ell * ell)
        + 32.0 * h * h * eta * eta * numerator / (ell * ell * ell)
    )
    if not all(math.isfinite(v) for v in (value, parameter, parameter2)):
        raise OverflowError("axial second parameter jet is outside binary64 range")
    return value, parameter, parameter2


def hierarchy_owned_lower_source_second_parameter_jet(
    hierarchy: Section5LowerHistoryPhiFourthMixedHierarchy,
    order: int,
    X: float,
    eta: float,
) -> LowerSourceSecondParameterJet:
    """Return hierarchy-owned first two eta derivatives of ``actualLowerSource``.

    ``order`` is the positive coefficient order being solved, so only strict
    lower coefficients ``0,...,order-1`` are consumed.  The existing
    ``hierarchy_owned_lower_source_parameter_jet`` owns the value and first eta
    derivative.  This function analytically assembles only the second derivative
    from the stronger hierarchy-owned jets.

    The strong hierarchy requirement is deliberate: if the genuine leading
    fourth/fifth mixed profile is absent, including at the axis, the call fails
    closed instead of replacing it by unrelated caller data.
    """

    if not isinstance(hierarchy, Section5LowerHistoryPhiFourthMixedHierarchy):
        raise TypeError(
            "hierarchy must be a Section5LowerHistoryPhiFourthMixedHierarchy"
        )
    order = _positive_order(order)
    X, eta = _point(X, eta)
    if len(hierarchy.sources) < order:
        raise ValueError(
            "requested source order is missing one or more strict-lower coefficients"
        )

    first = hierarchy_owned_lower_source_parameter_jet(
        hierarchy,
        order,
        X,
        eta,
    )

    phi = [hierarchy.phi_fourth_mixed_jet(j, X, eta) for j in range(order)]
    axial = [
        hierarchy.axial_fourth_mixed_jet(j, X, eta) for j in range(order)
    ]
    beta = [hierarchy.beta_fourth_mixed_jet(j, X, eta) for j in range(order)]

    angular_power = -1.0 - hierarchy.h
    axial_power = -0.5 - hierarchy.h
    angular_parameter2 = 0.0
    axial_parameter2 = 0.0
    pressure_parameter2 = 0.0

    for i in range(1, order):
        j = order - i
        lam_j = 2.0 * j * hierarchy.h

        angular_gradient = X * phi[j].radial + phi[j].value
        angular_gradient_parameter = (
            X * phi[j].radial_parameter + phi[j].parameter
        )
        angular_gradient_parameter2 = (
            X * phi[j].radial_parameter2 + phi[j].parameter2
        )
        angular_parameter2 += (
            beta[i].parameter2 * angular_gradient
            + 2.0 * beta[i].parameter * angular_gradient_parameter
            + beta[i].value * angular_gradient_parameter2
        )

        angular_z = _axial_value_parameter2(
            hierarchy.h,
            angular_power + lam_j,
            X,
            eta,
            phi[j],
        )
        angular_parameter2 += (
            axial[i].parameter2 * angular_z[0]
            + 2.0 * axial[i].parameter * angular_z[1]
            + axial[i].value * angular_z[2]
        )

        axial_parameter2 += (
            beta[i].parameter2 * X * axial[j].radial
            + 2.0 * beta[i].parameter * X * axial[j].radial_parameter
            + beta[i].value * X * axial[j].radial_parameter2
        )
        axial_z = _axial_value_parameter2(
            hierarchy.h,
            axial_power + lam_j,
            X,
            eta,
            _profile_fourth_from_axial(axial[j]),
        )
        axial_parameter2 += (
            axial[i].parameter2 * axial_z[0]
            + 2.0 * axial[i].parameter * axial_z[1]
            + axial[i].value * axial_z[2]
        )

        pressure_parameter2 += (
            phi[i].parameter2 * phi[j].value
            + 2.0 * phi[i].parameter * phi[j].parameter
            + phi[i].value * phi[j].parameter2
        )

    angular_preceding = preceding_diffusion_second_parameter_jet(
        hierarchy.h,
        angular_power,
        order,
        X,
        eta,
        phi[order - 1],
    )
    angular_parameter2 -= angular_preceding.parameter2

    axial_preceding = preceding_diffusion_second_parameter_jet(
        hierarchy.h,
        axial_power,
        order,
        X,
        eta,
        _profile_fourth_from_axial(axial[order - 1]),
    )
    axial_parameter2 -= axial_preceding.parameter2

    omega = hierarchy_owned_omega_second_parameter_jet(
        hierarchy,
        order - 1,
        X,
        eta,
    )

    parameter2 = PositiveAxisSourceJet(
        angular=angular_parameter2,
        axial=axial_parameter2,
        pressure_product=pressure_parameter2,
        omega_quotient=omega.parameter2,
    )
    return LowerSourceSecondParameterJet(
        value=first.value,
        parameter=first.parameter,
        parameter2=parameter2,
    )
