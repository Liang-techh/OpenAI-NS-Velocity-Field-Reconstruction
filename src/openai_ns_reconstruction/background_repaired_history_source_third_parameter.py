"""Hierarchy-owned analytic third eta jet of Section 5 ``actualLowerSource``.

The landed second-parameter bridge owns the value, first eta derivative, and
second eta derivative of the complete strict-lower PositiveAxis source.  The
stronger hierarchy now also owns exactly the fifth/sixth mixed coefficient data
needed to differentiate the same source once more:

* repaired fifth-mixed ``phi`` jets;
* repaired sixth-mixed ``U`` jets and their fifth-mixed projections;
* fifth-mixed ``beta=V/X`` derived only through analytic Eq. (5.2);
* hierarchy-owned angular and axial third-eta preceding-diffusion rows; and
* hierarchy-owned third-eta ``Omega/X`` from Eq. (5.6).

This module composes those landed pieces.  Lower source rows remain delegated to
the existing second-parameter implementation; only ``partial_eta^3`` is new.
Production contains no finite differences, fitted coefficients, generic cutoff,
caller-supplied ``SourceJet`` derivative table, or sampled ``V/X`` division.

Genuine leading strong data and unrepaired compact-repair inputs remain upstream
contracts, so this is Stage-2 ``formal-structure`` infrastructure rather than a
paper-exact recursive coefficient or an all-order convergence claim.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import math

from .background_moment_repair_phi_fifth_mixed_jets import ProfileFifthMixedJet
from .background_positive_axis import PositiveAxisSourceJet
from .background_regular_flux_fourth_mixed_jets import AxialFifthMixedJet
from .background_preceding_diffusion_third_parameter_jet import (
    hierarchy_owned_angular_preceding_diffusion_third_parameter_jet,
)
from .background_repaired_history_axial_preceding_diffusion_third_parameter import (
    hierarchy_owned_axial_preceding_diffusion_third_parameter_jet,
)
from .background_repaired_history_omega_third_parameter import (
    hierarchy_owned_omega_third_parameter_jet,
)
from .background_repaired_history_sixth_mixed import (
    Section5LowerHistorySixthMixedHierarchy,
)
from .background_repaired_history_source_second_parameter import (
    LowerSourceSecondParameterJet,
    hierarchy_owned_lower_source_second_parameter_jet,
)


@dataclass(frozen=True)
class LowerSourceThirdParameterJet:
    """Value and first three analytic eta derivatives of ``actualLowerSource``."""

    value: PositiveAxisSourceJet
    parameter: PositiveAxisSourceJet
    parameter2: PositiveAxisSourceJet
    parameter3: PositiveAxisSourceJet

    def __post_init__(self) -> None:
        for name in ("value", "parameter", "parameter2", "parameter3"):
            if not isinstance(getattr(self, name), PositiveAxisSourceJet):
                raise TypeError(f"{name} must be a PositiveAxisSourceJet")

    def second(self) -> LowerSourceSecondParameterJet:
        """Project exactly onto the landed lower source second-eta jet."""
        return LowerSourceSecondParameterJet(
            value=self.value,
            parameter=self.parameter,
            parameter2=self.parameter2,
        )


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


def _profile_fifth_from_axial(jet: AxialFifthMixedJet) -> ProfileFifthMixedJet:
    if not isinstance(jet, AxialFifthMixedJet):
        raise TypeError("jet must be an AxialFifthMixedJet")
    return ProfileFifthMixedJet(
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
        radial2_parameter3=jet.radial2_parameter3,
        radial_parameter4=jet.radial_parameter4,
        parameter5=jet.parameter5,
    )


def _axial_value_parameter3(
    h: float,
    power: float,
    X: float,
    eta: float,
    jet: ProfileFifthMixedJet,
) -> tuple[float, float, float, float]:
    """Return eta derivatives 0..3 of the displayed ``Z_power`` operator."""

    if not isinstance(jet, ProfileFifthMixedJet):
        raise TypeError("jet must be a ProfileFifthMixedJet")
    d = 1.0 - eta * eta
    ell = 1.0 - 2.0 * h * eta * eta
    if ell <= 0.0:
        raise ValueError("similarity factor L must be positive")

    n0 = (
        2.0 * eta * power * jet.value
        + d * jet.parameter
        - 2.0 * eta * X * jet.radial
    )
    n1 = (
        2.0 * power * jet.value
        + 2.0 * eta * (power - 1.0) * jet.parameter
        + d * jet.parameter2
        - 2.0 * X * jet.radial
        - 2.0 * eta * X * jet.radial_parameter
    )
    n2 = (
        (4.0 * power - 2.0) * jet.parameter
        + 2.0 * eta * (power - 2.0) * jet.parameter2
        + d * jet.parameter3
        - 4.0 * X * jet.radial_parameter
        - 2.0 * eta * X * jet.radial_parameter2
    )
    n3 = (
        6.0 * (power - 1.0) * jet.parameter2
        + 2.0 * eta * (power - 3.0) * jet.parameter3
        + d * jet.parameter4
        - 6.0 * X * jet.radial_parameter2
        - 2.0 * eta * X * jet.radial_parameter3
    )

    s0 = 1.0 / ell
    s1 = 4.0 * h * eta / ell**2
    s2 = 4.0 * h / ell**2 + 32.0 * h * h * eta * eta / ell**3
    s3 = 96.0 * h * h * eta / ell**3 + 384.0 * h**3 * eta**3 / ell**4
    out = (
        n0 * s0,
        n1 * s0 + n0 * s1,
        n2 * s0 + 2.0 * n1 * s1 + n0 * s2,
        n3 * s0 + 3.0 * n2 * s1 + 3.0 * n1 * s2 + n0 * s3,
    )
    if not all(math.isfinite(value) for value in out):
        raise OverflowError("axial third parameter jet is outside binary64 range")
    return out


def _product_parameter3(
    a0: float,
    a1: float,
    a2: float,
    a3: float,
    b0: float,
    b1: float,
    b2: float,
    b3: float,
) -> float:
    return math.fsum((a3 * b0, 3.0 * a2 * b1, 3.0 * a1 * b2, a0 * b3))


def hierarchy_owned_lower_source_third_parameter_jet(
    hierarchy: Section5LowerHistorySixthMixedHierarchy,
    order: int,
    X: float,
    eta: float,
) -> LowerSourceThirdParameterJet:
    """Return value through ``partial_eta^3 actualLowerSource`` from one hierarchy.

    ``order`` is the positive coefficient order being solved.  Only strict-lower
    coefficient orders ``0,...,order-1`` are queried.  The function first calls
    the landed second-eta bridge, so all lower rows retain exactly their existing
    semantics.  It then differentiates only the strict-lower convolution rows by
    the ordinary third-order product rule and obtains all three nonlocal/source
    seams from their hierarchy-owned analytic bridges.
    """

    if not isinstance(hierarchy, Section5LowerHistorySixthMixedHierarchy):
        raise TypeError("hierarchy must be a Section5LowerHistorySixthMixedHierarchy")
    order = _positive_order(order)
    X, eta = _point(X, eta)
    if len(hierarchy.sources) < order:
        raise ValueError(
            "requested source order is missing one or more strict-lower coefficients"
        )

    lower = hierarchy_owned_lower_source_second_parameter_jet(
        hierarchy,
        order,
        X,
        eta,
    )

    phi = [hierarchy.phi_fifth_mixed_jet(j, X, eta) for j in range(order)]
    axial = [hierarchy.axial_fifth_mixed_jet(j, X, eta) for j in range(order)]
    beta = [hierarchy.beta_fifth_mixed_jet(j, X, eta) for j in range(order)]

    angular_power = -1.0 - hierarchy.h
    axial_power = -0.5 - hierarchy.h
    angular_parameter3 = 0.0
    axial_parameter3 = 0.0
    pressure_parameter3 = 0.0

    for i in range(1, order):
        j = order - i
        lam_j = 2.0 * j * hierarchy.h

        angular_gradient = (
            X * phi[j].radial + phi[j].value,
            X * phi[j].radial_parameter + phi[j].parameter,
            X * phi[j].radial_parameter2 + phi[j].parameter2,
            X * phi[j].radial_parameter3 + phi[j].parameter3,
        )
        angular_parameter3 += _product_parameter3(
            beta[i].value,
            beta[i].parameter,
            beta[i].parameter2,
            beta[i].parameter3,
            *angular_gradient,
        )

        angular_z = _axial_value_parameter3(
            hierarchy.h,
            angular_power + lam_j,
            X,
            eta,
            phi[j],
        )
        angular_parameter3 += _product_parameter3(
            axial[i].value,
            axial[i].parameter,
            axial[i].parameter2,
            axial[i].parameter3,
            *angular_z,
        )

        axial_parameter3 += _product_parameter3(
            beta[i].value,
            beta[i].parameter,
            beta[i].parameter2,
            beta[i].parameter3,
            X * axial[j].radial,
            X * axial[j].radial_parameter,
            X * axial[j].radial_parameter2,
            X * axial[j].radial_parameter3,
        )
        axial_z = _axial_value_parameter3(
            hierarchy.h,
            axial_power + lam_j,
            X,
            eta,
            _profile_fifth_from_axial(axial[j]),
        )
        axial_parameter3 += _product_parameter3(
            axial[i].value,
            axial[i].parameter,
            axial[i].parameter2,
            axial[i].parameter3,
            *axial_z,
        )

        pressure_parameter3 += _product_parameter3(
            phi[i].value,
            phi[i].parameter,
            phi[i].parameter2,
            phi[i].parameter3,
            phi[j].value,
            phi[j].parameter,
            phi[j].parameter2,
            phi[j].parameter3,
        )

    angular_preceding = hierarchy_owned_angular_preceding_diffusion_third_parameter_jet(
        hierarchy,
        order,
        X,
        eta,
    )
    axial_preceding = hierarchy_owned_axial_preceding_diffusion_third_parameter_jet(
        hierarchy,
        order,
        X,
        eta,
    )
    omega = hierarchy_owned_omega_third_parameter_jet(
        hierarchy,
        order - 1,
        X,
        eta,
    )

    parameter3 = PositiveAxisSourceJet(
        angular=angular_parameter3 - angular_preceding.parameter3,
        axial=axial_parameter3 - axial_preceding.parameter3,
        pressure_product=pressure_parameter3,
        omega_quotient=omega.parameter3,
    )
    return LowerSourceThirdParameterJet(
        value=lower.value,
        parameter=lower.parameter,
        parameter2=lower.parameter2,
        parameter3=parameter3,
    )
