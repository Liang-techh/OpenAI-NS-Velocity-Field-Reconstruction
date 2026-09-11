"""Exact PositiveAxis lower-history source from already-constructed coefficients.

This module is a narrow executable transcription of the pinned
``PositiveAxisSystem.actualLowerSource`` / ``PositiveAxisExistence.lowerHistoryData``
bridge.  It removes one remaining arbitrary-input layer from Stage 2: once the
lower coefficient profiles are known through second ``(X, eta)`` jets, the
``BaseJet`` and ``SourceJet`` used by Eq. (5.7) are determined by the paper's
strictly-lower-order convolutions, preceding axial diffusions, and Eq. (5.6)
regular ``Omega/X`` quotient.

The lower jets themselves remain inputs.  Therefore this module is
``formal-structure`` infrastructure and must not be used to promote the
reconstruction to paper-exact before the recursive hierarchy is materialized.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import NamedTuple, Sequence
import math

from .background_positive_axis import (
    PositiveAxisBaseJet,
    PositiveAxisSourceJet,
    RadialParameterJet,
)
from .background_recurrence import RegularFluxJet, omega_over_x_eq_5_6
from .coordinates import validate_h


@dataclass(frozen=True)
class ProfileSecondJet:
    """Second ``(X, eta)`` jet of one scalar coefficient profile."""

    value: float
    radial: float
    radial2: float
    parameter: float
    radial_parameter: float
    parameter2: float

    def __post_init__(self) -> None:
        for name in (
            "value",
            "radial",
            "radial2",
            "parameter",
            "radial_parameter",
            "parameter2",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)

    def first_jet(self) -> RadialParameterJet:
        return RadialParameterJet(
            self.value, self.radial, self.radial2, self.parameter
        )

    def regular_flux_jet(self) -> RegularFluxJet:
        """Interpret the jet as the regular factor ``v`` in ``V=X v``."""

        return RegularFluxJet(
            value=self.value,
            dX=self.radial,
            dEta=self.parameter,
            dXX=self.radial2,
            dXdEta=self.radial_parameter,
            dEtaEta=self.parameter2,
        )


class PositiveAxisPointData(NamedTuple):
    base: PositiveAxisBaseJet
    source: PositiveAxisSourceJet


def _point(X: float, eta: float) -> tuple[float, float]:
    X = float(X)
    eta = float(eta)
    if not math.isfinite(X) or X < 0.0:
        raise ValueError("X must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return X, eta


def _positive_order(order: int) -> int:
    if isinstance(order, bool) or not isinstance(order, Integral) or int(order) < 1:
        raise ValueError("order must be a positive integer")
    return int(order)


def _history(
    values: Sequence[ProfileSecondJet], order: int, name: str
) -> tuple[ProfileSecondJet, ...]:
    if len(values) < order:
        raise ValueError(f"{name} must contain orders 0 through {order - 1}")
    out = tuple(values[j] for j in range(order))
    if any(not isinstance(jet, ProfileSecondJet) for jet in out):
        raise TypeError(f"{name} entries must be ProfileSecondJet instances")
    return out


def _scalars(h: float, eta: float) -> tuple[float, float, float]:
    h = validate_h(h)
    eta = float(eta)
    d = 1.0 - eta * eta
    ell = 1.0 - 2.0 * h * eta * eta
    D = 0.5 - h
    if ell <= 0.0:
        raise ValueError("ell(h,eta) must be positive")
    return D, d, ell


def _axial_value(
    h: float, power: float, X: float, eta: float, jet: ProfileSecondJet
) -> float:
    _D, d, ell = _scalars(h, eta)
    return (
        2.0 * eta * power * jet.value
        + d * jet.parameter
        - 2.0 * eta * X * jet.radial
    ) / ell


def preceding_diffusion_from_second_jet(
    h: float,
    power: float,
    order: int,
    X: float,
    eta: float,
    jet: ProfileSecondJet,
) -> float:
    """Evaluate pinned ``precedingDiffusion`` at positive order ``order``.

    The Lean definition is

    ``Z_(b-D) (Z_b F_(n-1))`` with ``b=power+lambda_(n-1)``.

    This routine differentiates the first ``Z`` analytically, so no numerical
    parameter/radial finite difference and no sampled surrogate is used in the
    production path.
    """

    h = validate_h(h)
    order = _positive_order(order)
    X, eta = _point(X, eta)
    if not isinstance(jet, ProfileSecondJet):
        raise TypeError("jet must be ProfileSecondJet")
    power = float(power)
    if not math.isfinite(power):
        raise ValueError("power must be finite")

    D, d, ell = _scalars(h, eta)
    b = power + 2.0 * (order - 1) * h

    numerator = (
        2.0 * eta * b * jet.value
        + d * jet.parameter
        - 2.0 * eta * X * jet.radial
    )
    numerator_X = (
        2.0 * eta * (b - 1.0) * jet.radial
        + d * jet.radial_parameter
        - 2.0 * eta * X * jet.radial2
    )
    numerator_eta = (
        2.0 * b * jet.value
        + 2.0 * eta * (b - 1.0) * jet.parameter
        + d * jet.parameter2
        - 2.0 * X * jet.radial
        - 2.0 * eta * X * jet.radial_parameter
    )

    first = numerator / ell
    first_X = numerator_X / ell
    first_eta = numerator_eta / ell + 4.0 * h * eta * numerator / (ell * ell)
    value = (
        2.0 * eta * (b - D) * first
        + d * first_eta
        - 2.0 * eta * X * first_X
    ) / ell
    if not math.isfinite(value):
        raise OverflowError("preceding diffusion is outside binary64 range")
    return value


def positive_axis_point_data_from_lower_history(
    h: float,
    order: int,
    X: float,
    eta: float,
    phi_jets: Sequence[ProfileSecondJet],
    axial_jets: Sequence[ProfileSecondJet],
    beta_jets: Sequence[ProfileSecondJet],
) -> PositiveAxisPointData:
    """Build exact ``BaseJet`` and ``actualLowerSource`` at one point.

    ``phi_jets``, ``axial_jets`` and ``beta_jets`` contain only orders strictly
    below ``order``.  ``beta_jets[j]`` is the regular quotient ``V_j/X`` from
    Eq. (5.2).  Consequently the previous ``Omega/X`` source is obtained from
    the already-landed exact Eq. (5.6) regular primitive rather than supplied as
    an unrelated scalar.
    """

    h = validate_h(h)
    order = _positive_order(order)
    X, eta = _point(X, eta)
    phi = _history(phi_jets, order, "phi_jets")
    axial = _history(axial_jets, order, "axial_jets")
    beta = _history(beta_jets, order, "beta_jets")

    angular_power = -1.0 - h
    axial_power = -0.5 - h

    angular = 0.0
    axial_source = 0.0
    pressure_product = 0.0
    for i in range(1, order):
        j = order - i
        lam_j = 2.0 * j * h
        angular += beta[i].value * (X * phi[j].radial + phi[j].value)
        angular += axial[i].value * _axial_value(
            h, angular_power + lam_j, X, eta, phi[j]
        )
        axial_source += beta[i].value * X * axial[j].radial
        axial_source += axial[i].value * _axial_value(
            h, axial_power + lam_j, X, eta, axial[j]
        )
        pressure_product += phi[i].value * phi[j].value

    angular -= preceding_diffusion_from_second_jet(
        h, angular_power, order, X, eta, phi[order - 1]
    )
    axial_source -= preceding_diffusion_from_second_jet(
        h, axial_power, order, X, eta, axial[order - 1]
    )

    omega_quotient = omega_over_x_eq_5_6(
        order - 1,
        X,
        eta,
        [jet.regular_flux_jet() for jet in beta],
        [jet.value for jet in axial],
        h=h,
    )

    base = PositiveAxisBaseJet(
        phi=phi[0].first_jet(),
        axial=axial[0].first_jet(),
        beta=beta[0].value,
    )
    source = PositiveAxisSourceJet(
        angular=angular,
        axial=axial_source,
        pressure_product=pressure_product,
        omega_quotient=omega_quotient,
    )
    return PositiveAxisPointData(base=base, source=source)
