"""Analytic second eta derivative of the regular Eq. (5.6) ``Omega_k / X`` row.

The landed :mod:`background_omega_parameter_jet` module owns the exact value and
first eta derivative of the regular quotient ``Omega_k/X``.  A second eta
forcing jet needs one more derivative of that same row.  Differentiating Eq.
(5.6) again requires a fourth mixed jet of ``beta_j = V_j/X`` only in the
eta-bearing entries that can actually occur:

``beta_XXetaeta``, ``beta_Xetaetaeta`` and ``beta_etaetaetaeta``.

This module makes that debt explicit and evaluates ``partial_eta^2(Omega_k/X)``
analytically.  Value and first derivative are delegated to the already-landed
parameter-jet implementation so this file cannot silently fork their semantics.
No finite difference is used in production.

The repaired Section-5 hierarchy does not yet own these fourth-mixed beta jets;
Eq. (5.2) would require a stronger U layer.  This is therefore fail-closed
``formal-structure`` infrastructure, not a paper-exact coefficient claim.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Sequence
import math

from .background_lower_history_source import ProfileSecondJet
from .background_omega_parameter_jet import (
    OmegaParameterJet,
    RegularFluxThirdMixedJet,
    omega_over_x_parameter_jet_eq_5_6,
)
from .coordinates import validate_h


@dataclass(frozen=True)
class RegularFluxFourthMixedJet:
    """Third mixed beta jet plus the eta-bearing total-order-four entries."""

    value: float
    radial: float
    radial2: float
    parameter: float
    radial_parameter: float
    parameter2: float
    radial2_parameter: float
    radial_parameter2: float
    parameter3: float
    radial2_parameter2: float
    radial_parameter3: float
    parameter4: float

    def __post_init__(self) -> None:
        for name in (
            "value",
            "radial",
            "radial2",
            "parameter",
            "radial_parameter",
            "parameter2",
            "radial2_parameter",
            "radial_parameter2",
            "parameter3",
            "radial2_parameter2",
            "radial_parameter3",
            "parameter4",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)

    def third(self) -> RegularFluxThirdMixedJet:
        return RegularFluxThirdMixedJet(
            value=self.value,
            radial=self.radial,
            radial2=self.radial2,
            parameter=self.parameter,
            radial_parameter=self.radial_parameter,
            parameter2=self.parameter2,
            radial2_parameter=self.radial2_parameter,
            radial_parameter2=self.radial_parameter2,
            parameter3=self.parameter3,
        )


@dataclass(frozen=True)
class OmegaSecondParameterJet:
    """Value and first two analytic eta derivatives of ``Omega_k/X``."""

    value: float
    parameter: float
    parameter2: float

    def __post_init__(self) -> None:
        for name in ("value", "parameter", "parameter2"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)


def _order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError("order must be a nonnegative integer")
    return int(value)


def _history(
    values: Sequence[object], order: int, expected_type: type, name: str
) -> tuple[object, ...]:
    if len(values) < order + 1:
        raise ValueError(f"{name} must contain orders 0 through {order}")
    out = tuple(values[j] for j in range(order + 1))
    if any(not isinstance(value, expected_type) for value in out):
        raise TypeError(f"{name} entries must be {expected_type.__name__} instances")
    return out


def _geometry(h: float, eta: float) -> tuple[float, float, float]:
    h = validate_h(h)
    D = 0.5 - h
    d = 1.0 - eta * eta
    ell = 1.0 - 2.0 * h * eta * eta
    if ell <= 0.0:
        raise ValueError("similarity factor L must be positive")
    return D, d, ell


def _quotient_parameter(
    numerator: float,
    numerator_eta: float,
    *,
    h: float,
    eta: float,
    ell: float,
) -> float:
    return numerator_eta / ell + 4.0 * h * eta * numerator / (ell * ell)


def _quotient_parameter2(
    numerator: float,
    numerator_eta: float,
    numerator_eta2: float,
    *,
    h: float,
    eta: float,
    ell: float,
) -> float:
    return (
        numerator_eta2 / ell
        + 8.0 * h * eta * numerator_eta / (ell * ell)
        + 4.0 * h * numerator / (ell * ell)
        + 32.0 * h * h * eta * eta * numerator / (ell * ell * ell)
    )


def _quotient_parameter3(
    numerator: float,
    numerator_eta: float,
    numerator_eta2: float,
    numerator_eta3: float,
    *,
    h: float,
    eta: float,
    ell: float,
) -> float:
    return (
        numerator_eta3 / ell
        + 12.0 * h * eta * numerator_eta2 / (ell * ell)
        + 12.0 * h * numerator_eta / (ell * ell)
        + 96.0 * h * h * eta * eta * numerator_eta / (ell * ell * ell)
        + 96.0 * h * h * eta * numerator / (ell * ell * ell)
        + 384.0 * h * h * h * eta * eta * eta * numerator / (ell**4)
    )


def _regular_T_parameter2(
    jet: RegularFluxFourthMixedJet,
    X: float,
    eta: float,
    *,
    b: float,
    h: float,
) -> float:
    D, _d, ell = _geometry(h, eta)
    numerator = (
        (1.0 - b) * jet.value
        + D * eta * jet.parameter
        + X * jet.radial
    )
    numerator_eta = (
        (1.0 - b) * jet.parameter
        + D * (jet.parameter + eta * jet.parameter2)
        + X * jet.radial_parameter
    )
    numerator_eta2 = (
        (1.0 - b) * jet.parameter2
        + D * (2.0 * jet.parameter2 + eta * jet.parameter3)
        + X * jet.radial_parameter2
    )
    return _quotient_parameter2(
        numerator,
        numerator_eta,
        numerator_eta2,
        h=h,
        eta=eta,
        ell=ell,
    )


def _regular_Z_value_parameter2(
    jet: RegularFluxFourthMixedJet,
    X: float,
    eta: float,
    *,
    b: float,
    h: float,
) -> tuple[float, float, float]:
    _D, d, ell = _geometry(h, eta)
    numerator = (
        2.0 * eta * (b - 1.0) * jet.value
        + d * jet.parameter
        - 2.0 * eta * X * jet.radial
    )
    numerator_eta = (
        2.0 * (b - 1.0) * jet.value
        + 2.0 * eta * (b - 2.0) * jet.parameter
        + d * jet.parameter2
        - 2.0 * X * jet.radial
        - 2.0 * eta * X * jet.radial_parameter
    )
    numerator_eta2 = (
        (4.0 * b - 6.0) * jet.parameter
        + 2.0 * eta * (b - 3.0) * jet.parameter2
        + d * jet.parameter3
        - 4.0 * X * jet.radial_parameter
        - 2.0 * eta * X * jet.radial_parameter2
    )
    return (
        numerator / ell,
        _quotient_parameter(
            numerator, numerator_eta, h=h, eta=eta, ell=ell
        ),
        _quotient_parameter2(
            numerator,
            numerator_eta,
            numerator_eta2,
            h=h,
            eta=eta,
            ell=ell,
        ),
    )


def _regular_Z2_parameter2(
    jet: RegularFluxFourthMixedJet,
    X: float,
    eta: float,
    *,
    b: float,
    h: float,
) -> float:
    """Return ``partial_eta^2[Z_(b-D) Z_b(X beta) / X]`` analytically."""

    D, d, ell = _geometry(h, eta)
    numerator = (
        2.0 * eta * (b - 1.0) * jet.value
        + d * jet.parameter
        - 2.0 * eta * X * jet.radial
    )
    numerator_X = (
        2.0 * eta * (b - 2.0) * jet.radial
        + d * jet.radial_parameter
        - 2.0 * eta * X * jet.radial2
    )
    numerator_eta = (
        2.0 * (b - 1.0) * jet.value
        + 2.0 * eta * (b - 2.0) * jet.parameter
        + d * jet.parameter2
        - 2.0 * X * jet.radial
        - 2.0 * eta * X * jet.radial_parameter
    )
    numerator_X_eta = (
        2.0 * (b - 2.0) * jet.radial
        + 2.0 * eta * (b - 3.0) * jet.radial_parameter
        + d * jet.radial_parameter2
        - 2.0 * X * jet.radial2
        - 2.0 * eta * X * jet.radial2_parameter
    )
    numerator_eta2 = (
        (4.0 * b - 6.0) * jet.parameter
        + 2.0 * eta * (b - 3.0) * jet.parameter2
        + d * jet.parameter3
        - 4.0 * X * jet.radial_parameter
        - 2.0 * eta * X * jet.radial_parameter2
    )
    numerator_X_eta2 = (
        (4.0 * b - 10.0) * jet.radial_parameter
        + 2.0 * eta * (b - 4.0) * jet.radial_parameter2
        + d * jet.radial_parameter3
        - 4.0 * X * jet.radial2_parameter
        - 2.0 * eta * X * jet.radial2_parameter2
    )
    numerator_eta3 = (
        6.0 * (b - 2.0) * jet.parameter2
        + 2.0 * eta * (b - 4.0) * jet.parameter3
        + d * jet.parameter4
        - 6.0 * X * jet.radial_parameter2
        - 2.0 * eta * X * jet.radial_parameter3
    )

    first = numerator / ell
    first_X = numerator_X / ell
    first_eta = _quotient_parameter(
        numerator, numerator_eta, h=h, eta=eta, ell=ell
    )
    first_X_eta = _quotient_parameter(
        numerator_X, numerator_X_eta, h=h, eta=eta, ell=ell
    )
    first_eta2 = _quotient_parameter2(
        numerator,
        numerator_eta,
        numerator_eta2,
        h=h,
        eta=eta,
        ell=ell,
    )
    first_X_eta2 = _quotient_parameter2(
        numerator_X,
        numerator_X_eta,
        numerator_X_eta2,
        h=h,
        eta=eta,
        ell=ell,
    )
    first_eta3 = _quotient_parameter3(
        numerator,
        numerator_eta,
        numerator_eta2,
        numerator_eta3,
        h=h,
        eta=eta,
        ell=ell,
    )

    outer_b = b - D - 1.0
    outer_numerator = (
        2.0 * eta * outer_b * first
        + d * first_eta
        - 2.0 * eta * X * first_X
    )
    outer_numerator_eta = (
        2.0 * outer_b * first
        + 2.0 * eta * (outer_b - 1.0) * first_eta
        + d * first_eta2
        - 2.0 * X * first_X
        - 2.0 * eta * X * first_X_eta
    )
    outer_numerator_eta2 = (
        (4.0 * outer_b - 2.0) * first_eta
        + 2.0 * eta * (outer_b - 2.0) * first_eta2
        + d * first_eta3
        - 4.0 * X * first_X_eta
        - 2.0 * eta * X * first_X_eta2
    )
    return _quotient_parameter2(
        outer_numerator,
        outer_numerator_eta,
        outer_numerator_eta2,
        h=h,
        eta=eta,
        ell=ell,
    )


def omega_over_x_second_parameter_jet_eq_5_6(
    order: int,
    X: float,
    eta: float,
    radial_flux_jets: Sequence[RegularFluxFourthMixedJet],
    axial_jets: Sequence[ProfileSecondJet],
    *,
    h: float,
) -> OmegaSecondParameterJet:
    """Return ``(Omega/X, partial_eta Omega/X, partial_eta^2 Omega/X)``.

    The existing first-parameter implementation remains authoritative for the
    first two components.  The new second derivative is the direct analytic
    second eta derivative of the same five Eq. (5.6) rows.
    """

    order = _order(order)
    beta = _history(
        radial_flux_jets, order, RegularFluxFourthMixedJet, "radial_flux_jets"
    )
    axial = _history(axial_jets, order, ProfileSecondJet, "axial_jets")
    first_result: OmegaParameterJet = omega_over_x_parameter_jet_eq_5_6(
        order,
        X,
        eta,
        [jet.third() for jet in beta],
        axial,
        h=h,
    )

    X = float(X)
    eta = float(eta)
    h = validate_h(h)
    lam = lambda n: 2.0 * n * h

    parameter2 = _regular_T_parameter2(
        beta[order], X, eta, b=lam(order), h=h
    )

    for i in range(order + 1):
        j = order - i
        right = 0.5 * beta[j].value + X * beta[j].radial
        right_eta = 0.5 * beta[j].parameter + X * beta[j].radial_parameter
        right_eta2 = (
            0.5 * beta[j].parameter2 + X * beta[j].radial_parameter2
        )
        parameter2 += (
            beta[i].parameter2 * right
            + 2.0 * beta[i].parameter * right_eta
            + beta[i].value * right_eta2
        )

        z_value, z_parameter, z_parameter2 = _regular_Z_value_parameter2(
            beta[j], X, eta, b=lam(j), h=h
        )
        parameter2 += (
            axial[i].parameter2 * z_value
            + 2.0 * axial[i].parameter * z_parameter
            + axial[i].value * z_parameter2
        )

    parameter2 += -2.0 * (
        2.0 * beta[order].radial_parameter2
        + X * beta[order].radial2_parameter2
    )
    if order >= 1:
        parameter2 -= _regular_Z2_parameter2(
            beta[order - 1], X, eta, b=lam(order - 1), h=h
        )

    if not math.isfinite(parameter2):
        raise OverflowError(
            "Eq. (5.6) omega second parameter jet is outside binary64 range"
        )
    return OmegaSecondParameterJet(
        value=first_result.value,
        parameter=first_result.parameter,
        parameter2=parameter2,
    )
