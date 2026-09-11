"""Analytic eta derivative of the regular Eq. (5.6) ``Omega_k / X`` row.

The landed lower-history bridge evaluates ``Omega_k / X`` from second jets of
``beta_j = V_j / X``.  The next Lemma-5.1 Picard application needs an analytic
``partial_eta f_n``.  Differentiating the forcing in turn requires
``partial_eta(Omega_{n-1}/X)``.

This module isolates that derivative exactly.  It accepts the minimal third
mixed jet of each regular flux factor beta and ordinary second jets of U,
then differentiates every term of Eq. (5.6) analytically, including the nested
shifted axial-viscosity term.  No finite difference is used in production.

The current Section-5 hierarchy does not yet own these third mixed beta jets:
Eq. (5.2) would require one more derivative level of U to construct them.
Therefore this is fail-closed solver infrastructure, not a paper-exact
coefficient or a claim that ``partial_eta f_n`` is already hierarchy-owned.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Sequence
import math

from .background_lower_history_source import ProfileSecondJet
from .background_recurrence import RegularFluxJet, omega_over_x_eq_5_6
from .coordinates import validate_h


@dataclass(frozen=True)
class RegularFluxThirdMixedJet:
    """Minimal beta-jet needed for ``partial_eta(Omega_k/X)``.

    The first six fields are the existing second jet of ``beta=V/X``.  The
    final three are exactly the derivatives introduced when eta differentiates
    the radial-viscosity and nested ``Z^[2]`` terms of Eq. (5.6).
    """

    value: float
    radial: float
    radial2: float
    parameter: float
    radial_parameter: float
    parameter2: float
    radial2_parameter: float
    radial_parameter2: float
    parameter3: float

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
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)

    def second(self) -> RegularFluxJet:
        return RegularFluxJet(
            value=self.value,
            dX=self.radial,
            dEta=self.parameter,
            dXX=self.radial2,
            dXdEta=self.radial_parameter,
            dEtaEta=self.parameter2,
        )


@dataclass(frozen=True)
class OmegaParameterJet:
    """Value and analytic eta derivative of ``Omega_k/X``."""

    value: float
    parameter: float

    def __post_init__(self) -> None:
        for name in ("value", "parameter"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)


def _order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError("order must be a nonnegative integer")
    return int(value)


def _point(X: float, eta: float) -> tuple[float, float]:
    X = float(X)
    eta = float(eta)
    if not math.isfinite(X) or X < 0.0:
        raise ValueError("X must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return X, eta


def _geometry(h: float, eta: float) -> tuple[float, float, float]:
    h = validate_h(h)
    d = 1.0 - eta * eta
    ell = 1.0 - 2.0 * h * eta * eta
    D = 0.5 - h
    if ell <= 0.0:
        raise ValueError("similarity factor L must be positive")
    return D, d, ell


def _history(
    values: Sequence[object], order: int, expected_type: type, name: str
) -> tuple[object, ...]:
    if len(values) < order + 1:
        raise ValueError(f"{name} must contain orders 0 through {order}")
    out = tuple(values[j] for j in range(order + 1))
    if any(not isinstance(value, expected_type) for value in out):
        raise TypeError(f"{name} entries must be {expected_type.__name__} instances")
    return out


def _quotient_eta(
    numerator: float,
    numerator_eta: float,
    *,
    h: float,
    eta: float,
    ell: float,
) -> float:
    return numerator_eta / ell + 4.0 * h * eta * numerator / (ell * ell)


def _regular_T_parameter(
    jet: RegularFluxThirdMixedJet,
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
    return _quotient_eta(
        numerator, numerator_eta, h=h, eta=eta, ell=ell
    )


def _regular_Z_value_parameter(
    jet: RegularFluxThirdMixedJet,
    X: float,
    eta: float,
    *,
    b: float,
    h: float,
) -> tuple[float, float]:
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
    return (
        numerator / ell,
        _quotient_eta(
            numerator, numerator_eta, h=h, eta=eta, ell=ell
        ),
    )


def _regular_Z2_parameter(
    jet: RegularFluxThirdMixedJet,
    X: float,
    eta: float,
    *,
    b: float,
    h: float,
) -> float:
    """Return ``partial_eta[Z_(b-D) Z_b(X beta) / X]`` analytically."""

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

    first = numerator / ell
    first_X = numerator_X / ell
    first_eta = _quotient_eta(
        numerator, numerator_eta, h=h, eta=eta, ell=ell
    )
    first_X_eta = _quotient_eta(
        numerator_X, numerator_X_eta, h=h, eta=eta, ell=ell
    )
    first_eta2 = (
        numerator_eta2 / ell
        + 8.0 * h * eta * numerator_eta / (ell * ell)
        + 4.0 * h * numerator / (ell * ell)
        + 32.0 * h * h * eta * eta * numerator / (ell * ell * ell)
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
    return _quotient_eta(
        outer_numerator,
        outer_numerator_eta,
        h=h,
        eta=eta,
        ell=ell,
    )


def omega_over_x_parameter_jet_eq_5_6(
    order: int,
    X: float,
    eta: float,
    radial_flux_jets: Sequence[RegularFluxThirdMixedJet],
    axial_jets: Sequence[ProfileSecondJet],
    *,
    h: float,
) -> OmegaParameterJet:
    """Return ``(Omega_order/X, partial_eta(Omega_order/X))`` from Eq. (5.6).

    The value is delegated to the already-landed regular Eq. (5.6) primitive.
    The parameter derivative is a direct product/quotient-rule differentiation
    of the same five rows.  In particular, the shifted axial-viscosity term is
    differentiated through both nested ``Z`` operators rather than sampled at
    neighboring eta values.
    """

    order = _order(order)
    X, eta = _point(X, eta)
    h = validate_h(h)
    beta = _history(
        radial_flux_jets, order, RegularFluxThirdMixedJet, "radial_flux_jets"
    )
    axial = _history(axial_jets, order, ProfileSecondJet, "axial_jets")
    lam = lambda n: 2.0 * n * h

    value = omega_over_x_eq_5_6(
        order,
        X,
        eta,
        [jet.second() for jet in beta],
        [jet.value for jet in axial],
        h=h,
    )

    parameter = _regular_T_parameter(
        beta[order], X, eta, b=lam(order), h=h
    )

    for i in range(order + 1):
        j = order - i
        right = 0.5 * beta[j].value + X * beta[j].radial
        right_eta = (
            0.5 * beta[j].parameter + X * beta[j].radial_parameter
        )
        parameter += beta[i].parameter * right + beta[i].value * right_eta

        z_value, z_parameter = _regular_Z_value_parameter(
            beta[j], X, eta, b=lam(j), h=h
        )
        parameter += (
            axial[i].parameter * z_value + axial[i].value * z_parameter
        )

    parameter += -2.0 * (
        2.0 * beta[order].radial_parameter
        + X * beta[order].radial2_parameter
    )
    if order >= 1:
        parameter -= _regular_Z2_parameter(
            beta[order - 1], X, eta, b=lam(order - 1), h=h
        )

    if not math.isfinite(parameter):
        raise OverflowError("Eq. (5.6) omega parameter jet is outside binary64 range")
    return OmegaParameterJet(value=value, parameter=parameter)
