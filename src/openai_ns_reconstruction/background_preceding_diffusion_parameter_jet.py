"""Analytic eta jet for the Section 5 preceding-diffusion source term.

The exact strict-lower PositiveAxis source contains

    precedingDiffusion = Z_(b-D) (Z_b F_(n-1)).

The landed value path in :mod:`background_lower_history_source` evaluates this
quantity from a second ``(X, eta)`` jet.  Differentiating it once in ``eta``
requires exactly three additional mixed derivatives:
``F_XXeta``, ``F_Xetaeta`` and ``F_etaetaeta``.  This module makes that
requirement explicit and evaluates both the value and its analytic eta
derivative without finite differences in production.

The upstream coefficient jet is still an input.  In particular, the repaired
``phi`` hierarchy does not yet own this stronger derivative layer, so this is
``formal-structure`` solver infrastructure rather than a paper-exact source or
a complete ``partial_eta f_n`` implementation.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import math

from .background_lower_history_source import (
    ProfileSecondJet,
    preceding_diffusion_from_second_jet,
)
from .coordinates import validate_h


@dataclass(frozen=True)
class ProfileThirdMixedJet:
    """Second scalar jet plus the mixed total-order-three eta derivatives."""

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

    def second(self) -> ProfileSecondJet:
        return ProfileSecondJet(
            value=self.value,
            radial=self.radial,
            radial2=self.radial2,
            parameter=self.parameter,
            radial_parameter=self.radial_parameter,
            parameter2=self.parameter2,
        )


@dataclass(frozen=True)
class PrecedingDiffusionParameterJet:
    """Value and analytic ``eta`` derivative of ``precedingDiffusion``."""

    value: float
    parameter: float

    def __post_init__(self) -> None:
        for name in ("value", "parameter"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)


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


def preceding_diffusion_parameter_jet(
    h: float,
    power: float,
    order: int,
    X: float,
    eta: float,
    jet: ProfileThirdMixedJet,
) -> PrecedingDiffusionParameterJet:
    """Return ``(precedingDiffusion, partial_eta precedingDiffusion)``.

    This analytically differentiates the exact landed expression
    ``Z_(b-D)(Z_b F)`` with ``b=power+lambda_(order-1)``.  The value is also
    cross-linked to the existing second-jet implementation so future changes to
    the two paths cannot silently disagree.
    """

    h = validate_h(h)
    order = _positive_order(order)
    X, eta = _point(X, eta)
    if not isinstance(jet, ProfileThirdMixedJet):
        raise TypeError("jet must be ProfileThirdMixedJet")
    power = float(power)
    if not math.isfinite(power):
        raise ValueError("power must be finite")

    D = 0.5 - h
    d = 1.0 - eta * eta
    ell = 1.0 - 2.0 * h * eta * eta
    if ell <= 0.0:
        raise ValueError("ell(h,eta) must be positive")
    b = power + 2.0 * (order - 1) * h

    N = (
        2.0 * eta * b * jet.value
        + d * jet.parameter
        - 2.0 * eta * X * jet.radial
    )
    N_X = (
        2.0 * eta * (b - 1.0) * jet.radial
        + d * jet.radial_parameter
        - 2.0 * eta * X * jet.radial2
    )
    N_eta = (
        2.0 * b * jet.value
        + 2.0 * eta * (b - 1.0) * jet.parameter
        + d * jet.parameter2
        - 2.0 * X * jet.radial
        - 2.0 * eta * X * jet.radial_parameter
    )
    N_X_eta = (
        2.0 * (b - 1.0) * jet.radial
        + 2.0 * eta * (b - 2.0) * jet.radial_parameter
        + d * jet.radial_parameter2
        - 2.0 * X * jet.radial2
        - 2.0 * eta * X * jet.radial2_parameter
    )
    N_eta_eta = (
        (4.0 * b - 2.0) * jet.parameter
        + 2.0 * eta * (b - 2.0) * jet.parameter2
        + d * jet.parameter3
        - 4.0 * X * jet.radial_parameter
        - 2.0 * eta * X * jet.radial_parameter2
    )

    first = N / ell
    first_X = N_X / ell
    first_eta = N_eta / ell + 4.0 * h * eta * N / (ell * ell)
    first_X_eta = N_X_eta / ell + 4.0 * h * eta * N_X / (ell * ell)
    first_eta_eta = (
        N_eta_eta / ell
        + 8.0 * h * eta * N_eta / (ell * ell)
        + 4.0 * h * N / (ell * ell)
        + 32.0 * h * h * eta * eta * N / (ell * ell * ell)
    )

    numerator = (
        2.0 * eta * (b - D) * first
        + d * first_eta
        - 2.0 * eta * X * first_X
    )
    numerator_eta = (
        2.0 * (b - D) * first
        + 2.0 * eta * (b - D) * first_eta
        - 2.0 * eta * first_eta
        + d * first_eta_eta
        - 2.0 * X * first_X
        - 2.0 * eta * X * first_X_eta
    )
    parameter = (
        numerator_eta / ell
        + 4.0 * h * eta * numerator / (ell * ell)
    )

    value = preceding_diffusion_from_second_jet(
        h, power, order, X, eta, jet.second()
    )
    if not math.isfinite(parameter):
        raise OverflowError(
            "preceding diffusion parameter derivative is outside binary64 range"
        )
    return PrecedingDiffusionParameterJet(value=value, parameter=parameter)
