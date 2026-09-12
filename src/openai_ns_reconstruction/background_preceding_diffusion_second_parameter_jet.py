"""Analytic second eta jet for the Section 5 preceding-diffusion term.

The strict-lower PositiveAxis source contains

    precedingDiffusion = Z_(b-D) (Z_b F_(n-1)).

The first eta jet is already implemented in
:mod:`background_preceding_diffusion_parameter_jet`.  Differentiating once
more requires exactly three additional fourth-order mixed entries beyond its
``ProfileThirdMixedJet``: ``F_XXetaeta``, ``F_Xetaetaeta`` and
``F_etaetaetaeta``.  This module exposes that stronger contract and evaluates
the second eta derivative analytically.

The stronger profile jet remains upstream data.  This module therefore adds
``formal-structure`` solver infrastructure only; it does not claim a
hierarchy-owned second forcing derivative or paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .background_preceding_diffusion_parameter_jet import (
    ProfileThirdMixedJet,
    preceding_diffusion_parameter_jet,
)
from .coordinates import validate_h


@dataclass(frozen=True)
class ProfileFourthMixedJet:
    """Third mixed scalar jet plus the eta-bearing total-order-four entries."""

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

    def third(self) -> ProfileThirdMixedJet:
        return ProfileThirdMixedJet(
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
class PrecedingDiffusionSecondParameterJet:
    """Value, first eta derivative and second eta derivative."""

    value: float
    parameter: float
    parameter2: float

    def __post_init__(self) -> None:
        for name in ("value", "parameter", "parameter2"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)


def preceding_diffusion_second_parameter_jet(
    h: float,
    power: float,
    order: int,
    X: float,
    eta: float,
    jet: ProfileFourthMixedJet,
) -> PrecedingDiffusionSecondParameterJet:
    """Return the first two analytic eta derivatives of preceding diffusion.

    The landed first-parameter implementation owns the value and first
    derivative.  This routine projects the stronger jet into that path and
    computes only the new second derivative analytically, preventing a silent
    fork of the existing value/first-derivative semantics.
    """

    if not isinstance(jet, ProfileFourthMixedJet):
        raise TypeError("jet must be ProfileFourthMixedJet")

    first_result = preceding_diffusion_parameter_jet(
        h, power, order, X, eta, jet.third()
    )
    h = validate_h(h)
    power = float(power)
    order = int(order)
    X = float(X)
    eta = float(eta)

    D = 0.5 - h
    d = 1.0 - eta * eta
    ell = 1.0 - 2.0 * h * eta * eta
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
    N_X_eta_eta = (
        (4.0 * b - 6.0) * jet.radial_parameter
        + 2.0 * eta * (b - 3.0) * jet.radial_parameter2
        + d * jet.radial_parameter3
        - 4.0 * X * jet.radial2_parameter
        - 2.0 * eta * X * jet.radial2_parameter2
    )
    N_eta_eta_eta = (
        6.0 * (b - 1.0) * jet.parameter2
        + 2.0 * eta * (b - 3.0) * jet.parameter3
        + d * jet.parameter4
        - 6.0 * X * jet.radial_parameter2
        - 2.0 * eta * X * jet.radial_parameter3
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
    first_X_eta_eta = (
        N_X_eta_eta / ell
        + 8.0 * h * eta * N_X_eta / (ell * ell)
        + 4.0 * h * N_X / (ell * ell)
        + 32.0 * h * h * eta * eta * N_X / (ell * ell * ell)
    )
    first_eta_eta_eta = (
        N_eta_eta_eta / ell
        + 12.0 * h * eta * N_eta_eta / (ell * ell)
        + 12.0 * h * N_eta / (ell * ell)
        + 96.0 * h * h * eta * eta * N_eta / (ell * ell * ell)
        + 96.0 * h * h * eta * N / (ell * ell * ell)
        + 384.0 * h * h * h * eta * eta * eta * N / (ell**4)
    )

    c = b - D
    numerator = (
        2.0 * eta * c * first
        + d * first_eta
        - 2.0 * eta * X * first_X
    )
    numerator_eta = (
        2.0 * c * first
        + 2.0 * eta * (c - 1.0) * first_eta
        + d * first_eta_eta
        - 2.0 * X * first_X
        - 2.0 * eta * X * first_X_eta
    )
    numerator_eta_eta = (
        (4.0 * c - 2.0) * first_eta
        + 2.0 * eta * (c - 2.0) * first_eta_eta
        + d * first_eta_eta_eta
        - 4.0 * X * first_X_eta
        - 2.0 * eta * X * first_X_eta_eta
    )
    parameter2 = (
        numerator_eta_eta / ell
        + 8.0 * h * eta * numerator_eta / (ell * ell)
        + 4.0 * h * numerator / (ell * ell)
        + 32.0 * h * h * eta * eta * numerator / (ell * ell * ell)
    )
    if not math.isfinite(parameter2):
        raise OverflowError(
            "preceding diffusion second parameter derivative is outside binary64 range"
        )

    return PrecedingDiffusionSecondParameterJet(
        value=first_result.value,
        parameter=first_result.parameter,
        parameter2=parameter2,
    )
