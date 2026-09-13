"""Analytic Eq. (5.2) fifth-mixed jets of the regular flux ``beta=V/X``.

The landed fourth-mixed adapter stops at
``beta_XXetaeta``, ``beta_Xetaetaeta`` and ``beta_etaetaetaeta``.
The third-eta Eq. (5.6) ``Omega/X`` row needs one additional eta derivative:
``beta_XXetaetaeta``, ``beta_Xetaetaetaeta`` and ``beta_etaetaetaetaeta``.

Eq. (5.2) closes those rows from a sixth-mixed axial ``U`` jet.  The lower
fourth-mixed beta jet remains authoritative; only the three new total-order-five
beta fields are assembled here.  No sampled ``V/X`` division, finite difference,
generic cutoff, or fitted coefficient is used in production.

The sixth-mixed U provider is an upstream contract.  This module is therefore
fail-closed Stage-2 ``formal-structure`` infrastructure, not a paper-exact
coefficient or an all-order claim.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import math

from .background_omega_second_parameter_jet import RegularFluxFourthMixedJet
from .background_regular_flux_fourth_mixed_jets import (
    AxialFifthMixedJet,
    regular_flux_fourth_mixed_jet_eq_5_2,
)
from .coordinates import validate_h
from .quadrature import unit_rule


@dataclass(frozen=True)
class AxialSixthMixedJet:
    """Minimal U jet needed for the Eq. (5.2) beta fifth-mixed jet."""

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
    radial2_parameter3: float
    radial_parameter4: float
    parameter5: float
    radial2_parameter4: float
    radial_parameter5: float
    parameter6: float

    def __post_init__(self) -> None:
        for name in (
            "value", "radial", "radial2", "parameter", "radial_parameter",
            "parameter2", "radial2_parameter", "radial_parameter2", "parameter3",
            "radial2_parameter2", "radial_parameter3", "parameter4",
            "radial2_parameter3", "radial_parameter4", "parameter5",
            "radial2_parameter4", "radial_parameter5", "parameter6",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)

    def fifth(self) -> AxialFifthMixedJet:
        return AxialFifthMixedJet(
            value=self.value, radial=self.radial, radial2=self.radial2,
            parameter=self.parameter, radial_parameter=self.radial_parameter,
            parameter2=self.parameter2, radial2_parameter=self.radial2_parameter,
            radial_parameter2=self.radial_parameter2, parameter3=self.parameter3,
            radial2_parameter2=self.radial2_parameter2,
            radial_parameter3=self.radial_parameter3, parameter4=self.parameter4,
            radial2_parameter3=self.radial2_parameter3,
            radial_parameter4=self.radial_parameter4, parameter5=self.parameter5,
        )


@dataclass(frozen=True)
class RegularFluxFifthMixedJet:
    """Fourth-mixed beta jet plus the eta-bearing total-order-five rows."""

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
    radial2_parameter3: float
    radial_parameter4: float
    parameter5: float

    def __post_init__(self) -> None:
        for name in (
            "value", "radial", "radial2", "parameter", "radial_parameter",
            "parameter2", "radial2_parameter", "radial_parameter2", "parameter3",
            "radial2_parameter2", "radial_parameter3", "parameter4",
            "radial2_parameter3", "radial_parameter4", "parameter5",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)

    def fourth(self) -> RegularFluxFourthMixedJet:
        return RegularFluxFourthMixedJet(
            value=self.value, radial=self.radial, radial2=self.radial2,
            parameter=self.parameter, radial_parameter=self.radial_parameter,
            parameter2=self.parameter2, radial2_parameter=self.radial2_parameter,
            radial_parameter2=self.radial_parameter2, parameter3=self.parameter3,
            radial2_parameter2=self.radial2_parameter2,
            radial_parameter3=self.radial_parameter3, parameter4=self.parameter4,
        )


AxialSixthMixedJetProvider = Callable[[float, float], AxialSixthMixedJet]


def _provider_value(provider: AxialSixthMixedJetProvider, X: float, eta: float) -> AxialSixthMixedJet:
    value = provider(float(X), float(eta))
    if not isinstance(value, AxialSixthMixedJet):
        raise TypeError("axial_jet_provider must return AxialSixthMixedJet")
    return value


def regular_flux_fifth_mixed_jet_eq_5_2(
    h: float, order: int, X: float, eta: float,
    axial_jet_provider: AxialSixthMixedJetProvider, *, quadrature_points: int = 32,
) -> RegularFluxFifthMixedJet:
    """Return the minimal fifth-mixed jet of ``beta_n=V_n/X`` from Eq. (5.2)."""

    h = validate_h(h)
    if isinstance(order, bool) or not isinstance(order, int) or order < 0:
        raise ValueError("order must be a nonnegative integer")
    X = float(X)
    eta = float(eta)
    if not math.isfinite(X) or X < 0.0:
        raise ValueError("X must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    if not callable(axial_jet_provider):
        raise TypeError("axial_jet_provider must be callable")

    def fifth_provider(x: float, e: float) -> AxialFifthMixedJet:
        return _provider_value(axial_jet_provider, x, e).fifth()

    fourth = regular_flux_fourth_mixed_jet_eq_5_2(
        h, order, X, eta, fifth_provider, quadrature_points=quadrature_points,
    )
    nodes, weights = unit_rule(quadrature_points)
    center = _provider_value(axial_jet_provider, X, eta)
    sums = {name: 0.0 for name in (
        "AXXEE", "AXXEEE", "AXXEEEE", "AXEEE", "AXEEEE", "AXEEEEE",
        "AEEEE", "AEEEEE", "AEEEEEE",
    )}
    for node, weight in zip(nodes, weights):
        s = float(node)
        w = float(weight)
        jet = _provider_value(axial_jet_provider, s * X, eta)
        sums["AXXEE"] += w * s * s * jet.radial2_parameter2
        sums["AXXEEE"] += w * s * s * jet.radial2_parameter3
        sums["AXXEEEE"] += w * s * s * jet.radial2_parameter4
        sums["AXEEE"] += w * s * jet.radial_parameter3
        sums["AXEEEE"] += w * s * jet.radial_parameter4
        sums["AXEEEEE"] += w * s * jet.radial_parameter5
        sums["AEEEE"] += w * jet.parameter4
        sums["AEEEEE"] += w * jet.parameter5
        sums["AEEEEEE"] += w * jet.parameter6
    if not all(math.isfinite(value) for value in sums.values()):
        raise OverflowError("radial-average derivative is outside binary64 range")

    coefficient = 0.5 - h + 2.0 * order * h
    d = 1.0 - eta * eta
    ell = 1.0 - 2.0 * h * eta * eta
    if ell <= 0.0:
        raise ValueError("L(h,eta) must be positive")

    NXXEEE = (
        2.0 * eta * center.radial2_parameter3 + 6.0 * center.radial2_parameter2
        - d * sums["AXXEEEE"]
        + 2.0 * eta * (3.0 - coefficient) * sums["AXXEEE"]
        + 6.0 * (1.0 - coefficient) * sums["AXXEE"]
    )
    NXEEEE = (
        2.0 * eta * center.radial_parameter4 + 8.0 * center.radial_parameter3
        - d * sums["AXEEEEE"]
        + 2.0 * eta * (4.0 - coefficient) * sums["AXEEEE"]
        + 4.0 * (3.0 - 2.0 * coefficient) * sums["AXEEE"]
    )
    NEEEEE = (
        2.0 * eta * center.parameter5 + 10.0 * center.parameter4
        - d * sums["AEEEEEE"]
        + 2.0 * eta * (5.0 - coefficient) * sums["AEEEEE"]
        + 10.0 * (2.0 - coefficient) * sums["AEEEE"]
    )

    ell_eta = -4.0 * h * eta
    ell_eta2 = -4.0 * h
    radial2_parameter3 = (
        NXXEEE - 3.0 * ell_eta * fourth.radial2_parameter2
        - 3.0 * ell_eta2 * fourth.radial2_parameter
    ) / ell
    radial_parameter4 = (
        NXEEEE - 4.0 * ell_eta * fourth.radial_parameter3
        - 6.0 * ell_eta2 * fourth.radial_parameter2
    ) / ell
    parameter5 = (
        NEEEEE - 5.0 * ell_eta * fourth.parameter4
        - 10.0 * ell_eta2 * fourth.parameter3
    ) / ell
    if not all(math.isfinite(value) for value in (radial2_parameter3, radial_parameter4, parameter5)):
        raise OverflowError("regular-flux fifth-mixed jet is outside binary64 range")

    return RegularFluxFifthMixedJet(
        value=fourth.value, radial=fourth.radial, radial2=fourth.radial2,
        parameter=fourth.parameter, radial_parameter=fourth.radial_parameter,
        parameter2=fourth.parameter2, radial2_parameter=fourth.radial2_parameter,
        radial_parameter2=fourth.radial_parameter2, parameter3=fourth.parameter3,
        radial2_parameter2=fourth.radial2_parameter2,
        radial_parameter3=fourth.radial_parameter3, parameter4=fourth.parameter4,
        radial2_parameter3=radial2_parameter3, radial_parameter4=radial_parameter4,
        parameter5=parameter5,
    )
