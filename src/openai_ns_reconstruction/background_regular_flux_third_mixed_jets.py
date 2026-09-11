"""Analytic Eq. (5.2) third-mixed jets of the regular flux ``beta=V/X``.

The landed Eq. (5.6) parameter-jet row needs exactly
``beta_XXeta``, ``beta_Xetaeta`` and ``beta_etaetaeta``.  The existing
Eq. (5.2) adapter stops at the ordinary second jet because those three
quantities require one additional derivative layer of the axial coefficient.

This module makes that derivative debt explicit.  It augments the landed
``AxialThirdMixedJet`` by exactly
``U_XXetaeta``, ``U_Xetaetaeta`` and ``U_etaetaetaeta`` and analytically
differentiates Eq. (5.2) to produce the three required beta derivatives.
No sampled ``V/X`` division and no finite-difference derivative is used in
production.

The fourth-mixed U jet is still an upstream input; the Lemma 5.2 repair and the
Section-5 hierarchy do not yet own it.  This is therefore fail-closed
``formal-structure`` solver infrastructure, not a paper-exact coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Callable
import math

from .background_omega_parameter_jet import RegularFluxThirdMixedJet
from .background_regular_flux_second_jets import (
    AxialThirdMixedJet,
    regular_flux_second_jet_eq_5_2,
)
from .coordinates import validate_h
from .quadrature import unit_rule


@dataclass(frozen=True)
class AxialFourthMixedJet:
    """Minimal U jet needed for the Eq. (5.2) beta third-mixed jet.

    The first nine fields are the landed :class:`AxialThirdMixedJet`.  Because
    Eq. (5.2) contains ``A[partial_eta U]``, one further eta differentiation of
    the beta second jet introduces exactly the final three total-order-four
    derivatives.  No pure ``U_XXX`` or ``U_XXXX`` derivative is needed.
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

    def third(self) -> AxialThirdMixedJet:
        """Project to the already-landed Eq. (5.2) second-jet input."""

        return AxialThirdMixedJet(
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


AxialFourthMixedJetProvider = Callable[[float, float], AxialFourthMixedJet]


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


def _provider_value(
    provider: AxialFourthMixedJetProvider, X: float, eta: float
) -> AxialFourthMixedJet:
    value = provider(X, eta)
    if not isinstance(value, AxialFourthMixedJet):
        raise TypeError("axial_jet_provider must return AxialFourthMixedJet")
    return value


def regular_flux_third_mixed_jet_eq_5_2(
    h: float,
    order: int,
    X: float,
    eta: float,
    axial_jet_provider: AxialFourthMixedJetProvider,
    *,
    quadrature_points: int = 32,
) -> RegularFluxThirdMixedJet:
    """Return the minimal third-mixed jet of ``beta_n=V_n/X`` from Eq. (5.2).

    The ordinary second jet is delegated to the landed analytic adapter.  The
    three additional derivatives are obtained by differentiating its underlying
    Eq. (5.2) numerator exactly.  The only new radial-average derivatives are

    ``A_XXetaeta``, ``A_Xetaetaeta`` and ``A_etaetaetaeta``;

    their chain-rule weights are respectively ``s^2``, ``s`` and ``1``.
    Hence the axis ``X=0`` is evaluated directly and never through division by
    a sampled ``V``.
    """

    h = validate_h(h)
    order = _order(order)
    X, eta = _point(X, eta)
    if not callable(axial_jet_provider):
        raise TypeError("axial_jet_provider must be callable")

    def third_provider(x: float, e: float) -> AxialThirdMixedJet:
        return _provider_value(axial_jet_provider, x, e).third()

    second = regular_flux_second_jet_eq_5_2(
        h,
        order,
        X,
        eta,
        third_provider,
        quadrature_points=quadrature_points,
    )

    nodes, weights = unit_rule(quadrature_points)
    center = _provider_value(axial_jet_provider, X, eta)
    sums = {
        "AXX": 0.0,
        "AXXE": 0.0,
        "AXXEE": 0.0,
        "AXE": 0.0,
        "AXEE": 0.0,
        "AXEEE": 0.0,
        "AEE": 0.0,
        "AEEE": 0.0,
        "AEEEE": 0.0,
    }
    for node, weight in zip(nodes, weights):
        s = float(node)
        w = float(weight)
        jet = _provider_value(axial_jet_provider, s * X, eta)
        sums["AXX"] += w * s * s * jet.radial2
        sums["AXXE"] += w * s * s * jet.radial2_parameter
        sums["AXXEE"] += w * s * s * jet.radial2_parameter2
        sums["AXE"] += w * s * jet.radial_parameter
        sums["AXEE"] += w * s * jet.radial_parameter2
        sums["AXEEE"] += w * s * jet.radial_parameter3
        sums["AEE"] += w * jet.parameter2
        sums["AEEE"] += w * jet.parameter3
        sums["AEEEE"] += w * jet.parameter4
    if not all(math.isfinite(value) for value in sums.values()):
        raise OverflowError("radial-average derivative is outside binary64 range")

    lam = 2.0 * order * h
    coefficient = 0.5 - h + lam
    d = 1.0 - eta * eta
    ell = 1.0 - 2.0 * h * eta * eta
    if ell <= 0.0:
        raise ValueError("L(h,eta) must be positive")

    # Third derivatives of the Eq. (5.2) numerator N.  These formulas are the
    # direct product-rule derivatives of
    # N = 2 eta U - 2 coefficient eta A[U] - d A[U_eta].
    NXXE = (
        2.0 * center.radial2
        + 2.0 * eta * center.radial2_parameter
        - 2.0 * coefficient * sums["AXX"]
        + 2.0 * eta * (1.0 - coefficient) * sums["AXXE"]
        - d * sums["AXXEE"]
    )
    NXEE = (
        4.0 * center.radial_parameter
        + 2.0 * eta * center.radial_parameter2
        + (2.0 - 4.0 * coefficient) * sums["AXE"]
        + 2.0 * eta * (2.0 - coefficient) * sums["AXEE"]
        - d * sums["AXEEE"]
    )
    NEEE = (
        6.0 * center.parameter2
        + 2.0 * eta * center.parameter3
        + 6.0 * (1.0 - coefficient) * sums["AEE"]
        + 2.0 * eta * (3.0 - coefficient) * sums["AEEE"]
        - d * sums["AEEEE"]
    )

    ell_eta = -4.0 * h * eta
    ell_eta2 = -4.0 * h
    ell2 = ell * ell
    ell3 = ell2 * ell
    ell4 = ell3 * ell

    # Recover the lower numerator derivatives from N = L beta and the landed
    # beta second jet.  This avoids duplicating the already-tested Eq. (5.2)
    # value/first/second derivative algebra.
    N = ell * second.value
    NE = ell * second.parameter + ell_eta * second.value
    NEE = (
        ell * second.parameter2
        + 2.0 * ell_eta * second.parameter
        + ell_eta2 * second.value
    )
    NX = ell * second.radial
    NXE = ell * second.radial_parameter + ell_eta * second.radial
    NXX = ell * second.radial2

    radial2_parameter = NXXE / ell - NXX * ell_eta / ell2
    radial_parameter2 = (
        NXEE / ell
        - 2.0 * NXE * ell_eta / ell2
        - NX * ell_eta2 / ell2
        + 2.0 * NX * ell_eta * ell_eta / ell3
    )
    parameter3 = (
        NEEE / ell
        - 3.0 * NEE * ell_eta / ell2
        - 3.0 * NE * ell_eta2 / ell2
        + 6.0 * NE * ell_eta * ell_eta / ell3
        + 6.0 * N * ell_eta * ell_eta2 / ell3
        - 6.0 * N * ell_eta * ell_eta * ell_eta / ell4
    )

    values = (radial2_parameter, radial_parameter2, parameter3)
    if not all(math.isfinite(value) for value in values):
        raise OverflowError("regular-flux third-mixed jet is outside binary64 range")
    return RegularFluxThirdMixedJet(
        value=second.value,
        radial=second.radial,
        radial2=second.radial2,
        parameter=second.parameter,
        radial_parameter=second.radial_parameter,
        parameter2=second.parameter2,
        radial2_parameter=radial2_parameter,
        radial_parameter2=radial_parameter2,
        parameter3=parameter3,
    )
