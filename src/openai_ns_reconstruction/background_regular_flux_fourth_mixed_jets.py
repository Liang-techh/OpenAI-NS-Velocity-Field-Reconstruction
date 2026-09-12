"""Analytic Eq. (5.2) fourth-mixed jets of the regular flux ``beta=V/X``.

The landed second-eta Eq. (5.6) ``Omega/X`` row consumes
``beta_XXetaeta``, ``beta_Xetaetaeta`` and ``beta_etaetaetaeta``.  The
existing Eq. (5.2) adapter stops one derivative earlier.  Because Eq. (5.2)
contains the radial average of ``partial_eta U``, those three beta entries need
exactly one additional axial layer:

``U_XXetaetaeta``, ``U_Xetaetaetaeta`` and ``U_etaetaetaetaeta``.

This module makes that derivative debt explicit.  The already-tested third
mixed beta jet remains authoritative for every lower field; only the final
three fourth-mixed entries are differentiated here.  No sampled ``V/X``
division and no finite-difference derivative is used in production.

The fifth-mixed U jet is an upstream contract, not a paper-derived coefficient.
This is therefore fail-closed ``formal-structure`` infrastructure and must not
be labelled paper-exact.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Callable
import math

from .background_omega_second_parameter_jet import RegularFluxFourthMixedJet
from .background_regular_flux_third_mixed_jets import (
    AxialFourthMixedJet,
    regular_flux_third_mixed_jet_eq_5_2,
)
from .coordinates import validate_h
from .quadrature import unit_rule


@dataclass(frozen=True)
class AxialFifthMixedJet:
    """Minimal U jet needed for the Eq. (5.2) beta fourth-mixed jet."""

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
            "radial2_parameter3",
            "radial_parameter4",
            "parameter5",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)

    def fourth(self) -> AxialFourthMixedJet:
        """Project to the already-landed Eq. (5.2) third-mixed input."""

        return AxialFourthMixedJet(
            value=self.value,
            radial=self.radial,
            radial2=self.radial2,
            parameter=self.parameter,
            radial_parameter=self.radial_parameter,
            parameter2=self.parameter2,
            radial2_parameter=self.radial2_parameter,
            radial_parameter2=self.radial_parameter2,
            parameter3=self.parameter3,
            radial2_parameter2=self.radial2_parameter2,
            radial_parameter3=self.radial_parameter3,
            parameter4=self.parameter4,
        )


AxialFifthMixedJetProvider = Callable[[float, float], AxialFifthMixedJet]


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
    provider: AxialFifthMixedJetProvider, X: float, eta: float
) -> AxialFifthMixedJet:
    value = provider(X, eta)
    if not isinstance(value, AxialFifthMixedJet):
        raise TypeError("axial_jet_provider must return AxialFifthMixedJet")
    return value


def regular_flux_fourth_mixed_jet_eq_5_2(
    h: float,
    order: int,
    X: float,
    eta: float,
    axial_jet_provider: AxialFifthMixedJetProvider,
    *,
    quadrature_points: int = 32,
) -> RegularFluxFourthMixedJet:
    """Return the minimal fourth-mixed jet of ``beta_n=V_n/X`` from Eq. (5.2).

    The first nine fields are delegated to
    :func:`regular_flux_third_mixed_jet_eq_5_2`.  The final three entries are
    obtained from ``N = L beta`` with

    ``N = 2 eta U - 2 c eta A[U] - (1-eta^2) A[U_eta]``

    and ``c = 1/2-h+lambda_n``.  The only newly required radial-average
    derivatives are ``A_XXetaetaeta``, ``A_Xetaetaetaeta`` and
    ``A_etaetaetaetaeta``.  Their chain-rule radial weights are ``s^2``, ``s``
    and ``1`` respectively, so the axis is evaluated directly without dividing
    by a sampled radial flux.
    """

    h = validate_h(h)
    order = _order(order)
    X, eta = _point(X, eta)
    if not callable(axial_jet_provider):
        raise TypeError("axial_jet_provider must be callable")

    def fourth_provider(x: float, e: float) -> AxialFourthMixedJet:
        return _provider_value(axial_jet_provider, x, e).fourth()

    third = regular_flux_third_mixed_jet_eq_5_2(
        h,
        order,
        X,
        eta,
        fourth_provider,
        quadrature_points=quadrature_points,
    )

    nodes, weights = unit_rule(quadrature_points)
    center = _provider_value(axial_jet_provider, X, eta)
    sums = {
        "AXXE": 0.0,
        "AXXEE": 0.0,
        "AXXEEE": 0.0,
        "AXEE": 0.0,
        "AXEEE": 0.0,
        "AXEEEE": 0.0,
        "AEEE": 0.0,
        "AEEEE": 0.0,
        "AEEEEE": 0.0,
    }
    for node, weight in zip(nodes, weights):
        s = float(node)
        w = float(weight)
        jet = _provider_value(axial_jet_provider, s * X, eta)
        sums["AXXE"] += w * s * s * jet.radial2_parameter
        sums["AXXEE"] += w * s * s * jet.radial2_parameter2
        sums["AXXEEE"] += w * s * s * jet.radial2_parameter3
        sums["AXEE"] += w * s * jet.radial_parameter2
        sums["AXEEE"] += w * s * jet.radial_parameter3
        sums["AXEEEE"] += w * s * jet.radial_parameter4
        sums["AEEE"] += w * jet.parameter3
        sums["AEEEE"] += w * jet.parameter4
        sums["AEEEEE"] += w * jet.parameter5
    if not all(math.isfinite(value) for value in sums.values()):
        raise OverflowError("radial-average derivative is outside binary64 range")

    lam = 2.0 * order * h
    coefficient = 0.5 - h + lam
    d = 1.0 - eta * eta
    ell = 1.0 - 2.0 * h * eta * eta
    if ell <= 0.0:
        raise ValueError("L(h,eta) must be positive")

    # Direct eta derivatives of the Eq. (5.2) numerator.  For any radial
    # derivative level p,
    #
    #   d_eta^q N_p = 2 eta U_{p,q} + 2 q U_{p,q-1}
    #       - d A_{p,q+1} + 2 eta (q-c) A_{p,q}
    #       + q(q-1-2c) A_{p,q-1}.
    #
    # We need (p,q)=(2,2),(1,3),(0,4).
    NXXEE = (
        2.0 * eta * center.radial2_parameter2
        + 4.0 * center.radial2_parameter
        - d * sums["AXXEEE"]
        + 2.0 * eta * (2.0 - coefficient) * sums["AXXEE"]
        + (2.0 - 4.0 * coefficient) * sums["AXXE"]
    )
    NXEEE = (
        2.0 * eta * center.radial_parameter3
        + 6.0 * center.radial_parameter2
        - d * sums["AXEEEE"]
        + 2.0 * eta * (3.0 - coefficient) * sums["AXEEE"]
        + 6.0 * (1.0 - coefficient) * sums["AXEE"]
    )
    NEEEE = (
        2.0 * eta * center.parameter4
        + 8.0 * center.parameter3
        - d * sums["AEEEEE"]
        + 2.0 * eta * (4.0 - coefficient) * sums["AEEEE"]
        + 4.0 * (3.0 - 2.0 * coefficient) * sums["AEEE"]
    )

    ell_eta = -4.0 * h * eta
    ell_eta2 = -4.0 * h

    # Differentiate N=L beta.  L is quadratic in eta, so all higher L
    # derivatives vanish and the fourth mixed beta entries close on the landed
    # third-mixed jet plus the three numerator derivatives above.
    radial2_parameter2 = (
        NXXEE
        - 2.0 * ell_eta * third.radial2_parameter
        - ell_eta2 * third.radial2
    ) / ell
    radial_parameter3 = (
        NXEEE
        - 3.0 * ell_eta * third.radial_parameter2
        - 3.0 * ell_eta2 * third.radial_parameter
    ) / ell
    parameter4 = (
        NEEEE
        - 4.0 * ell_eta * third.parameter3
        - 6.0 * ell_eta2 * third.parameter2
    ) / ell

    values = (radial2_parameter2, radial_parameter3, parameter4)
    if not all(math.isfinite(value) for value in values):
        raise OverflowError("regular-flux fourth-mixed jet is outside binary64 range")

    return RegularFluxFourthMixedJet(
        value=third.value,
        radial=third.radial,
        radial2=third.radial2,
        parameter=third.parameter,
        radial_parameter=third.radial_parameter,
        parameter2=third.parameter2,
        radial2_parameter=third.radial2_parameter,
        radial_parameter2=third.radial_parameter2,
        parameter3=third.parameter3,
        radial2_parameter2=radial2_parameter2,
        radial_parameter3=radial_parameter3,
        parameter4=parameter4,
    )
