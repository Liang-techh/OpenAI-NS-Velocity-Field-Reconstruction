"""Analytic Eq. (5.2) regular-flux second jets from axial-profile jets.

The strict-lower PositiveAxis source consumes second ``(X, eta)`` jets of
``beta_n = V_n / X``.  Eq. (5.2) determines beta_n from the axial coefficient
``U_n`` by

    beta_n = [2 eta U_n - 2 (D + lambda_n) eta A[U_n]
              - d A[partial_eta U_n]] / L,

where ``A[f](X)=integral_0^1 f(s X) ds``, ``D=1/2-h``, ``d=1-eta^2`` and
``L=1-2h eta^2``.

Differentiating beta_n twice requires the mixed third derivatives
``U_XXeta``, ``U_Xetaeta`` and ``U_etaetaeta``.  This module makes that
requirement explicit and performs the differentiation analytically.  Radial
averages are evaluated with the repository's cached Gauss-Legendre rule; no
sampled division by X and no finite-difference derivative occurs in production.

The upstream recursively solved/repaired U_n jet provider is still an input, so
this is ``formal-structure`` infrastructure rather than a paper-exact hierarchy.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Callable
import math

from .background_lower_history_source import ProfileSecondJet
from .coordinates import validate_h
from .quadrature import unit_rule


@dataclass(frozen=True)
class AxialThirdMixedJet:
    """Minimal U-jet needed for the full second jet of beta=V/X.

    The first six fields are the ordinary second ``(X,eta)`` jet.  Eq. (5.2)
    contains ``A[partial_eta U]``; differentiating that term twice adds exactly
    the three total-order-three mixed derivatives below.  No pure ``U_XXX`` is
    needed.
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


AxialThirdMixedJetProvider = Callable[[float, float], AxialThirdMixedJet]


def _point(X: float, eta: float) -> tuple[float, float]:
    X = float(X)
    eta = float(eta)
    if not math.isfinite(X) or X < 0.0:
        raise ValueError("X must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return X, eta


def _order(order: int) -> int:
    if isinstance(order, bool) or not isinstance(order, Integral) or int(order) < 0:
        raise ValueError("order must be a nonnegative integer")
    return int(order)


def _provider_value(
    provider: AxialThirdMixedJetProvider, X: float, eta: float
) -> AxialThirdMixedJet:
    value = provider(X, eta)
    if not isinstance(value, AxialThirdMixedJet):
        raise TypeError("axial_jet_provider must return AxialThirdMixedJet")
    return value


def regular_flux_second_jet_eq_5_2(
    h: float,
    order: int,
    X: float,
    eta: float,
    axial_jet_provider: AxialThirdMixedJetProvider,
    *,
    quadrature_points: int = 32,
) -> ProfileSecondJet:
    """Return the analytic second jet of ``beta_n=V_n/X`` from Eq. (5.2).

    ``order`` selects ``lambda_n=2 n h``.  The radial average derivatives use

        partial_X^a partial_eta^b A[U](X,eta)
          = integral_0^1 s^a partial_X^a partial_eta^b U(sX,eta) ds,

    so the axis ``X=0`` is evaluated directly and never through ``V/X``.
    """

    h = validate_h(h)
    order = _order(order)
    X, eta = _point(X, eta)
    if not callable(axial_jet_provider):
        raise TypeError("axial_jet_provider must be callable")

    nodes, weights = unit_rule(quadrature_points)
    center = _provider_value(axial_jet_provider, X, eta)

    # A_ab = partial_X^a partial_eta^b A[U].  Powers of s are the exact
    # chain-rule weights induced by U(sX,eta).
    sums = {
        "A": 0.0,
        "AX": 0.0,
        "AXX": 0.0,
        "AE": 0.0,
        "AXE": 0.0,
        "AEE": 0.0,
        "AXXE": 0.0,
        "AXEE": 0.0,
        "AEEE": 0.0,
    }
    for node, weight in zip(nodes, weights):
        s = float(node)
        w = float(weight)
        jet = _provider_value(axial_jet_provider, s * X, eta)
        sums["A"] += w * jet.value
        sums["AX"] += w * s * jet.radial
        sums["AXX"] += w * s * s * jet.radial2
        sums["AE"] += w * jet.parameter
        sums["AXE"] += w * s * jet.radial_parameter
        sums["AEE"] += w * jet.parameter2
        sums["AXXE"] += w * s * s * jet.radial2_parameter
        sums["AXEE"] += w * s * jet.radial_parameter2
        sums["AEEE"] += w * jet.parameter3

    if not all(math.isfinite(value) for value in sums.values()):
        raise OverflowError("radial-average derivative is outside binary64 range")

    lam = 2.0 * order * h
    coefficient = 0.5 - h + lam
    d = 1.0 - eta * eta
    ell = 1.0 - 2.0 * h * eta * eta
    if ell <= 0.0:
        raise ValueError("L(h,eta) must be positive")

    A = sums["A"]
    AX = sums["AX"]
    AXX = sums["AXX"]
    AE = sums["AE"]
    AXE = sums["AXE"]
    AEE = sums["AEE"]
    AXXE = sums["AXXE"]
    AXEE = sums["AXEE"]
    AEEE = sums["AEEE"]

    # N is the numerator of Eq. (5.2); all derivatives below are exact
    # product-rule differentiations before the final quotient by L.
    N = (
        2.0 * eta * center.value
        - 2.0 * coefficient * eta * A
        - d * AE
    )
    NX = (
        2.0 * eta * center.radial
        - 2.0 * coefficient * eta * AX
        - d * AXE
    )
    NXX = (
        2.0 * eta * center.radial2
        - 2.0 * coefficient * eta * AXX
        - d * AXXE
    )
    NE = (
        2.0 * center.value
        + 2.0 * eta * center.parameter
        - 2.0 * coefficient * A
        + 2.0 * eta * (1.0 - coefficient) * AE
        - d * AEE
    )
    NXE = (
        2.0 * center.radial
        + 2.0 * eta * center.radial_parameter
        - 2.0 * coefficient * AX
        + 2.0 * eta * (1.0 - coefficient) * AXE
        - d * AXEE
    )
    NEE = (
        4.0 * center.parameter
        + 2.0 * eta * center.parameter2
        + (2.0 - 4.0 * coefficient) * AE
        + 2.0 * eta * (2.0 - coefficient) * AEE
        - d * AEEE
    )

    ell_eta = -4.0 * h * eta
    ell_eta2 = -4.0 * h
    ell2 = ell * ell
    ell3 = ell2 * ell

    value = N / ell
    radial = NX / ell
    radial2 = NXX / ell
    parameter = NE / ell - N * ell_eta / ell2
    radial_parameter = NXE / ell - NX * ell_eta / ell2
    parameter2 = (
        NEE / ell
        - 2.0 * NE * ell_eta / ell2
        - N * ell_eta2 / ell2
        + 2.0 * N * ell_eta * ell_eta / ell3
    )

    values = (value, radial, radial2, parameter, radial_parameter, parameter2)
    if not all(math.isfinite(v) for v in values):
        raise OverflowError("regular-flux second jet is outside binary64 range")
    return ProfileSecondJet(
        value=value,
        radial=radial,
        radial2=radial2,
        parameter=parameter,
        radial_parameter=radial_parameter,
        parameter2=parameter2,
    )
