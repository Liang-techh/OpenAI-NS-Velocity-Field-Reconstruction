"""Exact real positive-order Eq. (5.7) coefficients from the pinned Lean system.

This module ports the displayed ``A0``, ``A1`` and forcing formulas from
``NavierStokes/PositiveAxisSystem.lean``.  It removes the arbitrary-matrix gap
in the landed Eq. (5.7) Picard solver: callers provide structured leading/base
jets and strictly lower-order source jets, and this module constructs the
six-component system used by Lemma 5.1.

The structured jets are still inputs.  Until Issue #1 materializes the actual
leading profile and the recursive lower-order history, outputs remain
formal-structure and must not be labelled paper-exact.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Callable, NamedTuple
import math

import numpy as np

from .coordinates import validate_h
from .background_inner_solver import first_picard_term_eq_5_7


@dataclass(frozen=True)
class RadialParameterJet:
    """One ``PositiveAxisSystem.Jet`` in the squared-radius variable X."""

    value: float
    radial: float
    radial2: float
    parameter: float

    def __post_init__(self) -> None:
        for name in ("value", "radial", "radial2", "parameter"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class PositiveAxisBaseJet:
    """Leading/base data entering the exact positive-order matrices."""

    phi: RadialParameterJet
    axial: RadialParameterJet
    beta: float

    def __post_init__(self) -> None:
        if not isinstance(self.phi, RadialParameterJet) or not isinstance(
            self.axial, RadialParameterJet
        ):
            raise TypeError("phi and axial must be RadialParameterJet instances")
        beta = float(self.beta)
        if not math.isfinite(beta):
            raise ValueError("beta must be finite")
        object.__setattr__(self, "beta", beta)


@dataclass(frozen=True)
class PositiveAxisSourceJet:
    """Strictly lower-order source data in ``PositiveAxisSystem.SourceJet``."""

    angular: float
    axial: float
    pressure_product: float
    omega_quotient: float

    def __post_init__(self) -> None:
        for name in ("angular", "axial", "pressure_product", "omega_quotient"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)


BaseJetProvider = Callable[[float, float], PositiveAxisBaseJet]
SourceJetProvider = Callable[[float, float], PositiveAxisSourceJet]


class PositiveAxisFields(NamedTuple):
    A0: Callable[[float, float], np.ndarray]
    A1: Callable[[float, float], np.ndarray]
    forcing: Callable[[float, float], np.ndarray]


def _xi(value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError("xi must be finite and nonnegative")
    return value


def _eta(value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or abs(value) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return value


def _positive_order(order: int) -> int:
    if isinstance(order, bool) or not isinstance(order, Integral) or int(order) < 1:
        raise ValueError("order must be a positive integer")
    return int(order)


def _positive_C(value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("C must be finite and positive")
    return value


def _base(value: PositiveAxisBaseJet) -> PositiveAxisBaseJet:
    if not isinstance(value, PositiveAxisBaseJet):
        raise TypeError("base provider must return PositiveAxisBaseJet")
    return value


def _source(value: PositiveAxisSourceJet) -> PositiveAxisSourceJet:
    if not isinstance(value, PositiveAxisSourceJet):
        raise TypeError("source provider must return PositiveAxisSourceJet")
    return value


def _scalars(h: float, eta: float) -> tuple[float, float, float, float]:
    h = validate_h(h)
    eta = _eta(eta)
    a = 0.5 + h
    d = 0.5 - h
    edge = 1.0 - eta * eta
    ell = 1.0 - 2.0 * h * eta * eta
    if ell <= 0.0:  # impossible on the validated real geometry domain, but fail closed.
        raise ValueError("ell(h,eta) must be positive")
    return a, d, edge, ell


def pressure_source(C: float, source: PositiveAxisSourceJet) -> float:
    """Pinned ``pressureSource = C^-2 pressureProduct - omegaQuotient/2``."""
    C = _positive_C(C)
    source = _source(source)
    return source.pressure_product / (C * C) - 0.5 * source.omega_quotient


def positive_axis_A0(
    h: float,
    order: int,
    C: float,
    xi: float,
    eta: float,
    base: PositiveAxisBaseJet,
) -> np.ndarray:
    """Return the exact displayed ``PositiveAxisSystem.A0`` at lambda_n=2nh."""
    h = validate_h(h)
    order = _positive_order(order)
    C = _positive_C(C)
    xi, eta, base = _xi(xi), _eta(eta), _base(base)
    lam = 2.0 * order * h
    a, d, _edge, ell = _scalars(h, eta)
    angular_power = -a - 0.5
    axial_power = -a
    X = xi * xi
    inv_c2 = 1.0 / (C * C)

    M = 1.0 - 2.0 * eta * base.axial.value
    R = M / ell + base.beta
    Gphi = X * base.phi.radial + base.phi.value
    Gu = X * base.axial.radial

    def axial_value(power: float, jet: RadialParameterJet) -> float:
        return (
            2.0 * eta * power * jet.value
            + (1.0 - eta * eta) * jet.parameter
            - 2.0 * eta * X * jet.radial
        ) / ell

    out = np.zeros((6, 6), dtype=float)
    out[0, 4] = 1.0
    out[1, 5] = 1.0
    out[2, 5] = -1.0
    out[3, 0] = 4.0 * xi * inv_c2 * base.phi.value
    out[4, 0] = 2.0 * (base.beta - (angular_power + lam) * M / ell)
    out[4, 1] = 2.0 * (
        axial_value(angular_power, base.phi)
        + 2.0 * eta * (a - lam) * Gphi / ell
    )
    out[4, 2] = -4.0 * eta * (d + lam) * Gphi / ell
    out[4, 4] = xi * R
    out[5, 0] = -8.0 * eta * X * inv_c2 * base.phi.value / ell
    out[5, 1] = 2.0 * (
        -(axial_power + lam) * M / ell
        + axial_value(axial_power, base.axial)
        + 2.0 * eta * (a - lam) * Gu / ell
    )
    out[5, 2] = -4.0 * eta * (d + lam) * Gu / ell
    out[5, 3] = 4.0 * eta * (-2.0 * a + lam) / ell
    out[5, 5] = xi * R
    if not np.all(np.isfinite(out)):
        raise OverflowError("positive-axis A0 is outside binary64 range")
    return out


def positive_axis_A1(
    h: float,
    xi: float,
    eta: float,
    base: PositiveAxisBaseJet,
) -> np.ndarray:
    """Return the exact sparse ``PositiveAxisSystem.A1`` matrix."""
    h = validate_h(h)
    xi, eta, base = _xi(xi), _eta(eta), _base(base)
    _a, d, edge, ell = _scalars(h, eta)
    X = xi * xi
    H = d * eta + edge * base.axial.value
    Gphi = X * base.phi.radial + base.phi.value
    Gu = X * base.axial.radial

    out = np.zeros((6, 6), dtype=float)
    out[4, 0] = 2.0 * H / ell
    out[4, 1] = -2.0 * edge * Gphi / ell
    out[4, 2] = -2.0 * edge * Gphi / ell
    out[5, 1] = 2.0 * (H - edge * Gu) / ell
    out[5, 2] = -2.0 * edge * Gu / ell
    out[5, 3] = 2.0 * edge / ell
    if not np.all(np.isfinite(out)):
        raise OverflowError("positive-axis A1 is outside binary64 range")
    return out


def positive_axis_forcing(
    h: float,
    C: float,
    xi: float,
    eta: float,
    source: PositiveAxisSourceJet,
) -> np.ndarray:
    """Return the exact six-vector ``PositiveAxisSystem.forcing``."""
    h = validate_h(h)
    C = _positive_C(C)
    xi, eta, source = _xi(xi), _eta(eta), _source(source)
    _a, _d, _edge, ell = _scalars(h, eta)
    psrc = pressure_source(C, source)
    X = xi * xi
    out = np.array(
        [
            0.0,
            0.0,
            0.0,
            2.0 * xi * psrc,
            2.0 * source.angular,
            2.0 * source.axial - 4.0 * eta * X * psrc / ell,
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(out)):
        raise OverflowError("positive-axis forcing is outside binary64 range")
    return out


def positive_axis_eq_5_7_fields(
    h: float,
    order: int,
    C: float,
    base_provider: BaseJetProvider,
    source_provider: SourceJetProvider,
) -> PositiveAxisFields:
    """Build Eq. (5.7) A0/A1/f_n callbacks from structured theorem data.

    This adapter is the bridge into ``background_inner_solver.picard_map_eq_5_7``.
    The formulas are pinned-paper/Lean exact, but the provider data are not
    certified here; caller data therefore remain formal-structure.
    """
    h = validate_h(h)
    order = _positive_order(order)
    C = _positive_C(C)
    if not callable(base_provider) or not callable(source_provider):
        raise TypeError("base_provider and source_provider must be callable")

    def A0(xi: float, eta: float) -> np.ndarray:
        return positive_axis_A0(h, order, C, xi, eta, _base(base_provider(xi, eta)))

    def A1(xi: float, eta: float) -> np.ndarray:
        return positive_axis_A1(h, xi, eta, _base(base_provider(xi, eta)))

    def forcing(xi: float, eta: float) -> np.ndarray:
        return positive_axis_forcing(h, C, xi, eta, _source(source_provider(xi, eta)))

    return PositiveAxisFields(A0=A0, A1=A1, forcing=forcing)


def first_positive_order_term_from_jets(
    order: int,
    xi: float,
    eta: float,
    *,
    h: float,
    C: float,
    base_provider: BaseJetProvider,
    source_provider: SourceJetProvider,
    quadrature_points: int = 32,
) -> np.ndarray:
    """Evaluate the genuine first Lemma-5.1 term ``G f_n`` for exact f_n.

    Unlike ``first_picard_term_eq_5_7`` with an arbitrary forcing callback,
    this function first constructs the displayed ``f_n`` from the structured
    lower-order ``SourceJet``.  It is still not the converged coefficient and
    is not paper-exact until the providers come from the materialized hierarchy.
    """
    fields = positive_axis_eq_5_7_fields(h, order, C, base_provider, source_provider)
    return first_picard_term_eq_5_7(
        order,
        xi,
        eta,
        fields.forcing,
        quadrature_points=quadrature_points,
    )
