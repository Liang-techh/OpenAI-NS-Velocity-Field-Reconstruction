"""Section-5 PositiveAxis solver bridge from strict lower coefficient history.

The landed :mod:`background_lower_history_source` module already evaluates the
pinned ``actualLowerSource`` / ``lowerHistoryData`` formulas at one ``(X,eta)``
point from second jets of orders strictly below ``n``.  The landed
:mod:`background_positive_axis` module already turns a ``BaseJet`` and
``SourceJet`` into the exact Eq. (5.7) matrices/forcing and the singular
Volterra inverse ``G``.

This module connects those two layers over the *whole radial integration
interval*.  A caller supplies lower-history second-jet functions; the adapter
maps the Eq. (5.7) radius ``xi`` to the paper's squared-radius variable
``X=xi^2``, reconstructs the strict lower history at every quadrature point,
and exposes the genuine theorem-shaped ``A0/A1/f_n`` fields.  In particular,
``first_positive_order_term_from_lower_history`` evaluates the actual
``G f_n`` implied by those lower-history jets, rather than accepting an
unrelated forcing callback.

The jet functions themselves remain upstream inputs.  Until they come from the
materialized, recursively solved and Lemma-5.2-repaired hierarchy (and the
leading order from Issue #1), this bridge is ``formal-structure`` only and must
not be labelled paper-exact.
"""
from __future__ import annotations

from numbers import Integral
from typing import Callable, NamedTuple
import math

import numpy as np

from .background_inner_solver import first_picard_term_eq_5_7
from .background_lower_history_source import (
    ProfileSecondJet,
    positive_axis_point_data_from_lower_history,
)
from .background_positive_axis import (
    PositiveAxisBaseJet,
    PositiveAxisFields,
    PositiveAxisSourceJet,
    positive_axis_eq_5_7_fields,
)
from .coordinates import validate_h


ProfileSecondJetProvider = Callable[[int, float, float], ProfileSecondJet]


class LowerHistoryPositiveAxisProviders(NamedTuple):
    """Structured Eq. (5.7) base/source providers derived from lower history."""

    base: Callable[[float, float], PositiveAxisBaseJet]
    source: Callable[[float, float], PositiveAxisSourceJet]


def _positive_order(order: int) -> int:
    if isinstance(order, bool) or not isinstance(order, Integral) or int(order) < 1:
        raise ValueError("order must be a positive integer")
    return int(order)


def _xi_eta(xi: float, eta: float) -> tuple[float, float]:
    xi = float(xi)
    eta = float(eta)
    if not math.isfinite(xi) or xi < 0.0:
        raise ValueError("xi must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return xi, eta


def _provider(value: object, name: str) -> ProfileSecondJetProvider:
    if not callable(value):
        raise TypeError(f"{name} must be callable")
    return value  # type: ignore[return-value]


def _jet(
    provider: ProfileSecondJetProvider,
    order: int,
    X: float,
    eta: float,
    name: str,
) -> ProfileSecondJet:
    value = provider(order, X, eta)
    if not isinstance(value, ProfileSecondJet):
        raise TypeError(f"{name} must return ProfileSecondJet")
    return value


def lower_history_positive_axis_providers(
    h: float,
    order: int,
    phi_provider: ProfileSecondJetProvider,
    axial_provider: ProfileSecondJetProvider,
    beta_provider: ProfileSecondJetProvider,
) -> LowerHistoryPositiveAxisProviders:
    """Build exact Eq. (5.7) providers from strict lower-history jet functions.

    Providers are queried only for indices ``0,...,order-1``.  At an Eq. (5.7)
    radial point ``xi`` the supplied second jets are evaluated at ``X=xi^2``;
    no interpolation, finite differencing, or sampled radial surrogate is
    introduced by this bridge.
    """

    h = validate_h(h)
    order = _positive_order(order)
    phi_provider = _provider(phi_provider, "phi_provider")
    axial_provider = _provider(axial_provider, "axial_provider")
    beta_provider = _provider(beta_provider, "beta_provider")

    def point_data(xi: float, eta: float):
        xi, eta = _xi_eta(xi, eta)
        X = xi * xi
        phi = [
            _jet(phi_provider, j, X, eta, "phi_provider")
            for j in range(order)
        ]
        axial = [
            _jet(axial_provider, j, X, eta, "axial_provider")
            for j in range(order)
        ]
        beta = [
            _jet(beta_provider, j, X, eta, "beta_provider")
            for j in range(order)
        ]
        return positive_axis_point_data_from_lower_history(
            h, order, X, eta, phi, axial, beta
        )

    def base(xi: float, eta: float) -> PositiveAxisBaseJet:
        return point_data(xi, eta).base

    def source(xi: float, eta: float) -> PositiveAxisSourceJet:
        return point_data(xi, eta).source

    return LowerHistoryPositiveAxisProviders(base=base, source=source)


def positive_axis_eq_5_7_fields_from_lower_history(
    h: float,
    order: int,
    C: float,
    phi_provider: ProfileSecondJetProvider,
    axial_provider: ProfileSecondJetProvider,
    beta_provider: ProfileSecondJetProvider,
) -> PositiveAxisFields:
    """Return the pinned Eq. (5.7) ``A0/A1/f_n`` from strict lower history.

    This removes the manual ``BaseJet``/``SourceJet`` layer for a recursively
    available lower hierarchy.  It does not claim that the supplied lower jets
    are the manuscript hierarchy or that the resulting Picard series converges.
    """

    providers = lower_history_positive_axis_providers(
        h, order, phi_provider, axial_provider, beta_provider
    )
    return positive_axis_eq_5_7_fields(
        h,
        order,
        C,
        providers.base,
        providers.source,
    )


def first_positive_order_term_from_lower_history(
    order: int,
    xi: float,
    eta: float,
    *,
    h: float,
    C: float,
    phi_provider: ProfileSecondJetProvider,
    axial_provider: ProfileSecondJetProvider,
    beta_provider: ProfileSecondJetProvider,
    quadrature_points: int = 32,
) -> np.ndarray:
    """Evaluate the genuine first Lemma-5.1 term ``G f_n`` from lower history.

    Unlike the generic Eq. (5.7) primitive, the forcing is reconstructed at
    every radial quadrature point from the pinned strict-lower convolutions,
    preceding diffusion, and regular Eq. (5.6) ``Omega/X`` source.  This is one
    real recursive solver step *conditional on the supplied lower hierarchy*;
    it is not the converged coefficient and is not paper-exact by itself.
    """

    order = _positive_order(order)
    xi, eta = _xi_eta(xi, eta)
    fields = positive_axis_eq_5_7_fields_from_lower_history(
        h,
        order,
        C,
        phi_provider,
        axial_provider,
        beta_provider,
    )
    return first_picard_term_eq_5_7(
        order,
        xi,
        eta,
        fields.forcing,
        quadrature_points=quadrature_points,
    )
