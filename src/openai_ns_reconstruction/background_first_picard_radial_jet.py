"""Positive-radius radial jet of the hierarchy-owned first Eq. (5.7) Picard term.

For one positive Section-5 order, the landed repaired-history bridge determines the
exact displayed Eq. (5.7) forcing ``f_n`` from strict lower coefficients and the
first Lemma-5.1 term is

    W^(0)_n = G f_n,

where component ``i`` of the singular inverse is

    W_i(xi,eta) = integral_0^xi (s/xi)^c_i f_i(s,eta) ds,
    c = (0,0,2,0,3,1).

On the positive axis ``xi > 0`` this identity can be differentiated without a
finite-difference approximation:

    d_xi W_i = f_i(xi,eta) - c_i W_i / xi.

This module exposes that radial derivative for the *hierarchy-owned* strict-lower
forcing.  It is a small step toward making solved Eq. (5.7) iterates usable as
coefficient jets.  It deliberately does not manufacture the missing eta derivative
needed for the next Picard map, and it does not infer an axis derivative at xi=0
from point data alone.  The result therefore remains formal-structure solver
infrastructure and never upgrades the reconstruction to paper-exact status.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .background_inner_solver import EQ_5_7_SINGULAR_WEIGHTS
from .background_repaired_history import Section5LowerHistoryJetHierarchy


@dataclass(frozen=True)
class FirstPicardRadialJet:
    """Value and exact positive-radius radial derivative of ``G f_n``.

    ``ode_residual`` records

        d_xi W + xi^-1 diag(c) W - f_n,

    using the same hierarchy-owned forcing.  It is a numerical consistency
    diagnostic, not an independent proof; regression tests differentiate the
    integral representation independently.
    """

    value: np.ndarray
    d_xi: np.ndarray
    forcing: np.ndarray
    ode_residual: np.ndarray

    def __post_init__(self) -> None:
        for name in ("value", "d_xi", "forcing", "ode_residual"):
            array = np.asarray(getattr(self, name), dtype=float)
            if array.shape != (6,) or not np.all(np.isfinite(array)):
                raise ValueError(f"{name} must be a finite vector of shape (6,)")
            array = np.array(array, dtype=float, copy=True)
            array.setflags(write=False)
            object.__setattr__(self, name, array)

    @property
    def max_abs_residual(self) -> float:
        return float(np.max(np.abs(self.ode_residual)))


def first_picard_radial_jet_from_hierarchy(
    hierarchy: Section5LowerHistoryJetHierarchy,
    order: int,
    xi: float,
    eta: float,
    *,
    quadrature_points: int | None = None,
) -> FirstPicardRadialJet:
    """Return ``(G f_n, d_xi G f_n)`` from owned strict lower history.

    Only ``xi > 0`` is admitted.  At the axis, differentiability of the singular
    integral requires a neighborhood regularity statement about ``f_n``; a
    single point value is not enough, so this adapter fails closed instead of
    silently asserting an axis derivative.

    The coefficient matrices ``A0`` and ``A1`` do not enter the first term.
    They start contributing at the next Picard application.  Consequently this
    routine is a genuine derivative of the paper's first term but is not the
    converged order-``n`` coefficient.
    """

    if not isinstance(hierarchy, Section5LowerHistoryJetHierarchy):
        raise TypeError("hierarchy must be a Section5LowerHistoryJetHierarchy")
    xi = float(xi)
    eta = float(eta)
    if not math.isfinite(xi) or xi <= 0.0:
        raise ValueError("xi must be finite and strictly positive")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")

    fields = hierarchy.positive_axis_fields(order)
    value = hierarchy.first_picard_term(
        order,
        xi,
        eta,
        quadrature_points=quadrature_points,
    )
    forcing = np.asarray(fields.forcing(xi, eta), dtype=float)
    if forcing.shape != (6,) or not np.all(np.isfinite(forcing)):
        raise ValueError("hierarchy forcing must be a finite vector of shape (6,)")

    weights = np.asarray(EQ_5_7_SINGULAR_WEIGHTS, dtype=float)
    d_xi = forcing - weights * value / xi
    residual = d_xi + weights * value / xi - forcing
    if not np.all(np.isfinite(d_xi)) or not np.all(np.isfinite(residual)):
        raise OverflowError("first Picard radial jet is outside binary64 range")

    return FirstPicardRadialJet(
        value=value,
        d_xi=d_xi,
        forcing=forcing,
        ode_residual=residual,
    )
