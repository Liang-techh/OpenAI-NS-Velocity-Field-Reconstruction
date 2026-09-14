"""Hierarchy-owned leading/base bridge for the Section 5 PositiveAxis matrices.

The displayed Eq. (5.7) matrices ``A0`` and ``A1`` depend on the order-zero
profiles only through ``PositiveAxisBaseJet``: the value/X/XX/eta rows of
``phi_0`` and ``U_0`` and the regular flux ``beta_0=V_0/X``.  Stage 2 already
owns substantially stronger order-zero data inside
``Section5LowerHistorySixthMixedHierarchy``.  This module projects exactly the
needed rows from that one authoritative hierarchy instead of accepting a
second caller-maintained base-jet table.

``beta_0`` is regenerated through the hierarchy's analytic Eq. (5.2) adapter;
no sampled ``V/X`` division, finite difference, fit, generic cutoff, or manual
repair jet is used.  The underlying order-zero strong profile is still an
Issue-1/upstream input, so these bridges remain Stage-2 ``formal-structure``
and do not make the coefficient or velocity paper-exact.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import math

import numpy as np

from .background_positive_axis import (
    PositiveAxisBaseJet,
    RadialParameterJet,
    positive_axis_A0,
    positive_axis_A1,
)
from .background_repaired_history_sixth_mixed import (
    Section5LowerHistorySixthMixedHierarchy,
)


def _positive_order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("order must be a positive integer")
    return int(value)


def _X_eta(X: float, eta: float) -> tuple[float, float]:
    X = float(X)
    eta = float(eta)
    if not math.isfinite(X) or X < 0.0:
        raise ValueError("X must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return X, eta


def _xi_eta(xi: float, eta: float) -> tuple[float, float]:
    xi = float(xi)
    eta = float(eta)
    if not math.isfinite(xi) or xi < 0.0:
        raise ValueError("xi must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return xi, eta


def hierarchy_owned_positive_axis_base_jet(
    hierarchy: Section5LowerHistorySixthMixedHierarchy,
    X: float,
    eta: float,
) -> PositiveAxisBaseJet:
    """Project the exact Eq. (5.7) order-zero base jet from one hierarchy.

    The bridge deliberately uses coefficient order zero.  Fifth-mixed ``phi``
    and sixth-mixed ``U`` are only projected downward; ``beta_0`` is obtained
    from the same ``U_0`` through the hierarchy-owned analytic Eq. (5.2)
    adapter.  An unrelated ``PositiveAxisBaseJet`` cannot be supplied.
    """

    if not isinstance(hierarchy, Section5LowerHistorySixthMixedHierarchy):
        raise TypeError("hierarchy must be a Section5LowerHistorySixthMixedHierarchy")
    X, eta = _X_eta(X, eta)

    phi = hierarchy.phi_fifth_mixed_jet(0, X, eta)
    axial = hierarchy.axial_sixth_mixed_jet(0, X, eta)
    beta = hierarchy.beta_fifth_mixed_jet(0, X, eta)

    return PositiveAxisBaseJet(
        phi=RadialParameterJet(
            value=phi.value,
            radial=phi.radial,
            radial2=phi.radial2,
            parameter=phi.parameter,
        ),
        axial=RadialParameterJet(
            value=axial.value,
            radial=axial.radial,
            radial2=axial.radial2,
            parameter=axial.parameter,
        ),
        beta=beta.value,
    )


def _matrix6(value: np.ndarray, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (6, 6) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite matrix of shape (6, 6)")
    out = out.copy()
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class HierarchyOwnedPositiveAxisMatrices:
    """Value rows of the exact displayed ``A0`` and ``A1`` matrices."""

    A0: np.ndarray
    A1: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "A0", _matrix6(self.A0, "A0"))
        object.__setattr__(self, "A1", _matrix6(self.A1, "A1"))


def hierarchy_owned_positive_axis_matrices(
    hierarchy: Section5LowerHistorySixthMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
) -> HierarchyOwnedPositiveAxisMatrices:
    """Build Eq. (5.7) ``A0/A1`` values from the hierarchy-owned base row."""

    if not isinstance(hierarchy, Section5LowerHistorySixthMixedHierarchy):
        raise TypeError("hierarchy must be a Section5LowerHistorySixthMixedHierarchy")
    order = _positive_order(order)
    xi, eta = _xi_eta(xi, eta)
    base = hierarchy_owned_positive_axis_base_jet(hierarchy, xi * xi, eta)
    return HierarchyOwnedPositiveAxisMatrices(
        A0=positive_axis_A0(hierarchy.h, order, hierarchy.C, xi, eta, base),
        A1=positive_axis_A1(hierarchy.h, xi, eta, base),
    )
