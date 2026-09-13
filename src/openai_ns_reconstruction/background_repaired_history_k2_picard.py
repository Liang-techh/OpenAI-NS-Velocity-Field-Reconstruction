"""Hierarchy-owned second nontrivial Picard value for Section 5 Eq. (5.7).

The landed Section-5 stack now owns both the value and analytic eta derivative
of the first nontrivial iterate

    W_n^(1) = G(A0 W_n^(0) + A1 partial_eta W_n^(0) + f_n).

This module uses that strong-history-owned pair as the previous iterate in the
same paper-derived Picard map and therefore evaluates

    W_n^(2) = G(A0 W_n^(1) + A1 partial_eta W_n^(1) + f_n).

No matrix, forcing, previous iterate, eta derivative, generic cutoff, fitted
coefficient family, or production finite difference can be supplied by the
caller.  Genuine order-zero fourth/fifth-mixed profile data remain an Issue-#1
dependency, so missing strong leading data fail closed.  This is still Stage-2
``formal-structure`` infrastructure: a second explicit Picard iterate is not a
proof of Picard convergence or a finalized paper-exact recursive coefficient.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from .background_inner_solver import picard_map_eq_5_7
from .background_repaired_history_forcing_second_parameter import (
    hierarchy_owned_positive_axis_forcing_second_parameter_jet,
)
from .background_repaired_history_k1_picard_parameter import (
    hierarchy_owned_k1_picard_parameter_jet_eq_5_7,
)
from .background_repaired_history_phi_fourth_mixed import (
    Section5LowerHistoryPhiFourthMixedHierarchy,
)


def hierarchy_owned_k2_picard_value_eq_5_7(
    hierarchy: Section5LowerHistoryPhiFourthMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
    *,
    quadrature_points: int = 32,
) -> np.ndarray:
    """Return the hierarchy-owned Eq. (5.7) Picard iterate ``W_n^(2)``.

    ``A0`` and ``A1`` come only from ``hierarchy.positive_axis_fields``.
    ``W_n^(1)`` and ``partial_eta W_n^(1)`` come only from the landed analytic
    strong-history k=1 parameter bridge.  ``f_n`` is read from the landed
    hierarchy-owned second-parameter forcing jet, so this public API never
    accepts a caller-maintained forcing table.

    The strong k=1 and forcing paths are explicitly preflighted at the axis
    before applying ``G``.  Thus the exact ``G(0)=0`` shortcut cannot hide
    missing Issue-#1 leading fourth/fifth-mixed data.
    """

    if not isinstance(hierarchy, Section5LowerHistoryPhiFourthMixedHierarchy):
        raise TypeError(
            "hierarchy must be a Section5LowerHistoryPhiFourthMixedHierarchy"
        )

    fields = hierarchy.positive_axis_fields(order)

    # Preserve the strong-history truth boundary even for xi == 0, where the
    # singular inverse would otherwise return before sampling its source.
    hierarchy_owned_k1_picard_parameter_jet_eq_5_7(
        hierarchy,
        order,
        0.0,
        eta,
        quadrature_points=quadrature_points,
    )
    hierarchy_owned_positive_axis_forcing_second_parameter_jet(
        hierarchy,
        order,
        0.0,
        eta,
    )
    fields.A0(0.0, eta)
    fields.A1(0.0, eta)

    @lru_cache(maxsize=None)
    def previous_jet(s: float, eta_value: float):
        return hierarchy_owned_k1_picard_parameter_jet_eq_5_7(
            hierarchy,
            order,
            float(s),
            float(eta_value),
            quadrature_points=quadrature_points,
        )

    @lru_cache(maxsize=None)
    def forcing_jet(s: float, eta_value: float):
        return hierarchy_owned_positive_axis_forcing_second_parameter_jet(
            hierarchy,
            order,
            float(s),
            float(eta_value),
        )

    def previous(s: float, eta_value: float):
        return previous_jet(float(s), float(eta_value)).value

    def previous_parameter(s: float, eta_value: float):
        return previous_jet(float(s), float(eta_value)).parameter

    def forcing(s: float, eta_value: float):
        return forcing_jet(float(s), float(eta_value)).value

    value = picard_map_eq_5_7(
        order,
        xi,
        eta,
        previous,
        previous_parameter,
        fields.A0,
        fields.A1,
        forcing,
        quadrature_points=quadrature_points,
    )
    if value.shape != (6,) or not np.all(np.isfinite(value)):
        raise OverflowError("Eq. (5.7) k=2 Picard value is outside binary64 range")
    return value
