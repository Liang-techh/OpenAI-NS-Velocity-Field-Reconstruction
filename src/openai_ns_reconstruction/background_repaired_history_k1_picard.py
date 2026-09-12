"""Hierarchy-owned nontrivial ``k=1`` Picard application for Section 5 Eq. (5.7).

PR #179 landed the strong-history-owned parameter jet

    W_n^(0) = G f_n,
    partial_eta W_n^(0) = G(partial_eta f_n).

The lower-history hierarchy already constructs the displayed Eq. (5.7)
``A0``/``A1`` matrices from its owned order-zero base jets.  This module joins
those landed pieces and evaluates the first nontrivial Picard application

    W_n^(1) = G(A0 W_n^(0) + A1 partial_eta W_n^(0) + f_n).

No matrix, forcing, previous iterate, or eta derivative can be supplied by the
caller.  Production uses the existing paper-derived matrix bridge and singular
inverse only; it introduces no finite difference, generic cutoff, or fitted
coefficient.

This is still Stage-2 ``formal-structure`` infrastructure.  In particular the
Issue-#1 leading mixed profile remains upstream, one Picard application is not
Picard convergence, and the returned value is not a paper-exact recursive
coefficient or velocity.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from .background_inner_solver import picard_map_eq_5_7
from .background_repaired_history_first_picard_parameter import (
    hierarchy_owned_first_picard_parameter_jet_eq_5_7,
)
from .background_repaired_history_forcing_parameter import (
    hierarchy_owned_positive_axis_forcing_parameter_jet,
)
from .background_repaired_history_phi_third_mixed import (
    Section5LowerHistoryPhiThirdMixedHierarchy,
)


def hierarchy_owned_k1_picard_value_eq_5_7(
    hierarchy: Section5LowerHistoryPhiThirdMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
    *,
    quadrature_points: int = 32,
) -> np.ndarray:
    """Return the hierarchy-owned Eq. (5.7) Picard iterate ``W_n^(1)``.

    ``A0`` and ``A1`` are obtained only from ``hierarchy.positive_axis_fields``.
    ``W_n^(0)`` and its analytic eta derivative are obtained only from the
    landed strong-history first-Picard parameter bridge, while ``f_n`` is
    obtained only from the landed hierarchy-owned forcing parameter bridge.

    The strong source path is explicitly preflighted at the axis before calling
    the singular inverse.  Therefore the exact ``G(0)=0`` shortcut cannot hide
    missing Issue-#1 fourth-mixed leading data or an incoherent strong history.
    """

    if not isinstance(hierarchy, Section5LowerHistoryPhiThirdMixedHierarchy):
        raise TypeError(
            "hierarchy must be a Section5LowerHistoryPhiThirdMixedHierarchy"
        )

    fields = hierarchy.positive_axis_fields(order)

    # Preserve the strong-history truth boundary even when xi == 0, where G
    # would otherwise return before evaluating the Picard right-hand side.
    hierarchy_owned_first_picard_parameter_jet_eq_5_7(
        hierarchy,
        order,
        0.0,
        eta,
        quadrature_points=quadrature_points,
    )
    hierarchy_owned_positive_axis_forcing_parameter_jet(
        hierarchy,
        order,
        0.0,
        eta,
    )
    fields.A0(0.0, eta)
    fields.A1(0.0, eta)

    @lru_cache(maxsize=None)
    def first_jet(s: float, eta_value: float):
        return hierarchy_owned_first_picard_parameter_jet_eq_5_7(
            hierarchy,
            order,
            s,
            eta_value,
            quadrature_points=quadrature_points,
        )

    @lru_cache(maxsize=None)
    def forcing_jet(s: float, eta_value: float):
        return hierarchy_owned_positive_axis_forcing_parameter_jet(
            hierarchy,
            order,
            s,
            eta_value,
        )

    def previous(s: float, eta_value: float):
        return first_jet(float(s), float(eta_value)).value

    def previous_parameter(s: float, eta_value: float):
        return first_jet(float(s), float(eta_value)).parameter

    def forcing(s: float, eta_value: float):
        return forcing_jet(float(s), float(eta_value)).value

    return picard_map_eq_5_7(
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
