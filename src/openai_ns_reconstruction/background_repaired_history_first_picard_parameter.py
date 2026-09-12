"""Hierarchy-owned first Picard parameter jet for Section 5 Eq. (5.7).

PR #174 landed a fail-closed analytic hierarchy-owned ``(f_n, partial_eta f_n)``.
PR #123 had already materialized the exact structural commutation

    partial_eta (G f_n) = G (partial_eta f_n)

for the singular diagonal inverse in Eq. (5.7).  This module composes those two
pieces without accepting an independent source or derivative callback.

The result is the genuine hierarchy-owned ``k=0`` Picard parameter jet
``(G f_n, partial_eta G f_n)`` whenever the strong repaired strict-lower
history exists.  It is still Stage-2 ``formal-structure`` infrastructure: it is
not the nontrivial ``k=1`` application, not a converged positive-order
coefficient, and not a paper-exact background/velocity.
"""
from __future__ import annotations

from .background_first_picard_parameter_jet import (
    FirstPicardParameterJet,
    first_picard_parameter_jet_eq_5_7,
)
from .background_repaired_history_forcing_parameter import (
    hierarchy_owned_positive_axis_forcing_parameter_jet,
)
from .background_repaired_history_phi_third_mixed import (
    Section5LowerHistoryPhiThirdMixedHierarchy,
)


def hierarchy_owned_first_picard_parameter_jet_eq_5_7(
    hierarchy: Section5LowerHistoryPhiThirdMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
    *,
    quadrature_points: int = 32,
) -> FirstPicardParameterJet:
    """Return hierarchy-owned ``(G f_n, partial_eta G f_n)``.

    The only source admitted by this interface is the strong repaired
    lower-history hierarchy.  At each radial quadrature point the value and
    analytic eta derivative are obtained from
    ``hierarchy_owned_positive_axis_forcing_parameter_jet`` and passed to the
    already-landed singular-inverse parameter-jet primitive.

    A zero-radius request is preflighted through the forcing bridge before the
    singular inverse returns its exact zero.  This preserves fail-closed
    behavior for missing strong order-zero data instead of allowing the axis
    shortcut to hide an incomplete Issue-#1 dependency.
    """

    if not isinstance(hierarchy, Section5LowerHistoryPhiThirdMixedHierarchy):
        raise TypeError(
            "hierarchy must be a Section5LowerHistoryPhiThirdMixedHierarchy"
        )

    # Preserve the hierarchy truth boundary even though G(0)=0 can otherwise
    # short-circuit without evaluating its source.
    hierarchy_owned_positive_axis_forcing_parameter_jet(
        hierarchy,
        order,
        0.0,
        eta,
    )

    def forcing_value(s: float, eta_value: float):
        return hierarchy_owned_positive_axis_forcing_parameter_jet(
            hierarchy,
            order,
            s,
            eta_value,
        ).value

    def forcing_parameter(s: float, eta_value: float):
        return hierarchy_owned_positive_axis_forcing_parameter_jet(
            hierarchy,
            order,
            s,
            eta_value,
        ).parameter

    return first_picard_parameter_jet_eq_5_7(
        order,
        xi,
        eta,
        forcing_value,
        forcing_parameter,
        quadrature_points=quadrature_points,
    )
