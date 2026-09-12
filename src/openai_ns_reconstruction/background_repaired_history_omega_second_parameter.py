"""Hierarchy-owned second eta jet for the regular Eq. (5.6) Omega/X row.

This bridge accepts only the strong Section-5 lower-history hierarchy whose
coefficient sources own fifth-mixed U data.  Every fourth-mixed beta jet is
derived on demand from that same U source through the landed analytic Eq. (5.2)
adapter, and ordinary axial data are lower projections of the same hierarchy.

Thus ``partial_eta^2(Omega_k/X)`` cannot be paired with caller-supplied beta
tables.  Missing strong data at any order, including the Issue-#1 leading
coefficient, fail closed.  This remains Stage-2 ``formal-structure`` solver
infrastructure and is not a paper-exact coefficient claim.
"""
from __future__ import annotations

from numbers import Integral

from .background_omega_second_parameter_jet import (
    OmegaSecondParameterJet,
    omega_over_x_second_parameter_jet_eq_5_6,
)
from .background_repaired_history_fifth_mixed import (
    Section5LowerHistoryFifthMixedHierarchy,
)


def hierarchy_owned_omega_second_parameter_jet(
    hierarchy: Section5LowerHistoryFifthMixedHierarchy,
    order: int,
    X: float,
    eta: float,
) -> OmegaSecondParameterJet:
    """Return hierarchy-owned value/first/second eta jets of ``Omega_order/X``."""

    if not isinstance(hierarchy, Section5LowerHistoryFifthMixedHierarchy):
        raise TypeError(
            "hierarchy must be a Section5LowerHistoryFifthMixedHierarchy"
        )
    if isinstance(order, bool) or not isinstance(order, Integral) or int(order) < 0:
        raise ValueError("order must be a nonnegative integer")
    order = int(order)
    if order >= len(hierarchy.sources):
        raise ValueError(
            "requested omega order is missing one or more required coefficient sources"
        )

    beta = [
        hierarchy.beta_fourth_mixed_jet(j, X, eta)
        for j in range(order + 1)
    ]
    axial = [
        hierarchy.axial_second_jet(j, X, eta)
        for j in range(order + 1)
    ]
    return omega_over_x_second_parameter_jet_eq_5_6(
        order,
        X,
        eta,
        beta,
        axial,
        h=hierarchy.h,
    )
