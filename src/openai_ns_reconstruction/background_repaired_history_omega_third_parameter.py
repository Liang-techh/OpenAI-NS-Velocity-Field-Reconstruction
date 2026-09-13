"""Hierarchy-owned third eta jet for the regular Eq. (5.6) ``Omega/X`` row.

Only the strong Section-5 lower-history hierarchy with authoritative repaired
sixth-mixed U data is accepted. Fifth-mixed beta jets are derived on demand from
that same hierarchy through analytic Eq. (5.2), and axial third eta data are
exact projections of the hierarchy-owned U jets. Missing strong history fails
closed. This remains Stage-2 ``formal-structure`` infrastructure.
"""
from __future__ import annotations

from numbers import Integral

from .background_omega_third_parameter_jet import (
    OmegaThirdParameterJet,
    omega_over_x_third_parameter_jet_eq_5_6,
)
from .background_repaired_history_sixth_mixed import (
    Section5LowerHistorySixthMixedHierarchy,
)


def hierarchy_owned_omega_third_parameter_jet(
    hierarchy: Section5LowerHistorySixthMixedHierarchy,
    order: int,
    X: float,
    eta: float,
) -> OmegaThirdParameterJet:
    """Return hierarchy-owned value through third eta derivative of ``Omega/X``."""
    if not isinstance(hierarchy, Section5LowerHistorySixthMixedHierarchy):
        raise TypeError("hierarchy must be a Section5LowerHistorySixthMixedHierarchy")
    if isinstance(order, bool) or not isinstance(order, Integral) or int(order) < 0:
        raise ValueError("order must be a nonnegative integer")
    order = int(order)
    if order >= len(hierarchy.sources):
        raise ValueError(
            "requested omega order is missing one or more required coefficient sources"
        )

    beta = [hierarchy.beta_fifth_mixed_jet(j, X, eta) for j in range(order + 1)]
    axial = [
        hierarchy.axial_sixth_mixed_jet(j, X, eta).fifth().fourth().third()
        for j in range(order + 1)
    ]
    return omega_over_x_third_parameter_jet_eq_5_6(
        order,
        X,
        eta,
        beta,
        axial,
        h=hierarchy.h,
    )
