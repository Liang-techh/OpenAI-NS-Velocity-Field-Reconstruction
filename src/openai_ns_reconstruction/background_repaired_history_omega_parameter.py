"""Hierarchy-owned analytic parameter jet for the regular Eq. (5.6) source.

The landed :mod:`background_repaired_history` hierarchy can own coherent
fourth-mixed axial jets and derive the corresponding third-mixed regular-flux
jets only through Eq. (5.2).  The landed
:func:`background_omega_parameter_jet.omega_over_x_parameter_jet_eq_5_6`
then gives the exact analytic eta derivative of ``Omega_k / X``.

This module is the narrow ownership bridge between those two pieces.  It accepts
only a ``Section5LowerHistoryJetHierarchy`` and an order ``k``; no independent
caller-supplied beta or axial derivative tables are accepted.  If any coefficient
``0,...,k`` lacks the stronger fourth-mixed U layer, the call fails closed.

Issue #1 still owns the real leading profile, so this is solver infrastructure /
``formal-structure`` only.  In particular it does not claim that the leading
fourth-mixed data are paper-exact, that ``partial_eta f_n`` is materialized, or
that the Eq. (5.7) Picard series has advanced beyond its first term.
"""
from __future__ import annotations

from numbers import Integral

from .background_omega_parameter_jet import (
    OmegaParameterJet,
    omega_over_x_parameter_jet_eq_5_6,
)
from .background_repaired_history import Section5LowerHistoryJetHierarchy


def hierarchy_owned_omega_parameter_jet(
    hierarchy: Section5LowerHistoryJetHierarchy,
    order: int,
    X: float,
    eta: float,
) -> OmegaParameterJet:
    """Return hierarchy-owned ``(Omega_order/X, partial_eta(Omega_order/X))``.

    Every regular-flux third-mixed jet is derived on demand from the same
    hierarchy-owned fourth-mixed U source through the analytic Eq. (5.2)
    adapter.  Ordinary axial second jets are projected from the same coherent
    U source.  Therefore the Eq. (5.6) parameter jet cannot be paired with an
    unrelated caller-provided beta derivative record.

    This routine deliberately requires all coefficient orders ``0,...,order``
    to be present and to own fourth-mixed U data.  That includes order zero;
    until Issue #1 exposes the real leading fourth-mixed profile, production use
    of this bridge must remain blocked rather than inventing that data.
    """

    if not isinstance(hierarchy, Section5LowerHistoryJetHierarchy):
        raise TypeError("hierarchy must be a Section5LowerHistoryJetHierarchy")
    if isinstance(order, bool) or not isinstance(order, Integral) or int(order) < 0:
        raise ValueError("order must be a nonnegative integer")
    order = int(order)
    if order >= len(hierarchy.sources):
        raise ValueError(
            "requested omega order is missing one or more required coefficient sources"
        )

    beta = [
        hierarchy.beta_third_mixed_jet(j, X, eta)
        for j in range(order + 1)
    ]
    axial = [
        hierarchy.axial_second_jet(j, X, eta)
        for j in range(order + 1)
    ]
    return omega_over_x_parameter_jet_eq_5_6(
        order,
        X,
        eta,
        beta,
        axial,
        h=hierarchy.h,
    )
