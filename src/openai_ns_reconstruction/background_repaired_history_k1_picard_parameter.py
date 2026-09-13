"""Hierarchy-owned analytic eta jet of the nontrivial Eq. (5.7) ``k=1`` iterate.

The landed Section-5 stack now owns every derivative needed to differentiate

    W_n^(1) = G(A0 W_n^(0) + A1 partial_eta W_n^(0) + f_n)

at fixed ``xi``.  In particular, the strong repaired hierarchy owns
``(f_n, partial_eta f_n, partial_eta^2 f_n)``; the first Picard primitive
commutes the eta derivatives through the eta-independent singular inverse; and
the Eq. (5.7) matrix bridge owns ``(A0, partial_eta A0)`` and
``(A1, partial_eta A1)``.

This module composes only those landed analytic layers.  Production accepts no
caller-supplied matrix, forcing, iterate derivative, finite-difference
surrogate, generic cutoff, or fitted coefficient family.  Genuine order-zero
fourth/fifth-mixed profile data remain an Issue-#1 dependency, so missing strong
leading data fail closed.  The result is still Stage-2 ``formal-structure``:
one differentiated Picard iterate is not Picard convergence or a paper-exact
recursive coefficient/velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Sequence

import numpy as np

from .background_first_picard_second_parameter_jet import (
    first_picard_second_parameter_jet_eq_5_7,
)
from .background_inner_solver import singular_inverse_eq_5_7
from .background_repaired_history_forcing_second_parameter import (
    hierarchy_owned_positive_axis_forcing_second_parameter_jet,
)
from .background_repaired_history_k1_picard import (
    hierarchy_owned_k1_picard_value_eq_5_7,
)
from .background_repaired_history_matrix_parameter import (
    hierarchy_owned_positive_axis_matrix_parameter_jets,
)
from .background_repaired_history_phi_fourth_mixed import (
    Section5LowerHistoryPhiFourthMixedHierarchy,
)


def _vector6(value: Sequence[float], name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (6,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite vector of shape (6,)")
    out = out.copy()
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class K1PicardParameterJet:
    """Value and analytic eta derivative of the hierarchy-owned ``W_n^(1)``."""

    value: np.ndarray
    parameter: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", _vector6(self.value, "value"))
        object.__setattr__(self, "parameter", _vector6(self.parameter, "parameter"))


def hierarchy_owned_k1_picard_parameter_jet_eq_5_7(
    hierarchy: Section5LowerHistoryPhiFourthMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
    *,
    quadrature_points: int = 32,
) -> K1PicardParameterJet:
    """Return hierarchy-owned ``(W_n^(1), partial_eta W_n^(1))``.

    Writing ``W0 = W_n^(0)``, the derivative of the Eq. (5.7) Picard right-hand
    side is assembled analytically as

    ``A0_eta W0 + A0 W0_eta + A1_eta W0_eta + A1 W0_etaeta + f_eta``.

    The singular inverse ``G`` is independent of ``eta``, so the parameter
    component is the same canonical ``G`` applied to that derivative.  All
    ingredients are obtained from the supplied strong repaired hierarchy:

    - matrix values/derivatives from the hierarchy-owned matrix parameter jet;
    - ``W0``, ``W0_eta``, ``W0_etaeta`` from the first-Picard second-parameter
      primitive driven only by hierarchy-owned forcing jets;
    - ``f_n`` derivatives from the hierarchy-owned forcing second-parameter
      bridge.

    The value is delegated unchanged to the already-landed hierarchy-owned
    ``k=1`` value path.  An explicit axis preflight evaluates the new strong
    forcing/matrix dependencies before ``G(0)=0`` can short-circuit, preserving
    the Issue-#1 fail-closed boundary.
    """

    if not isinstance(hierarchy, Section5LowerHistoryPhiFourthMixedHierarchy):
        raise TypeError(
            "hierarchy must be a Section5LowerHistoryPhiFourthMixedHierarchy"
        )

    @lru_cache(maxsize=None)
    def forcing_jet(s: float, eta_value: float):
        return hierarchy_owned_positive_axis_forcing_second_parameter_jet(
            hierarchy,
            order,
            float(s),
            float(eta_value),
        )

    @lru_cache(maxsize=None)
    def matrix_jets(s: float, eta_value: float):
        return hierarchy_owned_positive_axis_matrix_parameter_jets(
            hierarchy,
            order,
            float(s),
            float(eta_value),
        )

    # Preserve the stronger truth boundary at the axis.  The singular inverse
    # returns exact zero there without evaluating its source, so the new strong
    # dependencies must be checked explicitly first.
    forcing_jet(0.0, float(eta))
    matrix_jets(0.0, float(eta))

    value = hierarchy_owned_k1_picard_value_eq_5_7(
        hierarchy,
        order,
        xi,
        eta,
        quadrature_points=quadrature_points,
    )

    def forcing_value(s: float, eta_value: float):
        return forcing_jet(float(s), float(eta_value)).value

    def forcing_parameter(s: float, eta_value: float):
        return forcing_jet(float(s), float(eta_value)).parameter

    def forcing_second_parameter(s: float, eta_value: float):
        return forcing_jet(float(s), float(eta_value)).parameter2

    @lru_cache(maxsize=None)
    def first_picard_jet(s: float, eta_value: float):
        return first_picard_second_parameter_jet_eq_5_7(
            order,
            float(s),
            float(eta_value),
            forcing_value,
            forcing_parameter,
            forcing_second_parameter,
            quadrature_points=quadrature_points,
        )

    def rhs_parameter(s: float, eta_value: float) -> np.ndarray:
        s = float(s)
        eta_value = float(eta_value)
        matrices = matrix_jets(s, eta_value)
        previous = first_picard_jet(s, eta_value)
        forcing = forcing_jet(s, eta_value)
        derivative = (
            matrices.A0.parameter @ previous.value
            + matrices.A0.value @ previous.parameter
            + matrices.A1.parameter @ previous.parameter
            + matrices.A1.value @ previous.second_parameter
            + forcing.parameter
        )
        if derivative.shape != (6,) or not np.all(np.isfinite(derivative)):
            raise OverflowError("Eq. (5.7) k=1 parameter right-hand side is nonfinite")
        return derivative

    parameter = singular_inverse_eq_5_7(
        xi,
        eta,
        rhs_parameter,
        quadrature_points=quadrature_points,
    )
    if not np.all(np.isfinite(parameter)):
        raise OverflowError("Eq. (5.7) k=1 parameter jet is outside binary64 range")

    return K1PicardParameterJet(value=value, parameter=parameter)
