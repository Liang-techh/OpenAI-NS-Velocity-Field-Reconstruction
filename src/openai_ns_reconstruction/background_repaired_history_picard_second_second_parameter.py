"""Hierarchy-owned eta^2 jet of the second Section 5 Picard-series term.

For positive coefficient order ``n``, Lemma 5.1 writes

    W_n = sum_{k>=0} K^k G f_n,
    K W = G(A0 W + A1 partial_eta W).

The landed hierarchy path already owns the first term
``W_n^(0)=G f_n`` through ``partial_eta^3`` and owns ``A0/A1`` through
``partial_eta^2``.  This module closes the next finite derivative seam by
constructing

    W_n^(1) = K W_n^(0)

through ``partial_eta^2``.  Since the singular inverse ``G`` has an
eta-independent kernel and radial integration limits, eta derivatives commute
with ``G``.  The right-hand side is therefore differentiated only by exact
Leibniz rules.

No caller-supplied SourceJet, repair jet, forcing derivative table, matrix
jet, or independent normalization data are accepted.  Production uses no
finite differences, fitting, generic cutoff, or sampled ``V/X`` division.

This remains Stage-2 ``formal-structure`` infrastructure.  It is one finite
Picard-series term, not a converged positive-order coefficient, an all-order
hierarchy, or paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import math

import numpy as np

from .background_inner_solver import singular_inverse_eq_5_7
from .background_repaired_history_picard_first_third_parameter import (
    hierarchy_owned_first_picard_third_parameter_jet,
)
from .background_repaired_history_positive_axis_matrix_jets import (
    hierarchy_owned_positive_axis_matrix_second_parameter_jet,
)
from .background_repaired_history_sixth_mixed import (
    Section5LowerHistorySixthMixedHierarchy,
)


def _positive_order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("order must be a positive integer")
    return int(value)


def _point(xi: float, eta: float) -> tuple[float, float]:
    xi = float(xi)
    eta = float(eta)
    if not math.isfinite(xi) or xi < 0.0:
        raise ValueError("xi must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return xi, eta


def _vector6(value: np.ndarray, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (6,):
        raise ValueError(f"{name} must have shape (6,)")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be finite")
    out = out.copy()
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class PositiveAxisSecondPicardSecondParameterJet:
    """Value and first two eta derivatives of ``W_n^(1)=K W_n^(0)``."""

    value: np.ndarray
    parameter: np.ndarray
    parameter2: np.ndarray

    def __post_init__(self) -> None:
        for name in ("value", "parameter", "parameter2"):
            object.__setattr__(self, name, _vector6(getattr(self, name), name))


def hierarchy_owned_second_picard_second_parameter_jet(
    hierarchy: Section5LowerHistorySixthMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
) -> PositiveAxisSecondPicardSecondParameterJet:
    """Construct ``W_n^(1)=K W_n^(0)`` through ``partial_eta^2``.

    At each radial quadrature point this regenerates both the first Picard jet
    and the PositiveAxis matrix jet from the same strong hierarchy.  Writing

    ``R = A0 W0 + A1 W0_eta``, exact differentiation gives

    ``R_eta = A0_eta W0 + A0 W0_eta
              + A1_eta W0_eta + A1 W0_etaeta``

    and

    ``R_etaeta = A0_etaeta W0 + 2 A0_eta W0_eta + A0 W0_etaeta
                 + A1_etaeta W0_eta + 2 A1_eta W0_etaeta
                 + A1 W0_etaetaeta``.

    The three returned rows are ``G R``, ``G R_eta`` and ``G R_etaeta``.
    There is deliberately no additional ``f_n`` term: ``G f_n`` is already
    the separate ``W_n^(0)`` term of the Picard series.
    """

    if not isinstance(hierarchy, Section5LowerHistorySixthMixedHierarchy):
        raise TypeError("hierarchy must be a Section5LowerHistorySixthMixedHierarchy")
    order = _positive_order(order)
    xi, eta = _point(xi, eta)
    quadrature_points = hierarchy.quadrature_points

    def rhs_row(s: float, eta_value: float, derivative: int) -> np.ndarray:
        first = hierarchy_owned_first_picard_third_parameter_jet(
            hierarchy,
            order,
            s,
            eta_value,
        )
        matrices = hierarchy_owned_positive_axis_matrix_second_parameter_jet(
            hierarchy,
            order,
            s,
            eta_value,
        )

        if derivative == 0:
            value = matrices.A0 @ first.value + matrices.A1 @ first.parameter
        elif derivative == 1:
            value = (
                matrices.A0_parameter @ first.value
                + matrices.A0 @ first.parameter
                + matrices.A1_parameter @ first.parameter
                + matrices.A1 @ first.parameter2
            )
        elif derivative == 2:
            value = (
                matrices.A0_parameter2 @ first.value
                + 2.0 * (matrices.A0_parameter @ first.parameter)
                + matrices.A0 @ first.parameter2
                + matrices.A1_parameter2 @ first.parameter
                + 2.0 * (matrices.A1_parameter @ first.parameter2)
                + matrices.A1 @ first.parameter3
            )
        else:
            raise ValueError("derivative must be 0, 1, or 2")
        return _vector6(value, "second Picard right-hand side")

    value = singular_inverse_eq_5_7(
        xi,
        eta,
        lambda s, e: rhs_row(s, e, 0),
        quadrature_points=quadrature_points,
    )
    parameter = singular_inverse_eq_5_7(
        xi,
        eta,
        lambda s, e: rhs_row(s, e, 1),
        quadrature_points=quadrature_points,
    )
    parameter2 = singular_inverse_eq_5_7(
        xi,
        eta,
        lambda s, e: rhs_row(s, e, 2),
        quadrature_points=quadrature_points,
    )
    return PositiveAxisSecondPicardSecondParameterJet(
        value=value,
        parameter=parameter,
        parameter2=parameter2,
    )
