"""Hierarchy-owned eta jet of the third Section 5 Picard-series term.

For positive coefficient order ``n``, Lemma 5.1 writes

    W_n = sum_{k>=0} K^k G f_n,
    K W = G(A0 W + A1 partial_eta W).

The landed hierarchy path already owns ``W_n^(1)=K W_n^(0)`` through
``partial_eta^2`` and owns ``A0/A1`` through ``partial_eta^2``.  This module
takes exactly one further Picard step and constructs

    W_n^(2) = K W_n^(1)

through ``partial_eta``.  Since the singular inverse ``G`` has an
eta-independent kernel and radial integration limits, eta derivatives commute
with ``G`` and the derivative row follows the exact first Leibniz rule.

No caller-supplied SourceJet, repair jet, forcing derivative table, matrix jet,
or independent normalization data are accepted.  Production uses no finite
differences, fitting, generic cutoff, or sampled ``V/X`` division.

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
from .background_repaired_history_picard_second_second_parameter import (
    hierarchy_owned_second_picard_second_parameter_jet,
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
class PositiveAxisThirdPicardFirstParameterJet:
    """Value and first eta derivative of ``W_n^(2)=K W_n^(1)``."""

    value: np.ndarray
    parameter: np.ndarray

    def __post_init__(self) -> None:
        for name in ("value", "parameter"):
            object.__setattr__(self, name, _vector6(getattr(self, name), name))


def hierarchy_owned_third_picard_first_parameter_jet(
    hierarchy: Section5LowerHistorySixthMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
) -> PositiveAxisThirdPicardFirstParameterJet:
    """Construct ``W_n^(2)=K W_n^(1)`` through ``partial_eta``.

    At each radial quadrature point this regenerates both the hierarchy-owned
    second Picard jet and the hierarchy-owned PositiveAxis matrix jet.  Writing

    ``R = A0 W1 + A1 W1_eta``,

    exact differentiation gives

    ``R_eta = A0_eta W1 + A0 W1_eta
              + A1_eta W1_eta + A1 W1_etaeta``.

    The returned rows are ``G R`` and ``G R_eta``.  There is deliberately no
    additional ``f_n`` term because ``G f_n`` is the separate ``W_n^(0)`` term.
    """

    if not isinstance(hierarchy, Section5LowerHistorySixthMixedHierarchy):
        raise TypeError("hierarchy must be a Section5LowerHistorySixthMixedHierarchy")
    order = _positive_order(order)
    xi, eta = _point(xi, eta)
    quadrature_points = hierarchy.quadrature_points

    jet_cache: dict[tuple[float, float], tuple[object, object]] = {}

    def upstream(s: float, eta_value: float) -> tuple[object, object]:
        key = (float(s), float(eta_value))
        cached = jet_cache.get(key)
        if cached is None:
            cached = (
                hierarchy_owned_second_picard_second_parameter_jet(
                    hierarchy,
                    order,
                    s,
                    eta_value,
                ),
                hierarchy_owned_positive_axis_matrix_second_parameter_jet(
                    hierarchy,
                    order,
                    s,
                    eta_value,
                ),
            )
            jet_cache[key] = cached
        return cached

    def rhs_row(s: float, eta_value: float, derivative: int) -> np.ndarray:
        second, matrices = upstream(s, eta_value)
        if derivative == 0:
            value = matrices.A0 @ second.value + matrices.A1 @ second.parameter
        elif derivative == 1:
            value = (
                matrices.A0_parameter @ second.value
                + matrices.A0 @ second.parameter
                + matrices.A1_parameter @ second.parameter
                + matrices.A1 @ second.parameter2
            )
        else:
            raise ValueError("derivative must be 0 or 1")
        return _vector6(value, "third Picard right-hand side")

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
    return PositiveAxisThirdPicardFirstParameterJet(
        value=value,
        parameter=parameter,
    )
