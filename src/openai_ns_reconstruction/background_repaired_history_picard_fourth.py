"""Hierarchy-owned value of the fourth Section 5 Picard-series term.

For positive coefficient order ``n``, Lemma 5.1 writes

    W_n = sum_{k>=0} K^k G f_n,
    K W = G(A0 W + A1 partial_eta W).

The landed hierarchy path already owns ``W_n^(2)=K W_n^(1)`` through
``partial_eta``.  This module takes exactly one further Picard step and
constructs

    W_n^(3) = K W_n^(2)

at value level.  This is the last row of the finite derivative triangle
currently supported by the hierarchy-owned eta jets

    W_n^(0): eta^3 -> W_n^(1): eta^2 -> W_n^(2): eta -> W_n^(3): value.

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
from .background_repaired_history_picard_third_first_parameter import (
    hierarchy_owned_third_picard_first_parameter_jet,
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
class PositiveAxisFourthPicardValue:
    """Value of ``W_n^(3)=K W_n^(2)``."""

    value: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", _vector6(self.value, "value"))


def hierarchy_owned_fourth_picard_value(
    hierarchy: Section5LowerHistorySixthMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
) -> PositiveAxisFourthPicardValue:
    """Construct the value row ``W_n^(3)=K W_n^(2)``.

    At every radial quadrature point this regenerates ``W_n^(2)`` and its
    first eta derivative from the same strong hierarchy and regenerates the
    PositiveAxis matrices from coefficient-order-zero hierarchy data.  Thus

    ``W_n^(3) = G(A0 W_n^(2) + A1 partial_eta W_n^(2))``.

    There is deliberately no additional ``f_n`` term because ``G f_n`` is the
    separate ``W_n^(0)`` term of the Lemma 5.1 Picard series.
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
                hierarchy_owned_third_picard_first_parameter_jet(
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

    def rhs(s: float, eta_value: float) -> np.ndarray:
        third, matrices = upstream(s, eta_value)
        return _vector6(
            matrices.A0 @ third.value + matrices.A1 @ third.parameter,
            "fourth Picard right-hand side",
        )

    value = singular_inverse_eq_5_7(
        xi,
        eta,
        rhs,
        quadrature_points=quadrature_points,
    )
    return PositiveAxisFourthPicardValue(value=value)
