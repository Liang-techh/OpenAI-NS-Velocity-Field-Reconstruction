"""Hierarchy-owned third eta jet of the first Section 5 Picard image.

For positive coefficient order ``n``, Lemma 5.1 writes the Eq. (5.7) solution as

    W_n = sum_{k>=0} K^k G f_n,

with the first image ``W_n^(0)=G f_n`` and

    (G g)_i(xi,eta) = integral_0^xi (s/xi)^c_i g_i(s,eta) ds,
    c = (0,0,2,0,3,1).

The kernel and radial integration limits are independent of ``eta``. Therefore
``partial_eta^r`` commutes with ``G`` for the finite analytic jets used here:

    partial_eta^r W_n^(0) = G(partial_eta^r f_n),  0 <= r <= 3.

The forcing jet is obtained only from
``Section5LowerHistorySixthMixedHierarchy`` through the hierarchy-owned
third-parameter forcing bridge. No caller-supplied SourceJet, repair jet,
forcing derivative table, or independent normalization constant is accepted.

This is Stage-2 formal-structure infrastructure. The radial quadrature is the
same numerical Gauss-Legendre realization already used by the landed Eq. (5.7)
primitive; this module does not claim a converged Picard coefficient, an
all-order hierarchy, or paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import math

import numpy as np

from .background_inner_solver import (
    first_picard_term_eq_5_7,
    singular_inverse_eq_5_7,
)
from .background_repaired_history_forcing_third_parameter import (
    hierarchy_owned_positive_axis_forcing_third_parameter_jet,
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
class PositiveAxisFirstPicardThirdParameterJet:
    """Value and first three analytic eta derivatives of ``W_n^(0)=G f_n``."""

    value: np.ndarray
    parameter: np.ndarray
    parameter2: np.ndarray
    parameter3: np.ndarray

    def __post_init__(self) -> None:
        for name in ("value", "parameter", "parameter2", "parameter3"):
            object.__setattr__(self, name, _vector6(getattr(self, name), name))


def hierarchy_owned_first_picard_third_parameter_jet(
    hierarchy: Section5LowerHistorySixthMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
) -> PositiveAxisFirstPicardThirdParameterJet:
    """Construct ``W_n^(0)=G f_n`` through ``partial_eta^3``.

    Every forcing row is regenerated from the same strong hierarchy at the
    quadrature point. The value row deliberately delegates to the landed
    ``first_picard_term_eq_5_7`` primitive. Higher rows use the same singular
    inverse on the corresponding hierarchy-owned analytic forcing derivative;
    this is exactly differentiation under the Eq. (5.7) radial integral because
    ``G`` has no eta-dependent kernel or integration endpoint.
    """

    if not isinstance(hierarchy, Section5LowerHistorySixthMixedHierarchy):
        raise TypeError("hierarchy must be a Section5LowerHistorySixthMixedHierarchy")
    order = _positive_order(order)
    xi, eta = _point(xi, eta)
    quadrature_points = hierarchy.quadrature_points

    def forcing_row(s: float, eta_value: float, derivative: int) -> np.ndarray:
        jet = hierarchy_owned_positive_axis_forcing_third_parameter_jet(
            hierarchy,
            order,
            s,
            eta_value,
        )
        return (jet.value, jet.parameter, jet.parameter2, jet.parameter3)[derivative]

    value = first_picard_term_eq_5_7(
        order,
        xi,
        eta,
        lambda s, e: forcing_row(s, e, 0),
        quadrature_points=quadrature_points,
    )
    parameter = singular_inverse_eq_5_7(
        xi,
        eta,
        lambda s, e: forcing_row(s, e, 1),
        quadrature_points=quadrature_points,
    )
    parameter2 = singular_inverse_eq_5_7(
        xi,
        eta,
        lambda s, e: forcing_row(s, e, 2),
        quadrature_points=quadrature_points,
    )
    parameter3 = singular_inverse_eq_5_7(
        xi,
        eta,
        lambda s, e: forcing_row(s, e, 3),
        quadrature_points=quadrature_points,
    )
    return PositiveAxisFirstPicardThirdParameterJet(
        value=value,
        parameter=parameter,
        parameter2=parameter2,
        parameter3=parameter3,
    )
