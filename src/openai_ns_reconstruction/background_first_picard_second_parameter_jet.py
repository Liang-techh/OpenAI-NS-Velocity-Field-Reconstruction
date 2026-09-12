"""Second analytic eta derivative of the first Section-5 Picard term.

For the singular inverse in Eq. (5.7),

    (G g)_i(xi, eta) = integral_0^xi (s/xi)^c_i g_i(s, eta) ds,

its kernel is independent of ``eta``.  Hence, whenever the source carries a
controlled second parameter derivative,

    partial_eta^2 (G g) = G(partial_eta^2 g).

This module extends the landed first-Picard parameter jet by exactly one eta
order.  The second source derivative remains an explicit analytic provider;
this file does not infer it from samples and does not claim that the repaired
history already owns the higher mixed jets needed to construct it.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .background_first_picard_parameter_jet import (
    first_picard_parameter_jet_eq_5_7,
)
from .background_inner_solver import VectorField6, singular_inverse_eq_5_7


def _vector6(value: Sequence[float], name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (6,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite vector of shape (6,)")
    out = out.copy()
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class FirstPicardSecondParameterJet:
    """Value, first eta derivative, and second eta derivative of ``G f_n``."""

    value: np.ndarray
    parameter: np.ndarray
    second_parameter: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", _vector6(self.value, "value"))
        object.__setattr__(self, "parameter", _vector6(self.parameter, "parameter"))
        object.__setattr__(
            self,
            "second_parameter",
            _vector6(self.second_parameter, "second_parameter"),
        )


def first_picard_second_parameter_jet_eq_5_7(
    order: int,
    xi: float,
    eta: float,
    lower_order_source: VectorField6,
    lower_order_source_parameter: VectorField6,
    lower_order_source_second_parameter: VectorField6,
    *,
    quadrature_points: int = 32,
) -> FirstPicardSecondParameterJet:
    """Return ``(G f_n, partial_eta G f_n, partial_eta^2 G f_n)``.

    The value and first derivative reuse the already-landed first-Picard
    parameter bridge.  The second derivative is obtained by applying the same
    singular inverse directly to an explicit analytic ``partial_eta^2 f_n``
    provider.  Production performs no eta finite difference or fitting.

    This is intentionally a solver primitive rather than a hierarchy-owned
    constructor.  A caller-supplied second derivative is not paper-exact
    provenance; missing higher repaired-history jets must remain visible.
    """

    if not callable(lower_order_source_second_parameter):
        raise TypeError("lower_order_source_second_parameter must be callable")

    first = first_picard_parameter_jet_eq_5_7(
        order,
        xi,
        eta,
        lower_order_source,
        lower_order_source_parameter,
        quadrature_points=quadrature_points,
    )
    second_parameter = singular_inverse_eq_5_7(
        xi,
        eta,
        lower_order_source_second_parameter,
        quadrature_points=quadrature_points,
    )
    if not np.all(np.isfinite(second_parameter)):
        raise OverflowError("first-Picard second parameter jet is outside binary64 range")

    return FirstPicardSecondParameterJet(
        value=first.value,
        parameter=first.parameter,
        second_parameter=second_parameter,
    )
