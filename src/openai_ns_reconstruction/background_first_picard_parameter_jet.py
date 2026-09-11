"""Analytic eta derivative of the first Lemma-5.1 Picard term.

For the singular inverse in Eq. (5.7),

    (G g)_i(xi, eta) = integral_0^xi (s/xi)^c_i g_i(s, eta) ds,

its kernel is independent of ``eta``.  Therefore, whenever the source has a
controlled parameter derivative,

    partial_eta (G g) = G (partial_eta g).

This module materializes that exact-structure commutation for the genuine
first Picard term ``W_n^(0)=G f_n``.  It deliberately requires an explicit
analytic ``partial_eta f_n`` provider.  The currently landed repaired-history
hierarchy only owns second jets of ``beta=V/X``; differentiating the Eq. (5.6)
``Omega/X`` source one more time needs additional higher jets.  Consequently
this is a fail-closed solver primitive, not a claim that the hierarchy already
owns ``partial_eta f_n`` and not a paper-exact positive-order coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np

from .background_inner_solver import (
    VectorField6,
    first_picard_term_eq_5_7,
    singular_inverse_eq_5_7,
)


def _vector6(value: Sequence[float], name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (6,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite vector of shape (6,)")
    out = out.copy()
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class FirstPicardParameterJet:
    """Value and analytic ``eta`` derivative of ``W_n^(0)=G f_n``."""

    value: np.ndarray
    parameter: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", _vector6(self.value, "value"))
        object.__setattr__(self, "parameter", _vector6(self.parameter, "parameter"))


def first_picard_parameter_jet_eq_5_7(
    order: int,
    xi: float,
    eta: float,
    lower_order_source: VectorField6,
    lower_order_source_parameter: VectorField6,
    *,
    quadrature_points: int = 32,
) -> FirstPicardParameterJet:
    """Return ``(G f_n, partial_eta G f_n)`` by analytic commutation.

    The value is evaluated by the already-landed genuine first Picard routine.
    The derivative is *not* estimated from neighboring eta samples: the caller
    must supply ``partial_eta f_n`` and production applies the same singular
    inverse directly to that derivative.

    This interface is intentionally not named ``...from_hierarchy``.  Until the
    repaired lower-history object owns the higher jets needed to differentiate
    all of Eq. (5.6), accepting a derivative callback is only solver
    infrastructure.  A sampled or fitted callback must not be promoted to a
    paper-exact coefficient.
    """

    if not callable(lower_order_source):
        raise TypeError("lower_order_source must be callable")
    if not callable(lower_order_source_parameter):
        raise TypeError("lower_order_source_parameter must be callable")

    value = first_picard_term_eq_5_7(
        order,
        xi,
        eta,
        lower_order_source,
        quadrature_points=quadrature_points,
    )
    parameter = singular_inverse_eq_5_7(
        xi,
        eta,
        lower_order_source_parameter,
        quadrature_points=quadrature_points,
    )
    if not np.all(np.isfinite(value)) or not np.all(np.isfinite(parameter)):
        raise OverflowError("first-Picard parameter jet is outside binary64 range")
    return FirstPicardParameterJet(value=value, parameter=parameter)
