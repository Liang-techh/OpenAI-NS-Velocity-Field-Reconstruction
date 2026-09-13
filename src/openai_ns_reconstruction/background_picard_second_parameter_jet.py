"""Exact second-eta differentiation primitive for the Eq. (5.7) Picard map.

For one positive-order Section-5 Picard update

    W^+ = G(A0 W + A1 partial_eta W + f),

exact differentiation through the eta-independent singular inverse ``G`` gives

    partial_eta^2 W^+ = G(
        A0_ee W + 2 A0_e W_e + A0 W_ee
        + A1_ee W_e + 2 A1_e W_ee + A1 W_eee
        + f_ee
    ).

The final ``A1 W_eee`` term is important: second-order matrix and forcing jets
alone do *not* determine the second eta derivative of the next Picard iterate.
This module makes that derivative debt explicit with a structural provider
interface rather than filling it with finite differences or fitted data.

This is reusable solver infrastructure only.  The providers may be backed by
analytic test functions or by a future hierarchy-owned stronger jet, but a
caller-supplied provider is not thereby paper-exact.  The module contains no
cutoff construction, numerical eta differentiation, or coefficient fitting.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from numbers import Integral
from typing import Callable, Protocol, Sequence
import math

import numpy as np

from .background_inner_solver import singular_inverse_eq_5_7


class MatrixSecondParameterJetLike(Protocol):
    """Structural value/first/second eta jet for one 6x6 matrix."""

    value: Sequence[Sequence[float]]
    parameter: Sequence[Sequence[float]]
    parameter2: Sequence[Sequence[float]]


class MatrixSecondParameterJetsLike(Protocol):
    """Structural second-eta jets for the Eq. (5.7) pair ``(A0,A1)``."""

    A0: MatrixSecondParameterJetLike
    A1: MatrixSecondParameterJetLike


class VectorSecondParameterJetLike(Protocol):
    """Structural value/first/second eta jet for a six-vector forcing."""

    value: Sequence[float]
    parameter: Sequence[float]
    parameter2: Sequence[float]


class VectorThirdParameterJetLike(Protocol):
    """Structural value through third eta derivative for a six-vector iterate."""

    value: Sequence[float]
    parameter: Sequence[float]
    parameter2: Sequence[float]
    parameter3: Sequence[float]


MatrixJetProvider = Callable[[float, float], MatrixSecondParameterJetsLike]
PreviousJetProvider = Callable[[float, float], VectorThirdParameterJetLike]
ForcingJetProvider = Callable[[float, float], VectorSecondParameterJetLike]


def _vector6(value: Sequence[float], name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (6,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite vector of shape (6,)")
    out = out.copy()
    out.setflags(write=False)
    return out


def _matrix6(value: Sequence[Sequence[float]], name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (6, 6) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite matrix of shape (6, 6)")
    out = out.copy()
    out.setflags(write=False)
    return out


def _positive_order(order: int) -> int:
    if isinstance(order, bool) or not isinstance(order, Integral) or order < 1:
        raise ValueError("coefficient order must be a positive integer")
    return int(order)


def _coordinates(xi: float, eta: float) -> tuple[float, float]:
    xi = float(xi)
    eta = float(eta)
    if not math.isfinite(xi) or xi < 0.0:
        raise ValueError("xi must be finite and nonnegative")
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite with |eta| <= 1")
    return xi, eta


@dataclass(frozen=True)
class PicardSecondParameterJet:
    """Value and first two exact eta derivatives of one Picard update."""

    value: np.ndarray
    parameter: np.ndarray
    parameter2: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", _vector6(self.value, "value"))
        object.__setattr__(self, "parameter", _vector6(self.parameter, "parameter"))
        object.__setattr__(self, "parameter2", _vector6(self.parameter2, "parameter2"))


@dataclass(frozen=True)
class _RhsSecondParameterJet:
    value: np.ndarray
    parameter: np.ndarray
    parameter2: np.ndarray



def picard_second_parameter_jet_eq_5_7(
    order: int,
    xi: float,
    eta: float,
    matrix_jets: MatrixJetProvider,
    previous_jet: PreviousJetProvider,
    forcing_jet: ForcingJetProvider,
    *,
    quadrature_points: int = 32,
) -> PicardSecondParameterJet:
    """Apply one Eq. (5.7) Picard map and return its second eta jet.

    ``matrix_jets`` must provide ``A0`` and ``A1`` through two eta derivatives;
    ``forcing_jet`` must provide ``f`` through two eta derivatives; and
    ``previous_jet`` must provide the previous iterate through *three* eta
    derivatives.  Requiring the third previous-iterate derivative is the exact
    consequence of differentiating ``A1 partial_eta W`` twice.

    The three right-hand-side jets are assembled analytically and the canonical
    singular inverse ``G`` is applied to each.  No eta finite difference occurs
    in production.  Passing analytic caller data is permitted for solver tests,
    but this interface deliberately makes no claim that such data are the
    manuscript's materialized coefficient hierarchy.
    """

    _positive_order(order)
    xi, eta = _coordinates(xi, eta)
    for provider, name in (
        (matrix_jets, "matrix_jets"),
        (previous_jet, "previous_jet"),
        (forcing_jet, "forcing_jet"),
    ):
        if not callable(provider):
            raise TypeError(f"{name} must be callable")

    @lru_cache(maxsize=None)
    def rhs_jet(s: float, eta_value: float) -> _RhsSecondParameterJet:
        s = float(s)
        eta_value = float(eta_value)
        matrices = matrix_jets(s, eta_value)
        previous = previous_jet(s, eta_value)
        forcing = forcing_jet(s, eta_value)

        A0 = _matrix6(matrices.A0.value, "A0.value")
        A0e = _matrix6(matrices.A0.parameter, "A0.parameter")
        A0ee = _matrix6(matrices.A0.parameter2, "A0.parameter2")
        A1 = _matrix6(matrices.A1.value, "A1.value")
        A1e = _matrix6(matrices.A1.parameter, "A1.parameter")
        A1ee = _matrix6(matrices.A1.parameter2, "A1.parameter2")

        W = _vector6(previous.value, "previous.value")
        We = _vector6(previous.parameter, "previous.parameter")
        Wee = _vector6(previous.parameter2, "previous.parameter2")
        Weee = _vector6(previous.parameter3, "previous.parameter3")
        f = _vector6(forcing.value, "forcing.value")
        fe = _vector6(forcing.parameter, "forcing.parameter")
        fee = _vector6(forcing.parameter2, "forcing.parameter2")

        value = A0 @ W + A1 @ We + f
        parameter = A0e @ W + A0 @ We + A1e @ We + A1 @ Wee + fe
        parameter2 = (
            A0ee @ W
            + 2.0 * A0e @ We
            + A0 @ Wee
            + A1ee @ We
            + 2.0 * A1e @ Wee
            + A1 @ Weee
            + fee
        )
        for assembled, name in (
            (value, "Picard right-hand side"),
            (parameter, "Picard first-eta right-hand side"),
            (parameter2, "Picard second-eta right-hand side"),
        ):
            if not np.all(np.isfinite(assembled)):
                raise OverflowError(f"{name} is outside binary64 range")
        return _RhsSecondParameterJet(value, parameter, parameter2)

    # ``G`` returns exact zero at the axis without touching its source.  Run one
    # structural preflight so malformed analytic providers do not get hidden by
    # that short-circuit.
    rhs_jet(0.0, eta)

    value = singular_inverse_eq_5_7(
        xi,
        eta,
        lambda s, e: rhs_jet(float(s), float(e)).value,
        quadrature_points=quadrature_points,
    )
    parameter = singular_inverse_eq_5_7(
        xi,
        eta,
        lambda s, e: rhs_jet(float(s), float(e)).parameter,
        quadrature_points=quadrature_points,
    )
    parameter2 = singular_inverse_eq_5_7(
        xi,
        eta,
        lambda s, e: rhs_jet(float(s), float(e)).parameter2,
        quadrature_points=quadrature_points,
    )
    return PicardSecondParameterJet(value, parameter, parameter2)
