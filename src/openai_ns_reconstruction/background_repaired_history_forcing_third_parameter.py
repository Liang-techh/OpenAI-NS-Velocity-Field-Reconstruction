"""Hierarchy-owned third eta jet of the Section 5 PositiveAxis forcing ``f_n``.

The strict-lower source derivative chain is owned by
``Section5LowerHistorySixthMixedHierarchy`` through
``hierarchy_owned_lower_source_third_parameter_jet``.  This module performs the
remaining exact algebra in the displayed PositiveAxisSystem forcing vector:

    f_n = (0, 0, 0,
           2 xi p_n,
           2 source.angular,
           2 source.axial - 4 eta X p_n / ell),

where ``X=xi^2``, ``ell=1-2h eta^2`` and
``p_n=C^-2 source.pressure_product-source.omega_quotient/2``.

No caller-supplied SourceJet derivative table or normalization constant is
accepted: both the source rows and ``C`` come from the same hierarchy.  All eta
derivatives are analytic.  This is Stage-2 formal-structure infrastructure; it
does not materialize a Picard fixed point or a paper-exact coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import math

import numpy as np

from .background_positive_axis import positive_axis_forcing, pressure_source
from .background_repaired_history_sixth_mixed import (
    Section5LowerHistorySixthMixedHierarchy,
)
from .background_repaired_history_source_third_parameter import (
    hierarchy_owned_lower_source_third_parameter_jet,
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


def _forcing_row(value: np.ndarray, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (6,):
        raise ValueError(f"{name} must have shape (6,)")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be finite")
    out = out.copy()
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class PositiveAxisForcingThirdParameterJet:
    """Value and first three analytic eta derivatives of ``f_n``."""

    value: np.ndarray
    parameter: np.ndarray
    parameter2: np.ndarray
    parameter3: np.ndarray

    def __post_init__(self) -> None:
        for name in ("value", "parameter", "parameter2", "parameter3"):
            object.__setattr__(self, name, _forcing_row(getattr(self, name), name))


def _eta_over_ell_jet(h: float, eta: float) -> tuple[float, float, float, float]:
    """Return eta derivatives 0..3 of ``eta/(1-2h eta^2)`` exactly."""

    ell = 1.0 - 2.0 * h * eta * eta
    if ell <= 0.0:
        raise ValueError("similarity factor ell must be positive")
    q0 = eta / ell
    q1 = (1.0 + 2.0 * h * eta * eta) / ell**2
    q2 = 4.0 * h * eta * (3.0 + 2.0 * h * eta * eta) / ell**3
    q3 = 12.0 * h / ell**2 + 192.0 * h * h * eta * eta / ell**4
    out = (q0, q1, q2, q3)
    if not all(math.isfinite(value) for value in out):
        raise OverflowError("eta/ell third parameter jet is outside binary64 range")
    return out


def _product_rows(
    q: tuple[float, float, float, float],
    p: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    """Return derivatives 0..3 of the scalar product ``q*p``."""

    q0, q1, q2, q3 = q
    p0, p1, p2, p3 = p
    return (
        q0 * p0,
        q1 * p0 + q0 * p1,
        q2 * p0 + 2.0 * q1 * p1 + q0 * p2,
        q3 * p0 + 3.0 * q2 * p1 + 3.0 * q1 * p2 + q0 * p3,
    )


def hierarchy_owned_positive_axis_forcing_third_parameter_jet(
    hierarchy: Section5LowerHistorySixthMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
) -> PositiveAxisForcingThirdParameterJet:
    """Construct ``f_n`` through ``partial_eta^3`` from one strong hierarchy.

    The function consumes no caller-maintained source derivative table and no
    independent ``C``.  It first obtains the complete strict-lower source jet
    from the hierarchy, delegates the value row exactly to the landed
    ``positive_axis_forcing`` implementation, and differentiates only the
    displayed algebraic eta dependence of the forcing vector.
    """

    if not isinstance(hierarchy, Section5LowerHistorySixthMixedHierarchy):
        raise TypeError("hierarchy must be a Section5LowerHistorySixthMixedHierarchy")
    order = _positive_order(order)
    xi, eta = _point(xi, eta)
    X = xi * xi

    source = hierarchy_owned_lower_source_third_parameter_jet(
        hierarchy,
        order,
        X,
        eta,
    )
    rows = (source.value, source.parameter, source.parameter2, source.parameter3)
    p = tuple(pressure_source(hierarchy.C, row) for row in rows)
    qp = _product_rows(_eta_over_ell_jet(hierarchy.h, eta), p)

    value = positive_axis_forcing(
        hierarchy.h,
        hierarchy.C,
        xi,
        eta,
        source.value,
    )

    derivatives: list[np.ndarray] = []
    for derivative, row in enumerate(rows[1:], start=1):
        out = np.zeros(6, dtype=float)
        out[3] = 2.0 * xi * p[derivative]
        out[4] = 2.0 * row.angular
        out[5] = 2.0 * row.axial - 4.0 * X * qp[derivative]
        if not np.all(np.isfinite(out)):
            raise OverflowError(
                f"positive-axis forcing eta derivative {derivative} is outside binary64 range"
            )
        derivatives.append(out)

    return PositiveAxisForcingThirdParameterJet(
        value=value,
        parameter=derivatives[0],
        parameter2=derivatives[1],
        parameter3=derivatives[2],
    )
