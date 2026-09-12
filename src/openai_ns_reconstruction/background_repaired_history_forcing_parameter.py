"""Hierarchy-owned analytic eta jet of the Section 5 Eq. (5.7) forcing.

The exact positive-axis forcing is already landed in
:mod:`background_positive_axis`, while the strong repaired history now owns an
analytic eta derivative of the complete strict-lower ``actualLowerSource``.
This module composes those two pieces and differentiates only the explicit
``eta/ell`` geometry in the sixth forcing row.  Production uses no finite
difference, sampled derivative table, generic cutoff, or fitted coefficient.

The genuine leading mixed profile still belongs to Issue #1.  Missing strong
order-zero data therefore fail closed in the lower-source bridge.  This is
Stage-2 ``formal-structure`` solver infrastructure: it materializes a
hierarchy-owned ``partial_eta f_n`` but does not by itself solve the next Picard
iterate or make the background/velocity paper-exact.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

import numpy as np

from .background_positive_axis import positive_axis_forcing, pressure_source
from .background_repaired_history_phi_third_mixed import (
    Section5LowerHistoryPhiThirdMixedHierarchy,
)
from .background_repaired_history_source_parameter import (
    hierarchy_owned_lower_source_parameter_jet,
)


def _positive_order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("order must be a positive integer")
    return int(value)


def _vector6(value: np.ndarray, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (6,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite vector of shape (6,)")
    out = out.copy()
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class PositiveAxisForcingParameterJet:
    """Value and analytic ``eta`` derivative of the exact forcing ``f_n``."""

    value: np.ndarray
    parameter: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", _vector6(self.value, "value"))
        object.__setattr__(self, "parameter", _vector6(self.parameter, "parameter"))


def hierarchy_owned_positive_axis_forcing_parameter_jet(
    hierarchy: Section5LowerHistoryPhiThirdMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
) -> PositiveAxisForcingParameterJet:
    """Return hierarchy-owned ``(f_n, partial_eta f_n)`` for Eq. (5.7).

    ``order`` is the positive recursive coefficient being solved and only its
    strict-lower history is consumed.  The forcing value is delegated to the
    already-landed exact ``positive_axis_forcing`` implementation.  The source
    derivative comes from ``hierarchy_owned_lower_source_parameter_jet`` and the
    remaining explicit ``eta/ell`` factor is differentiated analytically.

    In particular, no independently supplied source derivative is accepted by
    this interface.  If the strong repaired hierarchy is incomplete (most
    importantly at order zero while Issue #1 remains upstream), the lower-source
    bridge fails closed instead of manufacturing ``partial_eta f_n``.
    """

    if not isinstance(hierarchy, Section5LowerHistoryPhiThirdMixedHierarchy):
        raise TypeError(
            "hierarchy must be a Section5LowerHistoryPhiThirdMixedHierarchy"
        )
    order = _positive_order(order)
    xi = float(xi)
    eta = float(eta)

    source = hierarchy_owned_lower_source_parameter_jet(
        hierarchy,
        order,
        xi * xi,
        eta,
    )
    value = positive_axis_forcing(
        hierarchy.h,
        hierarchy.C,
        xi,
        eta,
        source.value,
    )

    psrc = pressure_source(hierarchy.C, source.value)
    psrc_parameter = pressure_source(hierarchy.C, source.parameter)
    X = xi * xi
    ell = 1.0 - 2.0 * hierarchy.h * eta * eta

    parameter = np.zeros(6, dtype=float)
    parameter[3] = 2.0 * xi * psrc_parameter
    parameter[4] = 2.0 * source.parameter.angular
    parameter[5] = (
        2.0 * source.parameter.axial
        - 4.0
        * X
        * (
            (psrc + eta * psrc_parameter) / ell
            + 4.0 * hierarchy.h * eta * eta * psrc / (ell * ell)
        )
    )
    if not np.all(np.isfinite(parameter)):
        raise OverflowError("positive-axis forcing parameter jet is outside binary64 range")
    return PositiveAxisForcingParameterJet(value=value, parameter=parameter)
