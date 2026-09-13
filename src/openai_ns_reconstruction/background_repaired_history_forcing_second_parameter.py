"""Hierarchy-owned analytic second eta jet of the Section 5 Eq. (5.7) forcing.

The landed forcing bridge already owns ``(f_n, partial_eta f_n)`` from the
strong repaired lower-history hierarchy.  The strict-lower source now also owns
``partial_eta^2 actualLowerSource``.  This module composes those two landed
layers and differentiates only the explicit ``eta/ell(h,eta)`` geometry in the
sixth forcing row one more time.

Production uses no finite difference, sampled derivative table, generic cutoff,
or fitted coefficient.  Genuine order-zero fourth/fifth-mixed profile data
remain an Issue-#1 input, so missing strong leading data fail closed.  This is
Stage-2 ``formal-structure`` infrastructure: it does not establish Picard
convergence, recursive coefficient materialization, or paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

import numpy as np

from .background_positive_axis import pressure_source
from .background_repaired_history_forcing_parameter import (
    hierarchy_owned_positive_axis_forcing_parameter_jet,
)
from .background_repaired_history_phi_fourth_mixed import (
    Section5LowerHistoryPhiFourthMixedHierarchy,
)
from .background_repaired_history_source_second_parameter import (
    hierarchy_owned_lower_source_second_parameter_jet,
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
class PositiveAxisForcingSecondParameterJet:
    """Value and first two analytic eta derivatives of the exact forcing ``f_n``."""

    value: np.ndarray
    parameter: np.ndarray
    parameter2: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", _vector6(self.value, "value"))
        object.__setattr__(self, "parameter", _vector6(self.parameter, "parameter"))
        object.__setattr__(self, "parameter2", _vector6(self.parameter2, "parameter2"))


def hierarchy_owned_positive_axis_forcing_second_parameter_jet(
    hierarchy: Section5LowerHistoryPhiFourthMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
) -> PositiveAxisForcingSecondParameterJet:
    """Return hierarchy-owned ``(f_n, partial_eta f_n, partial_eta^2 f_n)``.

    The already-landed first-parameter forcing bridge remains authoritative for
    ``value`` and ``parameter``.  Only ``parameter2`` is new here.  It uses the
    hierarchy-owned second-eta strict-lower source and the exact derivative of
    the final-row geometry ``g(eta)=eta/ell(h,eta)``:

    ``g'' = 12 h eta / ell^2 + 32 h^2 eta^3 / ell^3``.

    No independent source derivative may be supplied by the caller.  The strong
    fourth/fifth-mixed hierarchy requirement is deliberate and keeps Issue-#1
    leading-profile gaps fail closed, including on the axis.
    """

    if not isinstance(hierarchy, Section5LowerHistoryPhiFourthMixedHierarchy):
        raise TypeError(
            "hierarchy must be a Section5LowerHistoryPhiFourthMixedHierarchy"
        )
    order = _positive_order(order)
    xi = float(xi)
    eta = float(eta)

    first = hierarchy_owned_positive_axis_forcing_parameter_jet(
        hierarchy,
        order,
        xi,
        eta,
    )
    source = hierarchy_owned_lower_source_second_parameter_jet(
        hierarchy,
        order,
        xi * xi,
        eta,
    )

    psrc = pressure_source(hierarchy.C, source.value)
    psrc_parameter = pressure_source(hierarchy.C, source.parameter)
    psrc_parameter2 = pressure_source(hierarchy.C, source.parameter2)

    X = xi * xi
    ell = 1.0 - 2.0 * hierarchy.h * eta * eta
    geometry = eta / ell
    geometry_parameter = (
        1.0 / ell
        + 4.0 * hierarchy.h * eta * eta / (ell * ell)
    )
    geometry_parameter2 = (
        12.0 * hierarchy.h * eta / (ell * ell)
        + 32.0
        * hierarchy.h
        * hierarchy.h
        * eta
        * eta
        * eta
        / (ell * ell * ell)
    )

    parameter2 = np.zeros(6, dtype=float)
    parameter2[3] = 2.0 * xi * psrc_parameter2
    parameter2[4] = 2.0 * source.parameter2.angular
    parameter2[5] = (
        2.0 * source.parameter2.axial
        - 4.0
        * X
        * (
            geometry_parameter2 * psrc
            + 2.0 * geometry_parameter * psrc_parameter
            + geometry * psrc_parameter2
        )
    )
    if not np.all(np.isfinite(parameter2)):
        raise OverflowError(
            "positive-axis forcing second parameter jet is outside binary64 range"
        )

    return PositiveAxisForcingSecondParameterJet(
        value=first.value,
        parameter=first.parameter,
        parameter2=parameter2,
    )
