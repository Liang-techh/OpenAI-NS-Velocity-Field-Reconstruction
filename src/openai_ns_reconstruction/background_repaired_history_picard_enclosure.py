"""Hierarchy-bound conditional enclosure for the Section 5 Picard value.

For one positive coefficient order ``n``, the landed hierarchy path now owns the
finite derivative triangle

    W_n^(0) = G f_n,
    W_n^(1) = K W_n^(0),
    W_n^(2) = K W_n^(1),
    W_n^(3) = K W_n^(2),

through the derivatives needed to materialize the four value rows.  Lemma 5.1
and Eq. (5.8) separately provide an analytic majorant for the infinite Picard
tail once a theorem-backed ``C_n`` and complex-strip radii are available.

This module joins those two pieces without inventing the missing analytic
bounds.  The retained prefix is regenerated only from
``Section5LowerHistorySixthMixedHierarchy``.  The Eq. (5.8) inputs are kept in a
typed object that is bound to the *same hierarchy instance and recursive order*
and the resulting tail arithmetic is evaluated by the landed
``picard_tail_certificate`` routine.  A point can be enclosed only when it lies
inside the radial/parameter domain named by that bound object.

The important truth boundary is intentional: caller-supplied ``C_n``, ``rho``
and ``rho_prime`` are still conditional theorem hypotheses.  They are not
inferred from samples, fitted residuals, finite differences, or a generic
cutoff, and this module does not promote them to hierarchy-certified analytic
bounds.  Agent 7 / the all-order analytic-closure lane must eventually replace
those conditional inputs with theorem-owned bounds.  Until then the returned
object is a conditional norm-ball enclosure, not a paper-exact coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from numbers import Integral

import numpy as np

from .background_picard_bounds import PicardTailCertificate, picard_tail_certificate
from .background_repaired_history_picard_first_third_parameter import (
    hierarchy_owned_first_picard_third_parameter_jet,
)
from .background_repaired_history_picard_second_second_parameter import (
    hierarchy_owned_second_picard_second_parameter_jet,
)
from .background_repaired_history_picard_third_first_parameter import (
    hierarchy_owned_third_picard_first_parameter_jet,
)
from .background_repaired_history_picard_fourth import (
    hierarchy_owned_fourth_picard_value,
)
from .background_repaired_history_sixth_mixed import (
    Section5LowerHistorySixthMixedHierarchy,
)


_PREFIX_TERMS = 4


def _positive_order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("order must be a positive integer")
    return int(value)


def _positive_finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


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
class HierarchyBoundPicardMajorantInputs:
    """Conditional Eq. (5.8) inputs bound to one hierarchy/order pair.

    ``C_n``, ``radial_extent``, ``rho`` and ``rho_prime`` are *not* derived in
    this class.  They must come from an external analytic-bound argument.  The
    class merely prevents a tail computed for one runtime hierarchy/order from
    being silently reused for another and evaluates the theorem-side Eq. (5.8)
    arithmetic at the first currently unmaterialized term ``k=4``.
    """

    hierarchy: Section5LowerHistorySixthMixedHierarchy
    order: int
    C_n: float
    radial_extent: float
    rho: float
    rho_prime: float
    provenance: str
    tail_certificate: PicardTailCertificate = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.hierarchy, Section5LowerHistorySixthMixedHierarchy):
            raise TypeError(
                "hierarchy must be a Section5LowerHistorySixthMixedHierarchy"
            )
        order = _positive_order(self.order)
        if order > len(self.hierarchy.sources):
            raise ValueError(
                "requested recursive order is missing one or more strict lower coefficients"
            )
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty string")

        C_n = _positive_finite(self.C_n, "C_n")
        radial_extent = _positive_finite(self.radial_extent, "radial_extent")
        rho = _positive_finite(self.rho, "rho")
        rho_prime = float(self.rho_prime)
        if not math.isfinite(rho_prime) or not 0.0 < rho_prime < rho:
            raise ValueError("rho_prime must be finite with 0 < rho_prime < rho")

        certificate = picard_tail_certificate(
            _PREFIX_TERMS,
            C_n,
            radial_extent,
            rho,
            rho_prime,
        )
        if not math.isfinite(certificate.tail_upper_bound):
            raise ValueError("Eq. (5.8) tail bound must be finite for this enclosure")

        object.__setattr__(self, "order", order)
        object.__setattr__(self, "C_n", C_n)
        object.__setattr__(self, "radial_extent", radial_extent)
        object.__setattr__(self, "rho", rho)
        object.__setattr__(self, "rho_prime", rho_prime)
        object.__setattr__(self, "tail_certificate", certificate)


@dataclass(frozen=True)
class ConditionalPositiveAxisPicardValueEnclosure:
    """Four-term hierarchy-owned center plus a conditional Eq. (5.8) radius.

    The scalar radius is the theorem-side norm bound on the complete omitted
    tail ``sum_{k>=4} K^k G f_n``.  No componentwise interval interpretation is
    added here because this module intentionally does not redefine the norm used
    by the manuscript's analytic estimate.
    """

    order: int
    xi: float
    eta: float
    prefix_value: np.ndarray
    tail_upper_bound: float
    start_order: int
    bound_provenance: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "order", _positive_order(self.order))
        xi, eta = _point(self.xi, self.eta)
        object.__setattr__(self, "xi", xi)
        object.__setattr__(self, "eta", eta)
        object.__setattr__(
            self, "prefix_value", _vector6(self.prefix_value, "prefix_value")
        )
        tail = float(self.tail_upper_bound)
        if not math.isfinite(tail) or tail < 0.0:
            raise ValueError("tail_upper_bound must be finite and nonnegative")
        object.__setattr__(self, "tail_upper_bound", tail)
        if self.start_order != _PREFIX_TERMS:
            raise ValueError(f"start_order must equal {_PREFIX_TERMS}")
        if not isinstance(self.bound_provenance, str) or not self.bound_provenance.strip():
            raise ValueError("bound_provenance must be a nonempty string")

    def certifies_tolerance(self, tolerance: float) -> bool:
        """Return whether the conditional omitted-tail radius meets tolerance."""

        tolerance = _positive_finite(tolerance, "tolerance")
        return self.tail_upper_bound <= tolerance


def hierarchy_owned_picard_prefix4_value(
    hierarchy: Section5LowerHistorySixthMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
) -> np.ndarray:
    """Return the exact currently landed four-term Picard value prefix.

    This function performs no tail extrapolation.  Each term is regenerated by
    its hierarchy-owned bridge, so there is no caller SourceJet, repair jet,
    forcing derivative table, matrix derivative table, or independent
    normalization constant on this path.
    """

    if not isinstance(hierarchy, Section5LowerHistorySixthMixedHierarchy):
        raise TypeError("hierarchy must be a Section5LowerHistorySixthMixedHierarchy")
    order = _positive_order(order)
    if order > len(hierarchy.sources):
        raise ValueError(
            "requested recursive order is missing one or more strict lower coefficients"
        )
    xi, eta = _point(xi, eta)

    first = hierarchy_owned_first_picard_third_parameter_jet(
        hierarchy, order, xi, eta
    )
    second = hierarchy_owned_second_picard_second_parameter_jet(
        hierarchy, order, xi, eta
    )
    third = hierarchy_owned_third_picard_first_parameter_jet(
        hierarchy, order, xi, eta
    )
    fourth = hierarchy_owned_fourth_picard_value(hierarchy, order, xi, eta)
    return _vector6(
        first.value + second.value + third.value + fourth.value,
        "four-term Picard prefix",
    )


def hierarchy_owned_conditional_picard_value_enclosure(
    hierarchy: Section5LowerHistorySixthMixedHierarchy,
    order: int,
    xi: float,
    eta: float,
    majorant: HierarchyBoundPicardMajorantInputs,
) -> ConditionalPositiveAxisPicardValueEnclosure:
    """Bind the hierarchy-owned prefix to a same-hierarchy Eq. (5.8) tail.

    The function fails closed on cross-wired hierarchy/order evidence and on
    points outside the radial/parameter domain named by ``majorant``.  It never
    infers or weakens those analytic hypotheses.
    """

    if not isinstance(hierarchy, Section5LowerHistorySixthMixedHierarchy):
        raise TypeError("hierarchy must be a Section5LowerHistorySixthMixedHierarchy")
    if not isinstance(majorant, HierarchyBoundPicardMajorantInputs):
        raise TypeError("majorant must be HierarchyBoundPicardMajorantInputs")
    order = _positive_order(order)
    xi, eta = _point(xi, eta)

    if majorant.hierarchy is not hierarchy:
        raise ValueError("Picard majorant is bound to a different hierarchy instance")
    if majorant.order != order:
        raise ValueError("Picard majorant is bound to a different recursive order")
    if xi > majorant.radial_extent:
        raise ValueError("xi lies outside the radial extent certified by the majorant")
    if abs(eta) > majorant.rho_prime:
        raise ValueError("eta lies outside the smaller parameter strip")

    prefix = hierarchy_owned_picard_prefix4_value(hierarchy, order, xi, eta)
    certificate = majorant.tail_certificate
    return ConditionalPositiveAxisPicardValueEnclosure(
        order=order,
        xi=xi,
        eta=eta,
        prefix_value=prefix,
        tail_upper_bound=certificate.tail_upper_bound,
        start_order=certificate.start_order,
        bound_provenance=majorant.provenance,
    )
