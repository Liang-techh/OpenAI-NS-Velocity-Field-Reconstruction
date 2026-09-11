"""Independent finite-order consistency checks for Section 10 endpoint traces.

The landed Section 10 endpoint ladder records conditional, independently supplied
majorants for ``partial_t D^n R`` on the official late plateau.  The existing
trace enclosure then turns each majorant into a ball containing the hypothetical
endpoint limit.  This module adds a different check: two *actual supplied*
pre-endpoint full-jet evaluations must themselves obey the integrated Cauchy
bound implied by the majorant.

This is useful as a fail-closed regression gate once production Section 9
residual jets become machine-connected.  It can falsify a claimed majorant
without defining the expected answer from the same residual computation.

Passing the check does not prove the majorant on a spatial compact, locally
uniform convergence, endpoint-limit existence, or all-order smooth Borel
gluing.  It therefore remains formal-structure infrastructure and must not
promote ``paper_exact_velocity_available``.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import math

import numpy as np

from .endpoint_jets import validate_full_spacetime_jet
from .section10_endpoint_ladder import Section10EndpointMajorantLadder

PastFullSpacetimeJetFamily = Callable[[int, float, float, float, float], object]


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real number")
    try:
        out = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not math.isfinite(out):
        raise ValueError(f"{name} must be a finite real number")
    return out


def _within_bound(observed: float, bound: float) -> bool:
    """Float-safe comparison without weakening a materially violated bound."""
    if observed <= bound:
        return True
    return math.isclose(observed, bound, rel_tol=8e-15, abs_tol=0.0)


@dataclass(frozen=True)
class EndpointTraceDegreeConsistency:
    """One pointwise necessary consequence of a supplied endpoint majorant."""

    degree: int
    spatial_window: int
    point: tuple[float, float, float]
    t0: float
    t1: float
    endpoint: float
    observed_distance: float
    cauchy_upper_bound: float
    earlier_tail_radius: float
    later_tail_radius: float

    @property
    def compatible(self) -> bool:
        values = (
            *self.point,
            self.t0,
            self.t1,
            self.endpoint,
            self.observed_distance,
            self.cauchy_upper_bound,
            self.earlier_tail_radius,
            self.later_tail_radius,
        )
        return (
            self.degree >= 0
            and self.spatial_window >= 0
            and all(math.isfinite(value) for value in values)
            and self.t0 <= self.t1 < self.endpoint
            and self.observed_distance >= 0.0
            and self.cauchy_upper_bound >= 0.0
            and self.earlier_tail_radius >= 0.0
            and self.later_tail_radius >= 0.0
            and _within_bound(self.observed_distance, self.cauchy_upper_bound)
            and _within_bound(
                self.observed_distance,
                self.earlier_tail_radius + self.later_tail_radius,
            )
        )


@dataclass(frozen=True)
class Section10EndpointTraceConsistencyCertificate:
    """Finite-order pointwise compatibility with one endpoint-majorant ladder."""

    point: tuple[float, float, float]
    t0: float
    t1: float
    endpoint: float
    spatial_window: int
    rows: tuple[EndpointTraceDegreeConsistency, ...]
    actual_residual_jet_family_verified: bool = False
    uniform_spatial_majorant_verified: bool = False
    locally_uniform_endpoint_limits_verified: bool = False
    all_order_borel_smoothness_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_consistency_ready(self) -> bool:
        return (
            bool(self.rows)
            and tuple(row.degree for row in self.rows) == tuple(range(len(self.rows)))
            and all(row.compatible for row in self.rows)
            and all(row.point == self.point for row in self.rows)
            and all(row.t0 == self.t0 and row.t1 == self.t1 for row in self.rows)
            and all(row.endpoint == self.endpoint for row in self.rows)
            and all(row.spatial_window == self.spatial_window for row in self.rows)
            and not self.actual_residual_jet_family_verified
            and not self.uniform_spatial_majorant_verified
            and not self.locally_uniform_endpoint_limits_verified
            and not self.all_order_borel_smoothness_verified
            and not self.paper_exact_velocity_available
        )


def certify_section10_endpoint_trace_consistency(
    ladder: Section10EndpointMajorantLadder,
    past_full_jets: PastFullSpacetimeJetFamily,
    x: object,
    y: object,
    z: object,
    t0: object,
    t1: object,
) -> Section10EndpointTraceConsistencyCertificate:
    """Check two pre-endpoint jet samples against the ladder's FTC budgets.

    For every degree ``0..N`` this evaluates the dense full spacetime jet at the
    same spatial point and two times.  The observed Euclidean tensor distance is
    compared with the independently supplied integrated bound for
    ``partial_t D^n R``.  A violation raises instead of producing a certificate.

    This is deliberately only a *necessary* pointwise consistency test.  It does
    not infer a supremum bound over the compact spatial window from the sampled
    point, and it never fits or shrinks a majorant from the observed data.
    """
    if not isinstance(ladder, Section10EndpointMajorantLadder):
        raise ValueError("ladder must be a Section10EndpointMajorantLadder")
    if not ladder.certified:
        raise ValueError("endpoint majorant ladder must be certified")
    if not callable(past_full_jets):
        raise ValueError("past_full_jets must be callable")

    point = (_finite(x, "x"), _finite(y, "y"), _finite(z, "z"))
    cauchy = ladder.cauchy_certificates(t0, t1)
    start = cauchy[0].t0
    end = cauchy[0].t1
    tail0 = ladder.endpoint_tail_bounds(start)
    tail1 = ladder.endpoint_tail_bounds(end)

    rows: list[EndpointTraceDegreeConsistency] = []
    for degree, implication in enumerate(cauchy):
        earlier = validate_full_spacetime_jet(
            past_full_jets(degree, point[0], point[1], point[2], start),
            degree,
        )
        later = validate_full_spacetime_jet(
            past_full_jets(degree, point[0], point[1], point[2], end),
            degree,
        )
        observed = float(np.linalg.norm((later - earlier).reshape(-1)))
        if not math.isfinite(observed):
            raise ArithmeticError("endpoint-trace sample distance is not finite")

        row = EndpointTraceDegreeConsistency(
            degree=degree,
            spatial_window=ladder.spatial_window,
            point=point,
            t0=start,
            t1=end,
            endpoint=ladder.endpoint,
            observed_distance=observed,
            cauchy_upper_bound=float(implication.interval_upper_bound),
            earlier_tail_radius=float(tail0[degree]),
            later_tail_radius=float(tail1[degree]),
        )
        if not row.compatible:
            raise ValueError(
                "sampled full residual jets violate the supplied endpoint "
                f"majorant at derivative degree {degree}"
            )
        rows.append(row)

    result = Section10EndpointTraceConsistencyCertificate(
        point=point,
        t0=start,
        t1=end,
        endpoint=ladder.endpoint,
        spatial_window=ladder.spatial_window,
        rows=tuple(rows),
    )
    if not result.formal_consistency_ready:
        raise ArithmeticError("Section 10 endpoint-trace consistency invariant failed")
    return result
