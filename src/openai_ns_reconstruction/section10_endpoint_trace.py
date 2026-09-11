"""Finite-order endpoint-trace enclosures for the Section 10 residual.

``section10_endpoint_ladder.py`` turns independently proved late-plateau
majorants for ``partial_t D^n R`` into locally uniform endpoint-tail budgets.
This module uses those budgets without inventing an endpoint jet: at any
pre-endpoint time in the common certified regime, the actual value
``D^n R(t, x)`` is a center for a closed ball that must contain its endpoint
limit, conditional on the same majorant hypothesis.

The dense tensor convention is the one used by :mod:`endpoint_jets`: an
order-``n`` full spacetime jet has shape ``(3,) + (4,)*n``.  Distances here use
the Euclidean norm of that dense coordinate tensor.  Therefore a ladder fed
into this adapter must have been proved in the same norm (or a stronger norm
that dominates it).

This is still finite-order, conditional infrastructure.  It does not derive
majorants from the actual Section 9 residual, prove that a queried point lies
in the compact window indexed by ``spatial_window``, construct an exact limit,
or prove all-order Borel smoothness.  It must not promote
``paper_exact_velocity_available``.
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


@dataclass(frozen=True)
class EndpointJetEnclosure:
    """One conditional closed ball containing an endpoint residual jet.

    ``certified`` records only the arithmetic/interface consequence of the
    supplied majorant ladder.  In particular it is not evidence that the
    majorant hypothesis has itself been proved for the production residual.
    """

    degree: int
    spatial_window: int
    point: tuple[float, float, float]
    sample_time: float
    endpoint: float
    approximation: np.ndarray
    error_radius: float

    @property
    def certified(self) -> bool:
        expected = (3,) + (4,) * self.degree
        return (
            self.degree >= 0
            and self.spatial_window >= 0
            and self.approximation.shape == expected
            and np.all(np.isfinite(self.approximation))
            and all(math.isfinite(value) for value in self.point)
            and math.isfinite(self.sample_time)
            and math.isfinite(self.endpoint)
            and math.isfinite(self.error_radius)
            and self.sample_time < self.endpoint
            and self.error_radius >= 0.0
        )

    def distance_to(self, candidate: object) -> float:
        """Dense Euclidean distance from the center to a candidate endpoint jet."""
        tensor = validate_full_spacetime_jet(candidate, self.degree)
        distance = float(np.linalg.norm((tensor - self.approximation).reshape(-1)))
        if not math.isfinite(distance):
            raise ArithmeticError("endpoint-jet distance is not finite")
        return distance

    def contains(self, candidate: object) -> bool:
        """Check membership in the conditional endpoint-limit ball.

        Passing this check is only a necessary consistency check for a supplied
        candidate at this finite sample time; it does not identify the actual
        locally-uniform endpoint limit.
        """
        return self.distance_to(candidate) <= self.error_radius


def section10_endpoint_trace_enclosures(
    ladder: Section10EndpointMajorantLadder,
    past_full_jets: PastFullSpacetimeJetFamily,
    x: object,
    y: object,
    z: object,
    t: object,
) -> tuple[EndpointJetEnclosure, ...]:
    """Enclose endpoint jets ``D^n R(1-,x)`` for every degree ``0..N``.

    The center is evaluated from the caller's pre-endpoint full-jet family at
    one common time.  Radii come only from the already certified ladder tail
    budgets; no residual samples are used to fit or shrink those bounds.
    """
    if not isinstance(ladder, Section10EndpointMajorantLadder):
        raise ValueError("ladder must be a Section10EndpointMajorantLadder")
    if not ladder.certified:
        raise ValueError("endpoint majorant ladder must be certified")
    if not callable(past_full_jets):
        raise ValueError("past_full_jets must be callable")

    point = (_finite(x, "x"), _finite(y, "y"), _finite(z, "z"))
    time = _finite(t, "t")
    radii = ladder.endpoint_tail_bounds(time)

    out: list[EndpointJetEnclosure] = []
    for degree, radius in enumerate(radii):
        tensor = validate_full_spacetime_jet(
            past_full_jets(degree, point[0], point[1], point[2], time), degree
        ).copy()
        tensor.setflags(write=False)
        enclosure = EndpointJetEnclosure(
            degree=degree,
            spatial_window=ladder.spatial_window,
            point=point,
            sample_time=time,
            endpoint=ladder.endpoint,
            approximation=tensor,
            error_radius=float(radius),
        )
        if not enclosure.certified:
            raise ArithmeticError("Section 10 endpoint-jet enclosure invariant failed")
        out.append(enclosure)
    return tuple(out)
