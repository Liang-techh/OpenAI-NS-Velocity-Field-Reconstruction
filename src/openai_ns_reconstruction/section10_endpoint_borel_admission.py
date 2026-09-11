"""Fail-closed finite-order admission from past endpoint enclosures to Borel jets.

Section 10 has two already-landed conditional interfaces:

* :mod:`section10_endpoint_trace` encloses each candidate endpoint tensor in a
  closed ball centered at an actual caller-supplied pre-endpoint residual jet,
  with radius coming only from an independently supplied derivative majorant;
* :mod:`section10_borel_prefix_right_jets` proves exact right-jet algebra for a
  finite Taylor--Borel prefix once full endpoint tensors are supplied.

This module connects those interfaces without upgrading either hypothesis.  A
candidate full spacetime endpoint jet family must first lie inside every
finite-order trace enclosure at one certified pre-endpoint sample time.  Only
then is the *same frozen candidate family* admitted to the finite Borel-prefix
right-jet certificate.

Membership in one enclosure is a necessary consistency check, not a proof that
this candidate is the locally-uniform endpoint limit.  The majorant ladder is
still caller supplied, different candidates may fit the same nonzero balls,
and no infinite Borel-tail estimate is proved here.  In particular this bridge
must never promote ``paper_exact_velocity_available``.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .endpoint_borel import LocalScale
from .endpoint_jets import FullSpacetimeJetFamily, validate_full_spacetime_jet
from .section10_borel_prefix_right_jets import (
    Section10BorelPrefixRightJetsCertificate,
    certify_section10_borel_prefix_right_jets,
)
from .section10_endpoint_ladder import Section10EndpointMajorantLadder
from .section10_endpoint_trace import (
    PastFullSpacetimeJetFamily,
    section10_endpoint_trace_enclosures,
)


@dataclass(frozen=True)
class EndpointCandidateAdmissionRow:
    """Necessary finite-time consistency check for one endpoint derivative."""

    degree: int
    distance_to_center: float
    enclosure_radius: float

    def __post_init__(self) -> None:
        if isinstance(self.degree, bool) or not isinstance(self.degree, int) or self.degree < 0:
            raise ValueError("degree must be a nonnegative integer")
        distance = float(self.distance_to_center)
        radius = float(self.enclosure_radius)
        if not (math.isfinite(distance) and math.isfinite(radius)):
            raise ValueError("candidate distance and enclosure radius must be finite")
        if distance < 0.0 or radius < 0.0:
            raise ValueError("candidate distance and enclosure radius must be nonnegative")
        object.__setattr__(self, "distance_to_center", distance)
        object.__setattr__(self, "enclosure_radius", radius)

    @property
    def admitted(self) -> bool:
        """No numerical tolerance is allowed to enlarge the analytic budget."""

        return self.distance_to_center <= self.enclosure_radius


@dataclass(frozen=True)
class Section10EndpointBorelAdmissionCertificate:
    """Finite necessary left-trace check followed by exact prefix right jets.

    The negative flags are part of the truth contract.  Even when
    ``formal_candidate_ready`` is true, the supplied candidate has not been
    identified as the actual Section 9 endpoint limit and the infinite Borel
    smooth-extension theorem has not been established.
    """

    point: tuple[float, float, float]
    sample_time: float
    endpoint: float
    max_degree: int
    spatial_window: int
    rows: tuple[EndpointCandidateAdmissionRow, ...]
    right_prefix: Section10BorelPrefixRightJetsCertificate
    endpoint_trace_candidate_consistency_verified: bool = True
    finite_prefix_right_jet_algebra_verified: bool = True
    source_majorants_derived_from_actual_residual_verified: bool = False
    actual_section9_residual_limits_verified: bool = False
    endpoint_limit_uniqueness_verified: bool = False
    analytic_template_bounds_verified: bool = False
    infinite_borel_right_jets_verified: bool = False
    all_order_borel_smoothness_verified: bool = False
    smooth_compact_forcing_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_candidate_ready(self) -> bool:
        expected_degrees = tuple(range(self.max_degree + 1))
        return (
            self.endpoint == 1.0
            and tuple(row.degree for row in self.rows) == expected_degrees
            and all(row.admitted for row in self.rows)
            and self.right_prefix.point == self.point
            and self.right_prefix.endpoint == self.endpoint
            and self.right_prefix.max_degree == self.max_degree
            and self.right_prefix.formal_prefix_ready
            and self.endpoint_trace_candidate_consistency_verified
            and self.finite_prefix_right_jet_algebra_verified
            and not self.source_majorants_derived_from_actual_residual_verified
            and not self.actual_section9_residual_limits_verified
            and not self.endpoint_limit_uniqueness_verified
            and not self.analytic_template_bounds_verified
            and not self.infinite_borel_right_jets_verified
            and not self.all_order_borel_smoothness_verified
            and not self.smooth_compact_forcing_verified
            and not self.paper_exact_velocity_available
        )


def admit_section10_endpoint_candidate_for_borel_prefix(
    ladder: Section10EndpointMajorantLadder,
    past_full_jets: PastFullSpacetimeJetFamily,
    candidate_full_jets: FullSpacetimeJetFamily,
    local_scale: LocalScale,
    x: object,
    y: object,
    z: object,
    t: object,
) -> Section10EndpointBorelAdmissionCertificate:
    """Admit one finite endpoint candidate to the finite Borel-prefix path.

    The order is deliberate and fail closed:

    1. build the existing trace enclosures from the certified finite ladder;
    2. evaluate and freeze every candidate tensor through ``ladder.max_degree``;
    3. reject immediately if any tensor lies outside its enclosure;
    4. only then feed that exact frozen family to the finite right-jet algebra.

    Thus a rejected candidate cannot obtain an apparently successful Borel
    prefix certificate by changing values between the left and right checks.
    Passing remains only a necessary finite-order consistency statement.
    """

    if not isinstance(ladder, Section10EndpointMajorantLadder):
        raise ValueError("ladder must be a Section10EndpointMajorantLadder")
    if not ladder.certified:
        raise ValueError("endpoint majorant ladder must be certified")
    if not callable(past_full_jets):
        raise ValueError("past_full_jets must be callable")
    if not callable(candidate_full_jets):
        raise ValueError("candidate_full_jets must be callable")
    if not callable(local_scale):
        raise ValueError("local_scale must be callable")

    enclosures = section10_endpoint_trace_enclosures(
        ladder, past_full_jets, x, y, z, t
    )
    if len(enclosures) != ladder.max_degree + 1:
        raise ArithmeticError("endpoint enclosure ladder length mismatch")

    point = enclosures[0].point
    sample_time = enclosures[0].sample_time
    frozen_candidates: list[np.ndarray] = []
    rows: list[EndpointCandidateAdmissionRow] = []

    for enclosure in enclosures:
        tensor = validate_full_spacetime_jet(
            candidate_full_jets(enclosure.degree, *point), enclosure.degree
        ).copy()
        tensor.setflags(write=False)
        distance = enclosure.distance_to(tensor)
        row = EndpointCandidateAdmissionRow(
            degree=enclosure.degree,
            distance_to_center=distance,
            enclosure_radius=enclosure.error_radius,
        )
        if not row.admitted:
            raise ValueError(
                f"degree-{enclosure.degree} endpoint candidate lies outside the certified trace enclosure"
            )
        frozen_candidates.append(tensor)
        rows.append(row)

    frozen_tuple = tuple(frozen_candidates)

    def frozen_full_jets(degree: int, qx: float, qy: float, qz: float) -> np.ndarray:
        if isinstance(degree, bool) or not isinstance(degree, int):
            raise ValueError("degree must be an integer in the frozen candidate prefix")
        if degree < 0 or degree >= len(frozen_tuple):
            raise ValueError("degree lies outside the frozen candidate prefix")
        query_point = (float(qx), float(qy), float(qz))
        if not all(math.isfinite(value) for value in query_point):
            raise ValueError("candidate query point must be finite")
        if query_point != point:
            raise ValueError("frozen endpoint candidate is certified only at the requested point")
        return frozen_tuple[degree].copy()

    right_prefix = certify_section10_borel_prefix_right_jets(
        frozen_full_jets,
        local_scale,
        *point,
        ladder.max_degree,
        endpoint=ladder.endpoint,
    )

    result = Section10EndpointBorelAdmissionCertificate(
        point=point,
        sample_time=sample_time,
        endpoint=ladder.endpoint,
        max_degree=ladder.max_degree,
        spatial_window=ladder.spatial_window,
        rows=tuple(rows),
        right_prefix=right_prefix,
    )
    if not result.formal_candidate_ready:
        raise ArithmeticError("Section 10 endpoint-to-Borel admission invariant failed")
    return result
