"""Bind a ladder of analytically derived Section 9 endpoint majorants to one residual source.

The single-degree Eq. (9.18) adapter derives an endpoint majorant from a theorem-facing
uniform residual envelope, while ``section9_endpoint_residual_source`` binds an already
assembled endpoint ladder to one stable residual-source revision.  This module connects
those two existing layers without adding another free majorant path.

Every endpoint degree ``0..N`` must be represented by one
``Section9UniformResidualEnvelopeWitness``.  The witnesses must share one source id,
source revision, Section 9 stage, and spatial window, and each row is passed through the
existing analytic Eq. (9.18) adapter.  The resulting majorants are then assembled into
the existing finite ladder record and finally admitted by the existing residual-source
gate, which rechecks evidence kinds/provenance against the frozen source witness.

This is still only ``formal-structure`` plumbing.  In particular, the uniform envelope
witnesses remain theorem inputs: this module does not prove that they were derived from
the actual Eq. (9.21) residual, does not construct endpoint limits, and does not promote
``paper_exact_velocity_available``.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Iterable

from .section10_endpoint_ladder import Section10EndpointMajorantLadder
from .section9_endpoint_ladder_bridge import Section9EndpointLadderBridgeRecord
from .section9_endpoint_majorant_adapter import (
    Section9DerivedEndpointMajorantRecord,
    Section9UniformResidualEnvelopeWitness,
    derive_section9_endpoint_majorant_from_uniform_envelope,
)
from .section9_endpoint_residual_source import (
    Section9ResidualEndpointSourceBinding,
    Section9ResidualEndpointSourceWitness,
)


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


@dataclass(frozen=True)
class Section9DerivedMajorantSourceBindingRecord:
    """Finite Eq. (9.18)-derived ladder bound to one residual-source revision."""

    source_id: str
    source_revision: str
    stage: int
    spatial_window: int
    max_endpoint_degree: int
    derived_records: tuple[Section9DerivedEndpointMajorantRecord, ...]
    ladder_bridge: Section9EndpointLadderBridgeRecord
    source_binding: Section9ResidualEndpointSourceBinding
    status: str = "formal-structure"
    source_majorants_derived_from_actual_residual_verified: bool = False
    actual_section9_sequence_verified: bool = False
    endpoint_limits_constructed: bool = False
    all_order_borel_smoothness_verified: bool = False
    smooth_compact_forcing_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def endpoint_degrees(self) -> tuple[int, ...]:
        return tuple(row.endpoint_derivative_degree for row in self.derived_records)

    @property
    def formal_source_binding_ready(self) -> bool:
        expected = tuple(range(self.max_endpoint_degree + 1))
        return (
            self.endpoint_degrees == expected
            and self.ladder_bridge.formal_bridge_ready
            and self.source_binding.formal_source_binding_ready
            and self.source_binding.source_key == self.source_key
            and self.ladder_bridge.stage == self.stage
            and self.ladder_bridge.spatial_window == self.spatial_window
            and self.ladder_bridge.max_endpoint_degree == self.max_endpoint_degree
            and all(row.formal_adapter_ready for row in self.derived_records)
            and all(row.source_key == self.source_key for row in self.derived_records)
            and not self.source_majorants_derived_from_actual_residual_verified
            and not self.actual_section9_sequence_verified
            and not self.endpoint_limits_constructed
            and not self.all_order_borel_smoothness_verified
            and not self.smooth_compact_forcing_verified
            and not self.paper_exact_velocity_available
        )


def bind_section9_uniform_envelope_ladder_to_residual_source(
    witnesses: Iterable[Section9UniformResidualEnvelopeWitness],
    source: Section9ResidualEndpointSourceWitness,
    *,
    max_endpoint_degree: int,
) -> Section9DerivedMajorantSourceBindingRecord:
    """Derive degrees ``0..N`` from Eq. (9.18) envelopes and bind one source.

    The function intentionally does not accept preassembled majorants.  Each row must
    start as the theorem-facing uniform Eq. (9.18) envelope used by the landed analytic
    adapter.  After derivation, the existing ladder and source-binding classes are used
    directly rather than reimplementing their checks.
    """

    max_endpoint_degree = _natural(max_endpoint_degree, "max_endpoint_degree")
    if not isinstance(source, Section9ResidualEndpointSourceWitness):
        raise TypeError("source must be a Section9ResidualEndpointSourceWitness")

    rows = tuple(witnesses)
    if not rows:
        raise ValueError("uniform residual-envelope ladder must be nonempty")
    if not all(isinstance(row, Section9UniformResidualEnvelopeWitness) for row in rows):
        raise TypeError(
            "every ladder entry must be a Section9UniformResidualEnvelopeWitness"
        )

    by_degree: dict[int, Section9UniformResidualEnvelopeWitness] = {}
    for row in rows:
        degree = row.endpoint_derivative_degree
        if degree in by_degree:
            raise ValueError(f"duplicate endpoint derivative degree n={degree}")
        by_degree[degree] = row

    required = set(range(max_endpoint_degree + 1))
    supplied = set(by_degree)
    if supplied != required:
        missing = sorted(required - supplied)
        extra = sorted(supplied - required)
        raise ValueError(
            "uniform residual-envelope ladder must cover exactly degrees "
            "0..max_endpoint_degree; "
            f"missing={missing}, extra={extra}"
        )

    ordered = tuple(by_degree[n] for n in range(max_endpoint_degree + 1))
    first = ordered[0]
    for row in ordered:
        if row.source_key != first.source_key:
            raise ValueError("all residual-envelope witnesses must use one source revision")
        if row.stage != first.stage:
            raise ValueError("all residual-envelope witnesses must use one Section 9 stage")
        if row.spatial_window != first.spatial_window:
            raise ValueError("all residual-envelope witnesses must use one spatial window")

    if first.source_key != source.source_key:
        raise ValueError("residual-envelope ladder and source witness must use one source revision")
    if first.stage != source.stage:
        raise ValueError("residual-envelope ladder and source witness must use one Section 9 stage")
    if first.spatial_window != source.spatial_window:
        raise ValueError("residual-envelope ladder and source witness must use one spatial window")

    derived = tuple(
        derive_section9_endpoint_majorant_from_uniform_envelope(row) for row in ordered
    )
    bridge_rows = tuple(row.bridge for row in derived)
    endpoint_ladder = Section10EndpointMajorantLadder(
        row.majorant for row in bridge_rows
    )
    ladder_bridge = Section9EndpointLadderBridgeRecord(
        stage=first.stage,
        spatial_window=first.spatial_window,
        max_endpoint_degree=max_endpoint_degree,
        bridge_records=bridge_rows,
        endpoint_ladder=endpoint_ladder,
    )
    if not ladder_bridge.formal_bridge_ready:
        raise ArithmeticError("derived Section 9 endpoint ladder bridge invariant failed")

    source_binding = Section9ResidualEndpointSourceBinding(
        bridge=ladder_bridge,
        source=source,
    )
    result = Section9DerivedMajorantSourceBindingRecord(
        source_id=first.source_id,
        source_revision=first.source_revision,
        stage=first.stage,
        spatial_window=first.spatial_window,
        max_endpoint_degree=max_endpoint_degree,
        derived_records=derived,
        ladder_bridge=ladder_bridge,
        source_binding=source_binding,
    )
    if not result.formal_source_binding_ready:
        raise ArithmeticError("derived majorant residual-source binding invariant failed")
    return result
