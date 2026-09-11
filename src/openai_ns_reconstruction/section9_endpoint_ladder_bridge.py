"""Finite Section 9 uniform-derivative ladder admission to Section 10.

PR #100 deliberately admits one endpoint jet only after a separately justified
uniform physical-time estimate for ``partial_t D^n R`` on one compact spatial
window.  PR #104, by contrast, records finite powers of the Section 9
``q``-flat remainder.  Those are different theorem inputs: a finite family of
``E(q) <= C_N q^N`` statements does not by itself produce a uniform
``(1-t)`` majorant on the Section 10 late-time compact.

This module batches the already-landed single-degree bridge without weakening
that distinction.  Every degree ``0..N`` must have its own genuine
:class:`Section9UniformEndpointDerivativeWitness`, all witnesses must refer to
one Section 9 stage and one spatial window, and every row is independently
admitted through :func:`admit_section9_uniform_endpoint_witness`.  Only then is
the resulting finite family packaged as a
:class:`Section10EndpointMajorantLadder`.

The result is finite-order formal infrastructure.  It neither constructs the
actual Section 9 correction sequence nor proves the supplied uniform estimates,
all-order endpoint limits, Borel gluing, smooth forcing, bounded energy, or
blow-up.  ``paper_exact_velocity_available`` therefore remains false.
"""
from __future__ import annotations

from dataclasses import dataclass
import operator
from typing import Iterable

from .section10_endpoint_ladder import Section10EndpointMajorantLadder
from .section9_endpoint_bridge import (
    Section9EndpointBridgeRecord,
    Section9UniformEndpointDerivativeWitness,
    admit_section9_uniform_endpoint_witness,
)
from .section9_flat_remainder import Section9FlatRemainderLadderCertificate


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a nonnegative integer")
    try:
        out = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be a nonnegative integer") from exc
    if out < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(out)


@dataclass(frozen=True)
class Section9EndpointLadderBridgeRecord:
    """Finite coherent admission record from Section 9 to Section 10."""

    stage: int
    spatial_window: int
    max_endpoint_degree: int
    bridge_records: tuple[Section9EndpointBridgeRecord, ...]
    endpoint_ladder: Section10EndpointMajorantLadder
    status: str = "formal-structure"
    uniform_witness_ladder_supplied: bool = True
    section9_flat_remainder_ladder_sufficient: bool = False
    source_theorems_machine_verified: bool = False
    actual_section9_sequence_verified: bool = False
    endpoint_limits_constructed: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def endpoint_degrees(self) -> tuple[int, ...]:
        return tuple(row.endpoint_derivative_degree for row in self.bridge_records)

    @property
    def formal_bridge_ready(self) -> bool:
        expected = tuple(range(self.max_endpoint_degree + 1))
        return (
            self.endpoint_degrees == expected
            and all(row.stage == self.stage for row in self.bridge_records)
            and all(row.spatial_window == self.spatial_window for row in self.bridge_records)
            and all(row.formal_bridge_ready for row in self.bridge_records)
            and self.endpoint_ladder.certified
            and self.endpoint_ladder.max_degree == self.max_endpoint_degree
            and self.endpoint_ladder.spatial_window == self.spatial_window
            and not self.section9_flat_remainder_ladder_sufficient
            and not self.source_theorems_machine_verified
            and not self.actual_section9_sequence_verified
            and not self.endpoint_limits_constructed
            and not self.paper_exact_velocity_available
        )


def admit_section9_uniform_endpoint_ladder(
    witnesses: Iterable[Section9UniformEndpointDerivativeWitness],
    *,
    max_endpoint_degree: int,
) -> Section9EndpointLadderBridgeRecord:
    """Admit exactly the finite endpoint degrees ``0..max_endpoint_degree``.

    Each witness is first sent through the single-degree PR #100 bridge, so the
    derivative-order shift, integrable exponent, official ``t>=3/4`` plateau,
    endpoint ``T=1``, evidence provenance, and Section 10 localization gate are
    all rechecked rather than duplicated here.

    A :class:`Section9FlatRemainderLadderCertificate` is rejected explicitly.
    Even a coherent finite ``q``-power ladder is not a physical-time uniform
    derivative witness and cannot bypass the missing Section 9 common-domain /
    support theorem.
    """

    max_endpoint_degree = _natural(max_endpoint_degree, "max_endpoint_degree")
    if isinstance(witnesses, Section9FlatRemainderLadderCertificate):
        raise TypeError(
            "a finite q-flat remainder ladder cannot be promoted to a physical-time endpoint ladder"
        )

    rows = tuple(witnesses)
    if not rows:
        raise ValueError("uniform endpoint witness ladder must be nonempty")
    if not all(isinstance(row, Section9UniformEndpointDerivativeWitness) for row in rows):
        raise TypeError(
            "every ladder entry must be a Section9UniformEndpointDerivativeWitness"
        )

    by_degree: dict[int, Section9UniformEndpointDerivativeWitness] = {}
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
            "uniform endpoint witness ladder must cover exactly degrees 0..max_endpoint_degree; "
            f"missing={missing}, extra={extra}"
        )

    ordered = tuple(by_degree[n] for n in range(max_endpoint_degree + 1))
    first = ordered[0]
    for row in ordered[1:]:
        if row.stage != first.stage:
            raise ValueError("all endpoint witnesses must use one Section 9 stage")
        if row.spatial_window != first.spatial_window:
            raise ValueError("all endpoint witnesses must use one spatial window")

    records = tuple(admit_section9_uniform_endpoint_witness(row) for row in ordered)
    endpoint_ladder = Section10EndpointMajorantLadder(row.majorant for row in records)
    record = Section9EndpointLadderBridgeRecord(
        stage=first.stage,
        spatial_window=first.spatial_window,
        max_endpoint_degree=max_endpoint_degree,
        bridge_records=records,
        endpoint_ladder=endpoint_ladder,
    )
    if not record.formal_bridge_ready:
        raise ArithmeticError("Section 9 endpoint ladder bridge invariant failed")
    return record
