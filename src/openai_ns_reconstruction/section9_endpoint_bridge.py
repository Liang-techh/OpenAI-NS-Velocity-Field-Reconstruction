"""Fail-closed admission from a uniform Section 9 derivative witness to Section 10.

Section 10 needs a locally uniform limit of every residual jet ``D^n R`` as
``t -> 1-``.  :mod:`endpoint_limit_majorant` records one sufficient condition:
an integrable bound for ``partial_t D^n R`` on a compact spatial window.  The
new Section 9 stage certificate, however, is deliberately only *pointwise* in
``q`` and therefore cannot establish that hypothesis by itself.

This module makes the missing interface explicit.  A caller must provide a
separately justified **uniform** witness for the order-``n+1`` Section 9
residual norm on the whole late-time window and spatial compact.  Only then is
that witness admitted as an :class:`EndpointPowerLawMajorant`, and the result
is immediately passed through the official Section 10 ``t >= 3/4`` plateau
gate.

No ``q``-to-``1-t`` comparison is inferred here, no pointwise certificate is
promoted to a uniform theorem, and no coefficient/exponent is fitted from
samples.  Until the genuine Section 9 correction sequence supplies the uniform
witnesses, this remains ``formal-structure`` and
``paper_exact_velocity_available`` stays false.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import operator

from .endpoint_limit_majorant import EndpointPowerLawMajorant
from .section9_stage_certificate import CertifiedBoundDatum, Section9PointwiseBoundCertificate
from .time_localization import (
    LATE_START,
    SECTION10_ENDPOINT,
    EndpointLocalizationTransferCertificate,
    section10_endpoint_localization_transfer,
)


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
class Section9UniformEndpointDerivativeWitness:
    """Theorem-level input needed to bridge one residual jet to Section 10.

    Semantically, ``coefficient_bound`` and ``singularity_exponent`` assert the
    independently justified estimate

    ``sup_{x in K} ||partial_t D^n R_j(t,x)||
       <= C (1-t)^(-alpha)``

    for every ``valid_from <= t < 1`` on the stated spatial compact ``K``.
    Because ``partial_t D^n`` has total spacetime derivative order ``n+1``, the
    Section 9 derivative order is required to be exactly one larger than the
    endpoint jet degree.  This catches an otherwise easy off-by-one promotion.

    ``CertifiedBoundDatum`` classifies the evidence behind ``C``.  A
    ``certified-numerical`` datum is admissible only when it is a rigorous
    *uniform interval enclosure* for the whole stated region, not a grid
    maximum or fitted sample.  The class cannot prove that semantic assertion;
    the provenance string is therefore part of the fail-closed audit trail.
    """

    endpoint_derivative_degree: int
    section9_derivative_order: int
    stage: int
    spatial_window: int
    coefficient_bound: CertifiedBoundDatum
    singularity_exponent: float
    valid_from: float = LATE_START
    endpoint: float = SECTION10_ENDPOINT
    status: str = "formal-structure"
    actual_section9_sequence_verified: bool = False
    paper_exact_velocity_available: bool = False

    def __post_init__(self) -> None:
        degree = _natural(self.endpoint_derivative_degree, "endpoint derivative degree")
        section9_order = _natural(self.section9_derivative_order, "Section 9 derivative order")
        stage = _natural(self.stage, "Section 9 stage")
        window = _natural(self.spatial_window, "spatial window")
        if section9_order != degree + 1:
            raise ValueError(
                "Section 9 derivative order must equal endpoint derivative degree + 1"
            )
        if not isinstance(self.coefficient_bound, CertifiedBoundDatum):
            raise TypeError("coefficient_bound must be a CertifiedBoundDatum")
        exponent = _finite(self.singularity_exponent, "singularity exponent")
        if not 0.0 <= exponent < 1.0:
            raise ValueError("singularity exponent must satisfy 0 <= alpha < 1")
        valid_from = _finite(self.valid_from, "valid_from")
        endpoint = _finite(self.endpoint, "endpoint")
        if endpoint != SECTION10_ENDPOINT:
            raise ValueError("Section 9 endpoint bridge requires endpoint=1")
        if valid_from < LATE_START:
            raise ValueError("uniform witness must start in the official t>=3/4 plateau")
        if valid_from >= endpoint:
            raise ValueError("valid_from must be strictly before endpoint")

        object.__setattr__(self, "endpoint_derivative_degree", degree)
        object.__setattr__(self, "section9_derivative_order", section9_order)
        object.__setattr__(self, "stage", stage)
        object.__setattr__(self, "spatial_window", window)
        object.__setattr__(self, "singularity_exponent", exponent)
        object.__setattr__(self, "valid_from", valid_from)
        object.__setattr__(self, "endpoint", endpoint)


@dataclass(frozen=True)
class Section9EndpointBridgeRecord:
    """Structural admission record; it does not prove the supplied witness."""

    endpoint_derivative_degree: int
    section9_derivative_order: int
    stage: int
    spatial_window: int
    evidence_kind: str
    evidence_provenance: str
    majorant: EndpointPowerLawMajorant
    localization_transfer: EndpointLocalizationTransferCertificate
    status: str = "formal-structure"
    uniform_bound_witness_supplied: bool = True
    source_theorem_machine_verified: bool = False
    actual_section9_sequence_verified: bool = False
    endpoint_limit_constructed: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_bridge_ready(self) -> bool:
        return (
            self.section9_derivative_order == self.endpoint_derivative_degree + 1
            and self.majorant.derivative_degree == self.endpoint_derivative_degree
            and self.majorant.spatial_window == self.spatial_window
            and self.localization_transfer.certified
            and bool(self.evidence_provenance.strip())
            and not self.paper_exact_velocity_available
        )


def admit_section9_uniform_endpoint_witness(
    witness: Section9UniformEndpointDerivativeWitness,
) -> Section9EndpointBridgeRecord:
    """Admit one genuinely uniform Section 9 witness to the Section 10 ladder.

    A :class:`Section9PointwiseBoundCertificate` is rejected explicitly.  Such
    a certificate can support an audit at one ``q`` but cannot establish the
    all-``t``/all-``x`` hypothesis required by the endpoint Cauchy argument.
    """

    if isinstance(witness, Section9PointwiseBoundCertificate):
        raise TypeError(
            "a pointwise Section 9 certificate cannot be promoted to an endpoint majorant"
        )
    if not isinstance(witness, Section9UniformEndpointDerivativeWitness):
        raise TypeError("witness must be a Section9UniformEndpointDerivativeWitness")

    majorant = EndpointPowerLawMajorant(
        derivative_degree=witness.endpoint_derivative_degree,
        spatial_window=witness.spatial_window,
        coefficient=witness.coefficient_bound.upper_bound,
        singularity_exponent=witness.singularity_exponent,
        valid_from=witness.valid_from,
        endpoint=witness.endpoint,
    )
    transfer = section10_endpoint_localization_transfer(majorant)
    record = Section9EndpointBridgeRecord(
        endpoint_derivative_degree=witness.endpoint_derivative_degree,
        section9_derivative_order=witness.section9_derivative_order,
        stage=witness.stage,
        spatial_window=witness.spatial_window,
        evidence_kind=witness.coefficient_bound.kind,
        evidence_provenance=witness.coefficient_bound.provenance,
        majorant=majorant,
        localization_transfer=transfer,
    )
    if not record.formal_bridge_ready:
        raise ArithmeticError("Section 9 to Section 10 bridge invariant failed")
    return record
