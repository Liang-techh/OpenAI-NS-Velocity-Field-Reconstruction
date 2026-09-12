"""Content-address the Section 9 Eq. (9.18) -> Eq. (9.21) source bridge.

The existing endpoint-majorant source binding freezes a semantic
``(source_id, source_revision)`` together with stage/window/evidence metadata.
Those labels are necessary but do not, by themselves, prevent an integration
caller from accidentally reusing the same revision string after replacing the
concrete residual artifact/provider.

This module adds one deliberately narrow provenance gate.  Every theorem-facing
uniform Eq. (9.18) envelope and the corresponding residual-source witness are
wrapped with the same canonical SHA-256 digest of the concrete residual artifact
being referenced.  Only after exact digest agreement is established does this
module delegate to the already-landed analytic-majorant/source binding path.

The digest is audit metadata, not a proof that the artifact is the manuscript's
actual Eq. (9.21) residual.  This module does not derive any residual estimate,
construct endpoint limits, prove all-order smoothness, or promote any paper-exact
truth flag.  ``paper_exact_velocity_available`` therefore remains false.
"""
from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Iterable

from .section9_endpoint_majorant_adapter import Section9UniformResidualEnvelopeWitness
from .section9_endpoint_majorant_source_binding import (
    Section9DerivedMajorantSourceBindingRecord,
    bind_section9_uniform_envelope_ladder_to_residual_source,
)
from .section9_endpoint_residual_source import Section9ResidualEndpointSourceWitness


_CANONICAL_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


def _canonical_sha256(value: object, name: str) -> str:
    if not isinstance(value, str) or _CANONICAL_SHA256.fullmatch(value) is None:
        raise ValueError(f"{name} must be a canonical lowercase 64-hex SHA-256 digest")
    return value


@dataclass(frozen=True)
class Section9ContentAddressedUniformEnvelopeWitness:
    """Pair one Eq. (9.18) envelope with its concrete residual-artifact digest."""

    envelope: Section9UniformResidualEnvelopeWitness
    residual_artifact_sha256: str

    def __post_init__(self) -> None:
        if not isinstance(self.envelope, Section9UniformResidualEnvelopeWitness):
            raise TypeError("envelope must be a Section9UniformResidualEnvelopeWitness")
        object.__setattr__(
            self,
            "residual_artifact_sha256",
            _canonical_sha256(
                self.residual_artifact_sha256, "residual_artifact_sha256"
            ),
        )


@dataclass(frozen=True)
class Section9ContentAddressedResidualSourceWitness:
    """Pair the frozen residual-source witness with that same artifact digest."""

    source: Section9ResidualEndpointSourceWitness
    residual_artifact_sha256: str

    def __post_init__(self) -> None:
        if not isinstance(self.source, Section9ResidualEndpointSourceWitness):
            raise TypeError("source must be a Section9ResidualEndpointSourceWitness")
        object.__setattr__(
            self,
            "residual_artifact_sha256",
            _canonical_sha256(
                self.residual_artifact_sha256, "residual_artifact_sha256"
            ),
        )


@dataclass(frozen=True)
class Section9ContentAddressedMajorantSourceBindingRecord:
    """Existing derived-majorant binding plus one common residual-artifact digest."""

    residual_artifact_sha256: str
    binding: Section9DerivedMajorantSourceBindingRecord
    status: str = "formal-structure"
    source_majorants_derived_from_actual_residual_verified: bool = False
    actual_section9_sequence_verified: bool = False
    endpoint_limits_constructed: bool = False
    all_order_borel_smoothness_verified: bool = False
    smooth_compact_forcing_verified: bool = False
    paper_exact_velocity_available: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "residual_artifact_sha256",
            _canonical_sha256(
                self.residual_artifact_sha256, "residual_artifact_sha256"
            ),
        )
        if not isinstance(self.binding, Section9DerivedMajorantSourceBindingRecord):
            raise TypeError(
                "binding must be a Section9DerivedMajorantSourceBindingRecord"
            )
        if not self.formal_content_addressed_binding_ready:
            raise ArithmeticError("content-addressed Section 9 source binding invariant failed")

    @property
    def source_key(self) -> tuple[str, str]:
        return self.binding.source_key

    @property
    def formal_content_addressed_binding_ready(self) -> bool:
        return (
            self.binding.formal_source_binding_ready
            and not self.source_majorants_derived_from_actual_residual_verified
            and not self.actual_section9_sequence_verified
            and not self.endpoint_limits_constructed
            and not self.all_order_borel_smoothness_verified
            and not self.smooth_compact_forcing_verified
            and not self.paper_exact_velocity_available
        )


def bind_content_addressed_section9_uniform_envelopes_to_residual_source(
    witnesses: Iterable[Section9ContentAddressedUniformEnvelopeWitness],
    source: Section9ContentAddressedResidualSourceWitness,
    *,
    max_endpoint_degree: int,
) -> Section9ContentAddressedMajorantSourceBindingRecord:
    """Require one artifact digest, then run the existing Eq. (9.18) source gate.

    This function intentionally delegates all degree/stage/window/evidence checks
    to ``bind_section9_uniform_envelope_ladder_to_residual_source``.  Its sole new
    responsibility is fail-closed content identity: every wrapped envelope and
    the wrapped source must name exactly the same canonical SHA-256 digest.
    """

    if not isinstance(source, Section9ContentAddressedResidualSourceWitness):
        raise TypeError("source must be a Section9ContentAddressedResidualSourceWitness")

    rows = tuple(witnesses)
    if not rows:
        raise ValueError("content-addressed uniform residual-envelope ladder must be nonempty")
    if not all(
        isinstance(row, Section9ContentAddressedUniformEnvelopeWitness) for row in rows
    ):
        raise TypeError(
            "every ladder entry must be a Section9ContentAddressedUniformEnvelopeWitness"
        )

    digests = {row.residual_artifact_sha256 for row in rows}
    if len(digests) != 1:
        raise ValueError(
            "all residual-envelope witnesses must use one content-addressed residual artifact fingerprint"
        )
    residual_artifact_sha256 = next(iter(digests))
    if residual_artifact_sha256 != source.residual_artifact_sha256:
        raise ValueError(
            "residual_artifact_sha256 mismatch between envelope ladder and source witness"
        )

    binding = bind_section9_uniform_envelope_ladder_to_residual_source(
        (row.envelope for row in rows),
        source.source,
        max_endpoint_degree=max_endpoint_degree,
    )
    result = Section9ContentAddressedMajorantSourceBindingRecord(
        residual_artifact_sha256=residual_artifact_sha256,
        binding=binding,
    )
    if not result.formal_content_addressed_binding_ready:
        raise ArithmeticError("content-addressed Section 9 source binding invariant failed")
    return result
