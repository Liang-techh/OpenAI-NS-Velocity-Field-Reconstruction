"""Bind a Section 9 endpoint-majorant ladder to one residual source revision.

The finite Section 9 -> Section 10 endpoint bridge already checks derivative
orders, the official ``t >= 3/4`` time-switch plateau, and a coherent spatial
window.  The endpoint-to-Borel admission, however, accepts its pre-endpoint
residual-jet evaluator separately from that ladder.  Without an explicit source
identity, a future caller could therefore combine theorem majorants from one
residual revision with dense residual jets from another revision and still pass
both interfaces independently.

This module closes only that provenance/identity gap.  A theorem-facing source
witness freezes a stable ``(source_id, source_revision)`` together with the exact
ordered evidence kinds and provenance strings of the already-admitted majorants
and the pre-endpoint full-jet provider.  A binding requires exact
stage/window/evidence-kind/provenance agreement before the Section 10
candidate-to-Borel-prefix admission may consume that provider.

No residual estimate is proved here.  The source witness remains theorem input,
finite candidate balls do not identify the actual endpoint limit, and the
all-order Borel smoothness argument is still absent.  Consequently this layer is
strictly ``formal-structure`` and ``paper_exact_velocity_available`` remains
false.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .endpoint_borel import LocalScale
from .endpoint_jets import FullSpacetimeJetFamily
from .section10_endpoint_borel_admission import (
    Section10EndpointBorelAdmissionCertificate,
    admit_section10_endpoint_candidate_for_borel_prefix,
)
from .section10_endpoint_trace import PastFullSpacetimeJetFamily
from .section9_endpoint_ladder_bridge import Section9EndpointLadderBridgeRecord


_ACCEPTED_SOURCE_EVIDENCE = frozenset({"analytic-theorem", "formal-theorem"})


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


@dataclass(frozen=True)
class Section9ResidualEndpointSourceWitness:
    """Stable identity of the residual source used by one endpoint ladder.

    ``majorant_evidence_kinds`` and ``majorant_evidence_provenance`` are
    deliberately ordered by endpoint degree ``0..N``.  The binding below
    compares both byte-for-byte with the evidence metadata carried by the
    already-admitted Section 9 ladder.  Matching provenance text alone is not
    enough: changing an evidence class (for example from ``paper-derived`` to
    ``certified-numerical``) while reusing the same provenance label must fail
    closed rather than look source-identical.

    The three boolean facts are theorem-facing assertions, not facts proved by
    this Python module.  They are separated so the audit trail distinguishes a
    residual-provider contract, applicability of the majorant theorems, and the
    identity of the dense pre-endpoint jet evaluator.
    """

    stage: int
    spatial_window: int
    source_id: str
    source_revision: str
    evidence_kind: str
    provenance: str
    majorant_evidence_kinds: tuple[str, ...]
    majorant_evidence_provenance: tuple[str, ...]
    past_full_jets: PastFullSpacetimeJetFamily
    residual_provider_contract_certified: bool
    majorants_apply_to_source_certified: bool
    past_jets_same_source_certified: bool

    def __post_init__(self) -> None:
        object.__setattr__(self, "stage", _natural(self.stage, "stage"))
        object.__setattr__(
            self, "spatial_window", _natural(self.spatial_window, "spatial_window")
        )
        object.__setattr__(self, "source_id", _nonempty_text(self.source_id, "source_id"))
        object.__setattr__(
            self,
            "source_revision",
            _nonempty_text(self.source_revision, "source_revision"),
        )
        object.__setattr__(
            self, "provenance", _nonempty_text(self.provenance, "source provenance")
        )
        if self.evidence_kind not in _ACCEPTED_SOURCE_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem or formal-theorem; sampled/fitted evidence is rejected"
            )
        if not isinstance(self.majorant_evidence_kinds, tuple) or not self.majorant_evidence_kinds:
            raise ValueError("majorant_evidence_kinds must be a nonempty tuple")
        normalized_kinds = tuple(
            _nonempty_text(value, f"majorant evidence kind degree {degree}")
            for degree, value in enumerate(self.majorant_evidence_kinds)
        )
        object.__setattr__(self, "majorant_evidence_kinds", normalized_kinds)
        if not isinstance(self.majorant_evidence_provenance, tuple) or not self.majorant_evidence_provenance:
            raise ValueError("majorant_evidence_provenance must be a nonempty tuple")
        normalized = tuple(
            _nonempty_text(value, f"majorant provenance degree {degree}")
            for degree, value in enumerate(self.majorant_evidence_provenance)
        )
        object.__setattr__(self, "majorant_evidence_provenance", normalized)
        if not callable(self.past_full_jets):
            raise TypeError("past_full_jets must be callable")
        for name in (
            "residual_provider_contract_certified",
            "majorants_apply_to_source_certified",
            "past_jets_same_source_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision


@dataclass(frozen=True)
class Section9ResidualEndpointSourceBinding:
    """Bind one finite endpoint ladder to one residual provider revision."""

    bridge: Section9EndpointLadderBridgeRecord
    source: Section9ResidualEndpointSourceWitness

    def __post_init__(self) -> None:
        if not isinstance(self.bridge, Section9EndpointLadderBridgeRecord):
            raise TypeError("bridge must be a Section9EndpointLadderBridgeRecord")
        if not isinstance(self.source, Section9ResidualEndpointSourceWitness):
            raise TypeError("source must be a Section9ResidualEndpointSourceWitness")
        failed = [name for name, ok in self.binding_checks().items() if not ok]
        if failed:
            raise ValueError(
                "uncertified Section 9 endpoint residual-source binding: "
                + ", ".join(failed)
            )

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source.source_key

    @property
    def expected_majorant_evidence_kinds(self) -> tuple[str, ...]:
        return tuple(row.evidence_kind for row in self.bridge.bridge_records)

    @property
    def expected_majorant_provenance(self) -> tuple[str, ...]:
        return tuple(row.evidence_provenance for row in self.bridge.bridge_records)

    def binding_checks(self) -> dict[str, bool]:
        return {
            "endpoint_ladder_bridge_intact": self.bridge.formal_bridge_ready,
            "stage_identity": self.source.stage == self.bridge.stage,
            "spatial_window_identity": self.source.spatial_window == self.bridge.spatial_window,
            "majorant_evidence_kind_identity": (
                self.source.majorant_evidence_kinds
                == tuple(row.evidence_kind for row in self.bridge.bridge_records)
            ),
            "majorant_evidence_kind_count_identity": (
                len(self.source.majorant_evidence_kinds)
                == self.bridge.max_endpoint_degree + 1
            ),
            "majorant_provenance_identity": (
                self.source.majorant_evidence_provenance
                == tuple(row.evidence_provenance for row in self.bridge.bridge_records)
            ),
            "majorant_degree_count_identity": (
                len(self.source.majorant_evidence_provenance)
                == self.bridge.max_endpoint_degree + 1
            ),
            "stable_source_identity_present": bool(
                self.source.source_id and self.source.source_revision
            ),
            "source_theorem_provenance_present": bool(self.source.provenance),
            "residual_provider_contract_certified": (
                self.source.residual_provider_contract_certified is True
            ),
            "majorants_apply_to_source_certified": (
                self.source.majorants_apply_to_source_certified is True
            ),
            "past_jets_same_source_certified": (
                self.source.past_jets_same_source_certified is True
            ),
            "theorem_evidence_not_sampled": (
                self.source.evidence_kind in _ACCEPTED_SOURCE_EVIDENCE
            ),
        }

    @property
    def formal_source_binding_ready(self) -> bool:
        return all(self.binding_checks().values())

    @property
    def source_majorants_derived_from_actual_residual_verified(self) -> bool:
        return False

    @property
    def actual_section9_sequence_verified(self) -> bool:
        return False

    @property
    def endpoint_limits_constructed(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False


@dataclass(frozen=True)
class Section9SourceBoundEndpointBorelCertificate:
    """Source-bound finite endpoint-candidate admission record.

    Passing means only that the same named residual revision supplied the
    theorem-facing majorant evidence metadata and the pre-endpoint dense-jet
    provider, and that the existing finite Section 10 consistency/Borel-prefix
    checks passed.  It deliberately does not upgrade any analytic truth flag.
    """

    binding: Section9ResidualEndpointSourceBinding
    admission: Section10EndpointBorelAdmissionCertificate
    source_majorants_derived_from_actual_residual_verified: bool = False
    actual_section9_residual_limits_verified: bool = False
    endpoint_limit_uniqueness_verified: bool = False
    infinite_borel_right_jets_verified: bool = False
    all_order_borel_smoothness_verified: bool = False
    smooth_compact_forcing_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def source_key(self) -> tuple[str, str]:
        return self.binding.source_key

    @property
    def formal_source_bound_candidate_ready(self) -> bool:
        return (
            self.binding.formal_source_binding_ready
            and self.admission.formal_candidate_ready
            and self.admission.spatial_window == self.binding.bridge.spatial_window
            and self.admission.max_degree == self.binding.bridge.max_endpoint_degree
            and not self.source_majorants_derived_from_actual_residual_verified
            and not self.actual_section9_residual_limits_verified
            and not self.endpoint_limit_uniqueness_verified
            and not self.infinite_borel_right_jets_verified
            and not self.all_order_borel_smoothness_verified
            and not self.smooth_compact_forcing_verified
            and not self.paper_exact_velocity_available
        )


def admit_section10_endpoint_candidate_from_bound_section9_source(
    binding: Section9ResidualEndpointSourceBinding,
    candidate_full_jets: FullSpacetimeJetFamily,
    local_scale: LocalScale,
    x: object,
    y: object,
    z: object,
    t: object,
) -> Section9SourceBoundEndpointBorelCertificate:
    """Run the finite endpoint admission using the source-bound past provider.

    There is intentionally no ``past_full_jets`` argument here.  Once the
    source binding is formed, the existing Section 10 admission can only consume
    the provider frozen into that same source/revision witness.
    """

    if not isinstance(binding, Section9ResidualEndpointSourceBinding):
        raise TypeError("binding must be a Section9ResidualEndpointSourceBinding")
    if not binding.formal_source_binding_ready:
        raise ValueError("Section 9 residual-source binding is not formally ready")
    admission = admit_section10_endpoint_candidate_for_borel_prefix(
        binding.bridge.endpoint_ladder,
        binding.source.past_full_jets,
        candidate_full_jets,
        local_scale,
        x,
        y,
        z,
        t,
    )
    result = Section9SourceBoundEndpointBorelCertificate(
        binding=binding,
        admission=admission,
    )
    if not result.formal_source_bound_candidate_ready:
        raise ArithmeticError("source-bound endpoint Borel admission invariant failed")
    return result
