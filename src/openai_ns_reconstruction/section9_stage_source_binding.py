"""Bind a materialized Section 9 stage to the exact source revision that produced it.

PR #343 checks the pinned ``ActualIterationLedger.StageMetadata`` inequalities
against an already materialized stage, but ``Section9MaterializedStage`` itself
predates source-revision metadata.  A syntactically valid revision on the
metadata packet is therefore not, by itself, evidence that the field callables
came from that revision.

This module closes only that identity gap.  It wraps one genuine materialized
stage with a producer artifact/revision identity and admits the existing exact
StageMetadata certificate only when the *same stage object*, source revision,
and producer artifact agree.  No field is synthesized, sampled, fitted, or
reconstructed here, and no finite-stage identity is promoted to summability or
Eq. (9.21).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Generic, TypeVar

from .section9_finite_prefix_fields import Section9MaterializedStage
from .section9_stage_metadata_gate import Section9StageMetadataAdmission


PointT = TypeVar("PointT")
SOURCE_BINDING_KIND = "machine-materialized-stage-source-binding"


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _git_sha(value: object) -> str:
    revision = _text(value, "source_revision")
    if len(revision) != 40 or any(ch not in "0123456789abcdef" for ch in revision):
        raise ValueError("source_revision must be an exact 40-character lowercase git SHA")
    return revision


@dataclass(frozen=True)
class Section9MaterializedStageSource(Generic[PointT]):
    """Producer identity for one already materialized stage field bundle."""

    stage: Section9MaterializedStage[PointT]
    source_revision: str
    producer_artifact_id: str
    provenance: str
    binding_kind: str = SOURCE_BINDING_KIND

    def __post_init__(self) -> None:
        if not isinstance(self.stage, Section9MaterializedStage):
            raise TypeError("stage must be a Section9MaterializedStage")
        object.__setattr__(self, "source_revision", _git_sha(self.source_revision))
        object.__setattr__(
            self, "producer_artifact_id", _text(self.producer_artifact_id, "producer_artifact_id")
        )
        object.__setattr__(self, "provenance", _text(self.provenance, "provenance"))
        if self.binding_kind != SOURCE_BINDING_KIND:
            raise ValueError(
                "binding_kind must identify a machine-materialized stage source binding; "
                "sampled/fitted/toy bindings are rejected"
            )


@dataclass(frozen=True)
class Section9SourceBoundStageMetadataAdmission(Generic[PointT]):
    """Require field and exact StageMetadata to originate from one producer revision."""

    source: Section9MaterializedStageSource[PointT]
    metadata_admission: Section9StageMetadataAdmission[PointT]

    def __post_init__(self) -> None:
        if not isinstance(self.source, Section9MaterializedStageSource):
            raise TypeError("source must be a Section9MaterializedStageSource")
        if not isinstance(self.metadata_admission, Section9StageMetadataAdmission):
            raise TypeError("metadata_admission must be Section9StageMetadataAdmission")

        # Object identity is intentional.  Two separately supplied callable field
        # bundles can compare similarly at selected points while belonging to
        # different producers; this gate certifies the exact admitted stage object.
        if self.metadata_admission.stage is not self.source.stage:
            raise ValueError("materialized stage object is cross-wired against metadata admission")

        metadata = self.metadata_admission.metadata
        if metadata.source_revision != self.source.source_revision:
            raise ValueError("source_revision is cross-wired between stage fields and metadata")
        if metadata.producer_artifact_id != self.source.producer_artifact_id:
            raise ValueError("producer_artifact_id is cross-wired between stage fields and metadata")

    @property
    def exact_stage_source_identity_bound(self) -> bool:
        return True

    @property
    def same_materialized_stage_object_consumed(self) -> bool:
        return True

    @property
    def source_revision_cross_checked(self) -> bool:
        return True

    @property
    def producer_artifact_cross_checked(self) -> bool:
        return True

    @property
    def field_values_or_norms_evaluated(self) -> bool:
        return False

    @property
    def support_estimates_consumed(self) -> bool:
        return False

    @property
    def actual_oscillatory_wave_defect_consumed(self) -> bool:
        return False

    @property
    def infinite_correction_sequence_certified(self) -> bool:
        return False

    @property
    def eq_9_21_summed_field_certified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
