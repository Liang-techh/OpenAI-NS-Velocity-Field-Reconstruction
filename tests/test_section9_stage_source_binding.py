from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_finite_prefix_fields import (
    MATERIALIZATION_KIND,
    Section9MaterializedStage,
)
from openai_ns_reconstruction.section9_increment_exponent_ledger import (
    PINNED_KAPPA,
    coordinate_A,
    mean_native,
    wave_native,
    wave_pressure_native,
)
from openai_ns_reconstruction.section9_stage_metadata_gate import (
    METADATA_KIND,
    Section9ExactStageMetadata,
    Section9StageMetadataAdmission,
)
from openai_ns_reconstruction.section9_stage_source_binding import (
    SOURCE_BINDING_KIND,
    Section9MaterializedStageSource,
    Section9SourceBoundStageMetadataAdmission,
)


SOURCE_REVISION = "a" * 40
PRODUCER = "agent3-section8-stage-7"


def _stage(index: int = 7) -> Section9MaterializedStage[int]:
    return Section9MaterializedStage(
        index=index,
        potential_family_id="finite-velocity-family",
        direct_family_id="finite-velocity-family",
        pressure_family_id="finite-pressure-family",
        stage_id=f"stage-{index}",
        curl_potential=lambda _: (0.0, 0.0, 0.0),
        direct_velocity=lambda _: (0.0, 0.0, 0.0),
        pressure=lambda _: 0.0,
        materialization_kind=MATERIALIZATION_KIND,
        provenance="test-only callable field fixture",
        field_materialization_certified=True,
    )


def _metadata(stage: Section9MaterializedStage[int], h: Fraction = Fraction(1, 4), **changes):
    index = stage.index
    values = dict(
        stage_index=index,
        stage_id=stage.stage_id,
        potential_family_id=stage.potential_family_id,
        direct_family_id=stage.direct_family_id,
        pressure_family_id=stage.pressure_family_id,
        source_revision=SOURCE_REVISION,
        producer_artifact_id=PRODUCER,
        provenance="test-only exact metadata fixture",
        wave_potential_alpha=wave_native(PINNED_KAPPA, index),
        wave_potential_shift=-h,
        mean_stream_alpha=mean_native(PINNED_KAPPA, index),
        direct_angular_alpha=mean_native(PINNED_KAPPA, index),
        wave_pressure_alpha=wave_pressure_native(PINNED_KAPPA, index),
        wave_pressure_shift=-2 * coordinate_A(h),
        mean_pressure_alpha=mean_native(PINNED_KAPPA, index),
        metadata_kind=METADATA_KIND,
    )
    values.update(changes)
    return Section9ExactStageMetadata(**values)


def _admission(stage: Section9MaterializedStage[int], **metadata_changes):
    h = Fraction(1, 4)
    return Section9StageMetadataAdmission(stage, _metadata(stage, h=h, **metadata_changes), h=h)


def _source(stage: Section9MaterializedStage[int], **changes):
    values = dict(
        stage=stage,
        source_revision=SOURCE_REVISION,
        producer_artifact_id=PRODUCER,
        provenance="test-only source-bound field fixture",
        binding_kind=SOURCE_BINDING_KIND,
    )
    values.update(changes)
    return Section9MaterializedStageSource(**values)


def test_source_revision_and_producer_are_bound_to_same_materialized_stage():
    stage = _stage()
    bound = Section9SourceBoundStageMetadataAdmission(_source(stage), _admission(stage))

    assert bound.exact_stage_source_identity_bound
    assert bound.same_materialized_stage_object_consumed
    assert bound.source_revision_cross_checked
    assert bound.producer_artifact_cross_checked


def test_source_revision_cross_wire_fails_closed():
    stage = _stage()
    with pytest.raises(ValueError, match="source_revision"):
        Section9SourceBoundStageMetadataAdmission(
            _source(stage, source_revision="b" * 40),
            _admission(stage),
        )


def test_producer_artifact_cross_wire_fails_closed():
    stage = _stage()
    with pytest.raises(ValueError, match="producer_artifact_id"):
        Section9SourceBoundStageMetadataAdmission(
            _source(stage, producer_artifact_id="another-producer"),
            _admission(stage),
        )


def test_equivalent_but_distinct_stage_object_cannot_certify_another_field_bundle():
    stage = _stage()
    other_stage = _stage()
    with pytest.raises(ValueError, match="stage object"):
        Section9SourceBoundStageMetadataAdmission(_source(stage), _admission(other_stage))


def test_source_binding_requires_exact_revision_and_machine_materialization_kind():
    stage = _stage()
    with pytest.raises(ValueError, match="40-character"):
        _source(stage, source_revision="not-a-revision")
    with pytest.raises(ValueError, match="sampled/fitted/toy"):
        _source(stage, binding_kind="sampled-stage-source")


def test_truth_boundary_stays_finite_and_nonconvergent():
    stage = _stage()
    bound = Section9SourceBoundStageMetadataAdmission(_source(stage), _admission(stage))

    assert not bound.field_values_or_norms_evaluated
    assert not bound.support_estimates_consumed
    assert not bound.actual_oscillatory_wave_defect_consumed
    assert not bound.infinite_correction_sequence_certified
    assert not bound.eq_9_21_summed_field_certified
    assert not bound.paper_exact_velocity_available
