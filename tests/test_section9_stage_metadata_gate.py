from decimal import Decimal
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


SOURCE_REVISION = "a" * 40


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


def _metadata(index: int = 7, h: Fraction = Fraction(1, 4), **changes):
    values = dict(
        stage_index=index,
        stage_id=f"stage-{index}",
        potential_family_id="finite-velocity-family",
        direct_family_id="finite-velocity-family",
        pressure_family_id="finite-pressure-family",
        source_revision=SOURCE_REVISION,
        producer_artifact_id=f"producer-stage-{index}",
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


def test_exact_threshold_metadata_replays_all_five_gain_inequalities():
    h = Fraction(1, 4)
    admission = Section9StageMetadataAdmission(_stage(), _metadata(h=h), h=h)

    assert admission.exact_stage_metadata_premises_verified
    assert admission.materialized_stage_identity_consumed
    assert admission.achieved_physical_exponents == admission.ledger.physical_component_exponents
    assert admission.gain_slacks == admission.ledger.component_slacks
    assert all(slack > 0 for slack in admission.gain_slacks.values())


def test_stronger_metadata_keeps_exact_nonnegative_gain_slack():
    h = Fraction(1, 5)
    md = _metadata(
        h=h,
        wave_potential_alpha=wave_native(PINNED_KAPPA, 7) + Fraction(1, 3),
        wave_potential_shift=-h + Fraction(1, 7),
        mean_stream_alpha=mean_native(PINNED_KAPPA, 7) + Fraction(2, 9),
        direct_angular_alpha=mean_native(PINNED_KAPPA, 7) + Fraction(3, 11),
        wave_pressure_alpha=wave_pressure_native(PINNED_KAPPA, 7) + Fraction(4, 13),
        wave_pressure_shift=-2 * coordinate_A(h) + Fraction(1, 17),
        mean_pressure_alpha=mean_native(PINNED_KAPPA, 7) + Fraction(5, 19),
    )
    admission = Section9StageMetadataAdmission(_stage(), md, h=h)
    assert all(slack > 0 for slack in admission.gain_slacks.values())
    assert admission.achieved_physical_exponents["wavePotential"] > admission.ledger.physical_component_exponents["wavePotential"]


@pytest.mark.parametrize(
    "field,required",
    [
        ("wave_potential_alpha", wave_native(PINNED_KAPPA, 7)),
        ("mean_stream_alpha", mean_native(PINNED_KAPPA, 7)),
        ("direct_angular_alpha", mean_native(PINNED_KAPPA, 7)),
        ("wave_pressure_alpha", wave_pressure_native(PINNED_KAPPA, 7)),
        ("mean_pressure_alpha", mean_native(PINNED_KAPPA, 7)),
    ],
)
def test_each_native_alpha_lower_bound_is_fail_closed(field, required):
    with pytest.raises(ValueError, match=field):
        Section9StageMetadataAdmission(
            _stage(),
            _metadata(**{field: required - Fraction(1, 1000)}),
            h=Fraction(1, 4),
        )


@pytest.mark.parametrize(
    "field,required",
    [
        ("wave_potential_shift", -Fraction(1, 4)),
        ("wave_pressure_shift", -2 * coordinate_A(Fraction(1, 4))),
    ],
)
def test_each_shift_lower_bound_is_fail_closed(field, required):
    with pytest.raises(ValueError, match=field):
        Section9StageMetadataAdmission(
            _stage(),
            _metadata(**{field: required - Fraction(1, 1000)}),
            h=Fraction(1, 4),
        )


def test_stage_identity_and_family_cross_wires_are_rejected():
    with pytest.raises(ValueError, match="stage_id"):
        Section9StageMetadataAdmission(
            _stage(), _metadata(stage_id="another-stage"), h=Fraction(1, 4)
        )
    with pytest.raises(ValueError, match="potential family"):
        Section9StageMetadataAdmission(
            _stage(), _metadata(potential_family_id="another-family"), h=Fraction(1, 4)
        )


def test_stage_zero_is_not_admitted_as_increment_metadata():
    with pytest.raises(ValueError, match="at least 1"):
        _metadata(index=0)


def test_metadata_requires_exact_source_revision_and_exact_arithmetic():
    with pytest.raises(ValueError, match="40-character"):
        _metadata(source_revision="not-a-commit")
    with pytest.raises(TypeError, match="exact integer/Fraction"):
        _metadata(wave_potential_alpha=0.9)
    with pytest.raises(TypeError, match="exact integer/Fraction"):
        Section9StageMetadataAdmission(_stage(), _metadata(), h=Decimal("0.25"))


def test_metadata_kind_and_pinned_kappa_are_fail_closed():
    with pytest.raises(ValueError, match="sampled/fitted/toy"):
        _metadata(metadata_kind="sampled-metadata")
    with pytest.raises(ValueError, match="1/100000"):
        Section9StageMetadataAdmission(
            _stage(), _metadata(), h=Fraction(1, 4), kappa=Fraction(1, 99999)
        )


def test_truth_boundary_stays_finite_and_nonconvergent():
    admission = Section9StageMetadataAdmission(_stage(), _metadata(), h=Fraction(1, 4))
    assert admission.finite_positive_stage_only
    assert not admission.field_values_or_norms_evaluated
    assert not admission.support_estimates_consumed
    assert not admission.actual_paper_stage_certified
    assert not admission.infinite_correction_sequence_certified
    assert not admission.eq_9_21_summed_field_certified
    assert not admission.paper_exact_velocity_available
