from dataclasses import replace

import pytest

from openai_ns_reconstruction.compact_mean_correction_formal import (
    PINNED_FORMAL_COMMIT,
    PINNED_RANK_ROWS_THEOREM,
    CompactMeanCorrectionAdmission,
    CompactMeanCorrectionWitness,
)


def _witness(**changes):
    base = CompactMeanCorrectionWitness(
        rank_data_id="lean:RankData:primary",
        context_id="lean:Context:canonical-base",
        state_id="lean:State:incoming-cycle",
        debt_id="lean:CorrectionState.debt:band-13:slow-point",
        slow_point_id="lean:slow-point:x",
        plane_point_id="lean:PressureStream.Plane:zero",
        angular_background_id="lean:rank-angular-model",
        axial_background_id="lean:rank-axial-model",
        application_id="lean-app:rank_rows_on_patch:band-13:x",
        producer_kind="lean-formal-export",
        provenance="synthetic theorem metadata for admission regression only",
        band_index=13,
        lambda_positive_certified=True,
        coefficient_nonzero_certified=True,
        inner_positive_certified=True,
        inner_lt_outer_certified=True,
        length_positive_certified=True,
        velocity_nonzero_certified=True,
        angular_background_agreement_on_support_certified=True,
        axial_background_agreement_on_support_certified=True,
        physical_five_rows_dependency_certified=True,
        theorem_application_certified=True,
    )
    return replace(base, **changes)


def test_admits_exact_formal_application_metadata():
    admission = CompactMeanCorrectionAdmission(_witness())

    assert admission.compact_five_row_theorem_admitted
    assert admission.actual_state_debt_identity_bound
    assert admission.status == "formal-structure"
    assert not admission.theorem_application_machine_replayed
    assert not admission.actual_debt_values_materialized
    assert not admission.actual_bump_values_materialized
    assert not admission.finite_head_signed_defect_link_verified
    assert not admission.compact_mean_correction_field_materialized
    assert not admission.paper_exact_velocity_available


@pytest.mark.parametrize("producer", ["sampled", "fitted", "numeric-scan"])
def test_rejects_nonformal_evidence(producer):
    with pytest.raises(ValueError, match="lean-formal-export"):
        _witness(producer_kind=producer)


def test_rejects_stale_formal_commit():
    assert PINNED_FORMAL_COMMIT
    with pytest.raises(ValueError, match="formal_commit"):
        _witness(formal_commit="deadbeef")


def test_rejects_wrong_theorem_symbol():
    assert PINNED_RANK_ROWS_THEOREM.endswith("rank_rows_on_patch")
    with pytest.raises(ValueError, match="theorem_symbol"):
        _witness(theorem_symbol="NavierStokes.CorrectionState.rank_five_rows")


def test_rejects_missing_five_row_dependency():
    with pytest.raises(ValueError, match="dependency_symbols"):
        _witness(dependency_symbols=(PINNED_RANK_ROWS_THEOREM,))


def test_rejects_negative_or_boolean_band_index():
    with pytest.raises(ValueError, match="nonnegative"):
        _witness(band_index=-1)
    with pytest.raises(TypeError, match="integer"):
        _witness(band_index=True)


@pytest.mark.parametrize(
    "field",
    [
        "lambda_positive_certified",
        "coefficient_nonzero_certified",
        "inner_positive_certified",
        "inner_lt_outer_certified",
        "length_positive_certified",
        "velocity_nonzero_certified",
        "angular_background_agreement_on_support_certified",
        "axial_background_agreement_on_support_certified",
        "physical_five_rows_dependency_certified",
        "theorem_application_certified",
    ],
)
def test_rejects_uncertified_formal_hypothesis(field):
    with pytest.raises(ValueError, match=field):
        _witness(**{field: False})


@pytest.mark.parametrize(
    "field",
    [
        "rank_data_id",
        "context_id",
        "state_id",
        "debt_id",
        "slow_point_id",
        "plane_point_id",
        "angular_background_id",
        "axial_background_id",
        "application_id",
        "provenance",
    ],
)
def test_rejects_empty_opaque_identity(field):
    with pytest.raises(ValueError, match=field):
        _witness(**{field: "   "})
