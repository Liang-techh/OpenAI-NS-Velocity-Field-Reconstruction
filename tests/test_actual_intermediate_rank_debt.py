from dataclasses import replace

import pytest

from openai_ns_reconstruction.actual_intermediate_rank_debt import (
    PINNED_AFTER_TEMPORAL_DEBT_THEOREM,
    PINNED_FORMAL_COMMIT,
    ActualIntermediateRankDebtAdmission,
    ActualIntermediateRankDebtWitness,
)
from openai_ns_reconstruction.compact_mean_correction_formal import (
    CompactMeanCorrectionAdmission,
    CompactMeanCorrectionWitness,
)


def _rank_witness(**changes):
    base = CompactMeanCorrectionWitness(
        rank_data_id="lean:ActualPrimary.rankData",
        context_id="lean:ActualPrimary.commonContext:19",
        state_id="lean:ActualIntermediateDebtBounds.postTemporal:cycle-x",
        debt_id="lean:CorrectionState.debt:postTemporal-x:band-13:slow-x",
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


def _intermediate(**changes):
    base = ActualIntermediateRankDebtWitness(
        B=19,
        N0=7,
        cycle_state_id="lean:CycleState:x",
        static_data_id="lean:CorrectionAnalyticStep.StaticData:D",
        cycle_invariant_id="lean:CycleAnalyticInvariant:H",
        particular_inputs_id="lean:ActualParticularMeanGain.Inputs:particular",
        step_data_id="lean:CorrectionAnalyticStep.StepData:d",
        sigma_repr="lean:sigma:1/5",
        context_id="lean:ActualPrimary.commonContext:19",
        post_temporal_state_id="lean:ActualIntermediateDebtBounds.postTemporal:cycle-x",
        debt_family_id="lean:CorrectionState.debt:postTemporal-x",
        band_index=13,
        slow_point_id="lean:slow-point:x",
        band_debt_id="lean:CorrectionState.debt:postTemporal-x:band-13:slow-x",
        application_id="lean-app:afterTemporal_debt_from_stepData:cycle-x",
        producer_kind="lean-formal-export",
        provenance="synthetic theorem metadata for admission regression only",
        same_B_N0_certified=True,
        same_cycle_state_certified=True,
        common_context_identity_certified=True,
        post_temporal_state_identity_certified=True,
        step_data_for_invariant_certified=True,
        sigma_lower_bound_certified=True,
        debt_family_identity_certified=True,
        band_debt_evaluation_identity_certified=True,
        theorem_application_certified=True,
    )
    return replace(base, **changes)


def _admission(intermediate=None, rank=None):
    return ActualIntermediateRankDebtAdmission(
        intermediate or _intermediate(),
        CompactMeanCorrectionAdmission(rank or _rank_witness()),
    )


def test_admits_same_actual_post_temporal_debt_instance():
    admission = _admission()

    assert admission.after_temporal_debt_theorem_admitted
    assert admission.actual_post_temporal_debt_bound_to_compact_rank
    assert admission.status == "formal-structure"
    assert not admission.theorem_application_machine_replayed
    assert not admission.actual_cycle_state_materialized
    assert not admission.actual_post_temporal_debt_values_materialized
    assert not admission.finite_head_signed_defect_link_verified
    assert not admission.compact_mean_correction_field_materialized
    assert not admission.paper_exact_velocity_available


@pytest.mark.parametrize("producer", ["sampled", "fitted", "numeric-scan"])
def test_rejects_nonformal_intermediate_evidence(producer):
    with pytest.raises(ValueError, match="lean-formal-export"):
        _intermediate(producer_kind=producer)


def test_rejects_stale_commit_or_theorem_symbol():
    assert PINNED_FORMAL_COMMIT
    assert PINNED_AFTER_TEMPORAL_DEBT_THEOREM.endswith("afterTemporal_debt_from_stepData")
    with pytest.raises(ValueError, match="formal_commit"):
        _intermediate(formal_commit="deadbeef")
    with pytest.raises(ValueError, match="theorem_symbol"):
        _intermediate(theorem_symbol="NavierStokes.ActualIntermediateDebtBounds.stage_debt_bounds")


@pytest.mark.parametrize(
    "field",
    [
        "same_B_N0_certified",
        "same_cycle_state_certified",
        "common_context_identity_certified",
        "post_temporal_state_identity_certified",
        "step_data_for_invariant_certified",
        "sigma_lower_bound_certified",
        "debt_family_identity_certified",
        "band_debt_evaluation_identity_certified",
        "theorem_application_certified",
    ],
)
def test_rejects_uncertified_intermediate_hypothesis(field):
    with pytest.raises(ValueError, match=field):
        _intermediate(**{field: False})


@pytest.mark.parametrize(
    "field",
    [
        "cycle_state_id",
        "static_data_id",
        "cycle_invariant_id",
        "particular_inputs_id",
        "step_data_id",
        "sigma_repr",
        "context_id",
        "post_temporal_state_id",
        "debt_family_id",
        "slow_point_id",
        "band_debt_id",
        "application_id",
        "provenance",
    ],
)
def test_rejects_empty_opaque_identity(field):
    with pytest.raises(ValueError, match=field):
        _intermediate(**{field: "   "})


def test_rejects_boolean_or_negative_naturals():
    with pytest.raises(ValueError, match="B"):
        _intermediate(B=True)
    with pytest.raises(ValueError, match="N0"):
        _intermediate(N0=-1)
    with pytest.raises(ValueError, match="band_index"):
        _intermediate(band_index=-1)


@pytest.mark.parametrize(
    ("intermediate_changes", "rank_changes", "message"),
    [
        ({"context_id": "lean:other-context"}, {}, "context"),
        ({"post_temporal_state_id": "lean:other-state"}, {}, "state"),
        ({"band_index": 14}, {}, "band"),
        ({"slow_point_id": "lean:other-point"}, {}, "slow point"),
        ({"band_debt_id": "lean:other-debt"}, {}, "debt"),
    ],
)
def test_rejects_cross_wired_compact_rank_application(
    intermediate_changes, rank_changes, message
):
    with pytest.raises(ValueError, match=message):
        _admission(_intermediate(**intermediate_changes), _rank_witness(**rank_changes))
