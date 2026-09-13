from dataclasses import replace

import pytest

from openai_ns_reconstruction.actual_intermediate_rank_debt import (
    ActualIntermediateRankDebtAdmission,
    ActualIntermediateRankDebtWitness,
)
from openai_ns_reconstruction.compact_mean_correction_formal import (
    CompactMeanCorrectionAdmission,
    CompactMeanCorrectionWitness,
)
from openai_ns_reconstruction.section8_local_rank_solve import (
    PINNED_DEBT_REMAINDERS_THEOREM,
    PINNED_PRESERVE_MASSES_THEOREM,
    PINNED_SOLVED_ROWS_THEOREM,
    Section8LocalRankSolveAdmission,
    Section8LocalRankSolveWitness,
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


def _intermediate_witness(**changes):
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


def _intermediate_admission(intermediate=None, rank=None):
    return ActualIntermediateRankDebtAdmission(
        intermediate or _intermediate_witness(),
        CompactMeanCorrectionAdmission(rank or _rank_witness()),
    )


def _witness(**changes):
    base = Section8LocalRankSolveWitness(
        B=19,
        N0=7,
        rank_data_id="lean:ActualPrimary.rankData",
        gauge_data_id="lean:ActualInitialization.gauge",
        local_domain_id="lean:ActualInitialization.slow-domain",
        context_id="lean:ActualPrimary.commonContext:19",
        input_state_id="lean:ActualIntermediateDebtBounds.postTemporal:cycle-x",
        debt_family_id="lean:CorrectionState.debt:postTemporal-x",
        band_index=13,
        slow_point_id="lean:slow-point:x",
        axial_direction_id="lean:ActualInitialization.axial",
        increment_state_id="lean:VariableGaugeMean.rankIncrementState:cycle-x",
        rank_stage_state_id="lean:VariableGaugeMean.rankStageState:cycle-x",
        remainder_family_id="lean:LocalRankDefect.remainders:cycle-x",
        application_id="lean-app:LocalRankDefect.solved_rows:cycle-x:band-13:x",
        producer_kind="lean-formal-export",
        provenance="synthetic theorem metadata for admission regression only",
        rank_geometry_certified=True,
        local_domain_open_certified=True,
        support_window_certified=True,
        local_operators_certified=True,
        base_smooth_certified=True,
        angular_base_slow_certified=True,
        axial_base_slow_certified=True,
        current_mean_local_triple_certified=True,
        covariance_local_shell_certified=True,
        five_rows_application_certified=True,
        solved_rows_application_certified=True,
        preserve_masses_application_certified=True,
        debt_eq_remainders_application_certified=True,
        support_application_certified=True,
    )
    return replace(base, **changes)


def _admission(intermediate=None, witness=None):
    return Section8LocalRankSolveAdmission(
        intermediate or _intermediate_admission(), witness or _witness()
    )


def test_admits_local_five_row_solve_without_materializing_field():
    admission = _admission()

    assert admission.local_five_row_solve_admitted
    assert admission.linear_debt_cancellation_admitted
    assert admission.required_masses_preserved_admitted
    assert admission.post_rank_debt_equals_remainders_admitted
    assert admission.moving_support_preserved_admitted
    assert admission.status == "formal-structure"
    assert not admission.theorem_chain_machine_replayed
    assert not admission.actual_rank_increment_materialized
    assert not admission.actual_rank_stage_state_materialized
    assert not admission.actual_remainder_values_materialized
    assert not admission.section9_iteration_input_materialized
    assert not admission.paper_exact_velocity_available


def test_pins_exact_local_rank_theorem_chain():
    assert PINNED_SOLVED_ROWS_THEOREM.endswith("RankGeometry.solved_rows")
    assert PINNED_PRESERVE_MASSES_THEOREM.endswith("RankGeometry.preserve_masses")
    assert PINNED_DEBT_REMAINDERS_THEOREM.endswith("RankGeometry.debt_eq_remainders")


@pytest.mark.parametrize("producer", ["sampled", "fitted", "numeric-scan", "formal-theorem"])
def test_rejects_nonexport_evidence(producer):
    with pytest.raises(ValueError, match="lean-formal-export"):
        _witness(producer_kind=producer)


@pytest.mark.parametrize(
    "field",
    [
        "rank_geometry_certified",
        "local_domain_open_certified",
        "support_window_certified",
        "local_operators_certified",
        "base_smooth_certified",
        "angular_base_slow_certified",
        "axial_base_slow_certified",
        "current_mean_local_triple_certified",
        "covariance_local_shell_certified",
        "five_rows_application_certified",
        "solved_rows_application_certified",
        "preserve_masses_application_certified",
        "debt_eq_remainders_application_certified",
        "support_application_certified",
    ],
)
def test_rejects_missing_formal_hypothesis(field):
    with pytest.raises(ValueError, match=field):
        _witness(**{field: False})


@pytest.mark.parametrize(
    "field",
    [
        "rank_data_id",
        "gauge_data_id",
        "local_domain_id",
        "context_id",
        "input_state_id",
        "debt_family_id",
        "slow_point_id",
        "axial_direction_id",
        "increment_state_id",
        "rank_stage_state_id",
        "remainder_family_id",
        "application_id",
        "provenance",
    ],
)
def test_rejects_empty_opaque_identity(field):
    with pytest.raises(ValueError, match=field):
        _witness(**{field: "  "})


def test_rejects_stale_commit_and_wrong_primary_symbol():
    with pytest.raises(ValueError, match="formal_commit"):
        _witness(formal_commit="deadbeef")
    with pytest.raises(ValueError, match="theorem_symbol"):
        _witness(theorem_symbol=PINNED_PRESERVE_MASSES_THEOREM)


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"B": 20}, "B/N0"),
        ({"N0": 8}, "B/N0"),
        ({"rank_data_id": "lean:other-rank-data"}, "rank data"),
        ({"context_id": "lean:other-context"}, "context"),
        ({"input_state_id": "lean:other-state"}, "input state"),
        ({"debt_family_id": "lean:other-debt-family"}, "debt family"),
        ({"band_index": 14}, "band"),
        ({"slow_point_id": "lean:other-point"}, "slow point"),
    ],
)
def test_rejects_cross_wired_actual_rank_solve(changes, message):
    with pytest.raises(ValueError, match=message):
        _admission(witness=_witness(**changes))


def test_rejects_boolean_or_negative_naturals():
    with pytest.raises(ValueError, match="B"):
        _witness(B=True)
    with pytest.raises(ValueError, match="N0"):
        _witness(N0=-1)
    with pytest.raises(ValueError, match="band_index"):
        _witness(band_index=-1)
