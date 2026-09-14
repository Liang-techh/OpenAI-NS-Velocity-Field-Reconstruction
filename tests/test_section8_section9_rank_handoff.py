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
    Section8LocalRankSolveAdmission,
    Section8LocalRankSolveWitness,
)
from openai_ns_reconstruction.section8_section9_rank_handoff import (
    PINNED_RANK_CLASS,
    PINNED_RANK_INPUT,
    Section8Section9RankHandoffAdmission,
    Section8Section9RankHandoffWitness,
)
from openai_ns_reconstruction.section9_actual_stage_estimates import (
    Section9ActualStageEstimatesAdmission,
    Section9ActualStageEstimatesWitness,
)


def _compact_rank():
    return CompactMeanCorrectionAdmission(
        CompactMeanCorrectionWitness(
            rank_data_id="lean:ActualPrimary.rankData",
            context_id="lean:ActualPrimary.commonContext:19",
            state_id="lean:postTemporal:cycle-x",
            debt_id="lean:debt:cycle-x:13:x",
            slow_point_id="lean:slow:x",
            plane_point_id="lean:plane:0",
            angular_background_id="lean:rank-angular",
            axial_background_id="lean:rank-axial",
            application_id="lean-app:rank-rows:13:x",
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
    )


def _intermediate():
    witness = ActualIntermediateRankDebtWitness(
        B=19,
        N0=7,
        cycle_state_id="lean:CycleState:x",
        static_data_id="lean:StaticData:D",
        cycle_invariant_id="lean:Invariant:H",
        particular_inputs_id="lean:particular:inputs",
        step_data_id="lean:RunData.step:5",
        sigma_repr="lean:sigma:5",
        context_id="lean:ActualPrimary.commonContext:19",
        post_temporal_state_id="lean:postTemporal:cycle-x",
        debt_family_id="lean:CorrectionState.debt:postTemporal-x",
        band_index=13,
        slow_point_id="lean:slow:x",
        band_debt_id="lean:debt:cycle-x:13:x",
        application_id="lean-app:afterTemporal_debt:5",
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
    return ActualIntermediateRankDebtAdmission(witness, _compact_rank())


def _local_rank():
    return Section8LocalRankSolveAdmission(
        _intermediate(),
        Section8LocalRankSolveWitness(
            B=19,
            N0=7,
            rank_data_id="lean:ActualPrimary.rankData",
            gauge_data_id="lean:ActualPrimary.commonGauge",
            local_domain_id="lean:slow-domain",
            context_id="lean:ActualPrimary.commonContext:19",
            input_state_id="lean:postTemporal:cycle-x",
            debt_family_id="lean:CorrectionState.debt:postTemporal-x",
            band_index=13,
            slow_point_id="lean:slow:x",
            axial_direction_id="lean:ActualInitialization.axial",
            increment_state_id="lean:rankIncrementState:cycle-x",
            rank_stage_state_id="lean:rankStageState:cycle-x",
            remainder_family_id="lean:LocalRankDefect.remainders:cycle-x",
            application_id="lean-app:solved_rows:cycle-x",
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
        ),
    )


def _stage_estimates():
    return Section9ActualStageEstimatesAdmission(
        Section9ActualStageEstimatesWitness(
            domain_id="lean:D#actual-stage",
            forcing_id="lean:f#actual-stage",
            q_id="lean:q#actual-stage",
            run_data_id="lean:RunData:B19:N7",
            physical_data_id="lean:physical#actual-stage",
            positive_bump_id="lean:bump#actual-stage",
            finite_velocity_family_id="lean:finiteVelocity",
            finite_pressure_family_id="lean:finitePressure",
            finite_forcing_family_id="lean:finiteForcing",
            mixed_physical_data_id="lean:mixedPhysicalData",
            stage_estimates_id="lean:StageEstimates#actual-stage",
            application_id="lean:stageEstimates_of_representations#application",
            producer_kind="lean-formal-export",
            provenance="synthetic theorem metadata for admission regression only",
            q_open_unit_interval_certified=True,
            input_normalized_certified=True,
            output_strip_controlled_certified=True,
            output_mean_balanced_certified=True,
            stress_controlled_certified=True,
            div_osc_controlled_certified=True,
            velocity_increments_controlled_certified=True,
            pressure_increments_controlled_certified=True,
            forcing_increments_controlled_certified=True,
            forcing_increment_margin_certified=True,
            ledger_velocity_identification_certified=True,
            ledger_pressure_identification_certified=True,
            ledger_forcing_identification_certified=True,
            theorem_application_certified=True,
        )
    )


def _witness(**changes):
    base = Section8Section9RankHandoffWitness(
        B=19,
        N0=7,
        stage_index=5,
        cycle_state_id="lean:CycleState:x",
        post_temporal_state_id="lean:postTemporal:cycle-x",
        debt_family_id="lean:CorrectionState.debt:postTemporal-x",
        rank_stage_state_id="lean:rankStageState:cycle-x",
        remainder_family_id="lean:LocalRankDefect.remainders:cycle-x",
        run_data_id="lean:RunData:B19:N7",
        rank_class_id="lean:RunData.rank_class:5",
        rank_input_id="lean:ActualStageEstimates.rankInput:5",
        stage_estimates_id="lean:StageEstimates#actual-stage",
        rank_class_application_id="lean-app:RunData.rank_class:5",
        rank_input_application_id="lean-app:rankInput:5",
        producer_kind="lean-formal-export",
        provenance="synthetic theorem metadata for admission regression only",
        same_cycle_state_certified=True,
        same_post_temporal_state_certified=True,
        same_debt_family_certified=True,
        rank_class_uses_after_temporal_debt_certified=True,
        rank_input_uses_rank_class_certified=True,
        local_rank_solve_uses_same_debt_certified=True,
        post_rank_remainder_identity_certified=True,
        same_run_data_certified=True,
        same_stage_estimates_certified=True,
        handoff_application_certified=True,
    )
    return replace(base, **changes)


def _admission(**changes):
    return Section8Section9RankHandoffAdmission(
        _local_rank(), _stage_estimates(), _witness(**changes)
    )


def test_binds_same_rank_debt_into_section9_without_truth_upgrade():
    admitted = _admission()
    assert admitted.same_after_temporal_debt_bound_to_section8_and_section9
    assert admitted.section9_rank_class_admitted
    assert admitted.section9_rank_input_identity_admitted
    assert admitted.post_rank_remainder_identity_preserved
    assert admitted.status == "formal-structure"
    assert not admitted.theorem_chain_machine_replayed
    assert not admitted.actual_rank_increment_materialized
    assert not admitted.actual_rank_stage_state_materialized
    assert not admitted.actual_section9_rank_input_materialized
    assert not admitted.section9_iteration_machine_materialized
    assert not admitted.paper_exact_velocity_available


def test_pins_actual_rank_class_and_rank_input_symbols():
    assert PINNED_RANK_CLASS.endswith("RunData.rank_class")
    assert PINNED_RANK_INPUT.endswith("ActualStageEstimates.rankInput")


@pytest.mark.parametrize("producer", ["sampled", "fitted", "numeric-scan", "formal-theorem"])
def test_rejects_nonexport_evidence(producer):
    with pytest.raises(ValueError, match="lean-formal-export"):
        _witness(producer_kind=producer)


@pytest.mark.parametrize(
    "field",
    [
        "same_cycle_state_certified",
        "same_post_temporal_state_certified",
        "same_debt_family_certified",
        "rank_class_uses_after_temporal_debt_certified",
        "rank_input_uses_rank_class_certified",
        "local_rank_solve_uses_same_debt_certified",
        "post_rank_remainder_identity_certified",
        "same_run_data_certified",
        "same_stage_estimates_certified",
        "handoff_application_certified",
    ],
)
def test_rejects_missing_formal_certificate(field):
    with pytest.raises(ValueError, match=field):
        _witness(**{field: False})


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"B": 20}, "B/N0"),
        ({"N0": 8}, "B/N0"),
        ({"cycle_state_id": "lean:other-cycle"}, "cycle state"),
        ({"post_temporal_state_id": "lean:other-postTemporal"}, "post-temporal"),
        ({"debt_family_id": "lean:other-debt"}, "debt family"),
        ({"rank_stage_state_id": "lean:other-rank-stage"}, "rank-stage"),
        ({"remainder_family_id": "lean:other-remainders"}, "remainder family"),
        ({"run_data_id": "lean:other-RunData"}, "RunData"),
        ({"stage_estimates_id": "lean:other-estimates"}, "StageEstimates"),
    ],
)
def test_rejects_cross_wired_handoff(changes, message):
    with pytest.raises(ValueError, match=message):
        _admission(**changes)


def test_rejects_stale_formal_source_and_bad_naturals():
    with pytest.raises(ValueError, match="formal_commit"):
        _witness(formal_commit="deadbeef")
    with pytest.raises(ValueError, match="stage_index"):
        _witness(stage_index=-1)
    with pytest.raises(ValueError, match="B"):
        _witness(B=True)
