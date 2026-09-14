from dataclasses import replace

import pytest

from openai_ns_reconstruction.section9_actual_local_residual_flatness import (
    PINNED_POTENTIAL_STAGES,
    PINNED_THEOREM,
    Section9ActualLocalResidualFlatnessWitness,
    admit_actual_local_residual_flatness,
)


def _valid() -> Section9ActualLocalResidualFlatnessWitness:
    return Section9ActualLocalResidualFlatnessWitness(
        schedule_id="lean:selected-schedule:actual-v1",
        application_id="lean:LocalResidualFlatness.selected_schedule:actual-v1",
        provenance="openai formal export at pinned commit",
        closed_actual_theorem_application_certified=True,
        one_schedule_outside_order_quantifiers_certified=True,
        selected_schedule_certified=True,
        three_cut_bounds_certified=True,
        all_residual_jet_rates_certified=True,
        actual_unlocalized_mixed_residual_certified=True,
        open_past_endpoint_filter_certified=True,
    )


def test_accepts_closed_actual_all_order_export_without_promoting_runtime_truth():
    admitted = admit_actual_local_residual_flatness(_valid())
    assert admitted.actual_selected_schedule_all_order_flatness_admitted is True
    assert admitted.status == "formal-structure"
    assert admitted.theorem_machine_replayed is False
    assert admitted.actual_fields_materialized is False
    assert admitted.official_late_plateau_uniform_bounds_verified is False
    assert admitted.support_exterior_all_order_zero_verified is False
    assert admitted.endpoint_residual_closure_verified is False
    assert admitted.forcing_artifact_ready is False
    assert admitted.paper_exact_velocity_available is False
    assert admitted.full_reconstruction is False


def test_rejects_finite_derivative_or_decay_frontier():
    with pytest.raises(ValueError, match="finite derivative"):
        replace(_valid(), max_derivative_order=64)
    with pytest.raises(ValueError, match="finite decay-order"):
        replace(_valid(), max_decay_order=32.0)


def test_rejects_sampled_or_manufactured_evidence():
    with pytest.raises(ValueError, match="sampled/fitted"):
        replace(_valid(), sampled_or_fitted_evidence=True)
    with pytest.raises(ValueError, match="manufactured residual"):
        replace(_valid(), manufactured_residual=True)


def test_rejects_formal_source_or_theorem_drift():
    with pytest.raises(ValueError, match="formal_commit"):
        replace(_valid(), formal_commit="deadbeef")
    with pytest.raises(ValueError, match="theorem_symbol"):
        replace(_valid(), theorem_symbol=PINNED_THEOREM + ".finite")


def test_rejects_actual_stage_family_cross_wiring():
    with pytest.raises(ValueError, match="potential_stage_family_id"):
        replace(_valid(), potential_stage_family_id=PINNED_POTENTIAL_STAGES + ".other")


def test_rejects_missing_closed_theorem_outputs():
    with pytest.raises(ValueError, match="one_schedule_outside_order_quantifiers_certified"):
        replace(_valid(), one_schedule_outside_order_quantifiers_certified=False)
    with pytest.raises(ValueError, match="all_residual_jet_rates_certified"):
        replace(_valid(), all_residual_jet_rates_certified=False)
