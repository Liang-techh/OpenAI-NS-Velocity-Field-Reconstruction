from dataclasses import replace

import pytest

from openai_ns_reconstruction.section10_exterior_residual_zero_germ import (
    PINNED_EXTERIOR_POINT_THEOREM,
    Section10ExteriorResidualZeroGermWitness,
    admit_exterior_residual_zero_germ,
)
from openai_ns_reconstruction.section10_paper_localization_spine import (
    Section10PaperLocalizationSpineWitness,
    admit_paper_localization_spine,
)


def _spine():
    witness = Section10PaperLocalizationSpineWitness(
        application_id="lean:PaperLocalization.local_theorem_with_compact_candidate:v1",
        schedule_id="lean:selected-schedule:paper-localization:v1",
        local_potential_id="lean:local:A:v1",
        local_direct_id="lean:local:D:v1",
        local_velocity_id="lean:local:v:v1",
        local_pressure_id="lean:local:q:v1",
        compact_velocity_id="lean:compact:u:v1",
        compact_pressure_id="lean:compact:p:v1",
        force_id="lean:compact:f:v1",
        compact_support_id="lean:compact:K:v1",
        plateau_open_set_id="lean:SpatialLocalization.plateau:v1",
        provenance="openai formal export at pinned commit",
        closed_theorem_application_certified=True,
        same_selected_schedule_certified=True,
        local_properties_certified=True,
        local_away_extensions_certified=True,
        local_all_order_residual_flatness_certified=True,
        local_exterior_residual_zero_certified=True,
        compact_candidate_properties_certified=True,
        compact_force_smooth_certified=True,
        compact_force_positive_time_support_certified=True,
        compact_spatial_support_certified=True,
        compact_divergence_free_certified=True,
        compact_exact_navier_stokes_certified=True,
        compact_uniform_finite_energy_certified=True,
        compact_speed_unbounded_certified=True,
        no_global_finite_energy_solution_certified=True,
        initial_rest_certified=True,
        plateau_open_neighborhood_contains_origin_certified=True,
        official_late_plateau_same_fields_certified=True,
    )
    return admit_paper_localization_spine(witness)


def _valid() -> Section10ExteriorResidualZeroGermWitness:
    return Section10ExteriorResidualZeroGermWitness(
        application_id="lean-proof-chain:exterior-residual-zero-germ:v1",
        schedule_id="lean:selected-schedule:paper-localization:v1",
        local_velocity_id="lean:local:v:v1",
        local_pressure_id="lean:local:q:v1",
        residual_id="lean:navierStokesResidual(local:v,local:q):v1",
        provenance="derived from pinned exterior zero-germ proof chain",
        same_schedule_as_paper_spine_certified=True,
        same_local_velocity_as_paper_spine_certified=True,
        same_local_pressure_as_paper_spine_certified=True,
        literal_navier_stokes_residual_certified=True,
        strict_exterior_conditions_certified=True,
        exterior_conditions_persist_on_neighborhood_certified=True,
        pointwise_residual_zero_on_persisting_neighborhood_certified=True,
        residual_zero_germ_certified=True,
        all_order_exterior_residual_jets_zero_certified=True,
    )


def test_accepts_same_spine_all_order_zero_germ_without_runtime_promotion():
    admitted = admit_exterior_residual_zero_germ(_spine(), _valid())
    assert admitted.formal_exterior_residual_zero_germ_admitted is True
    assert admitted.formal_all_order_exterior_residual_jets_zero_admitted is True
    assert admitted.status == "formal-structure"
    assert admitted.theorem_machine_replayed is False
    assert admitted.actual_fields_materialized is False
    assert admitted.support_exterior_zero_jets_runtime_verified is False
    assert admitted.residual_artifact_ready is False
    assert admitted.forcing_artifact_ready is False
    assert admitted.endpoint_residual_closure_verified is False
    assert admitted.compact_support_runtime_verified is False
    assert admitted.paper_exact_velocity_available is False
    assert admitted.full_reconstruction is False


def test_rejects_schedule_or_local_field_cross_wiring():
    with pytest.raises(ValueError, match="schedule_id"):
        admit_exterior_residual_zero_germ(_spine(), replace(_valid(), schedule_id="other"))
    with pytest.raises(ValueError, match="local_velocity_id"):
        admit_exterior_residual_zero_germ(
            _spine(), replace(_valid(), local_velocity_id="other")
        )
    with pytest.raises(ValueError, match="local_pressure_id"):
        admit_exterior_residual_zero_germ(
            _spine(), replace(_valid(), local_pressure_id="other")
        )


def test_rejects_finite_sampled_manufactured_or_f_equals_r_shortcuts():
    with pytest.raises(ValueError, match="finite derivative"):
        replace(_valid(), max_derivative_order=64)
    with pytest.raises(ValueError, match="sampled/fitted"):
        replace(_valid(), sampled_or_fitted_evidence=True)
    with pytest.raises(ValueError, match="manufactured residual"):
        replace(_valid(), manufactured_residual=True)
    with pytest.raises(ValueError, match="f=R"):
        replace(_valid(), force_defined_as_residual=True)


def test_rejects_formal_source_or_proof_symbol_drift():
    with pytest.raises(ValueError, match="formal_commit"):
        replace(_valid(), formal_commit="deadbeef")
    with pytest.raises(ValueError, match="exterior_point_theorem_symbol"):
        replace(
            _valid(),
            exterior_point_theorem_symbol=PINNED_EXTERIOR_POINT_THEOREM + ".approx",
        )


def test_rejects_missing_neighborhood_or_all_order_derivative_step():
    with pytest.raises(
        ValueError, match="exterior_conditions_persist_on_neighborhood_certified"
    ):
        replace(_valid(), exterior_conditions_persist_on_neighborhood_certified=False)
    with pytest.raises(
        ValueError, match="pointwise_residual_zero_on_persisting_neighborhood_certified"
    ):
        replace(
            _valid(),
            pointwise_residual_zero_on_persisting_neighborhood_certified=False,
        )
    with pytest.raises(ValueError, match="all_order_exterior_residual_jets_zero_certified"):
        replace(_valid(), all_order_exterior_residual_jets_zero_certified=False)
