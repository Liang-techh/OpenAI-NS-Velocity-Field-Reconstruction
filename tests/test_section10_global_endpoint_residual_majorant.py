from dataclasses import replace

import pytest

from openai_ns_reconstruction.section10_exterior_residual_zero_germ import (
    Section10ExteriorResidualZeroGermWitness,
    admit_exterior_residual_zero_germ,
)
from openai_ns_reconstruction.section10_global_endpoint_residual_majorant import (
    PINNED_RESIDUAL_FLATNESS_FIELD,
    Section10GlobalEndpointResidualMajorantWitness,
    admit_global_endpoint_residual_majorant,
)
from openai_ns_reconstruction.section10_paper_localization_spine import (
    Section10PaperLocalizationSpineWitness,
    admit_paper_localization_spine,
)


def _spine(application_id: str = "lean:paper-localization:v1"):
    witness = Section10PaperLocalizationSpineWitness(
        application_id=application_id,
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


def _zero_germ(spine):
    witness = Section10ExteriorResidualZeroGermWitness(
        application_id="lean-proof-chain:exterior-zero-germ:v1",
        schedule_id=spine.witness.schedule_id,
        local_velocity_id=spine.witness.local_velocity_id,
        local_pressure_id=spine.witness.local_pressure_id,
        residual_id="lean:navierStokesResidual(local:v,local:q):v1",
        provenance="pinned exterior zero-germ chain",
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
    return admit_exterior_residual_zero_germ(spine, witness)


def _valid() -> Section10GlobalEndpointResidualMajorantWitness:
    return Section10GlobalEndpointResidualMajorantWitness(
        application_id="lean-proof-chain:global-small-q-residual-majorant:v1",
        schedule_id="lean:selected-schedule:paper-localization:v1",
        local_velocity_id="lean:local:v:v1",
        local_pressure_id="lean:local:q:v1",
        residual_id="lean:navierStokesResidual(local:v,local:q):v1",
        provenance="inner residual_flatness at Xext plus exterior all-order zero germ",
        same_schedule_as_paper_spine_certified=True,
        same_local_velocity_as_paper_spine_certified=True,
        same_local_pressure_as_paper_spine_certified=True,
        same_residual_as_exterior_zero_germ_certified=True,
        inner_residual_flatness_at_outer_edge_all_orders_certified=True,
        qstar_positive_certified=True,
        outer_all_order_zero_jets_certified=True,
        inner_outer_radius_split_exhaustive_certified=True,
        global_delta_min_positive_certified=True,
        global_small_scale_majorant_all_orders_certified=True,
    )


def test_accepts_universal_inner_outer_majorant_without_runtime_promotion():
    spine = _spine()
    admitted = admit_global_endpoint_residual_majorant(spine, _zero_germ(spine), _valid())
    assert admitted.formal_global_small_scale_residual_majorant_admitted is True
    assert admitted.formal_all_order_residual_decay_family_admitted is True
    assert admitted.status == "formal-structure"
    assert admitted.actual_fields_materialized is False
    assert admitted.actual_majorant_constants_materialized is False
    assert admitted.official_full_late_plateau_majorant_verified is False
    assert admitted.residual_artifact_ready is False
    assert admitted.forcing_artifact_ready is False
    assert admitted.endpoint_residual_closure_verified is False
    assert admitted.paper_exact_velocity_available is False
    assert admitted.full_reconstruction is False


def test_rejects_schedule_field_or_residual_cross_wiring():
    spine = _spine()
    exterior = _zero_germ(spine)
    for field, value in (
        ("schedule_id", "other-schedule"),
        ("local_velocity_id", "other-velocity"),
        ("local_pressure_id", "other-pressure"),
        ("residual_id", "other-residual"),
    ):
        with pytest.raises(ValueError, match=field):
            admit_global_endpoint_residual_majorant(
                spine, exterior, replace(_valid(), **{field: value})
            )


def test_rejects_parent_spine_cross_wiring_even_when_field_ids_match():
    spine = _spine("lean:paper-localization:one")
    other_spine = _spine("lean:paper-localization:two")
    with pytest.raises(ValueError, match="same paper spine"):
        admit_global_endpoint_residual_majorant(
            spine, _zero_germ(other_spine), _valid()
        )


def test_rejects_finite_sampled_manufactured_or_f_equals_r_shortcuts():
    with pytest.raises(ValueError, match="finite derivative"):
        replace(_valid(), max_derivative_order=32)
    with pytest.raises(ValueError, match="finite decay"):
        replace(_valid(), max_decay_order=32)
    with pytest.raises(ValueError, match="sampled/fitted"):
        replace(_valid(), sampled_or_fitted_evidence=True)
    with pytest.raises(ValueError, match="manufactured residual"):
        replace(_valid(), manufactured_residual=True)
    with pytest.raises(ValueError, match="f=R"):
        replace(_valid(), force_defined_as_residual=True)


def test_rejects_missing_inner_outer_or_positive_delta_steps():
    for field in (
        "inner_residual_flatness_at_outer_edge_all_orders_certified",
        "outer_all_order_zero_jets_certified",
        "inner_outer_radius_split_exhaustive_certified",
        "global_delta_min_positive_certified",
        "global_small_scale_majorant_all_orders_certified",
    ):
        with pytest.raises(ValueError, match=field):
            replace(_valid(), **{field: False})


def test_rejects_formal_theorem_drift():
    with pytest.raises(ValueError, match="formal_commit"):
        replace(_valid(), formal_commit="deadbeef")
    with pytest.raises(ValueError, match="residual_flatness_field"):
        replace(
            _valid(), residual_flatness_field=PINNED_RESIDUAL_FLATNESS_FIELD + ".approx"
        )
