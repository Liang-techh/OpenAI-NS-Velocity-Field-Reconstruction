from dataclasses import replace
from fractions import Fraction

import pytest

from openai_ns_reconstruction.section10_paper_localization_spine import (
    PINNED_PLATEAU,
    PINNED_THEOREM,
    Section10PaperLocalizationSpineWitness,
    admit_paper_localization_spine,
)


def _valid() -> Section10PaperLocalizationSpineWitness:
    return Section10PaperLocalizationSpineWitness(
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


def test_accepts_closed_same_schedule_paper_spine_without_runtime_promotion():
    admitted = admit_paper_localization_spine(_valid())
    assert admitted.formal_paper_localization_spine_admitted is True
    assert admitted.official_late_plateau_field_agreement_formally_bound is True
    assert admitted.final_candidate_properties_formally_bound is True
    assert admitted.status == "formal-structure"
    assert admitted.theorem_machine_replayed is False
    assert admitted.actual_fields_materialized is False
    assert admitted.actual_section9_sequence_verified is False
    assert admitted.section9_all_order_endpoint_limits_verified is False
    assert admitted.residual_artifact_ready is False
    assert admitted.forcing_artifact_ready is False
    assert admitted.endpoint_residual_closure_verified is False
    assert admitted.compact_support_runtime_verified is False
    assert admitted.divergence_free_closure_verified is False
    assert admitted.finite_energy_closure_verified is False
    assert admitted.blow_up_closure_verified is False
    assert admitted.paper_exact_velocity_available is False
    assert admitted.full_reconstruction is False


def test_rejects_nonofficial_plateau_or_inexact_rational_aliases():
    with pytest.raises(ValueError, match="plateau_set_id"):
        replace(_valid(), plateau_set_id=PINNED_PLATEAU + ".generic")
    with pytest.raises(ValueError, match="late_start"):
        replace(_valid(), late_start=Fraction(2, 3))
    with pytest.raises(ValueError, match="late_start"):
        replace(_valid(), late_start=0.75)
    with pytest.raises(ValueError, match="viscosity"):
        replace(_valid(), viscosity=1.0)


def test_rejects_finite_sampled_manufactured_or_generic_substitutes():
    with pytest.raises(ValueError, match="finite derivative"):
        replace(_valid(), max_derivative_order=64)
    with pytest.raises(ValueError, match="sampled/fitted"):
        replace(_valid(), sampled_or_fitted_evidence=True)
    with pytest.raises(ValueError, match="manufactured residual"):
        replace(_valid(), manufactured_residual=True)
    with pytest.raises(ValueError, match="generic localization"):
        replace(_valid(), generic_localization_substitute=True)


def test_rejects_formal_source_or_theorem_drift():
    with pytest.raises(ValueError, match="formal_commit"):
        replace(_valid(), formal_commit="deadbeef")
    with pytest.raises(ValueError, match="theorem_symbol"):
        replace(_valid(), theorem_symbol=PINNED_THEOREM + ".approx")


def test_rejects_missing_same_schedule_or_all_order_local_evidence():
    with pytest.raises(ValueError, match="same_selected_schedule_certified"):
        replace(_valid(), same_selected_schedule_certified=False)
    with pytest.raises(ValueError, match="local_all_order_residual_flatness_certified"):
        replace(_valid(), local_all_order_residual_flatness_certified=False)
    with pytest.raises(ValueError, match="local_away_extensions_certified"):
        replace(_valid(), local_away_extensions_certified=False)


def test_rejects_missing_final_candidate_or_plateau_outputs():
    with pytest.raises(ValueError, match="compact_uniform_finite_energy_certified"):
        replace(_valid(), compact_uniform_finite_energy_certified=False)
    with pytest.raises(ValueError, match="compact_force_positive_time_support_certified"):
        replace(_valid(), compact_force_positive_time_support_certified=False)
    with pytest.raises(ValueError, match="official_late_plateau_same_fields_certified"):
        replace(_valid(), official_late_plateau_same_fields_certified=False)
