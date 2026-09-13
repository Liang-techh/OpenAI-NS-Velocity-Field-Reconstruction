import pytest

from openai_ns_reconstruction.section9_selected_estimates_binding import (
    PINNED_FORMAL_FILES,
    REQUIRED_FORMAL_SYMBOLS,
    Section9SelectedEstimatesBindingAdmission,
    Section9SelectedEstimatesBindingWitness,
)
from openai_ns_reconstruction.section9_selected_schedule_flatness import (
    Section9SelectedScheduleFlatnessAdmission,
    Section9SelectedScheduleFlatnessWitness,
)


def _schedule_witness(**changes):
    kwargs = dict(
        selected_estimates_id="lean:LocalResidualFlatness.selectedEstimates#closed",
        selected_qbig_id="lean:ActualCandidateConstruction.selectedQbig#closed",
        potential_stages_id="lean:ActualCandidateAssembly.selectedPotentialStages#closed",
        direct_stages_id="lean:ActualCandidateAssembly.selectedDirectStages#closed",
        pressure_stages_id="lean:ActualCandidateAssembly.selectedPressureStages#closed",
        selected_schedule_id="lean:selected_schedule#a",
        three_cut_bounds_id="lean:selected_schedule#ThreeCutBounds",
        all_residual_jet_rates_id="lean:selected_schedule#AllResidualJetRates",
        application_id="lean:LocalResidualFlatness.selected_schedule#application",
        producer_kind="lean-formal-export",
        provenance="synthetic formal-export metadata; not a replayed Lean theorem",
        selected_estimates_identification_certified=True,
        selected_qbig_identification_certified=True,
        selected_stage_families_identification_certified=True,
        schedule_witness_export_certified=True,
        three_cut_bounds_certified=True,
        all_residual_jet_rates_certified=True,
        theorem_application_certified=True,
    )
    kwargs.update(changes)
    return Section9SelectedScheduleFlatnessWitness(**kwargs)


def _binding_witness(**changes):
    kwargs = dict(
        selected_budget_id="lean:ActualCandidateConstruction.selectedBudget#closed",
        selected_threshold_id="lean:ActualCandidateConstruction.selectedThreshold#closed",
        selected_threshold_geometry_id=(
            "lean:ActualCandidateConstruction.selectedThreshold_geometry#closed"
        ),
        selected_qbig_id="lean:ActualCandidateConstruction.selectedQbig#closed",
        candidate_estimates_id="lean:ActualCandidateAssembly.estimates#closed",
        selected_estimates_id="lean:LocalResidualFlatness.selectedEstimates#closed",
        potential_stages_id="lean:ActualCandidateAssembly.selectedPotentialStages#closed",
        direct_stages_id="lean:ActualCandidateAssembly.selectedDirectStages#closed",
        pressure_stages_id="lean:ActualCandidateAssembly.selectedPressureStages#closed",
        candidate_estimates_application_id=(
            "lean:ActualCandidateAssembly.estimates#selected-specialization"
        ),
        binding_application_id=(
            "lean:LocalResidualFlatness.selectedEstimates#unfolded-binding"
        ),
        producer_kind="lean-formal-export",
        provenance="synthetic formal-export metadata; not a replayed Lean theorem",
        selected_parameters_identified_certified=True,
        selected_threshold_geometry_certified=True,
        selected_qbig_specialization_certified=True,
        selected_stage_families_specialization_certified=True,
        candidate_estimates_actual_constructor_certified=True,
        selected_estimates_eq_candidate_certified=True,
        binding_application_certified=True,
    )
    kwargs.update(changes)
    return Section9SelectedEstimatesBindingWitness(**kwargs)


def _admission(**binding_changes):
    return Section9SelectedEstimatesBindingAdmission(
        Section9SelectedScheduleFlatnessAdmission(_schedule_witness()),
        _binding_witness(**binding_changes),
    )


def test_closed_selected_estimates_binding_has_no_truth_upgrade():
    admitted = _admission()
    assert admitted.selected_schedule_bound_to_closed_actual_estimates is True
    assert admitted.selected_stage_families_bound_to_closed_actual_choice is True
    assert admitted.status == "formal-structure"
    assert admitted.closed_estimate_binding_machine_replayed is False
    assert admitted.selected_estimates_values_materialized is False
    assert admitted.selected_stage_fields_materialized is False
    assert admitted.selected_schedule_values_materialized is False
    assert admitted.section9_iteration_machine_materialized is False
    assert admitted.paper_exact_velocity_available is False


@pytest.mark.parametrize(
    "field,value",
    [
        ("selected_estimates_id", "lean:other#estimates"),
        ("selected_qbig_id", "lean:other#qbig"),
        ("potential_stages_id", "lean:other#potential"),
        ("direct_stages_id", "lean:other#direct"),
        ("pressure_stages_id", "lean:other#pressure"),
    ],
)
def test_cross_wired_selected_schedule_inputs_are_rejected(field, value):
    with pytest.raises(ValueError, match=field):
        _admission(**{field: value})


@pytest.mark.parametrize("producer", ["sampled", "fitted", "numeric-scan", "formal-theorem"])
def test_non_formal_export_producers_are_rejected(producer):
    with pytest.raises(ValueError, match="lean-formal-export"):
        _binding_witness(producer_kind=producer)


@pytest.mark.parametrize(
    "field",
    [
        "selected_parameters_identified_certified",
        "selected_threshold_geometry_certified",
        "selected_qbig_specialization_certified",
        "selected_stage_families_specialization_certified",
        "candidate_estimates_actual_constructor_certified",
        "selected_estimates_eq_candidate_certified",
        "binding_application_certified",
    ],
)
def test_missing_closed_specialization_certificate_fails_closed(field):
    with pytest.raises(ValueError, match="formal-export certified true"):
        _binding_witness(**{field: False})


@pytest.mark.parametrize(
    "field",
    [
        "selected_budget_id",
        "selected_threshold_id",
        "selected_threshold_geometry_id",
        "selected_qbig_id",
        "candidate_estimates_id",
        "selected_estimates_id",
        "potential_stages_id",
        "direct_stages_id",
        "pressure_stages_id",
        "candidate_estimates_application_id",
        "binding_application_id",
        "provenance",
    ],
)
def test_closed_formal_identities_must_be_nonempty(field):
    with pytest.raises(ValueError, match=f"nonempty {field}"):
        _binding_witness(**{field: " "})


def test_formal_source_and_symbol_drift_fail_closed():
    with pytest.raises(ValueError, match="formal_commit"):
        _binding_witness(formal_commit="0" * 40)
    with pytest.raises(ValueError, match="formal_files"):
        _binding_witness(formal_files=PINNED_FORMAL_FILES[:-1])
    with pytest.raises(ValueError, match="dependency_symbols"):
        _binding_witness(dependency_symbols=REQUIRED_FORMAL_SYMBOLS[:-1])


def test_wrong_composition_types_are_rejected():
    schedule = Section9SelectedScheduleFlatnessAdmission(_schedule_witness())
    witness = _binding_witness()
    with pytest.raises(TypeError, match="Section9SelectedScheduleFlatnessAdmission"):
        Section9SelectedEstimatesBindingAdmission(object(), witness)
    with pytest.raises(TypeError, match="Section9SelectedEstimatesBindingWitness"):
        Section9SelectedEstimatesBindingAdmission(schedule, object())
