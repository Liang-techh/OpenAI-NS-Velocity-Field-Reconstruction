import pytest

from openai_ns_reconstruction.section9_selected_schedule_flatness import (
    PINNED_SELECTED_SCHEDULE_THEOREM,
    REQUIRED_FORMAL_SYMBOLS,
    Section9SelectedScheduleFlatnessAdmission,
    Section9SelectedScheduleFlatnessWitness,
)


def _witness(**changes):
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


def test_selected_schedule_flatness_contract_has_no_truth_upgrade():
    admitted = Section9SelectedScheduleFlatnessAdmission(_witness())
    assert admitted.closed_selected_schedule_theorem_admitted is True
    assert admitted.all_quantitative_residual_jet_rates_admitted is True
    assert admitted.status == "formal-structure"
    assert admitted.selected_schedule_theorem_machine_replayed is False
    assert admitted.selected_schedule_values_materialized is False
    assert admitted.finite_stage_fields_materialized is False
    assert admitted.residual_jet_rate_values_materialized is False
    assert admitted.section9_iteration_machine_materialized is False
    assert admitted.paper_exact_velocity_available is False


@pytest.mark.parametrize("producer", ["sampled", "fitted", "numeric-scan", "formal-theorem"])
def test_non_formal_export_producers_are_rejected(producer):
    with pytest.raises(ValueError, match="lean-formal-export"):
        _witness(producer_kind=producer)


@pytest.mark.parametrize(
    "field",
    [
        "selected_estimates_identification_certified",
        "selected_qbig_identification_certified",
        "selected_stage_families_identification_certified",
        "schedule_witness_export_certified",
        "three_cut_bounds_certified",
        "all_residual_jet_rates_certified",
        "theorem_application_certified",
    ],
)
def test_missing_formal_hypothesis_or_conclusion_fails_closed(field):
    with pytest.raises(ValueError, match="formal-export certified true"):
        _witness(**{field: False})


@pytest.mark.parametrize(
    "field",
    [
        "selected_estimates_id",
        "selected_qbig_id",
        "potential_stages_id",
        "direct_stages_id",
        "pressure_stages_id",
        "selected_schedule_id",
        "three_cut_bounds_id",
        "all_residual_jet_rates_id",
        "application_id",
        "provenance",
    ],
)
def test_opaque_formal_identities_must_be_nonempty(field):
    with pytest.raises(ValueError, match=f"nonempty {field}"):
        _witness(**{field: " "})


def test_formal_source_and_symbol_drift_fail_closed():
    with pytest.raises(ValueError, match="formal_commit"):
        _witness(formal_commit="0" * 40)
    with pytest.raises(ValueError, match="formal_file"):
        _witness(formal_file="NavierStokes/Fake.lean")
    with pytest.raises(ValueError, match="theorem_symbol"):
        _witness(theorem_symbol="NavierStokes.fake")
    with pytest.raises(ValueError, match="dependency_symbols"):
        _witness(dependency_symbols=REQUIRED_FORMAL_SYMBOLS[:-1])
    assert _witness().theorem_symbol == PINNED_SELECTED_SCHEDULE_THEOREM


def test_admission_rejects_wrong_object_type():
    with pytest.raises(TypeError, match="Section9SelectedScheduleFlatnessWitness"):
        Section9SelectedScheduleFlatnessAdmission(object())
