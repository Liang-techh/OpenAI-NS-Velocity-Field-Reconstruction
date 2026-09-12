import pytest

from openai_ns_reconstruction.section9_actual_stage_estimates import (
    GAIN_DEFINITION,
    LOSS_DEFINITION,
    PINNED_STAGE_ESTIMATES_THEOREM,
    Q1_DEFINITION,
    REQUIRED_FORMAL_SYMBOLS,
    Section9ActualStageEstimatesAdmission,
    Section9ActualStageEstimatesWitness,
)


def _witness(**changes):
    kwargs = dict(
        domain_id="lean:D#actual-stage",
        forcing_id="lean:f#actual-stage",
        q_id="lean:q#actual-stage",
        run_data_id="lean:data#actual-stage",
        physical_data_id="lean:physical#actual-stage",
        positive_bump_id="lean:bump#actual-stage",
        finite_velocity_family_id="lean:ActualIterationLedger.finiteVelocity#D",
        finite_pressure_family_id="lean:ActualIterationLedger.finitePressure#D",
        finite_forcing_family_id="lean:ActualIterationLedger.finiteForcing#D",
        mixed_physical_data_id="lean:physicalData#stageEstimates",
        stage_estimates_id="lean:StageEstimates#actual-stage",
        application_id="lean:stageEstimates_of_representations#application",
        producer_kind="lean-formal-export",
        provenance="synthetic formal-export metadata; not a replayed Lean theorem",
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
    kwargs.update(changes)
    return Section9ActualStageEstimatesWitness(**kwargs)


def test_actual_stage_estimate_contract_has_no_truth_upgrade():
    admitted = Section9ActualStageEstimatesAdmission(_witness())
    assert admitted.literal_finite_stage_estimates_admitted is True
    assert admitted.canonical_ledger_families_admitted is True
    assert admitted.status == "formal-structure"
    assert admitted.stage_estimates_theorem_machine_replayed is False
    assert admitted.actual_run_data_materialized is False
    assert admitted.actual_physical_data_materialized is False
    assert admitted.finite_stage_fields_materialized is False
    assert admitted.stagewise_quantitative_bounds_materialized is False
    assert admitted.section9_iteration_machine_materialized is False
    assert admitted.paper_exact_velocity_available is False


@pytest.mark.parametrize("producer", ["sampled", "fitted", "numeric-scan", "formal-theorem"])
def test_non_formal_export_producers_are_rejected(producer):
    with pytest.raises(ValueError, match="lean-formal-export"):
        _witness(producer_kind=producer)


@pytest.mark.parametrize(
    "field",
    [
        "q_open_unit_interval_certified",
        "input_normalized_certified",
        "output_strip_controlled_certified",
        "output_mean_balanced_certified",
        "stress_controlled_certified",
        "div_osc_controlled_certified",
        "velocity_increments_controlled_certified",
        "pressure_increments_controlled_certified",
        "forcing_increments_controlled_certified",
        "forcing_increment_margin_certified",
        "ledger_velocity_identification_certified",
        "ledger_pressure_identification_certified",
        "ledger_forcing_identification_certified",
        "theorem_application_certified",
    ],
)
def test_missing_formal_hypothesis_fails_closed(field):
    with pytest.raises(ValueError, match="formal-export certified true"):
        _witness(**{field: False})


@pytest.mark.parametrize(
    "field",
    [
        "domain_id",
        "forcing_id",
        "q_id",
        "run_data_id",
        "physical_data_id",
        "positive_bump_id",
        "finite_velocity_family_id",
        "finite_pressure_family_id",
        "finite_forcing_family_id",
        "mixed_physical_data_id",
        "stage_estimates_id",
        "application_id",
        "provenance",
    ],
)
def test_opaque_formal_identities_must_be_nonempty(field):
    with pytest.raises(ValueError, match=f"nonempty {field}"):
        _witness(**{field: " "})


def test_schedule_definitions_fail_closed_on_drift():
    with pytest.raises(ValueError, match="q1_definition"):
        _witness(q1_definition="q")
    with pytest.raises(ValueError, match="gain_definition"):
        _witness(gain_definition="q ^ g")
    with pytest.raises(ValueError, match="loss_definition"):
        _witness(loss_definition="q ^ g")
    assert _witness().q1_definition == Q1_DEFINITION
    assert _witness().gain_definition == GAIN_DEFINITION
    assert _witness().loss_definition == LOSS_DEFINITION


def test_formal_source_and_symbol_drift_fail_closed():
    with pytest.raises(ValueError, match="formal_commit"):
        _witness(formal_commit="0" * 40)
    with pytest.raises(ValueError, match="formal_file"):
        _witness(formal_file="NavierStokes/Fake.lean")
    with pytest.raises(ValueError, match="theorem_symbol"):
        _witness(theorem_symbol="NavierStokes.fake")
    with pytest.raises(ValueError, match="dependency_symbols"):
        _witness(dependency_symbols=REQUIRED_FORMAL_SYMBOLS[:-1])
    assert _witness().theorem_symbol == PINNED_STAGE_ESTIMATES_THEOREM


def test_admission_rejects_wrong_object_type():
    with pytest.raises(TypeError, match="Section9ActualStageEstimatesWitness"):
        Section9ActualStageEstimatesAdmission(object())
