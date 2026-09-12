import pytest

from openai_ns_reconstruction.section9_residual_flatness import (
    PINNED_ALL_JETS_FLAT_THEOREM,
    REQUIRED_FORMAL_SYMBOLS,
    Section9ResidualFlatnessAdmission,
    Section9ResidualFlatnessWitness,
)


def _witness(**changes):
    kwargs = dict(
        domain_id="lean:U#section9",
        filter_id="lean:l#section9",
        scale_q_id="lean:q#section9",
        limit_velocity_id="lean:u#section9",
        limit_pressure_id="lean:p#section9",
        velocity_stage_family_id="lean:uStage#section9",
        pressure_stage_family_id="lean:pStage#section9",
        gain_schedule_id="lean:g#section9",
        background_loss_id="lean:Lbg#section9",
        tail_loss_id="lean:Ltail#section9",
        residual_loss_id="lean:Lres#section9",
        application_id="lean:allJetsFlat_residual_of_stages#application",
        producer_kind="lean-formal-export",
        provenance="synthetic formal-export metadata; not a replayed Lean theorem",
        domain_open_certified=True,
        filter_eventually_in_domain_certified=True,
        scale_eventually_unit_interval_certified=True,
        limit_velocity_smooth_certified=True,
        limit_pressure_smooth_certified=True,
        velocity_stages_smooth_certified=True,
        pressure_stages_smooth_certified=True,
        gain_tendsto_at_top_certified=True,
        background_stage_jet_rates_certified=True,
        velocity_tail_jet_rates_certified=True,
        pressure_tail_jet_rates_certified=True,
        stage_residual_jet_rates_certified=True,
        theorem_application_certified=True,
    )
    kwargs.update(changes)
    return Section9ResidualFlatnessWitness(**kwargs)


def test_formal_residual_flatness_contract_has_no_truth_upgrade():
    admitted = Section9ResidualFlatnessAdmission(_witness())
    assert admitted.all_jets_flat_residual_conclusion_admitted is True
    assert admitted.stage_convergence_contract_admitted is True
    assert admitted.status == "formal-structure"
    assert admitted.all_jets_flat_theorem_machine_replayed is False
    assert admitted.section9_iteration_machine_materialized is False
    assert admitted.stage_fields_materialized is False
    assert admitted.stage_rate_bounds_materialized is False
    assert admitted.limit_velocity_materialized is False
    assert admitted.navier_stokes_residual_values_materialized is False
    assert admitted.paper_exact_velocity_available is False


@pytest.mark.parametrize("producer", ["sampled", "fitted", "numeric-scan", "formal-theorem"])
def test_non_formal_export_producers_are_rejected(producer):
    with pytest.raises(ValueError, match="lean-formal-export"):
        _witness(producer_kind=producer)


@pytest.mark.parametrize(
    "field",
    [
        "domain_open_certified",
        "filter_eventually_in_domain_certified",
        "scale_eventually_unit_interval_certified",
        "limit_velocity_smooth_certified",
        "limit_pressure_smooth_certified",
        "velocity_stages_smooth_certified",
        "pressure_stages_smooth_certified",
        "gain_tendsto_at_top_certified",
        "background_stage_jet_rates_certified",
        "velocity_tail_jet_rates_certified",
        "pressure_tail_jet_rates_certified",
        "stage_residual_jet_rates_certified",
        "theorem_application_certified",
    ],
)
def test_missing_theorem_hypothesis_fails_closed(field):
    with pytest.raises(ValueError, match="formal-export certified true"):
        _witness(**{field: False})


@pytest.mark.parametrize(
    "field",
    [
        "domain_id",
        "filter_id",
        "scale_q_id",
        "limit_velocity_id",
        "limit_pressure_id",
        "velocity_stage_family_id",
        "pressure_stage_family_id",
        "gain_schedule_id",
        "background_loss_id",
        "tail_loss_id",
        "residual_loss_id",
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
    assert _witness().theorem_symbol == PINNED_ALL_JETS_FLAT_THEOREM


def test_admission_rejects_wrong_object_type():
    with pytest.raises(TypeError, match="Section9ResidualFlatnessWitness"):
        Section9ResidualFlatnessAdmission(object())
