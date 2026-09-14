import math

import pytest

from openai_ns_reconstruction.section9_actual_stage_estimates import (
    Section9ActualStageEstimatesAdmission,
    Section9ActualStageEstimatesWitness,
)
from openai_ns_reconstruction.section9_finite_prefix_fields import (
    MATERIALIZATION_KIND,
    Section9FinitePrefixFields,
    Section9MaterializedStage,
)


def _admission() -> Section9ActualStageEstimatesAdmission:
    witness = Section9ActualStageEstimatesWitness(
        domain_id="domain",
        forcing_id="forcing",
        q_id="q",
        run_data_id="run-data",
        physical_data_id="physical-data",
        positive_bump_id="positive-bump",
        finite_velocity_family_id="finite-velocity-family",
        finite_pressure_family_id="finite-pressure-family",
        finite_forcing_family_id="finite-forcing-family",
        mixed_physical_data_id="mixed-physical-data",
        stage_estimates_id="stage-estimates",
        application_id="stage-estimates-application",
        producer_kind="lean-formal-export",
        provenance="test-only exact identity fixture",
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
    return Section9ActualStageEstimatesAdmission(witness)


def _stage(index: int, scale: float = 1.0) -> Section9MaterializedStage[float]:
    return Section9MaterializedStage(
        index=index,
        potential_family_id="finite-velocity-family",
        direct_family_id="finite-velocity-family",
        pressure_family_id="finite-pressure-family",
        stage_id=f"stage-{index}",
        curl_potential=lambda x, i=index, s=scale: (s * (i + 1) * x, s * 2.0, -s * i),
        direct_velocity=lambda x, i=index, s=scale: (s, -s * x, s * (i + 0.5)),
        pressure=lambda x, i=index, s=scale: s * (i + 1) * (x + 2.0),
        materialization_kind=MATERIALIZATION_KIND,
        provenance="test-only callable field fixture",
        field_materialization_certified=True,
    )


def test_literal_prefix_matches_independent_component_sum() -> None:
    stages = (_stage(0), _stage(1), _stage(2))
    prefix = Section9FinitePrefixFields.from_stages(_admission(), stages)
    x = 1.25

    expected_velocity = [0.0, 0.0, 0.0]
    expected_pressure = 0.0
    for stage in stages:
        curl = stage.curl_potential(x)
        direct = stage.direct_velocity(x)
        for i in range(3):
            expected_velocity[i] += curl[i] + direct[i]
        expected_pressure += stage.pressure(x)

    assert prefix.J == 2
    assert prefix.stage_count == 3
    assert prefix.velocity(x) == pytest.approx(tuple(expected_velocity))
    assert prefix.pressure(x) == pytest.approx(expected_pressure)
    assert prefix.exact_prefix_shape_materialized
    assert prefix.status == "formal-structure"
    assert not prefix.actual_paper_stage_fields_available
    assert not prefix.section9_infinite_iteration_closed
    assert not prefix.paper_exact_velocity_available


def test_prefix_never_fills_a_missing_stage_with_zero() -> None:
    admission = _admission()
    with pytest.raises(ValueError, match="contiguous"):
        Section9FinitePrefixFields(admission, (_stage(0), _stage(2)), J=1)

    with pytest.raises(ValueError, match="exactly the stages"):
        Section9FinitePrefixFields(admission, (_stage(0), _stage(1)), J=2)


def test_prefix_rejects_cross_wired_stage_families() -> None:
    bad = Section9MaterializedStage(
        index=0,
        potential_family_id="wrong-velocity-family",
        direct_family_id="finite-velocity-family",
        pressure_family_id="finite-pressure-family",
        stage_id="stage-0",
        curl_potential=lambda _: (0.0, 0.0, 0.0),
        direct_velocity=lambda _: (0.0, 0.0, 0.0),
        pressure=lambda _: 0.0,
        materialization_kind=MATERIALIZATION_KIND,
        provenance="test-only cross-wire fixture",
        field_materialization_certified=True,
    )
    with pytest.raises(ValueError, match="finite-velocity family"):
        Section9FinitePrefixFields.from_stages(_admission(), (bad,))


def test_stage_rejects_sampled_or_uncertified_payloads() -> None:
    kwargs = dict(
        index=0,
        potential_family_id="finite-velocity-family",
        direct_family_id="finite-velocity-family",
        pressure_family_id="finite-pressure-family",
        stage_id="stage-0",
        curl_potential=lambda _: (0.0, 0.0, 0.0),
        direct_velocity=lambda _: (0.0, 0.0, 0.0),
        pressure=lambda _: 0.0,
        provenance="test-only invalid fixture",
    )
    with pytest.raises(ValueError, match="sampled/fitted/toy"):
        Section9MaterializedStage(
            **kwargs,
            materialization_kind="sampled-field",
            field_materialization_certified=True,
        )
    with pytest.raises(ValueError, match="exactly true"):
        Section9MaterializedStage(
            **kwargs,
            materialization_kind=MATERIALIZATION_KIND,
            field_materialization_certified=False,
        )


def test_nonfinite_field_values_fail_closed_at_evaluation() -> None:
    stage = Section9MaterializedStage(
        index=0,
        potential_family_id="finite-velocity-family",
        direct_family_id="finite-velocity-family",
        pressure_family_id="finite-pressure-family",
        stage_id="stage-0",
        curl_potential=lambda _: (math.inf, 0.0, 0.0),
        direct_velocity=lambda _: (0.0, 0.0, 0.0),
        pressure=lambda _: 0.0,
        materialization_kind=MATERIALIZATION_KIND,
        provenance="test-only nonfinite fixture",
        field_materialization_certified=True,
    )
    prefix = Section9FinitePrefixFields.from_stages(_admission(), (stage,))
    with pytest.raises(ValueError, match="finite real"):
        prefix.velocity(0.0)


def test_duplicate_stage_identity_is_rejected() -> None:
    stage0 = _stage(0)
    stage1 = Section9MaterializedStage(
        index=1,
        potential_family_id="finite-velocity-family",
        direct_family_id="finite-velocity-family",
        pressure_family_id="finite-pressure-family",
        stage_id=stage0.stage_id,
        curl_potential=lambda _: (0.0, 0.0, 0.0),
        direct_velocity=lambda _: (0.0, 0.0, 0.0),
        pressure=lambda _: 0.0,
        materialization_kind=MATERIALIZATION_KIND,
        provenance="test-only duplicate-id fixture",
        field_materialization_certified=True,
    )
    with pytest.raises(ValueError, match="distinct stage_id"):
        Section9FinitePrefixFields.from_stages(_admission(), (stage0, stage1))
