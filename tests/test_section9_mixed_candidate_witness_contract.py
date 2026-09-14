from __future__ import annotations

from dataclasses import replace
from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_mixed_candidate_witness_contract import (
    PINNED_FORMAL_COMMIT,
    MixedCandidateWitnessFormalExport,
    admit_mixed_candidate_witness_export,
)


def _valid_export() -> MixedCandidateWitnessFormalExport:
    return MixedCandidateWitnessFormalExport(
        source_id="actual-section9-mixed-candidate",
        source_revision="rev-actual-001",
        field_bundle_id="A-B-P-actual-bundle",
        stage_estimates_id="actual-stage-estimates",
        potential_stage_family_id="A-stages",
        direct_stage_family_id="B-stages",
        pressure_stage_family_id="P-stages",
        provenance="pinned Lean theorem export",
        same_actual_stage_families_certified=True,
        stage_estimates_certified=True,
        shrinking_support_all_stages_certified=True,
        one_sided_extensions_all_stages_certified=True,
        one_common_selected_schedule_certified=True,
        selected_schedule_smooth_sums_certified=True,
        selected_schedule_vanishing_joint_jets_certified=True,
        away_extensions_all_three_sums_certified=True,
        smooth_force_constructed_certified=True,
        candidate_properties_certified=True,
        candidate_consequences_certified=True,
        derivative_h3_blowup_certified=True,
        force_rapid_decay_all_orders_certified=True,
        endpoint_force_jets_all_orders_certified=True,
    )


def _admit(export: MixedCandidateWitnessFormalExport):
    return admit_mixed_candidate_witness_export(
        export,
        expected_source_id="actual-section9-mixed-candidate",
        expected_source_revision="rev-actual-001",
        expected_field_bundle_id="A-B-P-actual-bundle",
        expected_potential_stage_family_id="A-stages",
        expected_direct_stage_family_id="B-stages",
        expected_pressure_stage_family_id="P-stages",
    )


def test_accepts_only_pinned_all_order_formal_witness_chain() -> None:
    certificate = _admit(_valid_export())
    assert certificate.formal_witness_chain_ready
    assert certificate.status == "formal-structure"
    assert certificate.export.formal_commit == PINNED_FORMAL_COMMIT
    assert certificate.export.endpoint == Fraction(1, 1)
    assert certificate.forcing_artifact_ready is False
    assert certificate.paper_exact_velocity_available is False
    assert certificate.full_reconstruction is False


def test_rejects_finite_derivative_frontier() -> None:
    with pytest.raises(ValueError, match="finite derivative-order frontiers"):
        replace(_valid_export(), max_derivative_order=12)


def test_rejects_sampled_or_fitted_evidence() -> None:
    with pytest.raises(ValueError, match="sampled/fitted"):
        replace(_valid_export(), sampled_or_fitted_evidence=True)


def test_rejects_unpinned_formal_commit() -> None:
    with pytest.raises(ValueError, match="formal commit mismatch"):
        replace(_valid_export(), formal_commit="deadbeef")


def test_rejects_missing_all_stage_extension_hypothesis() -> None:
    with pytest.raises(ValueError, match="one_sided_extensions_all_stages_certified"):
        replace(_valid_export(), one_sided_extensions_all_stages_certified=False)


def test_rejects_cross_wired_actual_stage_family() -> None:
    export = _valid_export()
    with pytest.raises(ValueError, match="stage-family identity mismatch"):
        admit_mixed_candidate_witness_export(
            export,
            expected_source_id="actual-section9-mixed-candidate",
            expected_source_revision="rev-actual-001",
            expected_field_bundle_id="A-B-P-actual-bundle",
            expected_potential_stage_family_id="A-stages",
            expected_direct_stage_family_id="other-B-stages",
            expected_pressure_stage_family_id="P-stages",
        )


def test_rejects_manufactured_residual_flag() -> None:
    with pytest.raises(ValueError, match="manufactured residual"):
        replace(_valid_export(), manufactured_residual=True)
