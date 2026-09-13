from dataclasses import replace
from fractions import Fraction

import pytest

from openai_ns_reconstruction.section10_mixed_periodic_candidate_force import (
    Section10MixedPeriodicCandidateForceAdmission,
    Section10MixedPeriodicCandidateForceWitness,
)
from openai_ns_reconstruction.section10_r3_compact_candidate import (
    PINNED_FORMAL_COMMIT,
    PINNED_FORMAL_FILE,
    PINNED_FORMAL_REPOSITORY,
    PINNED_R3_COMPACT_THEOREM,
    REQUIRED_FORMAL_SYMBOLS,
    Section10R3CompactCandidateAdmission,
    Section10R3CompactCandidateWitness,
)


def _periodic_admission() -> Section10MixedPeriodicCandidateForceAdmission:
    witness = Section10MixedPeriodicCandidateForceWitness(
        potential_field_id="actual-A",
        direct_field_id="actual-direct-v",
        pressure_field_id="actual-p",
        periodic_velocity_id="periodicVelocity(actual-A,actual-direct-v)",
        periodic_pressure_id="periodicPressure(actual-p)",
        activated_velocity_id="activatedVelocity(periodicVelocity)",
        activated_pressure_id="activatedPressure(periodicPressure)",
        boundary_limits_family_id="boundaryLimits(actual-A,actual-direct-v,actual-p)",
        candidate_force_id="CandidateFromLimits.force(actual-mixed-periodic-data)",
        application_id="lean-export:MixedPeriodicAssembly.exists_candidate_force:1",
        producer_kind="lean-formal-export",
        provenance="pinned formal theorem application export",
        endpoint_time=Fraction(1, 1),
        potential_open_past_smooth_certified=True,
        direct_open_past_smooth_certified=True,
        pressure_open_past_smooth_certified=True,
        localized_direct_divergence_zero_certified=True,
        original_residual_vanishing_joint_jets_certified=True,
        potential_away_extensions_certified=True,
        direct_away_extensions_certified=True,
        pressure_away_extensions_certified=True,
        axis_speed_tends_to_infinity_certified=True,
        periodic_velocity_identity_certified=True,
        periodic_pressure_identity_certified=True,
        activated_velocity_identity_certified=True,
        activated_pressure_identity_certified=True,
        boundary_limits_identity_certified=True,
        theorem_application_certified=True,
        candidate_properties_output_certified=True,
        force_global_contdiff_output_certified=True,
        force_endpoint_jets_equal_boundary_limits_certified=True,
    )
    return Section10MixedPeriodicCandidateForceAdmission(witness)


def _witness() -> Section10R3CompactCandidateWitness:
    periodic = _periodic_admission()
    return Section10R3CompactCandidateWitness(
        periodic_candidate=periodic,
        potential_field_id=periodic.witness.potential_field_id,
        direct_field_id=periodic.witness.direct_field_id,
        pressure_field_id=periodic.witness.pressure_field_id,
        periodic_candidate_velocity_id=periodic.witness.activated_velocity_id,
        periodic_candidate_pressure_id=periodic.witness.activated_pressure_id,
        periodic_candidate_force_id=periodic.witness.candidate_force_id,
        compact_velocity_id="R3CompactCandidate.velocity(actual-A,actual-direct-v)",
        compact_pressure_id="R3CompactCandidate.pressure(actual-p)",
        compact_force_id="R3CompactCandidate.compactForce(candidate-force)",
        support_cylinder_id="SpatialLocalization.supportCylinder",
        outer_force_support_id="R3CompactCandidate.outerSupport",
        application_id="lean-export:R3CompactCandidate.of_localized_fields:1",
        producer_kind="lean-formal-export",
        provenance="pinned compact whole-space theorem application export",
        candidate_properties_input_identity_certified=True,
        localized_velocity_identity_certified=True,
        localized_pressure_identity_certified=True,
        compact_force_identity_certified=True,
        support_cylinder_identity_certified=True,
        outer_force_support_identity_certified=True,
        velocity_supported_certified=True,
        pressure_supported_certified=True,
        velocity_locally_eq_periodic_certified=True,
        pressure_locally_eq_periodic_certified=True,
        theorem_application_certified=True,
        output_velocity_smooth_certified=True,
        output_pressure_smooth_certified=True,
        output_force_smooth_certified=True,
        output_velocity_compact_support_certified=True,
        output_pressure_compact_support_certified=True,
        output_force_compact_support_certified=True,
        output_zero_initial_velocity_certified=True,
        output_force_time_support_certified=True,
        output_divergence_free_certified=True,
        output_navier_stokes_certified=True,
        output_speed_unbounded_certified=True,
    )


def test_admits_compact_r3_theorem_without_runtime_truth_promotion() -> None:
    witness = _witness()
    admission = Section10R3CompactCandidateAdmission(witness)

    assert witness.formal_repository == PINNED_FORMAL_REPOSITORY
    assert witness.formal_commit == PINNED_FORMAL_COMMIT
    assert witness.formal_file == PINNED_FORMAL_FILE
    assert witness.theorem_symbol == PINNED_R3_COMPACT_THEOREM
    assert witness.dependency_symbols == REQUIRED_FORMAL_SYMBOLS

    assert admission.r3_compact_candidate_theorem_admitted
    assert admission.compact_support_theorem_output_admitted
    assert admission.divergence_free_theorem_output_admitted
    assert admission.navier_stokes_theorem_output_admitted
    assert admission.blow_up_theorem_output_admitted
    assert admission.status == "formal-structure"

    assert not admission.theorem_machine_replayed
    assert not admission.r3_compact_fields_materialized
    assert not admission.compact_force_field_materialized
    assert not admission.residual_artifact_ready
    assert not admission.forcing_artifact_ready
    assert not admission.compact_support_runtime_verified
    assert not admission.divergence_free_closure_verified
    assert not admission.blow_up_closure_verified
    assert not admission.finite_energy_closure_verified
    assert not admission.endpoint_residual_closure_verified
    assert not admission.paper_exact_velocity_available


def test_rejects_cross_wired_periodic_candidate_identities() -> None:
    with pytest.raises(ValueError, match="potential_field_id"):
        replace(_witness(), potential_field_id="other-A")

    with pytest.raises(ValueError, match="periodic_candidate_velocity_id"):
        replace(_witness(), periodic_candidate_velocity_id="wrong-velocity")

    with pytest.raises(ValueError, match="periodic_candidate_force_id"):
        replace(_witness(), periodic_candidate_force_id="residual-defined-force")


def test_rejects_nonformal_evidence_or_missing_theorem_certificates() -> None:
    with pytest.raises(ValueError, match="lean-formal-export"):
        replace(_witness(), producer_kind="numeric-scan")

    with pytest.raises(ValueError, match="velocity_locally_eq_periodic_certified"):
        replace(_witness(), velocity_locally_eq_periodic_certified=False)

    with pytest.raises(ValueError, match="output_divergence_free_certified"):
        replace(_witness(), output_divergence_free_certified=False)

    with pytest.raises(ValueError, match="output_speed_unbounded_certified"):
        replace(_witness(), output_speed_unbounded_certified=False)


def test_rejects_source_or_symbol_drift() -> None:
    with pytest.raises(ValueError, match="pinned source"):
        replace(_witness(), formal_commit="deadbeef")

    with pytest.raises(ValueError, match="R3CompactCandidate.lean"):
        replace(_witness(), formal_file="NavierStokes/Other.lean")

    with pytest.raises(ValueError, match="of_localized_fields"):
        replace(_witness(), theorem_symbol="some.other.theorem")

    with pytest.raises(ValueError, match="dependency_symbols"):
        replace(_witness(), dependency_symbols=REQUIRED_FORMAL_SYMBOLS[:-1])


def test_admission_rejects_wrong_witness_type() -> None:
    with pytest.raises(TypeError, match="Section10R3CompactCandidateWitness"):
        Section10R3CompactCandidateAdmission(object())
