from dataclasses import replace
from fractions import Fraction

import pytest

from openai_ns_reconstruction.section10_mixed_periodic_candidate_force import (
    PINNED_CANDIDATE_FORCE_THEOREM,
    PINNED_FORMAL_COMMIT,
    PINNED_FORMAL_FILE,
    PINNED_FORMAL_REPOSITORY,
    REQUIRED_FORMAL_SYMBOLS,
    Section10MixedPeriodicCandidateForceAdmission,
    Section10MixedPeriodicCandidateForceWitness,
)


def _witness() -> Section10MixedPeriodicCandidateForceWitness:
    return Section10MixedPeriodicCandidateForceWitness(
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


def test_admits_only_formal_structure_without_truth_promotion() -> None:
    witness = _witness()
    admission = Section10MixedPeriodicCandidateForceAdmission(witness)

    assert witness.formal_repository == PINNED_FORMAL_REPOSITORY
    assert witness.formal_commit == PINNED_FORMAL_COMMIT
    assert witness.formal_file == PINNED_FORMAL_FILE
    assert witness.theorem_symbol == PINNED_CANDIDATE_FORCE_THEOREM
    assert witness.dependency_symbols == REQUIRED_FORMAL_SYMBOLS

    assert admission.endpoint_time == Fraction(1, 1)
    assert admission.mixed_periodic_candidate_force_theorem_admitted
    assert admission.candidate_properties_theorem_output_admitted
    assert admission.smooth_force_theorem_output_admitted
    assert admission.force_endpoint_jet_theorem_output_admitted
    assert admission.status == "formal-structure"

    assert not admission.theorem_machine_replayed
    assert not admission.actual_input_fields_materialized
    assert not admission.residual_limit_values_materialized
    assert not admission.candidate_force_field_materialized
    assert not admission.actual_section9_sequence_verified
    assert not admission.section9_field_smooth_extension_through_t1_constructed
    assert not admission.residual_artifact_ready
    assert not admission.forcing_artifact_ready
    assert not admission.divergence_free_closure_verified
    assert not admission.blow_up_closure_verified
    assert not admission.finite_energy_closure_verified
    assert not admission.endpoint_residual_closure_verified
    assert not admission.paper_exact_velocity_available


def test_rejects_non_formal_producer_and_float_endpoint() -> None:
    with pytest.raises(ValueError, match="lean-formal-export"):
        replace(_witness(), producer_kind="numeric-scan")

    with pytest.raises(TypeError, match="Fraction"):
        replace(_witness(), endpoint_time=1.0)

    with pytest.raises(ValueError, match="paper endpoint"):
        replace(_witness(), endpoint_time=Fraction(3, 4))


def test_rejects_missing_theorem_hypothesis_or_output_certificate() -> None:
    with pytest.raises(ValueError, match="original_residual_vanishing_joint_jets_certified"):
        replace(_witness(), original_residual_vanishing_joint_jets_certified=False)

    with pytest.raises(ValueError, match="axis_speed_tends_to_infinity_certified"):
        replace(_witness(), axis_speed_tends_to_infinity_certified=False)

    with pytest.raises(ValueError, match="force_global_contdiff_output_certified"):
        replace(_witness(), force_global_contdiff_output_certified=False)

    with pytest.raises(ValueError, match="force_endpoint_jets_equal_boundary_limits_certified"):
        replace(_witness(), force_endpoint_jets_equal_boundary_limits_certified=False)


def test_rejects_source_symbol_or_identity_drift() -> None:
    with pytest.raises(ValueError, match="pinned source"):
        replace(_witness(), formal_commit="deadbeef")

    with pytest.raises(ValueError, match="MixedPeriodicAssembly.lean"):
        replace(_witness(), formal_file="NavierStokes/Other.lean")

    with pytest.raises(ValueError, match="exists_candidate_force"):
        replace(_witness(), theorem_symbol="some.other.theorem")

    with pytest.raises(ValueError, match="dependency_symbols"):
        replace(_witness(), dependency_symbols=REQUIRED_FORMAL_SYMBOLS[:-1])

    with pytest.raises(ValueError, match="candidate_force_id"):
        replace(_witness(), candidate_force_id="  ")


def test_admission_rejects_wrong_witness_type() -> None:
    with pytest.raises(TypeError, match="Section10MixedPeriodicCandidateForceWitness"):
        Section10MixedPeriodicCandidateForceAdmission(object())
