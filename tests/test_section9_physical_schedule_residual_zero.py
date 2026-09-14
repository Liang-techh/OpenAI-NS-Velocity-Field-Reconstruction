from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_physical_schedule_residual_zero import (
    FORMAL_COMMIT,
    FORMAL_REPOSITORY,
    FORMAL_SOURCE,
    THEOREM,
    Section9PhysicalScheduleResidualZeroWitness,
    admit_section9_physical_schedule_residual_zero,
)


BASE = dict(
    source_id="section9-actual-mixed-source",
    source_revision="rev-actual-001",
    potential_family_id="A-actual",
    direct_family_id="B-actual",
    pressure_family_id="P-actual",
    raw_potential_bounds_id="raw-A-actual",
    raw_direct_bounds_id="raw-B-actual",
    raw_pressure_bounds_id="raw-P-actual",
    background_rate_family_id="background-rates-actual",
    uncut_residual_rate_family_id="residual-rates-actual",
    selected_schedule_id="physical-schedule-actual",
    physical_q_id="physicalQ-h-actual",
    support_region_id="preterminal-support-actual",
    mixed_residual_id="mixed-residual-actual",
    provenance="pinned Lean export bound to one actual Section 9 source",
    h=Fraction(1, 10),
    qbig=Fraction(1, 8),
    lower_stage=7,
    support_open_certified=True,
    support_preterminal_certified=True,
    support_eventually_near_endpoint_origin_certified=True,
    all_stage_fields_smooth_on_physical_sublevel_certified=True,
    raw_potential_stage_bounds_certified=True,
    raw_direct_stage_bounds_certified=True,
    raw_pressure_stage_bounds_certified=True,
    gauge_nonnegative_at_zero_certified=True,
    gauge_positive_after_zero_certified=True,
    gauge_monotone_certified=True,
    gauge_tends_to_infinity_certified=True,
    background_jet_rate_all_J_m_certified=True,
    uncut_residual_jet_rate_all_J_m_certified=True,
    theorem_applied=True,
    selected_schedule_positive_certified=True,
    selected_schedule_doubling_certified=True,
    selected_schedule_strict_mono_certified=True,
    selected_schedule_tends_to_infinity_certified=True,
    selected_schedule_reciprocal_sublevel_certified=True,
    three_cut_bounds_certified=True,
    three_smooth_sums_certified=True,
    vanishing_joint_jets_certified=True,
    same_actual_field_source_certified=True,
)


def _admit(witness: Section9PhysicalScheduleResidualZeroWitness):
    return admit_section9_physical_schedule_residual_zero(
        witness,
        expected_source_id=BASE["source_id"],
        expected_source_revision=BASE["source_revision"],
        expected_potential_family_id=BASE["potential_family_id"],
        expected_direct_family_id=BASE["direct_family_id"],
        expected_pressure_family_id=BASE["pressure_family_id"],
        expected_raw_potential_bounds_id=BASE["raw_potential_bounds_id"],
        expected_raw_direct_bounds_id=BASE["raw_direct_bounds_id"],
        expected_raw_pressure_bounds_id=BASE["raw_pressure_bounds_id"],
        expected_mixed_residual_id=BASE["mixed_residual_id"],
    )


def test_admits_only_the_all_order_physical_schedule_theorem_export():
    cert = _admit(Section9PhysicalScheduleResidualZeroWitness(**BASE))

    assert cert.status == "formal-structure"
    assert cert.formal_physical_schedule_residual_zero_ready
    assert cert.common_physical_schedule_formally_selected
    assert cert.actual_mixed_residual_vanishing_joint_jets_formally_bound

    assert not cert.actual_section9_sequence_verified
    assert not cert.section9_all_order_endpoint_limits_verified
    assert not cert.residual_artifact_ready
    assert not cert.forcing_artifact_ready
    assert not cert.endpoint_residual_closure_verified
    assert not cert.paper_exact_velocity_available
    assert not cert.full_reconstruction


def test_rejects_nonexact_or_out_of_range_physical_h():
    with pytest.raises(ValueError, match="exact Fraction"):
        Section9PhysicalScheduleResidualZeroWitness(**{**BASE, "h": 0.1})

    with pytest.raises(ValueError, match="0 < h < 1/2"):
        Section9PhysicalScheduleResidualZeroWitness(**{**BASE, "h": Fraction(1, 2)})


def test_rejects_finite_stage_or_derivative_frontiers():
    with pytest.raises(ValueError, match="finite prefixes"):
        Section9PhysicalScheduleResidualZeroWitness(**{**BASE, "max_stage_index": 12})

    with pytest.raises(ValueError, match="finite ladders"):
        Section9PhysicalScheduleResidualZeroWitness(**{**BASE, "max_derivative_order": 8})

    with pytest.raises(ValueError, match="finite must be False"):
        Section9PhysicalScheduleResidualZeroWitness(**{**BASE, "finite": True})


def test_rejects_sampled_evidence_or_missing_universal_residual_rate():
    with pytest.raises(ValueError, match="lean-formal-export"):
        Section9PhysicalScheduleResidualZeroWitness(**{**BASE, "evidence_kind": "sampled"})

    with pytest.raises(ValueError, match="uncut_residual_jet_rate_all_J_m_certified"):
        Section9PhysicalScheduleResidualZeroWitness(
            **{**BASE, "uncut_residual_jet_rate_all_J_m_certified": False}
        )


def test_rejects_cross_wired_field_or_raw_bound_identities():
    witness = Section9PhysicalScheduleResidualZeroWitness(**BASE)

    with pytest.raises(ValueError, match="field identity mismatch"):
        admit_section9_physical_schedule_residual_zero(
            witness,
            expected_source_id=BASE["source_id"],
            expected_source_revision=BASE["source_revision"],
            expected_potential_family_id="A-from-another-source",
            expected_direct_family_id=BASE["direct_family_id"],
            expected_pressure_family_id=BASE["pressure_family_id"],
            expected_raw_potential_bounds_id=BASE["raw_potential_bounds_id"],
            expected_raw_direct_bounds_id=BASE["raw_direct_bounds_id"],
            expected_raw_pressure_bounds_id=BASE["raw_pressure_bounds_id"],
            expected_mixed_residual_id=BASE["mixed_residual_id"],
        )

    with pytest.raises(ValueError, match="raw-bound identity mismatch"):
        admit_section9_physical_schedule_residual_zero(
            witness,
            expected_source_id=BASE["source_id"],
            expected_source_revision=BASE["source_revision"],
            expected_potential_family_id=BASE["potential_family_id"],
            expected_direct_family_id=BASE["direct_family_id"],
            expected_pressure_family_id=BASE["pressure_family_id"],
            expected_raw_potential_bounds_id="raw-A-from-another-source",
            expected_raw_direct_bounds_id=BASE["raw_direct_bounds_id"],
            expected_raw_pressure_bounds_id=BASE["raw_pressure_bounds_id"],
            expected_mixed_residual_id=BASE["mixed_residual_id"],
        )


def test_rejects_formal_source_drift_and_unproved_theorem_output():
    assert FORMAL_REPOSITORY == "openai/NavierStokesAndEuler"
    assert FORMAL_COMMIT == "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
    assert FORMAL_SOURCE == "NavierStokes/MixedDiagonalResidual.lean"
    assert THEOREM.endswith("exists_physical_schedule_residual_zero")

    with pytest.raises(ValueError, match="pinned formal source value"):
        Section9PhysicalScheduleResidualZeroWitness(
            **{**BASE, "formal_commit": "not-the-pinned-commit"}
        )

    with pytest.raises(ValueError, match="vanishing_joint_jets_certified"):
        Section9PhysicalScheduleResidualZeroWitness(
            **{**BASE, "vanishing_joint_jets_certified": False}
        )
