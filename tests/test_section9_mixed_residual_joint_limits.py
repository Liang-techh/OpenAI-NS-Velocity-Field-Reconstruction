from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_mixed_residual_joint_limits import (
    BOUNDARY_LIMITS_LOCALLY_UNIFORM_THEOREM,
    FORMAL_COMMIT,
    RESIDUAL_IDENTITY_THEOREM,
    SECTION10_ENDPOINT_EXACT,
    Section9MixedResidualJointLimitsWitness,
    admit_section9_mixed_residual_joint_limits,
)


SOURCE_ID = "section9-selected-physical-mixed-residual"
SOURCE_REVISION = "lean-export-r1"
POTENTIAL_ID = "section9-potential-sum"
DIRECT_ID = "section9-direct-sum"
PRESSURE_ID = "section9-pressure-sum"
MIXED_RESIDUAL_ID = "MixedDiagonalResidual.residual[a,q,A,B,P]"
ORIGINAL_RESIDUAL_ID = "MixedPeriodicAssembly.originalResidual[ASum,BSum,PSum]"
CUT_RESIDUAL_ID = "MixedPeriodicAssembly.cutResidual[ASum,BSum,PSum]"


def _witness(**overrides):
    kwargs = dict(
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        scale_schedule_id="section9-paper-selected-scale-schedule",
        physical_q_id="PhysicalWaveSum.physicalQ(h)",
        potential_family_id=POTENTIAL_ID,
        direct_family_id=DIRECT_ID,
        pressure_family_id=PRESSURE_ID,
        mixed_residual_id=MIXED_RESIDUAL_ID,
        original_residual_id=ORIGINAL_RESIDUAL_ID,
        cut_residual_id=CUT_RESIDUAL_ID,
        boundary_limits_family_id="section10-cut-residual-boundary-limits-all-orders",
        provenance="pinned Lean theorem export; no finite/sampled surrogate",
        evidence_kind="lean-formal-export",
        endpoint=SECTION10_ENDPOINT_EXACT,
        stage_rate_quantifier="forall-natural-J-m",
        derivative_order_quantifier="forall-natural-derivative-orders",
        max_stage_index=None,
        max_derivative_order=None,
        finite=False,
        physical_h_range_certified=True,
        scale_schedule_tends_to_infinity_certified=True,
        preterminal_open_neighborhood_certified=True,
        all_stage_fields_smooth_certified=True,
        gauge_monotone_unbounded_certified=True,
        cut_stage_bounds_potential_certified=True,
        cut_stage_bounds_direct_certified=True,
        cut_stage_bounds_pressure_certified=True,
        uncut_velocity_jet_rate_all_J_m_certified=True,
        actual_residual_jet_rate_all_J_m_certified=True,
        residual_identity_theorem_applied=True,
        physical_vanishing_joint_jets_theorem_applied=True,
        potential_away_extensions_certified=True,
        direct_away_extensions_certified=True,
        pressure_away_extensions_certified=True,
        cut_residual_vanishing_joint_jets_theorem_applied=True,
        cut_residual_away_extensions_theorem_applied=True,
        boundary_limits_locally_uniform_theorem_applied=True,
        same_actual_field_source_certified=True,
    )
    kwargs.update(overrides)
    return Section9MixedResidualJointLimitsWitness(**kwargs)


def _admit(witness):
    return admit_section9_mixed_residual_joint_limits(
        witness,
        expected_source_id=SOURCE_ID,
        expected_source_revision=SOURCE_REVISION,
        expected_potential_family_id=POTENTIAL_ID,
        expected_direct_family_id=DIRECT_ID,
        expected_pressure_family_id=PRESSURE_ID,
        expected_mixed_residual_id=MIXED_RESIDUAL_ID,
        expected_original_residual_id=ORIGINAL_RESIDUAL_ID,
        expected_cut_residual_id=CUT_RESIDUAL_ID,
    )


def test_exact_pinned_all_order_theorem_chain_admits_without_runtime_promotion():
    certificate = _admit(_witness())

    assert certificate.formal_joint_limits_bridge_ready
    assert certificate.actual_mixed_residual_identity_formally_bound
    assert certificate.original_residual_vanishing_joint_jets_formally_bound
    assert certificate.cut_residual_boundary_limits_locally_uniform_formally_bound
    assert certificate.witness.endpoint == Fraction(1, 1)
    assert not certificate.actual_section9_sequence_verified
    assert not certificate.section9_all_order_endpoint_limits_verified
    assert not certificate.residual_artifact_ready
    assert not certificate.forcing_artifact_ready
    assert not certificate.endpoint_residual_closure_verified
    assert not certificate.paper_exact_velocity_available
    assert not certificate.full_reconstruction


def test_rejects_finite_stage_prefix_or_finite_derivative_ladder():
    with pytest.raises(ValueError, match="max_stage_index must be None"):
        _witness(max_stage_index=200)
    with pytest.raises(ValueError, match="max_derivative_order must be None"):
        _witness(max_derivative_order=40)
    with pytest.raises(ValueError, match="finite must be False"):
        _witness(finite=True)
    with pytest.raises(ValueError, match="universally quantify all natural J,m"):
        _witness(stage_rate_quantifier="J,m<=100")


def test_rejects_sampled_fitted_or_numeric_evidence():
    for kind in ("sampled", "fitted", "certified-numerical", "formal-theorem"):
        with pytest.raises(ValueError, match="must be lean-formal-export"):
            _witness(evidence_kind=kind)


def test_rejects_float_endpoint_and_pinned_formal_source_drift():
    with pytest.raises(ValueError, match="exact rational T=1"):
        _witness(endpoint=1.0)
    with pytest.raises(ValueError, match="formal_commit must equal pinned"):
        _witness(formal_commit="deadbeef")
    with pytest.raises(ValueError, match="residual_identity_theorem must equal pinned"):
        _witness(residual_identity_theorem="Some.Other.residual_identity")
    with pytest.raises(
        ValueError, match="boundary_limits_locally_uniform_theorem must equal pinned"
    ):
        _witness(
            boundary_limits_locally_uniform_theorem=(
                "Some.Other.boundaryLimits_locallyUniform"
            )
        )


def test_rejects_missing_universal_actual_residual_rate_or_away_extensions():
    with pytest.raises(ValueError, match="actual_residual_jet_rate_all_J_m_certified"):
        _witness(actual_residual_jet_rate_all_J_m_certified=False)
    with pytest.raises(ValueError, match="potential_away_extensions_certified"):
        _witness(potential_away_extensions_certified=False)
    with pytest.raises(ValueError, match="pressure_away_extensions_certified"):
        _witness(pressure_away_extensions_certified=False)


def test_rejects_missing_theorem_chain_application():
    with pytest.raises(ValueError, match="residual_identity_theorem_applied"):
        _witness(residual_identity_theorem_applied=False)
    with pytest.raises(ValueError, match="physical_vanishing_joint_jets_theorem_applied"):
        _witness(physical_vanishing_joint_jets_theorem_applied=False)
    with pytest.raises(ValueError, match="boundary_limits_locally_uniform_theorem_applied"):
        _witness(boundary_limits_locally_uniform_theorem_applied=False)


def test_rejects_cross_wired_source_fields_or_residual_identities():
    witness = _witness()

    with pytest.raises(ValueError, match="source identity mismatch"):
        admit_section9_mixed_residual_joint_limits(
            witness,
            expected_source_id=SOURCE_ID,
            expected_source_revision="other-revision",
            expected_potential_family_id=POTENTIAL_ID,
            expected_direct_family_id=DIRECT_ID,
            expected_pressure_family_id=PRESSURE_ID,
            expected_mixed_residual_id=MIXED_RESIDUAL_ID,
            expected_original_residual_id=ORIGINAL_RESIDUAL_ID,
            expected_cut_residual_id=CUT_RESIDUAL_ID,
        )

    with pytest.raises(ValueError, match="field identity mismatch"):
        admit_section9_mixed_residual_joint_limits(
            witness,
            expected_source_id=SOURCE_ID,
            expected_source_revision=SOURCE_REVISION,
            expected_potential_family_id="other-potential",
            expected_direct_family_id=DIRECT_ID,
            expected_pressure_family_id=PRESSURE_ID,
            expected_mixed_residual_id=MIXED_RESIDUAL_ID,
            expected_original_residual_id=ORIGINAL_RESIDUAL_ID,
            expected_cut_residual_id=CUT_RESIDUAL_ID,
        )

    with pytest.raises(ValueError, match="mixed residual identity mismatch"):
        admit_section9_mixed_residual_joint_limits(
            witness,
            expected_source_id=SOURCE_ID,
            expected_source_revision=SOURCE_REVISION,
            expected_potential_family_id=POTENTIAL_ID,
            expected_direct_family_id=DIRECT_ID,
            expected_pressure_family_id=PRESSURE_ID,
            expected_mixed_residual_id="manufactured-residual",
            expected_original_residual_id=ORIGINAL_RESIDUAL_ID,
            expected_cut_residual_id=CUT_RESIDUAL_ID,
        )


def test_constants_pin_expected_official_theorems():
    assert FORMAL_COMMIT == "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
    assert RESIDUAL_IDENTITY_THEOREM.endswith("residual_eq_originalResidual")
    assert BOUNDARY_LIMITS_LOCALLY_UNIFORM_THEOREM.endswith(
        "boundaryLimits_locallyUniform"
    )
