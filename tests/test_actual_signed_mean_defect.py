import pytest

from openai_ns_reconstruction.actual_signed_mean_defect import (
    PINNED_REQUESTED_CROSS_DEFECT_THEOREM,
    REQUIRED_DEFECT_FILES,
    REQUIRED_DEFECT_SYMBOLS,
    ActualSignedMeanDefectAdmission,
    ActualSignedMeanDefectWitness,
)
from test_actual_signed_canonical_scope import _scope


def _witness(scope, **changes):
    export = scope.canonical_export.export
    kwargs = dict(
        B=export.B,
        N0=export.N0,
        prepared_N=export.prepared_N,
        band=export.prepared_N,
        component=0,
        cycle_state_repr="lean:CycleState#mean-defect",
        point_repr="lean:Point#mean-defect",
        application_id="requested-cross-defect-application",
        evidence_kind="formal-theorem",
        provenance="synthetic theorem metadata; not a replayed Lean proof",
        active_labels_identity_certified=True,
        finite_head_band_certified=True,
        actual_strip_membership_certified=True,
        requested_stress_identity_certified=True,
        missing_weight_factor_identity_certified=True,
        theorem_application_certified=True,
    )
    kwargs.update(changes)
    return ActualSignedMeanDefectWitness(**kwargs)


def test_finite_head_mean_defect_binds_to_same_canonical_scope_without_truth_upgrade():
    scope = _scope()
    admitted = ActualSignedMeanDefectAdmission(scope, _witness(scope))
    assert admitted.requested_cross_defect_identity_admitted is True
    assert admitted.finite_head_band_admitted is True
    assert admitted.status == "formal-structure"
    assert admitted.requested_cross_defect_theorem_machine_replayed is False
    assert admitted.actual_mean_cross_values_materialized is False
    assert admitted.missing_weight_values_materialized is False
    assert admitted.finite_head_mean_debt_materialized is False
    assert admitted.compact_mean_correction_available is False
    assert admitted.paper_exact_velocity_available is False


def test_finite_head_boundary_is_exact_and_components_are_only_theta_or_axial():
    scope = _scope()
    prepared_n = scope.canonical_export.export.prepared_N
    for band in (0, prepared_n):
        ActualSignedMeanDefectAdmission(scope, _witness(scope, band=band))
    with pytest.raises(ValueError, match=r"n <= Prepared\.N"):
        _witness(scope, band=prepared_n + 1)
    for component in (0, 1):
        ActualSignedMeanDefectAdmission(scope, _witness(scope, component=component))
    with pytest.raises(ValueError, match="component must be 0"):
        _witness(scope, component=2)


def test_cross_wired_canonical_parameters_are_rejected():
    scope = _scope()
    export = scope.canonical_export.export
    with pytest.raises(ValueError, match="B/N0"):
        ActualSignedMeanDefectAdmission(scope, _witness(scope, B=export.B + 1))
    with pytest.raises(ValueError, match="Prepared N"):
        ActualSignedMeanDefectAdmission(
            scope,
            _witness(
                scope,
                prepared_N=export.prepared_N + 1,
                band=export.prepared_N,
            ),
        )


@pytest.mark.parametrize("evidence", ["sampled", "fitted", "numeric-scan"])
def test_non_theorem_mean_defect_evidence_is_rejected(evidence):
    scope = _scope()
    with pytest.raises(ValueError, match="formal-theorem"):
        _witness(scope, evidence_kind=evidence)


@pytest.mark.parametrize(
    "field",
    [
        "active_labels_identity_certified",
        "finite_head_band_certified",
        "actual_strip_membership_certified",
        "requested_stress_identity_certified",
        "missing_weight_factor_identity_certified",
        "theorem_application_certified",
    ],
)
def test_missing_theorem_facts_fail_closed(field):
    scope = _scope()
    with pytest.raises(ValueError, match="theorem-certified true"):
        _witness(scope, **{field: False})


def test_formal_symbol_and_file_drift_fail_closed():
    scope = _scope()
    with pytest.raises(ValueError, match="theorem_symbol"):
        _witness(scope, theorem_symbol="NavierStokes.fake")
    with pytest.raises(ValueError, match="dependency_symbols"):
        _witness(scope, dependency_symbols=REQUIRED_DEFECT_SYMBOLS[:-1])
    with pytest.raises(ValueError, match="dependency_files"):
        _witness(scope, dependency_files=REQUIRED_DEFECT_FILES[:-1])
    witness = _witness(scope)
    assert witness.theorem_symbol == PINNED_REQUESTED_CROSS_DEFECT_THEOREM


def test_formal_revision_drift_is_rejected():
    scope = _scope()
    with pytest.raises(ValueError, match="lean_commit"):
        _witness(scope, lean_commit="0" * 40)
