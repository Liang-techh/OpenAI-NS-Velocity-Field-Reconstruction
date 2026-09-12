import pytest

from openai_ns_reconstruction.actual_signed_mean_defect import (
    PINNED_DEFECT_THEOREM,
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
        band=0,
        component=0,
        cycle_context_repr="lean:Context#mean-defect",
        cycle_state_repr="lean:State#mean-defect",
        point_repr="lean:Point#mean-defect",
        mean_cross_repr="lean:meanBar(actualCross)#mean-defect",
        requested_stress_repr="lean:requestedStress#mean-defect",
        physical_scale_repr="lean:physicalScale#mean-defect",
        missing_weight_repr="lean:missingWeight#mean-defect",
        application_id="requested-cross-defect-application",
        producer_kind="lean-formal-export",
        provenance="synthetic formal-export metadata; not a replayed Lean proof",
        canonical_active_family_identity_certified=True,
        actual_strip_membership_certified=True,
        actual_cross_identity_certified=True,
        requested_stress_identity_certified=True,
        physical_scale_identity_certified=True,
        prepared_N_identity_certified=True,
        theorem_application_certified=True,
    )
    kwargs.update(changes)
    return ActualSignedMeanDefectWitness(**kwargs)


def test_structural_mean_defect_binds_to_canonical_scope_without_truth_upgrade():
    scope = _scope()
    admitted = ActualSignedMeanDefectAdmission(scope, _witness(scope))
    assert admitted.requested_cross_defect_identity_admitted is True
    assert admitted.status == "formal-structure"
    assert admitted.requested_cross_defect_theorem_machine_replayed is False
    assert admitted.actual_mean_cross_values_materialized is False
    assert admitted.missing_weight_values_materialized is False
    assert admitted.finite_head_mean_defect_values_materialized is False
    assert admitted.finite_head_mean_defect_solved is False
    assert admitted.compact_mean_correction_available is False
    assert admitted.paper_exact_velocity_available is False


def test_defect_identity_is_structural_and_does_not_require_tail_threshold():
    scope = _scope()
    prepared_n = scope.canonical_export.export.prepared_N
    ActualSignedMeanDefectAdmission(scope, _witness(scope, band=0))
    ActualSignedMeanDefectAdmission(scope, _witness(scope, band=prepared_n))
    ActualSignedMeanDefectAdmission(scope, _witness(scope, band=prepared_n + 1))
    for component in (0, 1):
        ActualSignedMeanDefectAdmission(scope, _witness(scope, component=component))
    with pytest.raises(ValueError, match="component must be 0"):
        _witness(scope, component=2)


def test_cross_wired_canonical_parameters_are_rejected():
    scope = _scope()
    with pytest.raises(ValueError, match="B/N0"):
        ActualSignedMeanDefectAdmission(
            scope, _witness(scope, B=scope.canonical_export.export.B + 1)
        )
    with pytest.raises(ValueError, match="Prepared N"):
        ActualSignedMeanDefectAdmission(
            scope,
            _witness(scope, prepared_N=scope.canonical_export.export.prepared_N + 1),
        )


@pytest.mark.parametrize("producer", ["sampled", "fitted", "numeric-scan", "formal-theorem"])
def test_non_export_mean_defect_evidence_is_rejected(producer):
    scope = _scope()
    with pytest.raises(ValueError, match="lean-formal-export"):
        _witness(scope, producer_kind=producer)


@pytest.mark.parametrize(
    "field",
    [
        "canonical_active_family_identity_certified",
        "actual_strip_membership_certified",
        "actual_cross_identity_certified",
        "requested_stress_identity_certified",
        "physical_scale_identity_certified",
        "prepared_N_identity_certified",
        "theorem_application_certified",
    ],
)
def test_missing_formal_export_fact_fails_closed(field):
    scope = _scope()
    with pytest.raises(ValueError, match="formal-export certified true"):
        _witness(scope, **{field: False})


def test_symbol_and_revision_drift_fail_closed():
    scope = _scope()
    with pytest.raises(ValueError, match="theorem_symbol"):
        _witness(scope, theorem_symbol="NavierStokes.fake")
    with pytest.raises(ValueError, match="dependency_symbols"):
        _witness(scope, dependency_symbols=REQUIRED_DEFECT_SYMBOLS[:-1])
    with pytest.raises(ValueError, match="lean_commit"):
        _witness(scope, lean_commit="0" * 40)
    assert _witness(scope).theorem_symbol == PINNED_DEFECT_THEOREM


@pytest.mark.parametrize(
    "field",
    [
        "cycle_context_repr",
        "cycle_state_repr",
        "point_repr",
        "mean_cross_repr",
        "requested_stress_repr",
        "physical_scale_repr",
        "missing_weight_repr",
        "application_id",
    ],
)
def test_opaque_formal_identities_must_be_nonempty(field):
    scope = _scope()
    with pytest.raises(ValueError, match=f"nonempty {field}"):
        _witness(scope, **{field: "  "})
