import pytest

from openai_ns_reconstruction.actual_signed_mean_tail import (
    PINNED_LITERAL_TAIL_THEOREM,
    REQUIRED_MEAN_SYMBOLS,
    ActualSignedMeanTailAdmission,
    ActualSignedMeanTailWitness,
)
from test_actual_signed_canonical_scope import _scope


def _witness(scope, **changes):
    export = scope.canonical_export.export
    kwargs = dict(
        B=export.B,
        N0=export.N0,
        prepared_N=export.prepared_N,
        band=export.prepared_N + 1,
        component=0,
        cycle_state_repr="lean:CycleState#mean-tail",
        point_repr="lean:Point#mean-tail",
        application_id="literal-requested-cross-tail-application",
        evidence_kind="formal-theorem",
        provenance="synthetic theorem metadata; not a replayed Lean proof",
        active_labels_identity_certified=True,
        tail_band_hypothesis_certified=True,
        actual_strip_membership_certified=True,
        literal_cycle_parameters_identity_certified=True,
        theorem_application_certified=True,
    )
    kwargs.update(changes)
    return ActualSignedMeanTailWitness(**kwargs)


def test_literal_tail_mean_cross_binds_to_same_canonical_scope_without_truth_upgrade():
    scope = _scope()
    admitted = ActualSignedMeanTailAdmission(scope, _witness(scope))
    assert admitted.requested_cross_tail_identity_admitted is True
    assert admitted.status == "formal-structure"
    assert admitted.literal_requested_cross_tail_theorem_machine_replayed is False
    assert admitted.actual_mean_cross_values_materialized is False
    assert admitted.finite_head_mean_defect_solved is False
    assert admitted.compact_mean_correction_available is False
    assert admitted.paper_exact_velocity_available is False


def test_tail_threshold_is_exact_and_components_are_only_theta_or_axial():
    scope = _scope()
    n0 = scope.canonical_export.export.prepared_N
    with pytest.raises(ValueError, match="Prepared.N \+ 1 <= n"):
        _witness(scope, band=n0)
    for component in (0, 1):
        ActualSignedMeanTailAdmission(scope, _witness(scope, component=component))
    with pytest.raises(ValueError, match="component must be 0"):
        _witness(scope, component=2)


def test_cross_wired_canonical_parameters_are_rejected():
    scope = _scope()
    with pytest.raises(ValueError, match="B/N0"):
        ActualSignedMeanTailAdmission(scope, _witness(scope, B=scope.canonical_export.export.B + 1))
    with pytest.raises(ValueError, match="Prepared N"):
        ActualSignedMeanTailAdmission(
            scope,
            _witness(scope, prepared_N=scope.canonical_export.export.prepared_N + 1,
                     band=scope.canonical_export.export.prepared_N + 2),
        )


@pytest.mark.parametrize("evidence", ["sampled", "fitted", "numeric-scan"])
def test_non_theorem_mean_tail_evidence_is_rejected(evidence):
    scope = _scope()
    with pytest.raises(ValueError, match="formal-theorem"):
        _witness(scope, evidence_kind=evidence)


def test_missing_theorem_fact_and_symbol_drift_fail_closed():
    scope = _scope()
    with pytest.raises(ValueError, match="theorem-certified true"):
        _witness(scope, active_labels_identity_certified=False)
    with pytest.raises(ValueError, match="theorem_symbol"):
        _witness(scope, theorem_symbol="NavierStokes.fake")
    with pytest.raises(ValueError, match="dependency_symbols"):
        _witness(scope, dependency_symbols=REQUIRED_MEAN_SYMBOLS[:-1])
    witness = _witness(scope)
    assert witness.theorem_symbol == PINNED_LITERAL_TAIL_THEOREM


def test_formal_revision_drift_is_rejected():
    scope = _scope()
    with pytest.raises(ValueError, match="lean_commit"):
        _witness(scope, lean_commit="0" * 40)
