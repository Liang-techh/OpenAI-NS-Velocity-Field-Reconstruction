import pytest

from openai_ns_reconstruction.actual_signed_common_curl import (
    PINNED_COMMON_CURL_THEOREM,
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_FILE,
    PINNED_LEAN_REPOSITORY,
    REQUIRED_DEPENDENCY_SYMBOLS,
    ActualSignedCommonCurlAdmission,
    ActualSignedCommonCurlWitness,
)


def _witness(**overrides):
    kwargs = dict(
        B=7,
        N0=11,
        harmonic=3,
        signed_label_repr="synthetic SignedLabel fixture",
        point_repr="synthetic FullPoint fixture",
        request_repr="synthetic request fixture",
        beta_repr="synthetic beta fixture",
        application_id="fixture-actual-signed-common-curl-001",
        evidence_kind="formal-theorem",
        provenance="synthetic theorem-application metadata; no field values",
        uniform_local_jets_hypothesis_certified=True,
        full_strip_domain_membership_certified=True,
        actual_signed_copy_identity_certified=True,
        actual_parameter_identity_certified=True,
        theorem_application_certified=True,
    )
    kwargs.update(overrides)
    return ActualSignedCommonCurlWitness(**kwargs)


def test_actual_signed_common_curl_admits_only_complete_pinned_theorem_contract():
    admission = ActualSignedCommonCurlAdmission(_witness())

    assert admission.application_key == (
        7,
        11,
        "synthetic SignedLabel fixture",
        3,
        "synthetic FullPoint fixture",
        "fixture-actual-signed-common-curl-001",
    )
    assert all(admission.admission_checks().values())
    assert admission.formal_common_curl_theorem_admitted is True
    assert admission.formal_divergence_zero_theorem_admitted is True


def test_actual_signed_common_curl_pins_formal_source_and_dependencies():
    witness = _witness()

    assert witness.lean_repository == PINNED_LEAN_REPOSITORY
    assert witness.lean_commit == PINNED_LEAN_COMMIT
    assert witness.lean_file == PINNED_LEAN_FILE
    assert witness.theorem == PINNED_COMMON_CURL_THEOREM
    assert witness.dependency_symbols == REQUIRED_DEPENDENCY_SYMBOLS

    for field in ("lean_repository", "lean_commit", "lean_file", "theorem"):
        with pytest.raises(ValueError, match="pinned formal source"):
            _witness(**{field: "stale-or-wrong-pin"})

    with pytest.raises(ValueError, match="dependency_symbols"):
        _witness(dependency_symbols=REQUIRED_DEPENDENCY_SYMBOLS[:-1])
    with pytest.raises(ValueError, match="dependency_symbols"):
        _witness(dependency_symbols=tuple(reversed(REQUIRED_DEPENDENCY_SYMBOLS)))


def test_actual_signed_common_curl_rejects_non_theorem_evidence_and_missing_facts():
    for evidence_kind in ("sampled", "fitted", "numeric-scan", ""):
        with pytest.raises(ValueError, match="evidence_kind"):
            _witness(evidence_kind=evidence_kind)

    for field in (
        "uniform_local_jets_hypothesis_certified",
        "full_strip_domain_membership_certified",
        "actual_signed_copy_identity_certified",
        "actual_parameter_identity_certified",
        "theorem_application_certified",
    ):
        with pytest.raises(ValueError, match=field):
            _witness(**{field: False})


def test_actual_signed_common_curl_keeps_lean_indices_opaque_and_fail_closed():
    with pytest.raises(ValueError, match="B must be a nonnegative integer"):
        _witness(B=-1)
    with pytest.raises(ValueError, match="N0 must be a nonnegative integer"):
        _witness(N0=-1)
    with pytest.raises(ValueError, match="harmonic must be a nonnegative integer"):
        _witness(harmonic=-1)
    with pytest.raises(ValueError, match="nonempty signed_label_repr"):
        _witness(signed_label_repr=" ")
    with pytest.raises(ValueError, match="nonempty point_repr"):
        _witness(point_repr="")

    zero_harmonic = ActualSignedCommonCurlAdmission(_witness(harmonic=0))
    assert zero_harmonic.witness.harmonic == 0


def test_actual_signed_common_curl_does_not_promote_unreplayed_theorem_to_velocity():
    admission = ActualSignedCommonCurlAdmission(_witness())

    assert admission.status == "formal-structure"
    assert admission.theorem_application_machine_replayed is False
    assert admission.actual_signed_wave_values_materialized is False
    assert admission.canonical_scope_link_verified is False
    assert admission.amplitude_ode_inputs_verified is False
    assert admission.paper_exact_divergence_free_wave_available is False
    assert admission.paper_exact_velocity_available is False
