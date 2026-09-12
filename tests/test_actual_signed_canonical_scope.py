from dataclasses import replace

import pytest

from openai_ns_reconstruction.actual_signed_canonical_scope import (
    LINK_PRODUCER_KIND,
    REQUIRED_LINK_SYMBOLS,
    CanonicalActualSignedLabelLinkWitness,
    CanonicalSignedInvariantCommonCurlBinding,
    LargeBandCanonicalSignedCommonCurlScopeAdmission,
)
from openai_ns_reconstruction.actual_signed_common_curl import (
    ActualSignedCommonCurlAdmission,
    ActualSignedCommonCurlWitness,
)
from openai_ns_reconstruction.actual_signed_common_curl_invariant import (
    ActualSignedInvariantCommonCurlAdmission,
    ActualSignedInvariantJetsWitness,
)
from openai_ns_reconstruction.phase_large_band_canonical_application_export import (
    LargeBandCanonicalApplicationExportAdmission,
)
from test_phase_large_band_canonical_application_export import _coverage, _export


def _canonical_export():
    coverage = _coverage()
    return LargeBandCanonicalApplicationExportAdmission(coverage, _export(coverage))


def _actual_application(*, B, N0, signed_label_repr, suffix):
    request = f"lean:LocalSignedRequest.fullRequest#{suffix}"
    beta = f"lean:sigma#{suffix}"
    curl = ActualSignedCommonCurlAdmission(
        ActualSignedCommonCurlWitness(
            B=B,
            N0=N0,
            harmonic=11,
            signed_label_repr=signed_label_repr,
            point_repr=f"lean:FullPoint#{suffix}",
            request_repr=request,
            beta_repr=beta,
            application_id=f"common-curl-{suffix}",
            evidence_kind="formal-theorem",
            provenance="synthetic theorem metadata; not a Lean replay",
            uniform_local_jets_hypothesis_certified=True,
            full_strip_domain_membership_certified=True,
            actual_signed_copy_identity_certified=True,
            actual_parameter_identity_certified=True,
            theorem_application_certified=True,
        )
    )
    jets = ActualSignedInvariantJetsWitness(
        B=B,
        N0=N0,
        request_repr=request,
        beta_repr=beta,
        cycle_state_repr=f"lean:CycleState#{suffix}",
        geometry_repr=f"lean:SignedMeanGain.Geometry#{suffix}",
        context_repr=f"lean:CorrectionState.Context#{suffix}",
        primary_repr=f"lean:primary#{suffix}",
        pressure_repr=f"lean:P#{suffix}",
        label_carrier_repr=f"lean:labelCarrier#{suffix}",
        application_id=f"full-request-jets-{suffix}",
        evidence_kind="formal-theorem",
        provenance="synthetic theorem metadata; not a Lean replay",
        cycle_analytic_invariant_certified=True,
        full_request_identity_certified=True,
        geometry_strip_is_actual_strip_certified=True,
        phase_cell_instantiation_certified=True,
        beta_is_invariant_sigma_certified=True,
        theorem_application_certified=True,
    )
    return ActualSignedInvariantCommonCurlAdmission(jets=jets, common_curl=curl)


def _binding(canonical_export, label, *, suffix=None, **link_changes):
    ell, a, sigma = label
    suffix = suffix or f"{ell}-{a[0]}-{a[1]}-{a[2]}-{sigma}"
    signed_label_repr = link_changes.pop(
        "signed_label_repr", f"lean:SignedLabel#{suffix}"
    )
    link = CanonicalActualSignedLabelLinkWitness(
        ell=ell,
        a=a,
        sigma=sigma,
        signed_label_repr=signed_label_repr,
        scope_sha256=canonical_export.export.scope_sha256,
        prepared_N=canonical_export.export.prepared_N,
        application_id=f"canonical-signed-link-{suffix}",
        producer_kind=LINK_PRODUCER_KIND,
        provenance="synthetic formal-export metadata; not a replayed Lean artifact",
        canonical_family_label_identity_certified=True,
        signed_label_export_identity_certified=True,
        canonical_phase_cell_identity_certified=True,
        same_prepared_application_certified=True,
        formal_export_application_certified=True,
        **link_changes,
    )
    actual = _actual_application(
        B=canonical_export.export.B,
        N0=canonical_export.export.N0,
        signed_label_repr=signed_label_repr,
        suffix=suffix,
    )
    return CanonicalSignedInvariantCommonCurlBinding(link=link, actual_application=actual)


def _scope(canonical_export=None):
    canonical_export = canonical_export or _canonical_export()
    bindings = tuple(
        _binding(canonical_export, label)
        for label in canonical_export.export.signed_labels
    )
    return LargeBandCanonicalSignedCommonCurlScopeAdmission(canonical_export, bindings)


def test_exact_canonical_signed_scope_binds_to_actual_common_curl_without_truth_upgrade():
    admitted = _scope()
    assert admitted.canonical_to_actual_signed_scope_machine_checked is True
    assert admitted.covered_canonical_labels == frozenset(
        admitted.canonical_export.export.signed_labels
    )
    assert admitted.status == "formal-structure"
    assert admitted.actual_signed_label_export_machine_verified is False
    assert admitted.canonical_prepared_application_machine_replayed is False
    assert admitted.actual_cycle_invariant_machine_verified is False
    assert admitted.common_curl_theorem_machine_replayed is False
    assert admitted.actual_signed_wave_values_materialized is False
    assert admitted.paper_exact_divergence_free_wave_available is False
    assert admitted.paper_exact_velocity_available is False


def test_missing_or_duplicate_canonical_label_is_rejected():
    canonical = _canonical_export()
    bindings = tuple(_binding(canonical, label) for label in canonical.export.signed_labels)
    with pytest.raises(ValueError, match="cover the canonical signed label scope exactly"):
        LargeBandCanonicalSignedCommonCurlScopeAdmission(canonical, bindings[:-1])
    with pytest.raises(ValueError, match="must not be duplicated"):
        LargeBandCanonicalSignedCommonCurlScopeAdmission(
            canonical, bindings[:-1] + (bindings[0],)
        )


def test_cross_wired_scope_digest_and_prepared_N_are_rejected():
    canonical = _canonical_export()
    label = canonical.export.signed_labels[0]
    wrong_scope = _binding(canonical, label, scope_sha256="1" * 64)
    rest = tuple(_binding(canonical, x) for x in canonical.export.signed_labels[1:])
    with pytest.raises(ValueError, match="scope digest"):
        LargeBandCanonicalSignedCommonCurlScopeAdmission(canonical, (wrong_scope,) + rest)

    wrong_n = _binding(canonical, label, prepared_N=canonical.export.prepared_N + 1)
    with pytest.raises(ValueError, match="Prepared N"):
        LargeBandCanonicalSignedCommonCurlScopeAdmission(canonical, (wrong_n,) + rest)


def test_cross_wired_B_N0_actual_application_is_rejected():
    canonical = _canonical_export()
    label = canonical.export.signed_labels[0]
    good = _binding(canonical, label)
    wrong_actual = _actual_application(
        B=canonical.export.B + 1,
        N0=canonical.export.N0,
        signed_label_repr=good.link.signed_label_repr,
        suffix="wrong-B",
    )
    wrong = CanonicalSignedInvariantCommonCurlBinding(good.link, wrong_actual)
    rest = tuple(_binding(canonical, x) for x in canonical.export.signed_labels[1:])
    with pytest.raises(ValueError, match="B/N0"):
        LargeBandCanonicalSignedCommonCurlScopeAdmission(canonical, (wrong,) + rest)


def test_opaque_signed_label_must_match_actual_application_and_be_unique():
    canonical = _canonical_export()
    label0, label1, *tail = canonical.export.signed_labels
    first = _binding(canonical, label0, signed_label_repr="lean:SignedLabel#shared")
    second = _binding(canonical, label1, signed_label_repr="lean:SignedLabel#shared")
    rest = tuple(_binding(canonical, x) for x in tail)
    with pytest.raises(ValueError, match="cannot represent multiple canonical labels"):
        LargeBandCanonicalSignedCommonCurlScopeAdmission(canonical, (first, second) + rest)

    link = first.link
    mismatched_actual = _actual_application(
        B=canonical.export.B,
        N0=canonical.export.N0,
        signed_label_repr="lean:SignedLabel#different",
        suffix="different",
    )
    with pytest.raises(ValueError, match="opaque SignedLabel identity"):
        CanonicalSignedInvariantCommonCurlBinding(link, mismatched_actual)


@pytest.mark.parametrize("producer", ["sampled", "fitted", "numeric-scan"])
def test_non_formal_export_label_mapping_is_rejected(producer):
    canonical = _canonical_export()
    label = canonical.export.signed_labels[0]
    with pytest.raises(ValueError, match="lean-formal-export"):
        _binding(canonical, label, producer_kind=producer)


def test_missing_export_fact_or_symbol_drift_is_rejected():
    canonical = _canonical_export()
    label = canonical.export.signed_labels[0]
    with pytest.raises(ValueError, match="formal-export certified true"):
        _binding(canonical, label, same_prepared_application_certified=False)
    with pytest.raises(ValueError, match="dependency_symbols"):
        _binding(canonical, label, dependency_symbols=REQUIRED_LINK_SYMBOLS[:-1])
