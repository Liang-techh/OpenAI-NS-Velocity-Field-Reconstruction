from dataclasses import replace

import pytest

from openai_ns_reconstruction.actual_signed_common_curl import (
    ActualSignedCommonCurlAdmission,
    ActualSignedCommonCurlWitness,
)
from openai_ns_reconstruction.actual_signed_common_curl_invariant import (
    ActualSignedInvariantCommonCurlAdmission,
    ActualSignedInvariantJetsWitness,
    PINNED_JETS_THEOREM,
)


def _curl_witness(**changes):
    base = ActualSignedCommonCurlWitness(
        B=8,
        N0=16,
        harmonic=11,
        signed_label_repr="lean:SignedLabel#fixture",
        point_repr="lean:FullPoint#fixture",
        request_repr="lean:LocalSignedRequest.fullRequest#fixture",
        beta_repr="lean:sigma#fixture",
        application_id="common-curl-fixture",
        evidence_kind="formal-theorem",
        provenance="synthetic theorem metadata; not a Lean replay",
        uniform_local_jets_hypothesis_certified=True,
        full_strip_domain_membership_certified=True,
        actual_signed_copy_identity_certified=True,
        actual_parameter_identity_certified=True,
        theorem_application_certified=True,
    )
    return replace(base, **changes)


def _jets_witness(**changes):
    base = ActualSignedInvariantJetsWitness(
        B=8,
        N0=16,
        request_repr="lean:LocalSignedRequest.fullRequest#fixture",
        beta_repr="lean:sigma#fixture",
        cycle_state_repr="lean:CycleState#fixture",
        geometry_repr="lean:SignedMeanGain.Geometry#fixture",
        context_repr="lean:CorrectionState.Context#fixture",
        primary_repr="lean:primary#fixture",
        pressure_repr="lean:P#fixture",
        label_carrier_repr="lean:labelCarrier#fixture",
        application_id="full-request-jets-fixture",
        evidence_kind="formal-theorem",
        provenance="synthetic theorem metadata; not a Lean replay",
        cycle_analytic_invariant_certified=True,
        full_request_identity_certified=True,
        geometry_strip_is_actual_strip_certified=True,
        phase_cell_instantiation_certified=True,
        beta_is_invariant_sigma_certified=True,
        theorem_application_certified=True,
    )
    return replace(base, **changes)


def _admit(jets=None, curl=None):
    return ActualSignedInvariantCommonCurlAdmission(
        jets=jets or _jets_witness(),
        common_curl=ActualSignedCommonCurlAdmission(curl or _curl_witness()),
    )


def test_invariant_jets_close_the_common_curl_hypothesis():
    admission = _admit()
    assert admission.invariant_supplies_common_curl_jets is True
    assert admission.formal_divergence_zero_theorem_admitted is True
    assert admission.status == "formal-structure"
    assert admission.full_request_jets_theorem_machine_replayed is False
    assert admission.actual_cycle_invariant_machine_verified is False
    assert admission.common_curl_theorem_machine_replayed is False
    assert admission.actual_signed_wave_values_materialized is False
    assert admission.paper_exact_divergence_free_wave_available is False
    assert admission.paper_exact_velocity_available is False


@pytest.mark.parametrize(
    ("jets_changes", "curl_changes"),
    [
        ({"B": 7}, {}),
        ({"N0": 15}, {}),
        ({"request_repr": "lean:different-request"}, {}),
        ({"beta_repr": "lean:different-sigma"}, {}),
    ],
)
def test_cross_wired_common_curl_application_is_rejected(jets_changes, curl_changes):
    with pytest.raises(ValueError, match="uncertified invariant-to-common-curl composition"):
        _admit(
            jets=_jets_witness(**jets_changes),
            curl=_curl_witness(**curl_changes),
        )


@pytest.mark.parametrize(
    "field",
    [
        "cycle_analytic_invariant_certified",
        "full_request_identity_certified",
        "geometry_strip_is_actual_strip_certified",
        "phase_cell_instantiation_certified",
        "beta_is_invariant_sigma_certified",
        "theorem_application_certified",
    ],
)
def test_missing_theorem_fact_is_rejected(field):
    with pytest.raises(ValueError, match="theorem-certified true"):
        _jets_witness(**{field: False})


@pytest.mark.parametrize("kind", ["sampled", "fitted", "numeric-scan"])
def test_non_theorem_jets_evidence_is_rejected(kind):
    with pytest.raises(ValueError, match="sampled/fitted/numeric-scan"):
        _jets_witness(evidence_kind=kind)


def test_stale_jets_theorem_pin_is_rejected():
    with pytest.raises(ValueError, match="pinned formal source"):
        _jets_witness(theorem="ActualSignedStageControls.fullRequest_jets_from_residuals")


def test_dependency_drift_is_rejected():
    with pytest.raises(ValueError, match="dependency_symbols"):
        _jets_witness(dependency_symbols=("LocalSignedRequest.fullRequest",))


def test_pinned_theorem_name_is_exact():
    assert PINNED_JETS_THEOREM == "ActualSignedStageControls.fullRequest_jets_from_invariant"
