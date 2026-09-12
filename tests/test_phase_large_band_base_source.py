from fractions import Fraction

import pytest

from openai_ns_reconstruction.phase_large_band_base_source import (
    LargeBandBaseSourceBinding,
    LargeBandBaseSourceWitness,
)
from openai_ns_reconstruction.phase_large_band_local_base import (
    AsymptoticSlowLabel,
    LargeBandLocalBaseAdmission,
    LargeBandLocalBaseWitness,
)
from openai_ns_reconstruction.phase_large_band_scale import LargeBandPhaseScaleCertificate


def _admission(*, ell=100000, a=(7, -3, 11), sigma=1, M=2.0):
    scale = LargeBandPhaseScaleCertificate(ell=ell, h=Fraction(1, 200), M=M)
    witness = LargeBandLocalBaseWitness(
        label=AsymptoticSlowLabel(ell=ell, a=a, sigma=sigma),
        M=M,
        diameter_times_S3=0.75,
        second_F0_over_M=0.8,
        second_G0_over_M=0.7,
        first_F0_over_M=0.6,
        value_F_error_over_M_epsilon2=0.5,
        first_F_error_over_M_epsilon2=0.4,
        first_G_error_over_M_epsilon2=0.3,
        radial_F_over_M=0.9,
        axial_F_over_M=0.8,
        axial_G_over_M=0.7,
        evidence_kind="formal-theorem",
        provenance="formal LocalBase fixture",
        convex_domain_certified=True,
        differentiability_certified=True,
        representative_membership_certified=True,
        local_base_suprema_certified=True,
    )
    return LargeBandLocalBaseAdmission(scale=scale, witness=witness)


def _source(
    *,
    ell=100000,
    a=(7, -3, 11),
    M=2.0,
    source_id="formal-base-provider",
    source_revision="fixture-r1",
    evidence_kind="formal-theorem",
    provenance="formal source-identity fixture",
    provider_contract_certified=True,
    sign_independent_source_certified=True,
    local_base_source_identity_certified=True,
):
    return LargeBandBaseSourceWitness(
        ell=ell,
        a=a,
        M=M,
        source_id=source_id,
        source_revision=source_revision,
        evidence_kind=evidence_kind,
        provenance=provenance,
        provider_contract_certified=provider_contract_certified,
        sign_independent_source_certified=sign_independent_source_certified,
        local_base_source_identity_certified=local_base_source_identity_certified,
    )


def test_binding_turns_one_sign_free_source_into_the_exact_label_pair():
    binding = LargeBandBaseSourceBinding(admission=_admission(sigma=1), source=_source())

    minus, plus = binding.sign_pair
    assert minus.sigma == -1
    assert plus.sigma == 1
    assert minus.box_key == plus.box_key == binding.box_key == (100000, (7, -3, 11))
    assert binding.signed_label(-1) == minus
    assert binding.signed_label(1) == plus
    assert binding.source_key == ("formal-base-provider", "fixture-r1")
    assert all(binding.binding_checks().values())


def test_binding_is_independent_of_which_sign_carried_the_localbase_witness():
    plus_binding = LargeBandBaseSourceBinding(admission=_admission(sigma=1), source=_source())
    minus_binding = LargeBandBaseSourceBinding(admission=_admission(sigma=-1), source=_source())

    assert plus_binding.box_key == minus_binding.box_key
    assert plus_binding.source_key == minus_binding.source_key
    assert plus_binding.sign_pair == minus_binding.sign_pair


def test_binding_rejects_box_band_or_M_mismatch():
    with pytest.raises(ValueError, match="box key"):
        LargeBandBaseSourceBinding(admission=_admission(), source=_source(a=(7, -3, 12)))

    with pytest.raises(ValueError, match="box key"):
        LargeBandBaseSourceBinding(admission=_admission(), source=_source(ell=100001))

    with pytest.raises(ValueError, match="base source M"):
        LargeBandBaseSourceBinding(admission=_admission(), source=_source(M=3.0))


def test_source_identity_fails_closed_on_sampled_or_blank_evidence():
    for kind in ("sampled", "fitted", "numeric-scan", ""):
        with pytest.raises(ValueError, match="evidence_kind"):
            _source(evidence_kind=kind)

    for field in ("source_id", "source_revision", "provenance"):
        kwargs = {field: "   "}
        with pytest.raises(ValueError, match="nonempty"):
            _source(**kwargs)


def test_source_requires_all_three_identity_theorem_facts():
    for field in (
        "provider_contract_certified",
        "sign_independent_source_certified",
        "local_base_source_identity_certified",
    ):
        kwargs = {field: False}
        with pytest.raises(ValueError, match=field):
            _source(**kwargs)


def test_signed_label_rejects_non_paper_signs():
    binding = LargeBandBaseSourceBinding(admission=_admission(), source=_source())
    for sigma in (0, 2, -2):
        with pytest.raises(ValueError, match="sigma"):
            binding.signed_label(sigma)


def test_binding_does_not_upgrade_paper_truth_status():
    binding = LargeBandBaseSourceBinding(admission=_admission(), source=_source())

    assert binding.status == "formal-structure"
    assert binding.actual_base_fields_verified is False
    assert binding.uniform_eq_7_9_to_7_11_verified is False
    assert binding.paper_exact_velocity_available is False
