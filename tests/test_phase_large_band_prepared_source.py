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
from openai_ns_reconstruction.phase_large_band_prepared_source import (
    PINNED_CONSTRUCTION_DEF,
    PINNED_EXISTS_PREPARED,
    PINNED_FAMILY_DEF,
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_FILE,
    PINNED_LEAN_REPOSITORY,
    PINNED_PHASES_DEF,
    PINNED_PREPARED_DECL,
    PINNED_PREPARED_DEF,
    LargeBandPreparedSourceBinding,
    PrimaryGeometryPreparedSourceWitness,
)
from openai_ns_reconstruction.phase_large_band_scale import LargeBandPhaseScaleCertificate


def _base_binding(*, sigma=1, source_id="formal-base-provider", source_revision="fixture-r1"):
    ell = 100000
    a = (7, -3, 11)
    M = 2.0
    scale = LargeBandPhaseScaleCertificate(ell=ell, h=Fraction(1, 200), M=M)
    local = LargeBandLocalBaseWitness(
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
        provenance="synthetic theorem fixture for LocalBase admission",
        convex_domain_certified=True,
        differentiability_certified=True,
        representative_membership_certified=True,
        local_base_suprema_certified=True,
    )
    admission = LargeBandLocalBaseAdmission(scale=scale, witness=local)
    source = LargeBandBaseSourceWitness(
        ell=ell,
        a=a,
        M=M,
        source_id=source_id,
        source_revision=source_revision,
        evidence_kind="formal-theorem",
        provenance="synthetic theorem fixture for source identity",
        provider_contract_certified=True,
        sign_independent_source_certified=True,
        local_base_source_identity_certified=True,
    )
    return LargeBandBaseSourceBinding(admission=admission, source=source)


def _prepared(**overrides):
    kwargs = dict(
        source_id="formal-base-provider",
        source_revision="fixture-r1",
        prepared_instance_id="prepared-fixture-001",
        evidence_kind="formal-theorem",
        provenance="synthetic theorem fixture for pinned Prepared identity",
        exists_prepared_application_certified=True,
        prepared_choice_identity_certified=True,
        family_base_source_identity_certified=True,
        phases_reuse_prepared_certified=True,
        both_signs_share_prepared_certified=True,
    )
    kwargs.update(overrides)
    return PrimaryGeometryPreparedSourceWitness(**kwargs)


def test_prepared_binding_preserves_source_revision_and_both_signs():
    binding = LargeBandPreparedSourceBinding(base_source=_base_binding(), prepared=_prepared())

    assert binding.source_key == ("formal-base-provider", "fixture-r1")
    assert binding.prepared_key == (
        "formal-base-provider",
        "fixture-r1",
        "prepared-fixture-001",
    )
    minus, plus = binding.sign_pair
    assert (minus.sigma, plus.sigma) == (-1, 1)
    assert minus.box_key == plus.box_key == (100000, (7, -3, 11))
    assert all(binding.binding_checks().values())


def test_prepared_binding_is_independent_of_localbase_sign():
    plus = LargeBandPreparedSourceBinding(base_source=_base_binding(sigma=1), prepared=_prepared())
    minus = LargeBandPreparedSourceBinding(base_source=_base_binding(sigma=-1), prepared=_prepared())

    assert plus.prepared_key == minus.prepared_key
    assert plus.sign_pair == minus.sign_pair


def test_prepared_binding_rejects_source_or_revision_drift():
    with pytest.raises(ValueError, match="source identity"):
        LargeBandPreparedSourceBinding(
            base_source=_base_binding(),
            prepared=_prepared(source_id="different-provider"),
        )

    with pytest.raises(ValueError, match="source identity"):
        LargeBandPreparedSourceBinding(
            base_source=_base_binding(),
            prepared=_prepared(source_revision="fixture-r2"),
        )


def test_prepared_witness_pins_exact_formal_source_symbols():
    assert _prepared().lean_repository == PINNED_LEAN_REPOSITORY
    assert _prepared().lean_commit == PINNED_LEAN_COMMIT
    assert _prepared().lean_file == PINNED_LEAN_FILE
    assert _prepared().prepared_decl == PINNED_PREPARED_DECL
    assert _prepared().exists_prepared_theorem == PINNED_EXISTS_PREPARED
    assert _prepared().prepared_def == PINNED_PREPARED_DEF
    assert _prepared().family_def == PINNED_FAMILY_DEF
    assert _prepared().construction_def == PINNED_CONSTRUCTION_DEF
    assert _prepared().phases_def == PINNED_PHASES_DEF

    for field in (
        "lean_repository",
        "lean_commit",
        "lean_file",
        "prepared_decl",
        "exists_prepared_theorem",
        "prepared_def",
        "family_def",
        "construction_def",
        "phases_def",
    ):
        with pytest.raises(ValueError, match="pinned formal source"):
            _prepared(**{field: "stale-or-wrong-pin"})


def test_prepared_witness_rejects_sampled_or_blank_provenance():
    for kind in ("sampled", "fitted", "numeric-scan", ""):
        with pytest.raises(ValueError, match="evidence_kind"):
            _prepared(evidence_kind=kind)

    for field in ("source_id", "source_revision", "prepared_instance_id", "provenance"):
        with pytest.raises(ValueError, match="nonempty"):
            _prepared(**{field: "   "})


def test_prepared_witness_requires_each_formal_identity_fact():
    for field in (
        "exists_prepared_application_certified",
        "prepared_choice_identity_certified",
        "family_base_source_identity_certified",
        "phases_reuse_prepared_certified",
        "both_signs_share_prepared_certified",
    ):
        with pytest.raises(ValueError, match=field):
            _prepared(**{field: False})


def test_prepared_bridge_does_not_upgrade_paper_truth_status():
    binding = LargeBandPreparedSourceBinding(base_source=_base_binding(), prepared=_prepared())

    assert binding.status == "formal-structure"
    assert binding.actual_prepared_application_machine_verified is False
    assert binding.actual_base_fields_verified is False
    assert binding.uniform_eq_7_9_to_7_11_verified is False
    assert binding.paper_exact_velocity_available is False
