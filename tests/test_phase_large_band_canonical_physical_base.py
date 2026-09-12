from fractions import Fraction

import pytest

from openai_ns_reconstruction.phase_large_band_base_source import (
    LargeBandBaseSourceBinding,
    LargeBandBaseSourceWitness,
)
from openai_ns_reconstruction.phase_large_band_canonical_physical_base import (
    PINNED_ACTIVE_RIGHT_DEF,
    PINNED_AXIAL_DEF,
    PINNED_CANONICAL_PHASES_DEF,
    PINNED_CANONICAL_PREPARED_DEF,
    PINNED_FINAL_SLOW_BASE_SCALES_DEF,
    PINNED_FREQUENCY_DEF,
    PINNED_FULL_TRUE_CONE,
    LargeBandCanonicalPhysicalBaseBinding,
    PrimaryGeometryCanonicalPreparedWitness,
)
from openai_ns_reconstruction.phase_large_band_local_base import (
    AsymptoticSlowLabel,
    LargeBandLocalBaseAdmission,
    LargeBandLocalBaseWitness,
)
from openai_ns_reconstruction.phase_large_band_physical_base import (
    LargeBandPhysicalBaseBinding,
    PrimaryGeometryPhysicalBaseWitness,
)
from openai_ns_reconstruction.phase_large_band_prepared_source import (
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_FILE,
    PINNED_LEAN_REPOSITORY,
    LargeBandPreparedSourceBinding,
    PrimaryGeometryPreparedSourceWitness,
)
from openai_ns_reconstruction.phase_large_band_scale import LargeBandPhaseScaleCertificate


def _physical_base(*, source_id="formal-base-provider", source_revision="fixture-r1", prepared_id="prepared-fixture-001"):
    ell = 100000
    a = (7, -3, 11)
    M = 2.0
    scale = LargeBandPhaseScaleCertificate(ell=ell, h=Fraction(1, 200), M=M)
    local = LargeBandLocalBaseWitness(
        label=AsymptoticSlowLabel(ell=ell, a=a, sigma=1),
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
    base = LargeBandBaseSourceBinding(
        admission=LargeBandLocalBaseAdmission(scale=scale, witness=local),
        source=LargeBandBaseSourceWitness(
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
        ),
    )
    prepared = LargeBandPreparedSourceBinding(
        base_source=base,
        prepared=PrimaryGeometryPreparedSourceWitness(
            source_id=source_id,
            source_revision=source_revision,
            prepared_instance_id=prepared_id,
            evidence_kind="formal-theorem",
            provenance="synthetic theorem fixture for pinned Prepared identity",
            exists_prepared_application_certified=True,
            prepared_choice_identity_certified=True,
            family_base_source_identity_certified=True,
            phases_reuse_prepared_certified=True,
            both_signs_share_prepared_certified=True,
        ),
    )
    witness = PrimaryGeometryPhysicalBaseWitness(
        ell=ell,
        a=a,
        source_id=source_id,
        source_revision=source_revision,
        prepared_instance_id=prepared_id,
        evidence_kind="formal-theorem",
        provenance="synthetic theorem fixture for physical-base theorem chain",
        local_base_same_frequency_axial_certified=True,
        construction_frequency_certified=True,
        construction_axial_certified=True,
        carrier_time_positive_certified=True,
        prepared_radius_positive_certified=True,
        frequency_eq_physical_certified=True,
        axial_eq_physical_certified=True,
        both_signs_share_physical_base_certified=True,
    )
    return LargeBandPhysicalBaseBinding(prepared_source=prepared, witness=witness)


def _canonical(**overrides):
    kwargs = dict(
        source_id="formal-base-provider",
        source_revision="fixture-r1",
        prepared_instance_id="prepared-fixture-001",
        B=4,
        r0=0.125,
        N0=17,
        evidence_kind="formal-theorem",
        provenance="synthetic theorem fixture for canonical Prepared identity",
        full_true_cone_witness_certified=True,
        canonical_prepared_application_certified=True,
        canonical_phases_application_certified=True,
        canonical_upper_twice_active_right_certified=True,
        frequency_uses_canonical_scales_certified=True,
        axial_uses_canonical_scales_certified=True,
        both_signs_share_canonical_phases_certified=True,
    )
    kwargs.update(overrides)
    return PrimaryGeometryCanonicalPreparedWitness(**kwargs)


def test_canonical_binding_connects_physical_base_to_canonical_prepared():
    binding = LargeBandCanonicalPhysicalBaseBinding(
        physical_base=_physical_base(), canonical=_canonical()
    )

    assert binding.box_key == (100000, (7, -3, 11))
    assert binding.source_key == ("formal-base-provider", "fixture-r1")
    assert binding.prepared_key == (
        "formal-base-provider",
        "fixture-r1",
        "prepared-fixture-001",
    )
    assert binding.canonical_parameter_key == (4, 0.125, 17)
    minus, plus = binding.sign_pair
    assert (minus.sigma, plus.sigma) == (-1, 1)
    assert binding.canonical_physical_base_identity_theorem_certified is True
    assert all(binding.binding_checks().values())


def test_canonical_witness_pins_exact_formal_symbols():
    witness = _canonical()
    assert witness.lean_repository == PINNED_LEAN_REPOSITORY
    assert witness.lean_commit == PINNED_LEAN_COMMIT
    assert witness.lean_file == PINNED_LEAN_FILE
    assert witness.canonical_prepared_def == PINNED_CANONICAL_PREPARED_DEF
    assert witness.canonical_phases_def == PINNED_CANONICAL_PHASES_DEF
    assert witness.active_right_def == PINNED_ACTIVE_RIGHT_DEF
    assert witness.frequency_def == PINNED_FREQUENCY_DEF
    assert witness.axial_def == PINNED_AXIAL_DEF
    assert witness.final_slow_base_scales_def == PINNED_FINAL_SLOW_BASE_SCALES_DEF
    assert witness.full_true_cone_decl == PINNED_FULL_TRUE_CONE

    for field in (
        "lean_repository",
        "lean_commit",
        "lean_file",
        "canonical_prepared_def",
        "canonical_phases_def",
        "active_right_def",
        "frequency_def",
        "axial_def",
        "final_slow_base_scales_def",
        "full_true_cone_decl",
    ):
        with pytest.raises(ValueError, match="pinned formal source"):
            _canonical(**{field: "stale-or-wrong-pin"})


def test_canonical_binding_rejects_source_or_prepared_drift():
    with pytest.raises(ValueError, match="source revision"):
        LargeBandCanonicalPhysicalBaseBinding(
            physical_base=_physical_base(), canonical=_canonical(source_revision="fixture-r2")
        )
    with pytest.raises(ValueError, match="Prepared identity"):
        LargeBandCanonicalPhysicalBaseBinding(
            physical_base=_physical_base(),
            canonical=_canonical(prepared_instance_id="different-prepared"),
        )


def test_canonical_witness_rejects_non_theorem_evidence_and_bad_parameters():
    for kind in ("sampled", "fitted", "numeric-scan", ""):
        with pytest.raises(ValueError, match="evidence_kind"):
            _canonical(evidence_kind=kind)

    for field in ("source_id", "source_revision", "prepared_instance_id", "provenance"):
        with pytest.raises(ValueError, match="nonempty"):
            _canonical(**{field: "   "})

    for field in ("B", "N0"):
        with pytest.raises(ValueError, match="nonnegative integer"):
            _canonical(**{field: -1})
    for bad_r0 in (0.0, -1.0, float("inf"), float("nan")):
        with pytest.raises(ValueError, match="r0"):
            _canonical(r0=bad_r0)


def test_canonical_witness_requires_each_formal_identity_fact():
    for field in (
        "full_true_cone_witness_certified",
        "canonical_prepared_application_certified",
        "canonical_phases_application_certified",
        "canonical_upper_twice_active_right_certified",
        "frequency_uses_canonical_scales_certified",
        "axial_uses_canonical_scales_certified",
        "both_signs_share_canonical_phases_certified",
    ):
        with pytest.raises(ValueError, match=field):
            _canonical(**{field: False})


def test_canonical_bridge_does_not_upgrade_paper_truth_status():
    binding = LargeBandCanonicalPhysicalBaseBinding(
        physical_base=_physical_base(), canonical=_canonical()
    )

    assert binding.status == "formal-structure"
    assert binding.actual_canonical_prepared_application_machine_verified is False
    assert binding.actual_physical_base_values_materialized is False
    assert binding.actual_base_fields_verified is False
    assert binding.uniform_eq_7_9_to_7_11_verified is False
    assert binding.paper_exact_velocity_available is False
