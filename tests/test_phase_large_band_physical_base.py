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
from openai_ns_reconstruction.phase_large_band_physical_base import (
    PINNED_AXIAL_EQ_PHYSICAL,
    PINNED_CARRIER_TIME_POSITIVE,
    PINNED_CELL_BAND,
    PINNED_CONSTRUCTION_AXIAL,
    PINNED_CONSTRUCTION_FREQUENCY,
    PINNED_FINAL_SLOW_BASE_VELOCITY,
    PINNED_FREQUENCY_EQ_PHYSICAL,
    PINNED_PREPARED_RADIUS_POS,
    LargeBandPhysicalBaseBinding,
    PrimaryGeometryPhysicalBaseWitness,
)
from openai_ns_reconstruction.phase_large_band_prepared_source import (
    LargeBandPreparedSourceBinding,
    PrimaryGeometryPreparedSourceWitness,
)
from openai_ns_reconstruction.phase_large_band_scale import LargeBandPhaseScaleCertificate


def _prepared_binding(*, source_id="formal-base-provider", source_revision="fixture-r1"):
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
    prepared = PrimaryGeometryPreparedSourceWitness(
        source_id=source_id,
        source_revision=source_revision,
        prepared_instance_id="prepared-fixture-001",
        evidence_kind="formal-theorem",
        provenance="synthetic theorem fixture for pinned Prepared identity",
        exists_prepared_application_certified=True,
        prepared_choice_identity_certified=True,
        family_base_source_identity_certified=True,
        phases_reuse_prepared_certified=True,
        both_signs_share_prepared_certified=True,
    )
    return LargeBandPreparedSourceBinding(base_source=base, prepared=prepared)


def _witness(**overrides):
    kwargs = dict(
        ell=100000,
        a=(7, -3, 11),
        source_id="formal-base-provider",
        source_revision="fixture-r1",
        prepared_instance_id="prepared-fixture-001",
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
    kwargs.update(overrides)
    return PrimaryGeometryPhysicalBaseWitness(**kwargs)


def test_physical_base_binding_composes_prepared_source_with_pinned_field_chain():
    binding = LargeBandPhysicalBaseBinding(
        prepared_source=_prepared_binding(), witness=_witness()
    )

    assert binding.box_key == (100000, (7, -3, 11))
    assert binding.source_key == ("formal-base-provider", "fixture-r1")
    assert binding.prepared_key == (
        "formal-base-provider",
        "fixture-r1",
        "prepared-fixture-001",
    )
    minus, plus = binding.sign_pair
    assert (minus.sigma, plus.sigma) == (-1, 1)
    assert binding.physical_base_identity_theorem_certified is True
    assert all(binding.binding_checks().values())


def test_physical_base_witness_pins_exact_formal_field_symbols():
    witness = _witness()
    assert witness.construction_frequency_theorem == PINNED_CONSTRUCTION_FREQUENCY
    assert witness.construction_axial_theorem == PINNED_CONSTRUCTION_AXIAL
    assert witness.carrier_time_positive_theorem == PINNED_CARRIER_TIME_POSITIVE
    assert witness.prepared_radius_pos_field == PINNED_PREPARED_RADIUS_POS
    assert witness.frequency_eq_physical_theorem == PINNED_FREQUENCY_EQ_PHYSICAL
    assert witness.axial_eq_physical_theorem == PINNED_AXIAL_EQ_PHYSICAL
    assert witness.physical_velocity_definition == PINNED_FINAL_SLOW_BASE_VELOCITY
    assert witness.cell_band_definition == PINNED_CELL_BAND

    for field in (
        "construction_frequency_theorem",
        "construction_axial_theorem",
        "carrier_time_positive_theorem",
        "prepared_radius_pos_field",
        "frequency_eq_physical_theorem",
        "axial_eq_physical_theorem",
        "physical_velocity_definition",
        "cell_band_definition",
    ):
        with pytest.raises(ValueError, match="pinned formal source"):
            _witness(**{field: "stale-or-wrong-pin"})


def test_physical_base_binding_rejects_box_source_or_prepared_drift():
    with pytest.raises(ValueError, match="box key"):
        LargeBandPhysicalBaseBinding(
            prepared_source=_prepared_binding(), witness=_witness(a=(8, -3, 11))
        )
    with pytest.raises(ValueError, match="source revision"):
        LargeBandPhysicalBaseBinding(
            prepared_source=_prepared_binding(), witness=_witness(source_revision="fixture-r2")
        )
    with pytest.raises(ValueError, match="Prepared identity"):
        LargeBandPhysicalBaseBinding(
            prepared_source=_prepared_binding(),
            witness=_witness(prepared_instance_id="different-prepared"),
        )


def test_physical_base_witness_rejects_non_theorem_evidence_and_missing_facts():
    for kind in ("sampled", "fitted", "numeric-scan", ""):
        with pytest.raises(ValueError, match="evidence_kind"):
            _witness(evidence_kind=kind)

    for field in (
        "local_base_same_frequency_axial_certified",
        "construction_frequency_certified",
        "construction_axial_certified",
        "carrier_time_positive_certified",
        "prepared_radius_positive_certified",
        "frequency_eq_physical_certified",
        "axial_eq_physical_certified",
        "both_signs_share_physical_base_certified",
    ):
        with pytest.raises(ValueError, match=field):
            _witness(**{field: False})


def test_physical_base_bridge_does_not_materialize_or_upgrade_truth_status():
    binding = LargeBandPhysicalBaseBinding(
        prepared_source=_prepared_binding(), witness=_witness()
    )

    assert binding.status == "formal-structure"
    assert binding.actual_physical_base_values_materialized is False
    assert binding.actual_base_fields_verified is False
    assert binding.uniform_eq_7_9_to_7_11_verified is False
    assert binding.paper_exact_velocity_available is False
