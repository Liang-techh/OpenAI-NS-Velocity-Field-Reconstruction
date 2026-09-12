from dataclasses import replace
from fractions import Fraction

import pytest

from openai_ns_reconstruction.phase_large_band_base_source import (
    LargeBandBaseSourceBinding,
    LargeBandBaseSourceWitness,
)
from openai_ns_reconstruction.phase_large_band_canonical_application_export import (
    EXPORT_PRODUCER_KIND,
    EXPORT_SCHEMA,
    REQUIRED_APPLICATION_SYMBOLS,
    CanonicalPreparedApplicationExport,
    LargeBandCanonicalApplicationExportAdmission,
    canonical_application_export_sha256,
    canonical_scope_sha256,
)
from openai_ns_reconstruction.phase_large_band_canonical_physical_base import (
    LargeBandCanonicalPhysicalBaseBinding,
    PrimaryGeometryCanonicalPreparedWitness,
)
from openai_ns_reconstruction.phase_large_band_canonical_scope_coverage import (
    LargeBandCanonicalPhysicalBaseScopeCoverage,
)
from openai_ns_reconstruction.phase_large_band_family_coverage import (
    LargeBandActiveFamilyScopeWitness,
)
from openai_ns_reconstruction.phase_large_band_local_base import (
    AsymptoticSlowLabel,
    LargeBandLocalBaseAdmission,
    LargeBandLocalBaseWitness,
)
from openai_ns_reconstruction.phase_large_band_phase_estimates import (
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_REPOSITORY,
)
from openai_ns_reconstruction.phase_large_band_physical_base import (
    LargeBandPhysicalBaseBinding,
    PrimaryGeometryPhysicalBaseWitness,
)
from openai_ns_reconstruction.phase_large_band_physical_base_coverage import (
    LargeBandPhysicalBaseScopeCoverage,
)
from openai_ns_reconstruction.phase_large_band_prepared_source import (
    LargeBandPreparedSourceBinding,
    PrimaryGeometryPreparedSourceWitness,
)
from openai_ns_reconstruction.phase_large_band_primary_geometry_scope import (
    PINNED_CELL_DOMAIN_DEFINITION,
    PINNED_CELL_INDEX_DEFINITION,
    PINNED_FAMILY_DEFINITION,
    LargeBandPrimaryGeometryScopeAdmission,
    PrimaryGeometryBandSliceWitness,
)
from openai_ns_reconstruction.phase_large_band_scale import LargeBandPhaseScaleCertificate


ELL = 100000
A0 = (7, -3, 11)
A1 = (8, -3, 11)
SOURCE_ID = "formal-base-provider"
SOURCE_REVISION = "fixture-r1"
PREPARED_ID = "canonical-prepared-fixture"
M = 2.0


def _labels():
    return tuple(
        AsymptoticSlowLabel(ell=ELL, a=a, sigma=sigma)
        for a in (A0, A1)
        for sigma in (-1, 1)
    )


def _scope_admission():
    scope = LargeBandActiveFamilyScopeWitness(
        scope_id="fixture-primary-geometry-band-slice",
        ell=ELL,
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        labels=_labels(),
        evidence_kind="formal-theorem",
        provenance="synthetic formal scope fixture; not manuscript data",
        finite_scope_certified=True,
        exact_active_signed_family_certified=True,
        carrier_slot_scope_certified=True,
        sign_duplication_certified=True,
    )
    witness = PrimaryGeometryBandSliceWitness(
        scope_id=scope.scope_id,
        ell=ELL,
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        unsigned_boxes=(A0, A1),
        lean_repository=PINNED_LEAN_REPOSITORY,
        lean_commit=PINNED_LEAN_COMMIT,
        cell_index_definition=PINNED_CELL_INDEX_DEFINITION,
        cell_domain_definition=PINNED_CELL_DOMAIN_DEFINITION,
        family_definition=PINNED_FAMILY_DEFINITION,
        evidence_kind="formal-theorem",
        provenance="synthetic theorem export fixture; not a CellIndex enumeration",
        complete_unsigned_box_manifest_certified=True,
        cell_index_band_slice_certified=True,
        band_and_cell_projection_certified=True,
        cell_domain_carrier_identity_certified=True,
        family_index_identity_certified=True,
        fin2_sign_duplication_certified=True,
    )
    return LargeBandPrimaryGeometryScopeAdmission(scope=scope, witness=witness)


def _physical_binding(a):
    scale = LargeBandPhaseScaleCertificate(ell=ELL, h=Fraction(1, 200), M=M)
    local = LargeBandLocalBaseWitness(
        label=AsymptoticSlowLabel(ell=ELL, a=a, sigma=1),
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
    base_source = LargeBandBaseSourceBinding(
        admission=LargeBandLocalBaseAdmission(scale=scale, witness=local),
        source=LargeBandBaseSourceWitness(
            ell=ELL,
            a=a,
            M=M,
            source_id=SOURCE_ID,
            source_revision=SOURCE_REVISION,
            evidence_kind="formal-theorem",
            provenance="synthetic theorem fixture for source identity",
            provider_contract_certified=True,
            sign_independent_source_certified=True,
            local_base_source_identity_certified=True,
        ),
    )
    prepared = LargeBandPreparedSourceBinding(
        base_source=base_source,
        prepared=PrimaryGeometryPreparedSourceWitness(
            source_id=SOURCE_ID,
            source_revision=SOURCE_REVISION,
            prepared_instance_id=PREPARED_ID,
            evidence_kind="formal-theorem",
            provenance="synthetic theorem fixture for one shared Prepared identity",
            exists_prepared_application_certified=True,
            prepared_choice_identity_certified=True,
            family_base_source_identity_certified=True,
            phases_reuse_prepared_certified=True,
            both_signs_share_prepared_certified=True,
        ),
    )
    return LargeBandPhysicalBaseBinding(
        prepared_source=prepared,
        witness=PrimaryGeometryPhysicalBaseWitness(
            ell=ELL,
            a=a,
            source_id=SOURCE_ID,
            source_revision=SOURCE_REVISION,
            prepared_instance_id=PREPARED_ID,
            evidence_kind="formal-theorem",
            provenance="synthetic theorem fixture for physical-base chain",
            local_base_same_frequency_axial_certified=True,
            construction_frequency_certified=True,
            construction_axial_certified=True,
            carrier_time_positive_certified=True,
            prepared_radius_positive_certified=True,
            frequency_eq_physical_certified=True,
            axial_eq_physical_certified=True,
            both_signs_share_physical_base_certified=True,
        ),
    )


def _canonical(physical):
    return LargeBandCanonicalPhysicalBaseBinding(
        physical_base=physical,
        canonical=PrimaryGeometryCanonicalPreparedWitness(
            source_id=SOURCE_ID,
            source_revision=SOURCE_REVISION,
            prepared_instance_id=PREPARED_ID,
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
        ),
    )


def _coverage():
    p0 = _physical_binding(A0)
    p1 = _physical_binding(A1)
    physical = LargeBandPhysicalBaseScopeCoverage(
        scope_admission=_scope_admission(),
        physical_base_bindings=(p0, p1),
    )
    return LargeBandCanonicalPhysicalBaseScopeCoverage(
        physical_base_coverage=physical,
        canonical_bindings=(_canonical(p1), _canonical(p0)),
    )


def _export(coverage, **overrides):
    B, r0, N0 = coverage.canonical_parameter_key
    source_id, source_revision, prepared_id = coverage.shared_prepared_key
    kwargs = dict(
        schema=EXPORT_SCHEMA,
        producer_kind=EXPORT_PRODUCER_KIND,
        formal_repository=PINNED_LEAN_REPOSITORY,
        formal_commit=PINNED_LEAN_COMMIT,
        source_id=source_id,
        source_revision=source_revision,
        prepared_instance_id=prepared_id,
        prepared_N=23,
        B=B,
        r0_hex=float(r0).hex(),
        N0=N0,
        unsigned_boxes=tuple(coverage.covered_box_keys),
        signed_labels=tuple(
            (label.ell, label.a, label.sigma)
            for label in coverage.covered_signed_labels
        ),
        scope_sha256=canonical_scope_sha256(coverage),
        application_symbols=REQUIRED_APPLICATION_SYMBOLS,
        provenance="synthetic external formal-export fixture; not a replayed Lean artifact",
        payload_sha256="0" * 64,
    )
    kwargs.update(overrides)
    raw = CanonicalPreparedApplicationExport(**kwargs)
    return replace(raw, payload_sha256=canonical_application_export_sha256(raw))


def test_exact_content_addressed_export_is_admitted_without_truth_upgrade():
    coverage = _coverage()
    export = _export(coverage)
    admitted = LargeBandCanonicalApplicationExportAdmission(coverage, export)

    assert admitted.export_integrity_machine_checked is True
    assert admitted.prepared_N == 23
    assert admitted.status == "formal-structure"
    assert admitted.actual_canonical_prepared_application_machine_verified is False
    assert admitted.actual_cell_index_enumeration_machine_verified is False
    assert admitted.actual_physical_base_values_materialized is False
    assert admitted.actual_base_fields_verified is False
    assert admitted.theorem_applications_machine_replayed is False
    assert admitted.uniform_eq_7_9_to_7_11_verified is False
    assert admitted.paper_exact_velocity_available is False


def test_payload_digest_detects_post_export_tampering():
    coverage = _coverage()
    export = _export(coverage)
    tampered = replace(export, prepared_N=export.prepared_N + 1)
    with pytest.raises(ValueError, match="payload SHA-256"):
        LargeBandCanonicalApplicationExportAdmission(coverage, tampered)


def test_scope_digest_prevents_cross_scope_relabeling():
    coverage = _coverage()
    export = _export(coverage, scope_sha256="1" * 64)
    with pytest.raises(ValueError, match="scope digest"):
        LargeBandCanonicalApplicationExportAdmission(coverage, export)


def test_export_must_cover_exact_unsigned_and_signed_scope():
    coverage = _coverage()
    export = _export(coverage)

    missing_box = _export(coverage, unsigned_boxes=export.unsigned_boxes[:-1])
    with pytest.raises(ValueError, match="unsigned box set"):
        LargeBandCanonicalApplicationExportAdmission(coverage, missing_box)

    missing_label = _export(coverage, signed_labels=export.signed_labels[:-1])
    with pytest.raises(ValueError, match="signed label set"):
        LargeBandCanonicalApplicationExportAdmission(coverage, missing_label)


def test_export_must_match_prepared_identity_and_canonical_parameters():
    coverage = _coverage()

    wrong_prepared = _export(coverage, prepared_instance_id="other-prepared")
    with pytest.raises(ValueError, match="Prepared identity"):
        LargeBandCanonicalApplicationExportAdmission(coverage, wrong_prepared)

    wrong_B = _export(coverage, B=5)
    with pytest.raises(ValueError, match="canonical .* parameters"):
        LargeBandCanonicalApplicationExportAdmission(coverage, wrong_B)


def test_export_requires_pinned_symbols_and_formal_producer_kind():
    coverage = _coverage()
    export = _export(coverage, application_symbols=REQUIRED_APPLICATION_SYMBOLS[:-1])
    with pytest.raises(ValueError, match="application symbols"):
        LargeBandCanonicalApplicationExportAdmission(coverage, export)

    with pytest.raises(ValueError, match="lean-formal-export"):
        CanonicalPreparedApplicationExport(
            schema=EXPORT_SCHEMA,
            producer_kind="numeric-scan",
            formal_repository=PINNED_LEAN_REPOSITORY,
            formal_commit=PINNED_LEAN_COMMIT,
            source_id=SOURCE_ID,
            source_revision=SOURCE_REVISION,
            prepared_instance_id=PREPARED_ID,
            prepared_N=23,
            B=4,
            r0_hex=(0.125).hex(),
            N0=17,
            unsigned_boxes=((ELL, A0),),
            signed_labels=((ELL, A0, -1), (ELL, A0, 1)),
            scope_sha256="0" * 64,
            application_symbols=REQUIRED_APPLICATION_SYMBOLS,
            provenance="rejected numeric fixture",
            payload_sha256="0" * 64,
        )


def test_prepared_N_is_recorded_but_not_misrepresented_as_replayed_proof():
    coverage = _coverage()
    admitted = LargeBandCanonicalApplicationExportAdmission(
        coverage,
        _export(coverage, prepared_N=0),
    )
    assert admitted.prepared_N == 0
    assert admitted.actual_canonical_prepared_application_machine_verified is False
    assert admitted.theorem_applications_machine_replayed is False
