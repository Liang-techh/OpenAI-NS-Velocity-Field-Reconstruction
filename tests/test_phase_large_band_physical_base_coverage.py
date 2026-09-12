from fractions import Fraction

import pytest

from openai_ns_reconstruction.phase_large_band_base_source import (
    LargeBandBaseSourceBinding,
    LargeBandBaseSourceWitness,
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
A2 = (9, -3, 11)
SOURCE_ID = "formal-base-provider"
SOURCE_REVISION = "fixture-r1"
M = 2.0


def _labels(boxes=(A0, A1)):
    return tuple(
        AsymptoticSlowLabel(ell=ELL, a=box, sigma=sigma)
        for box in boxes
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


def _physical_binding(
    a,
    *,
    source_id=SOURCE_ID,
    source_revision=SOURCE_REVISION,
    prepared_instance_id=None,
):
    if prepared_instance_id is None:
        prepared_instance_id = f"prepared-{a[0]}-{a[1]}-{a[2]}"

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
        base_source=base_source,
        prepared=PrimaryGeometryPreparedSourceWitness(
            source_id=source_id,
            source_revision=source_revision,
            prepared_instance_id=prepared_instance_id,
            evidence_kind="formal-theorem",
            provenance="synthetic theorem fixture for pinned Prepared identity",
            exists_prepared_application_certified=True,
            prepared_choice_identity_certified=True,
            family_base_source_identity_certified=True,
            phases_reuse_prepared_certified=True,
            both_signs_share_prepared_certified=True,
        ),
    )
    physical = PrimaryGeometryPhysicalBaseWitness(
        ell=ELL,
        a=a,
        source_id=source_id,
        source_revision=source_revision,
        prepared_instance_id=prepared_instance_id,
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
    return LargeBandPhysicalBaseBinding(prepared_source=prepared, witness=physical)


def test_exact_physical_bindings_cover_every_unsigned_box_and_both_signs():
    coverage = LargeBandPhysicalBaseScopeCoverage(
        scope_admission=_scope_admission(),
        physical_base_bindings=(_physical_binding(A1), _physical_binding(A0)),
    )

    assert coverage.physical_binding_count == 2
    assert coverage.covered_box_keys == frozenset({(ELL, A0), (ELL, A1)})
    assert coverage.covered_signed_labels == frozenset(_labels())
    assert coverage.declared_scope_physical_base_coverage_machine_checked is True


def test_missing_or_extra_physical_box_fails_exact_set_coverage():
    with pytest.raises(ValueError, match="physical-base box coverage must equal"):
        LargeBandPhysicalBaseScopeCoverage(
            scope_admission=_scope_admission(),
            physical_base_bindings=(_physical_binding(A0),),
        )

    with pytest.raises(ValueError, match="physical-base box coverage must equal"):
        LargeBandPhysicalBaseScopeCoverage(
            scope_admission=_scope_admission(),
            physical_base_bindings=(
                _physical_binding(A0),
                _physical_binding(A1),
                _physical_binding(A2),
            ),
        )


def test_duplicate_box_cannot_masquerade_as_scope_coverage():
    first = _physical_binding(A0)
    with pytest.raises(ValueError, match="exactly one binding per box"):
        LargeBandPhysicalBaseScopeCoverage(
            scope_admission=_scope_admission(),
            physical_base_bindings=(first, first, _physical_binding(A1)),
        )


def test_every_box_must_reuse_the_scope_source_revision():
    with pytest.raises(ValueError, match="scope source revision"):
        LargeBandPhysicalBaseScopeCoverage(
            scope_admission=_scope_admission(),
            physical_base_bindings=(
                _physical_binding(A0),
                _physical_binding(A1, source_revision="fixture-r2"),
            ),
        )


def test_coverage_requires_typed_nonempty_binding_tuple():
    with pytest.raises(TypeError, match="must be a tuple"):
        LargeBandPhysicalBaseScopeCoverage(
            scope_admission=_scope_admission(),
            physical_base_bindings=[_physical_binding(A0), _physical_binding(A1)],
        )

    with pytest.raises(ValueError, match="must be nonempty"):
        LargeBandPhysicalBaseScopeCoverage(
            scope_admission=_scope_admission(),
            physical_base_bindings=(),
        )


def test_exact_set_check_does_not_upgrade_the_upstream_truth_boundary():
    coverage = LargeBandPhysicalBaseScopeCoverage(
        scope_admission=_scope_admission(),
        physical_base_bindings=(_physical_binding(A0), _physical_binding(A1)),
    )

    assert coverage.status == "formal-structure"
    assert coverage.actual_cell_index_enumeration_machine_verified is False
    assert coverage.actual_physical_base_values_materialized is False
    assert coverage.actual_base_fields_verified is False
    assert coverage.uniform_eq_7_9_to_7_11_verified is False
    assert coverage.paper_exact_velocity_available is False
