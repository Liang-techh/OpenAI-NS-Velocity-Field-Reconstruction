from fractions import Fraction

import pytest

from openai_ns_reconstruction.phase_large_band_base_source import (
    LargeBandBaseSourceBinding,
    LargeBandBaseSourceWitness,
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


def _physical_binding(a, *, prepared_id=PREPARED_ID, M_value=M):
    scale = LargeBandPhaseScaleCertificate(ell=ELL, h=Fraction(1, 200), M=M_value)
    local = LargeBandLocalBaseWitness(
        label=AsymptoticSlowLabel(ell=ELL, a=a, sigma=1),
        M=M_value,
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
            M=M_value,
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
            prepared_instance_id=prepared_id,
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
            prepared_instance_id=prepared_id,
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


def _physical_coverage(*, p0=None, p1=None):
    if p0 is None:
        p0 = _physical_binding(A0)
    if p1 is None:
        p1 = _physical_binding(A1)
    return LargeBandPhysicalBaseScopeCoverage(
        scope_admission=_scope_admission(),
        physical_base_bindings=(p0, p1),
    )


def _canonical(physical, **overrides):
    kwargs = dict(
        source_id=physical.source_key[0],
        source_revision=physical.source_key[1],
        prepared_instance_id=physical.prepared_key[2],
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
    return LargeBandCanonicalPhysicalBaseBinding(
        physical_base=physical,
        canonical=PrimaryGeometryCanonicalPreparedWitness(**kwargs),
    )


def test_exact_scope_uses_one_shared_canonical_prepared_for_all_boxes():
    p0 = _physical_binding(A0)
    p1 = _physical_binding(A1)
    coverage = LargeBandCanonicalPhysicalBaseScopeCoverage(
        physical_base_coverage=_physical_coverage(p0=p0, p1=p1),
        canonical_bindings=(_canonical(p1), _canonical(p0)),
    )

    assert coverage.canonical_binding_count == 2
    assert coverage.covered_box_keys == frozenset({(ELL, A0), (ELL, A1)})
    assert coverage.covered_signed_labels == frozenset(_labels())
    assert coverage.shared_prepared_key == (
        SOURCE_ID,
        SOURCE_REVISION,
        PREPARED_ID,
    )
    assert coverage.canonical_parameter_key == (4, 0.125, 17)
    assert coverage.canonical_scope_machine_checked is True


def test_per_box_canonical_bindings_cannot_use_distinct_prepared_instances():
    p0 = _physical_binding(A0, prepared_id="canonical-a")
    p1 = _physical_binding(A1, prepared_id="canonical-b")
    with pytest.raises(ValueError, match="one canonical Prepared instance"):
        LargeBandCanonicalPhysicalBaseScopeCoverage(
            physical_base_coverage=_physical_coverage(p0=p0, p1=p1),
            canonical_bindings=(_canonical(p0), _canonical(p1)),
        )


def test_shared_prepared_requires_one_canonical_constructor_parameter_tuple():
    p0 = _physical_binding(A0)
    p1 = _physical_binding(A1)
    with pytest.raises(ValueError, match="same canonical"):
        LargeBandCanonicalPhysicalBaseScopeCoverage(
            physical_base_coverage=_physical_coverage(p0=p0, p1=p1),
            canonical_bindings=(_canonical(p0), _canonical(p1, B=5)),
        )


def test_canonical_box_coverage_is_exact_and_has_no_duplicates():
    p0 = _physical_binding(A0)
    p1 = _physical_binding(A1)
    physical = _physical_coverage(p0=p0, p1=p1)
    c0 = _canonical(p0)
    c1 = _canonical(p1)

    with pytest.raises(ValueError, match="canonical box coverage must equal"):
        LargeBandCanonicalPhysicalBaseScopeCoverage(
            physical_base_coverage=physical,
            canonical_bindings=(c0,),
        )
    with pytest.raises(ValueError, match="exactly one binding per box"):
        LargeBandCanonicalPhysicalBaseScopeCoverage(
            physical_base_coverage=physical,
            canonical_bindings=(c0, c0, c1),
        )


def test_canonical_binding_must_wrap_the_physical_binding_from_the_scope():
    p0 = _physical_binding(A0)
    p1 = _physical_binding(A1)
    alternate_p0 = _physical_binding(A0, M_value=3.0)
    with pytest.raises(ValueError, match="exact physical-base binding"):
        LargeBandCanonicalPhysicalBaseScopeCoverage(
            physical_base_coverage=_physical_coverage(p0=p0, p1=p1),
            canonical_bindings=(_canonical(alternate_p0), _canonical(p1)),
        )


def test_scope_composition_does_not_upgrade_paper_exact_truth():
    p0 = _physical_binding(A0)
    p1 = _physical_binding(A1)
    coverage = LargeBandCanonicalPhysicalBaseScopeCoverage(
        physical_base_coverage=_physical_coverage(p0=p0, p1=p1),
        canonical_bindings=(_canonical(p0), _canonical(p1)),
    )

    assert coverage.status == "formal-structure"
    assert coverage.actual_canonical_prepared_application_machine_verified is False
    assert coverage.actual_cell_index_enumeration_machine_verified is False
    assert coverage.actual_physical_base_values_materialized is False
    assert coverage.actual_base_fields_verified is False
    assert coverage.theorem_applications_machine_replayed is False
    assert coverage.uniform_eq_7_9_to_7_11_verified is False
    assert coverage.paper_exact_velocity_available is False
