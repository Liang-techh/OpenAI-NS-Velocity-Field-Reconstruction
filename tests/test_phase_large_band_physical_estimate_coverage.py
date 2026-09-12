from dataclasses import replace
from fractions import Fraction

import pytest

from openai_ns_reconstruction.phase_large_band_base_source import (
    LargeBandBaseSourceBinding,
    LargeBandBaseSourceWitness,
)
from openai_ns_reconstruction.phase_large_band_coordinate_errors import (
    PINNED_BASE_ESTIMATES_THEOREM,
    PINNED_COORDINATE_ERRORS_THEOREM,
    LargeBandCoordinateErrorsAdmission,
    LargeBandCoordinateErrorsWitness,
)
from openai_ns_reconstruction.phase_large_band_damping import (
    PINNED_DAMPING_ERROR_THEOREM,
    LargeBandDampingErrorAdmission,
    LargeBandDampingErrorWitness,
)
from openai_ns_reconstruction.phase_large_band_family_coverage import (
    LargeBandActiveFamilyScopeWitness,
    LargeBandFamilyEstimateCoverage,
)
from openai_ns_reconstruction.phase_large_band_family_inputs import (
    LargeBandFamilyInputAdmission,
    LargeBandFamilyInputWitness,
)
from openai_ns_reconstruction.phase_large_band_frame import LargeBandPhaseFrameCertificate
from openai_ns_reconstruction.phase_large_band_local_base import (
    AsymptoticSlowLabel,
    LargeBandLocalBaseAdmission,
    LargeBandLocalBaseWitness,
)
from openai_ns_reconstruction.phase_large_band_phase_estimates import (
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_REPOSITORY,
    PINNED_PHASE_ESTIMATES_THEOREM,
    LargeBandPhaseEstimatesAdmission,
    LargeBandPhaseEstimatesWitness,
)
from openai_ns_reconstruction.phase_large_band_physical_base import (
    LargeBandPhysicalBaseBinding,
    PrimaryGeometryPhysicalBaseWitness,
)
from openai_ns_reconstruction.phase_large_band_physical_base_coverage import (
    LargeBandPhysicalBaseScopeCoverage,
)
from openai_ns_reconstruction.phase_large_band_physical_estimate_coverage import (
    LargeBandPhysicalEstimateScopeCoverage,
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
A = (7, -3, 11)
M = 2.0
SOURCE_ID = "formal-base-provider"
SOURCE_REVISION = "fixture-r1"


def _labels():
    return tuple(AsymptoticSlowLabel(ell=ELL, a=A, sigma=sigma) for sigma in (-1, 1))


def _scope(**overrides):
    kwargs = dict(
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
    kwargs.update(overrides)
    return LargeBandActiveFamilyScopeWitness(**kwargs)


def _primary_scope_admission(scope=None):
    if scope is None:
        scope = _scope()
    witness = PrimaryGeometryBandSliceWitness(
        scope_id=scope.scope_id,
        ell=ELL,
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        unsigned_boxes=(A,),
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


def _local_admission(*, sigma, M_value=M):
    scale = LargeBandPhaseScaleCertificate(ell=ELL, h=Fraction(1, 200), M=M_value)
    local = LargeBandLocalBaseWitness(
        label=AsymptoticSlowLabel(ell=ELL, a=A, sigma=sigma),
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
    return LargeBandLocalBaseAdmission(scale=scale, witness=local)


def _base_source(*, sigma, M_value=M):
    admission = _local_admission(sigma=sigma, M_value=M_value)
    source = LargeBandBaseSourceWitness(
        ell=ELL,
        a=A,
        M=M_value,
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        evidence_kind="formal-theorem",
        provenance="synthetic theorem fixture for sign-free source identity",
        provider_contract_certified=True,
        sign_independent_source_certified=True,
        local_base_source_identity_certified=True,
    )
    return LargeBandBaseSourceBinding(admission=admission, source=source)


def _physical_coverage(*, scope=None, M_value=M):
    if scope is None:
        scope = _scope()
    base_source = _base_source(sigma=1, M_value=M_value)
    prepared_id = "prepared-fixture-a"
    prepared = LargeBandPreparedSourceBinding(
        base_source=base_source,
        prepared=PrimaryGeometryPreparedSourceWitness(
            source_id=SOURCE_ID,
            source_revision=SOURCE_REVISION,
            prepared_instance_id=prepared_id,
            evidence_kind="formal-theorem",
            provenance="synthetic theorem fixture for Prepared identity",
            exists_prepared_application_certified=True,
            prepared_choice_identity_certified=True,
            family_base_source_identity_certified=True,
            phases_reuse_prepared_certified=True,
            both_signs_share_prepared_certified=True,
        ),
    )
    physical = LargeBandPhysicalBaseBinding(
        prepared_source=prepared,
        witness=PrimaryGeometryPhysicalBaseWitness(
            ell=ELL,
            a=A,
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
    return LargeBandPhysicalBaseScopeCoverage(
        scope_admission=_primary_scope_admission(scope),
        physical_base_bindings=(physical,),
    )


def _family(*, sigma):
    binding = _base_source(sigma=sigma)
    frame = LargeBandPhaseFrameCertificate(
        scale=binding.admission.scale, u=1.0, B=1.0, viscosity=1.5
    )
    witness = LargeBandFamilyInputWitness(
        label=binding.admission.witness.label,
        M=M,
        u=1.0,
        B=1.0,
        viscosity=1.5,
        reference_radius=1.0,
        Tg=0.1,
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        evidence_kind="formal-theorem",
        provenance="synthetic FamilyData input fixture",
        family_index_source_identity_certified=True,
        carrier_membership_certified=True,
        slot_membership_certified=True,
        representative_data_identity_certified=True,
        reference_scale_identity_certified=True,
        viscosity_identity_certified=True,
    )
    return LargeBandFamilyInputAdmission(binding=binding, frame=frame, witness=witness)


def _bundle(*, sigma):
    family = _family(sigma=sigma)
    label = family.label
    phase = LargeBandPhaseEstimatesAdmission(
        family=family,
        witness=LargeBandPhaseEstimatesWitness(
            label=label,
            source_id=SOURCE_ID,
            source_revision=SOURCE_REVISION,
            lean_repository=PINNED_LEAN_REPOSITORY,
            lean_commit=PINNED_LEAN_COMMIT,
            theorem_name=PINNED_PHASE_ESTIMATES_THEOREM,
            evidence_kind="formal-theorem",
            provenance="synthetic pinned phase_estimates fixture",
            theorem_application_certified=True,
            actual_normal_identity_certified=True,
            actual_phase_velocity_identity_certified=True,
            uniform_carrier_slot_conclusion_certified=True,
        ),
    )
    coordinate = LargeBandCoordinateErrorsAdmission(
        phase=phase,
        witness=LargeBandCoordinateErrorsWitness(
            label=label,
            source_id=SOURCE_ID,
            source_revision=SOURCE_REVISION,
            lean_repository=PINNED_LEAN_REPOSITORY,
            lean_commit=PINNED_LEAN_COMMIT,
            base_estimates_theorem_name=PINNED_BASE_ESTIMATES_THEOREM,
            coordinate_errors_theorem_name=PINNED_COORDINATE_ERRORS_THEOREM,
            evidence_kind="formal-theorem",
            provenance="synthetic pinned coordinate fixture",
            base_estimates_application_certified=True,
            actual_base_value_identity_certified=True,
            actual_shear_vector_identity_certified=True,
            uniform_carrier_base_conclusion_certified=True,
            coordinate_errors_application_certified=True,
            actual_frame_errors_identity_certified=True,
            phase_estimates_dependency_certified=True,
            base_estimates_dependency_certified=True,
            uniform_carrier_slot_conclusion_certified=True,
        ),
    )
    damping = LargeBandDampingErrorAdmission(
        family=family,
        witness=LargeBandDampingErrorWitness(
            label=label,
            source_id=SOURCE_ID,
            source_revision=SOURCE_REVISION,
            lean_repository=PINNED_LEAN_REPOSITORY,
            lean_commit=PINNED_LEAN_COMMIT,
            theorem_name=PINNED_DAMPING_ERROR_THEOREM,
            evidence_kind="formal-theorem",
            provenance="synthetic pinned damping fixture",
            damping_error_application_certified=True,
            carrier_membership_certified=True,
            actual_damping_term_identity_certified=True,
            reference_damping_term_identity_certified=True,
            viscosity_damping_family_identity_certified=True,
            reference_viscosity_identity_certified=True,
            reference_radius_center_identity_certified=True,
            scale_identity_certified=True,
            family_input_dependency_certified=True,
        ),
    )
    return phase, coordinate, damping


def _estimate_coverage(*, scope=None):
    if scope is None:
        scope = _scope()
    minus = _bundle(sigma=-1)
    plus = _bundle(sigma=1)
    return LargeBandFamilyEstimateCoverage(
        scope=scope,
        phase_admissions=(minus[0], plus[0]),
        coordinate_admissions=(minus[1], plus[1]),
        damping_admissions=(minus[2], plus[2]),
    )


def _coverage():
    scope = _scope()
    return LargeBandPhysicalEstimateScopeCoverage(
        physical_base_coverage=_physical_coverage(scope=scope),
        family_estimate_coverage=_estimate_coverage(scope=scope),
    )


def test_same_scope_composes_physical_base_and_all_estimate_paths():
    coverage = _coverage()

    assert coverage.covered_labels == frozenset(_labels())
    assert coverage.covered_box_keys == frozenset({(ELL, A)})
    assert coverage.source_key == (SOURCE_ID, SOURCE_REVISION)
    assert coverage.physical_estimate_scope_machine_checked is True
    assert coverage.theorem_envelopes() == coverage.family_estimate_coverage.uniform_scope_bounds()
    for label in _labels():
        assert coverage.physical_binding_for_label(label).box_key == label.box_key


def test_distinct_scope_metadata_cannot_be_cross_wired_even_with_same_labels():
    physical_scope = _scope()
    estimate_scope = replace(
        physical_scope,
        scope_id="different-formal-scope",
        provenance="different synthetic theorem scope fixture",
    )
    with pytest.raises(ValueError, match="same exact active scope witness"):
        LargeBandPhysicalEstimateScopeCoverage(
            physical_base_coverage=_physical_coverage(scope=physical_scope),
            family_estimate_coverage=_estimate_coverage(scope=estimate_scope),
        )


def test_localbase_M_must_match_between_estimate_and_physical_prepared_binding():
    scope = _scope()
    with pytest.raises(ValueError, match="same LocalBase M"):
        LargeBandPhysicalEstimateScopeCoverage(
            physical_base_coverage=_physical_coverage(scope=scope, M_value=3.0),
            family_estimate_coverage=_estimate_coverage(scope=scope),
        )


def test_label_lookup_is_typed_and_exact():
    coverage = _coverage()
    with pytest.raises(TypeError, match="AsymptoticSlowLabel"):
        coverage.physical_binding_for_label((ELL, A, 1))

    missing = AsymptoticSlowLabel(ell=ELL, a=(99, 0, 0), sigma=1)
    with pytest.raises(ValueError, match="exactly one physical-base binding"):
        coverage.physical_binding_for_label(missing)


def test_composition_does_not_promote_paper_exact_truth():
    coverage = _coverage()

    assert coverage.status == "formal-structure"
    assert coverage.actual_cell_index_enumeration_machine_verified is False
    assert coverage.actual_physical_base_values_materialized is False
    assert coverage.actual_base_fields_verified is False
    assert coverage.theorem_applications_machine_replayed is False
    assert coverage.global_all_band_uniformity_verified is False
    assert coverage.uniform_eq_7_9_to_7_11_verified is False
    assert coverage.paper_exact_velocity_available is False
