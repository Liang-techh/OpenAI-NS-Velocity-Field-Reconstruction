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
from openai_ns_reconstruction.phase_large_band_scale import LargeBandPhaseScaleCertificate


ELL = 100000
A = (7, -3, 11)
M = 2.0
SOURCE_ID = "formal-base-provider"
SOURCE_REVISION = "fixture-r1"


def _family(*, sigma: int, family_provenance: str = "formal FamilyData input fixture"):
    scale = LargeBandPhaseScaleCertificate(ell=ELL, h=Fraction(1, 200), M=M)
    label = AsymptoticSlowLabel(ell=ELL, a=A, sigma=sigma)
    local = LargeBandLocalBaseWitness(
        label=label,
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
    local_admission = LargeBandLocalBaseAdmission(scale=scale, witness=local)
    source = LargeBandBaseSourceWitness(
        ell=ELL,
        a=A,
        M=M,
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        evidence_kind="formal-theorem",
        provenance="formal source fixture",
        provider_contract_certified=True,
        sign_independent_source_certified=True,
        local_base_source_identity_certified=True,
    )
    binding = LargeBandBaseSourceBinding(admission=local_admission, source=source)
    frame = LargeBandPhaseFrameCertificate(scale=scale, u=1.0, B=1.0, viscosity=1.5)
    family_witness = LargeBandFamilyInputWitness(
        label=label,
        M=M,
        u=1.0,
        B=1.0,
        viscosity=1.5,
        reference_radius=1.0,
        Tg=0.1,
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        evidence_kind="formal-theorem",
        provenance=family_provenance,
        family_index_source_identity_certified=True,
        carrier_membership_certified=True,
        slot_membership_certified=True,
        representative_data_identity_certified=True,
        reference_scale_identity_certified=True,
        viscosity_identity_certified=True,
    )
    return LargeBandFamilyInputAdmission(binding=binding, frame=frame, witness=family_witness)


def _bundle(*, sigma: int, family_provenance: str = "formal FamilyData input fixture"):
    family = _family(sigma=sigma, family_provenance=family_provenance)
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
            provenance="pinned phase_estimates fixture",
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
            provenance="pinned coordinate fixture",
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
            provenance="pinned damping fixture",
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


def _scope(**overrides):
    labels = (
        AsymptoticSlowLabel(ell=ELL, a=A, sigma=-1),
        AsymptoticSlowLabel(ell=ELL, a=A, sigma=1),
    )
    kwargs = dict(
        scope_id="fixture-band-active-family",
        ell=ELL,
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        labels=labels,
        evidence_kind="formal-theorem",
        provenance="formal exact finite active-family scope fixture",
        finite_scope_certified=True,
        exact_active_signed_family_certified=True,
        carrier_slot_scope_certified=True,
        sign_duplication_certified=True,
    )
    kwargs.update(overrides)
    return LargeBandActiveFamilyScopeWitness(**kwargs)


def _coverage():
    minus = _bundle(sigma=-1)
    plus = _bundle(sigma=1)
    return LargeBandFamilyEstimateCoverage(
        scope=_scope(),
        phase_admissions=(minus[0], plus[0]),
        coordinate_admissions=(minus[1], plus[1]),
        damping_admissions=(minus[2], plus[2]),
    )


def test_exact_sign_pair_scope_exposes_only_finite_scope_uniform_envelopes():
    coverage = _coverage()
    bounds = coverage.uniform_scope_bounds()

    minus = _bundle(sigma=-1)
    plus = _bundle(sigma=1)
    assert bounds["phase_normal_error"] == max(minus[0].normal_error_bound, plus[0].normal_error_bound)
    assert bounds["phase_velocity_error"] == max(minus[0].phase_velocity_bound, plus[0].phase_velocity_bound)
    assert bounds["coordinate_error"] == max(
        minus[1].coordinate_error_A_bound,
        minus[1].coordinate_error_B_bound,
        minus[1].coordinate_error_C_bound,
        plus[1].coordinate_error_A_bound,
        plus[1].coordinate_error_B_bound,
        plus[1].coordinate_error_C_bound,
    )
    assert bounds["damping_error"] == max(
        minus[2].damping_theorem_envelope,
        plus[2].damping_theorem_envelope,
    )
    assert coverage.scope_coverage_complete is True
    assert coverage.covered_labels == _scope().expected_labels


def test_each_theorem_path_must_cover_exactly_the_declared_scope():
    minus = _bundle(sigma=-1)
    plus = _bundle(sigma=1)
    kwargs = dict(
        scope=_scope(),
        phase_admissions=(minus[0], plus[0]),
        coordinate_admissions=(minus[1], plus[1]),
        damping_admissions=(minus[2], plus[2]),
    )

    for field, incomplete in (
        ("phase_admissions", (minus[0],)),
        ("coordinate_admissions", (minus[1],)),
        ("damping_admissions", (minus[2],)),
    ):
        changed = dict(kwargs)
        changed[field] = incomplete
        with pytest.raises(ValueError, match="coverage must equal the declared active scope"):
            LargeBandFamilyEstimateCoverage(**changed)


def test_duplicates_and_cross_path_family_drift_fail_closed():
    minus = _bundle(sigma=-1)
    plus = _bundle(sigma=1)
    with pytest.raises(ValueError, match="duplicate phase admission"):
        LargeBandFamilyEstimateCoverage(
            scope=_scope(),
            phase_admissions=(minus[0], minus[0], plus[0]),
            coordinate_admissions=(minus[1], plus[1]),
            damping_admissions=(minus[2], plus[2]),
        )

    alternate_plus = _bundle(sigma=1, family_provenance="different theorem-path fixture")
    with pytest.raises(ValueError, match="coordinate path must reuse"):
        LargeBandFamilyEstimateCoverage(
            scope=_scope(),
            phase_admissions=(minus[0], plus[0]),
            coordinate_admissions=(minus[1], alternate_plus[1]),
            damping_admissions=(minus[2], plus[2]),
        )


def test_scope_itself_requires_theorem_evidence_and_complete_sign_pairs():
    lone = (AsymptoticSlowLabel(ell=ELL, a=A, sigma=1),)
    with pytest.raises(ValueError, match="closed under the Section 6 sign pair"):
        _scope(labels=lone)

    for kind in ("sampled", "fitted", "numeric-scan", ""):
        with pytest.raises(ValueError, match="evidence_kind"):
            _scope(evidence_kind=kind)

    for field in (
        "finite_scope_certified",
        "exact_active_signed_family_certified",
        "carrier_slot_scope_certified",
        "sign_duplication_certified",
    ):
        with pytest.raises(ValueError, match=field):
            _scope(**{field: False})


def test_success_does_not_promote_global_or_paper_exact_truth():
    coverage = _coverage()

    assert coverage.status == "formal-structure"
    assert coverage.exact_active_family_scope_machine_verified is False
    assert coverage.theorem_witnesses_machine_verified is False
    assert coverage.global_all_band_uniformity_verified is False
    assert coverage.uniform_eq_7_9_to_7_11_verified is False
    assert coverage.paper_exact_velocity_available is False
