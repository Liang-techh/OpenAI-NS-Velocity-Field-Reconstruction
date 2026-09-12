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


def _phase(*, sigma=1):
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
        provenance="formal FamilyData input fixture",
        family_index_source_identity_certified=True,
        carrier_membership_certified=True,
        slot_membership_certified=True,
        representative_data_identity_certified=True,
        reference_scale_identity_certified=True,
        viscosity_identity_certified=True,
    )
    family = LargeBandFamilyInputAdmission(binding=binding, frame=frame, witness=family_witness)
    phase_witness = LargeBandPhaseEstimatesWitness(
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
    )
    return LargeBandPhaseEstimatesAdmission(family=family, witness=phase_witness)


def _witness(*, sigma=1, **overrides):
    kwargs = dict(
        label=AsymptoticSlowLabel(ell=ELL, a=A, sigma=sigma),
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        lean_repository=PINNED_LEAN_REPOSITORY,
        lean_commit=PINNED_LEAN_COMMIT,
        base_estimates_theorem_name=PINNED_BASE_ESTIMATES_THEOREM,
        coordinate_errors_theorem_name=PINNED_COORDINATE_ERRORS_THEOREM,
        evidence_kind="formal-theorem",
        provenance="pinned base_estimates -> coordinate_errors fixture",
        base_estimates_application_certified=True,
        actual_base_value_identity_certified=True,
        actual_shear_vector_identity_certified=True,
        uniform_carrier_base_conclusion_certified=True,
        coordinate_errors_application_certified=True,
        actual_frame_errors_identity_certified=True,
        phase_estimates_dependency_certified=True,
        base_estimates_dependency_certified=True,
        uniform_carrier_slot_conclusion_certified=True,
    )
    kwargs.update(overrides)
    return LargeBandCoordinateErrorsWitness(**kwargs)


def test_admission_exposes_the_two_pinned_theorem_envelopes():
    phase = _phase()
    admission = LargeBandCoordinateErrorsAdmission(phase=phase, witness=_witness())

    # Independent reconstruction of base_estimates' literal 16*M^2/D.scale(i)
    # with the admitted large-band identity D.scale(i)=S_*=ell^2.
    expected_base = 16.0 * M * M / float(ELL * ELL)
    assert admission.base_value_error_bound == pytest.approx(expected_base)
    assert admission.base_shear_vector_error_bound == pytest.approx(expected_base)

    # coordinate_errors has one common coordinateConstant(M,u)/D.scale(i)
    # envelope for errorA/errorB/errorC.  The scalar formula itself is already
    # independently regression-tested in phase_frame_bounds; this bridge must
    # inherit exactly that admitted envelope rather than fit three tolerances.
    expected_coordinate = phase.family.frame.frame_error_bound
    assert admission.coordinate_error_A_bound == expected_coordinate
    assert admission.coordinate_error_B_bound == expected_coordinate
    assert admission.coordinate_error_C_bound == expected_coordinate
    assert admission.coordinate_error_A_bound > 0.0
    assert all(admission.admission_checks().values())


def test_signed_family_and_source_revision_must_match_phase_admission():
    with pytest.raises(ValueError, match="signed_label_identity"):
        LargeBandCoordinateErrorsAdmission(phase=_phase(sigma=1), witness=_witness(sigma=-1))

    changed = _witness(source_revision="fixture-r2")
    with pytest.raises(ValueError, match="source_revision_identity"):
        LargeBandCoordinateErrorsAdmission(phase=_phase(), witness=changed)


def test_witness_is_pinned_to_exact_formal_source_and_both_theorems():
    with pytest.raises(ValueError, match="lean_repository"):
        _witness(lean_repository="example/other")
    with pytest.raises(ValueError, match="lean_commit"):
        _witness(lean_commit="deadbeef")
    with pytest.raises(ValueError, match="base_estimates_theorem_name"):
        _witness(base_estimates_theorem_name=PINNED_COORDINATE_ERRORS_THEOREM)
    with pytest.raises(ValueError, match="coordinate_errors_theorem_name"):
        _witness(coordinate_errors_theorem_name=PINNED_BASE_ESTIMATES_THEOREM)


def test_sampled_or_incomplete_theorem_claims_fail_closed():
    for kind in ("sampled", "fitted", "numeric-scan", ""):
        with pytest.raises(ValueError, match="evidence_kind"):
            _witness(evidence_kind=kind)

    for field in ("source_id", "source_revision", "provenance"):
        with pytest.raises(ValueError, match="nonempty"):
            _witness(**{field: "   "})

    flags = (
        "base_estimates_application_certified",
        "actual_base_value_identity_certified",
        "actual_shear_vector_identity_certified",
        "uniform_carrier_base_conclusion_certified",
        "coordinate_errors_application_certified",
        "actual_frame_errors_identity_certified",
        "phase_estimates_dependency_certified",
        "base_estimates_dependency_certified",
        "uniform_carrier_slot_conclusion_certified",
    )
    for field in flags:
        with pytest.raises(ValueError, match=field):
            _witness(**{field: False})


def test_success_does_not_promote_field_or_reconstruction_truth():
    admission = LargeBandCoordinateErrorsAdmission(phase=_phase(), witness=_witness())

    assert admission.status == "formal-structure"
    assert admission.actual_base_values_materialized is False
    assert admission.actual_shear_vectors_materialized is False
    assert admission.actual_coordinate_errors_materialized is False
    assert admission.base_estimates_machine_verified is False
    assert admission.coordinate_errors_machine_verified is False
    assert admission.uniform_eq_7_9_to_7_11_verified is False
    assert admission.paper_exact_velocity_available is False
