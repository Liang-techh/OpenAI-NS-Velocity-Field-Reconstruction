from fractions import Fraction
import math

import pytest

from openai_ns_reconstruction.phase_frame_bounds import base_phase_constant
from openai_ns_reconstruction.phase_large_band_base_source import (
    LargeBandBaseSourceBinding,
    LargeBandBaseSourceWitness,
)
from openai_ns_reconstruction.phase_large_band_damping import (
    PINNED_DAMPING_ERROR_THEOREM,
    LargeBandDampingErrorAdmission,
    LargeBandDampingErrorWitness,
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
)
from openai_ns_reconstruction.phase_large_band_scale import LargeBandPhaseScaleCertificate


ELL = 100000
A = (7, -3, 11)
M = 2.0
SOURCE_ID = "formal-base-provider"
SOURCE_REVISION = "fixture-r1"


def _scale():
    return LargeBandPhaseScaleCertificate(ell=ELL, h=Fraction(1, 200), M=M)


def _family(*, sigma=1):
    scale = _scale()
    local = LargeBandLocalBaseWitness(
        label=AsymptoticSlowLabel(ell=ELL, a=A, sigma=sigma),
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
    frame = LargeBandPhaseFrameCertificate(
        scale=scale,
        u=1.0,
        B=1.0,
        viscosity=1.5,
    )
    witness = LargeBandFamilyInputWitness(
        label=AsymptoticSlowLabel(ell=ELL, a=A, sigma=sigma),
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
    return LargeBandFamilyInputAdmission(binding=binding, frame=frame, witness=witness)


def _witness(**overrides):
    kwargs = dict(
        label=AsymptoticSlowLabel(ell=ELL, a=A, sigma=1),
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        lean_repository=PINNED_LEAN_REPOSITORY,
        lean_commit=PINNED_LEAN_COMMIT,
        theorem_name=PINNED_DAMPING_ERROR_THEOREM,
        evidence_kind="formal-theorem",
        provenance="formal damping_error theorem fixture",
        damping_error_application_certified=True,
        carrier_membership_certified=True,
        actual_damping_term_identity_certified=True,
        reference_damping_term_identity_certified=True,
        viscosity_damping_family_identity_certified=True,
        reference_viscosity_identity_certified=True,
        reference_radius_center_identity_certified=True,
        scale_identity_certified=True,
        family_input_dependency_certified=True,
    )
    kwargs.update(overrides)
    return LargeBandDampingErrorWitness(**kwargs)


def test_admission_exposes_exact_named_constant_over_large_band_scale():
    family = _family()
    admission = LargeBandDampingErrorAdmission(family=family, witness=_witness())

    # Independent closed-form reconstruction of
    # dampingConstant(M) = 4*M*(6*M+5)*phaseConstant(M).
    expected = 4.0 * M * (6.0 * M + 5.0) * base_phase_constant(M) / float(ELL * ELL)

    assert math.isclose(admission.damping_theorem_envelope, expected, rel_tol=1e-12)
    assert math.isclose(admission.constant_over_scale_envelope, expected, rel_tol=1e-12)
    assert all(admission.admission_checks().values())


def test_admission_rejects_signed_label_or_source_revision_drift():
    family = _family()

    wrong_sign = _witness(label=AsymptoticSlowLabel(ell=ELL, a=A, sigma=-1))
    with pytest.raises(ValueError, match="signed_label_identity"):
        LargeBandDampingErrorAdmission(family=family, witness=wrong_sign)

    wrong_source = _witness(source_revision="fixture-r2")
    with pytest.raises(ValueError, match="source_revision_identity"):
        LargeBandDampingErrorAdmission(family=family, witness=wrong_source)


def test_witness_pins_exact_formal_source_and_theorem():
    with pytest.raises(ValueError, match="lean_repository"):
        _witness(lean_repository="example/wrong")
    with pytest.raises(ValueError, match="lean_commit"):
        _witness(lean_commit="0" * 40)
    with pytest.raises(ValueError, match="theorem_name"):
        _witness(theorem_name="BasePhaseGeometry.FamilyData.coordinate_errors")


def test_witness_rejects_sampled_evidence_and_missing_damping_identities():
    for kind in ("sampled", "fitted", "numeric-scan", ""):
        with pytest.raises(ValueError, match="evidence_kind"):
            _witness(evidence_kind=kind)

    theorem_flags = (
        "damping_error_application_certified",
        "carrier_membership_certified",
        "actual_damping_term_identity_certified",
        "reference_damping_term_identity_certified",
        "viscosity_damping_family_identity_certified",
        "reference_viscosity_identity_certified",
        "reference_radius_center_identity_certified",
        "scale_identity_certified",
        "family_input_dependency_certified",
    )
    for field in theorem_flags:
        with pytest.raises(ValueError, match=field):
            _witness(**{field: False})


def test_either_sign_can_be_admitted_only_with_same_signed_family_identity():
    family = _family(sigma=-1)
    admitted = LargeBandDampingErrorAdmission(
        family=family,
        witness=_witness(label=AsymptoticSlowLabel(ell=ELL, a=A, sigma=-1)),
    )
    assert admitted.label.sigma == -1
    assert admitted.source_key == family.source_key


def test_successful_admission_does_not_upgrade_reconstruction_truth():
    admission = LargeBandDampingErrorAdmission(family=_family(), witness=_witness())

    assert admission.status == "formal-structure"
    assert admission.actual_damping_field_materialized is False
    assert admission.damping_error_machine_verified is False
    assert admission.all_active_family_damping_verified is False
    assert admission.uniform_eq_7_9_to_7_11_verified is False
    assert admission.paper_exact_velocity_available is False
