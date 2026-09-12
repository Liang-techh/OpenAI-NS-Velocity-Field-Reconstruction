from fractions import Fraction

import pytest

from openai_ns_reconstruction.phase_large_band_base_source import (
    LargeBandBaseSourceBinding,
    LargeBandBaseSourceWitness,
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
from openai_ns_reconstruction.phase_large_band_scale import LargeBandPhaseScaleCertificate


ELL = 100000
A = (7, -3, 11)
M = 2.0
SOURCE_ID = "formal-base-provider"
SOURCE_REVISION = "fixture-r1"


def _scale():
    return LargeBandPhaseScaleCertificate(ell=ELL, h=Fraction(1, 200), M=M)


def _binding(*, sigma=1):
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
    admission = LargeBandLocalBaseAdmission(scale=scale, witness=local)
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
    return LargeBandBaseSourceBinding(admission=admission, source=source)


def _frame():
    return LargeBandPhaseFrameCertificate(
        scale=_scale(),
        u=1.0,
        B=1.0,
        viscosity=1.5,
    )


def _witness(**overrides):
    kwargs = dict(
        label=AsymptoticSlowLabel(ell=ELL, a=A, sigma=1),
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
    kwargs.update(overrides)
    return LargeBandFamilyInputWitness(**kwargs)


def test_admission_binds_signed_family_index_to_same_source_and_scalar_geometry():
    binding = _binding(sigma=-1)
    frame = _frame()
    witness = _witness()

    admission = LargeBandFamilyInputAdmission(binding=binding, frame=frame, witness=witness)

    assert admission.label == binding.signed_label(1)
    assert admission.source_key == binding.source_key
    assert admission.family_phase_bound == frame.family_normal_delta
    assert admission.family_phase_bound > 0.0
    assert all(admission.admission_checks().values())


def test_either_sign_can_be_admitted_but_source_revision_cannot_change():
    binding = _binding(sigma=1)
    minus = _witness(label=AsymptoticSlowLabel(ell=ELL, a=A, sigma=-1))
    admitted = LargeBandFamilyInputAdmission(binding=binding, frame=_frame(), witness=minus)
    assert admitted.label == binding.signed_label(-1)

    changed_source = _witness(source_revision="fixture-r2")
    with pytest.raises(ValueError, match="source_revision_identity"):
        LargeBandFamilyInputAdmission(binding=binding, frame=_frame(), witness=changed_source)


def test_admission_rejects_label_M_or_scalar_parameter_mismatch():
    binding = _binding()
    frame = _frame()

    wrong_box = _witness(label=AsymptoticSlowLabel(ell=ELL, a=(7, -3, 12), sigma=1))
    with pytest.raises(ValueError, match="signed_label_identity"):
        LargeBandFamilyInputAdmission(binding=binding, frame=frame, witness=wrong_box)

    wrong_M = _witness(M=3.0)
    with pytest.raises(ValueError, match="M_identity"):
        LargeBandFamilyInputAdmission(binding=binding, frame=frame, witness=wrong_M)

    wrong_u = _witness(u=0.75)
    with pytest.raises(ValueError, match="u_identity"):
        LargeBandFamilyInputAdmission(binding=binding, frame=frame, witness=wrong_u)

    wrong_B = _witness(B=0.75)
    with pytest.raises(ValueError, match="B_identity"):
        LargeBandFamilyInputAdmission(binding=binding, frame=frame, witness=wrong_B)

    wrong_viscosity = _witness(viscosity=2.0)
    with pytest.raises(ValueError, match="viscosity_identity"):
        LargeBandFamilyInputAdmission(binding=binding, frame=frame, witness=wrong_viscosity)


def test_witness_checks_family_global_scalar_hypotheses_directly():
    with pytest.raises(ValueError, match="0<u<=M"):
        _witness(u=0.0)
    with pytest.raises(ValueError, match="0<B<=M"):
        _witness(B=0.0)
    with pytest.raises(ValueError, match="\[0,4\]"):
        _witness(viscosity=4.1)
    with pytest.raises(ValueError, match="reference_radius"):
        _witness(reference_radius=0.0)
    with pytest.raises(ValueError, match="1/\(2\*r0\)<=M"):
        _witness(reference_radius=0.2)
    with pytest.raises(ValueError, match="4\*r0\*Tg<=M"):
        _witness(Tg=0.6)


def test_witness_rejects_sampled_evidence_and_missing_theorem_facts():
    for kind in ("sampled", "fitted", "numeric-scan", ""):
        with pytest.raises(ValueError, match="evidence_kind"):
            _witness(evidence_kind=kind)

    for field in ("source_id", "source_revision", "provenance"):
        with pytest.raises(ValueError, match="nonempty"):
            _witness(**{field: "   "})

    theorem_flags = (
        "family_index_source_identity_certified",
        "carrier_membership_certified",
        "slot_membership_certified",
        "representative_data_identity_certified",
        "reference_scale_identity_certified",
        "viscosity_identity_certified",
    )
    for field in theorem_flags:
        with pytest.raises(ValueError, match=field):
            _witness(**{field: False})


def test_successful_admission_does_not_upgrade_reconstruction_truth():
    admission = LargeBandFamilyInputAdmission(
        binding=_binding(),
        frame=_frame(),
        witness=_witness(),
    )

    assert admission.status == "formal-structure"
    assert admission.actual_partition_verified is False
    assert admission.actual_base_fields_verified is False
    assert admission.family_phase_estimates_machine_verified is False
    assert admission.uniform_eq_7_9_to_7_11_verified is False
    assert admission.paper_exact_velocity_available is False
