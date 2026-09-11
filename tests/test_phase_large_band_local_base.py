from fractions import Fraction
import math

import pytest

from openai_ns_reconstruction.phase_large_band_local_base import (
    AsymptoticSlowLabel,
    LargeBandLocalBaseAdmission,
    LargeBandLocalBaseWitness,
)
from openai_ns_reconstruction.phase_large_band_scale import LargeBandPhaseScaleCertificate


def _witness(*, ell=100000, M=2.0, evidence_kind="analytic-theorem", provenance="Lemma-localbase"):
    return LargeBandLocalBaseWitness(
        label=AsymptoticSlowLabel(ell=ell, a=(7, -3, 11), sigma=1),
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
        evidence_kind=evidence_kind,
        provenance=provenance,
        convex_domain_certified=True,
        differentiability_certified=True,
        representative_membership_certified=True,
        local_base_suprema_certified=True,
    )


def test_large_band_admission_never_forms_underflowed_Q_or_epsilon():
    ell = 100000
    h = Fraction(1, 200)
    scale = LargeBandPhaseScaleCertificate(ell=ell, h=h, M=2.0)
    admission = LargeBandLocalBaseAdmission(scale=scale, witness=_witness(ell=ell))

    # Independent closed-form checks.  Q itself has long since underflowed in
    # binary64, but the theorem-facing band and epsilon^2 exponent remain exact.
    assert math.ldexp(1.0, -ell) == 0.0
    assert admission.S_star == ell**2
    assert admission.max_slow_box_diameter == Fraction(1, ell**6)
    assert admission.epsilon_squared_log2 == Fraction(-1000, 1)
    assert all(admission.admission_checks().values())


def test_normalized_localbase_ratios_match_pinned_threshold_form():
    witness = _witness()
    assert witness.second_F0_over_M == pytest.approx(0.8)
    assert witness.value_F_error_over_M_epsilon2 == pytest.approx(0.5)
    assert witness.first_G_error_over_M_epsilon2 == pytest.approx(0.3)
    assert witness.diameter_times_S3 == pytest.approx(0.75)
    assert witness.quantitative_hypotheses_admitted is True
    assert witness.qualitative_hypotheses_admitted is True


def test_asymptotic_label_keeps_sign_pair_without_binary64_band_cap():
    label = AsymptoticSlowLabel(ell=250000, a=(1, 2, 3), sigma=-1)
    opposite = label.opposite()
    assert label.box_key == opposite.box_key == (250000, (1, 2, 3))
    assert opposite.sigma == 1


def test_witness_rejects_sampled_fitted_or_missing_theorem_evidence():
    for kind in ("sampled", "fitted", "numeric-scan", ""):
        with pytest.raises(ValueError, match="evidence_kind"):
            _witness(evidence_kind=kind)

    with pytest.raises(ValueError, match="provenance"):
        _witness(provenance="   ")


def test_witness_fails_closed_on_threshold_or_qualitative_hypothesis_failure():
    base = dict(
        label=AsymptoticSlowLabel(ell=100000, a=(0, 0, 0), sigma=1),
        M=2.0,
        diameter_times_S3=1.0,
        second_F0_over_M=1.0,
        second_G0_over_M=1.0,
        first_F0_over_M=1.0,
        value_F_error_over_M_epsilon2=1.0,
        first_F_error_over_M_epsilon2=1.0,
        first_G_error_over_M_epsilon2=1.0,
        radial_F_over_M=1.0,
        axial_F_over_M=1.0,
        axial_G_over_M=1.0,
        evidence_kind="formal-theorem",
        provenance="pinned LocalBase theorem witness",
        convex_domain_certified=True,
        differentiability_certified=True,
        representative_membership_certified=True,
        local_base_suprema_certified=True,
    )

    too_large = dict(base)
    too_large["first_F_error_over_M_epsilon2"] = 1.0001
    with pytest.raises(ValueError, match=r"\[0,1\]"):
        LargeBandLocalBaseWitness(**too_large)

    too_wide = dict(base)
    too_wide["diameter_times_S3"] = 1.01
    with pytest.raises(ValueError, match=r"\[0,1\]"):
        LargeBandLocalBaseWitness(**too_wide)

    not_convex = dict(base)
    not_convex["convex_domain_certified"] = False
    with pytest.raises(ValueError, match="convex_domain_certified"):
        LargeBandLocalBaseWitness(**not_convex)


def test_admission_requires_same_band_and_same_M():
    scale = LargeBandPhaseScaleCertificate(ell=100000, h=Fraction(1, 200), M=2.0)

    with pytest.raises(ValueError, match="label band"):
        LargeBandLocalBaseAdmission(scale=scale, witness=_witness(ell=100001))

    with pytest.raises(ValueError, match="LocalBase M"):
        LargeBandLocalBaseAdmission(scale=scale, witness=_witness(M=3.0))


def test_admission_remains_formal_structure_and_does_not_upgrade_truth_status():
    scale = LargeBandPhaseScaleCertificate(ell=100000, h=0.005, M=2.0)
    admission = LargeBandLocalBaseAdmission(scale=scale, witness=_witness())

    assert admission.status == "formal-structure"
    assert admission.actual_base_fields_verified is False
    assert admission.uniform_eq_7_9_to_7_11_verified is False
    assert admission.paper_exact_velocity_available is False
