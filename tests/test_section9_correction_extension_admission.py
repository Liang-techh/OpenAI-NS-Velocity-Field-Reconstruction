from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_correction_extension_admission import (
    Section9CorrectionExtensionWitness,
    admit_section9_correction_stage,
)
from openai_ns_reconstruction.section9_local_sum_admission import Section9LocalSumHypothesis


def _hypothesis(**kwargs):
    data = dict(
        q_big=Fraction(1, 4),
        first_scale=8,
        kind="paper-derived",
        provenance="Proposition 9.9 + Lemma 5.4 infinite schedule",
    )
    data.update(kwargs)
    return Section9LocalSumHypothesis(**data)


def _witness(component, **kwargs):
    data = dict(
        stage=3,
        component=component,
        q_big=Fraction(1, 4),
        scale=40,
        kind="paper-derived",
        provenance=f"Eq. (9.21) stage 3 {component} smooth-extension theorem",
        common_domain_verified=True,
        schedule_membership_verified=True,
        cutoff_support_verified=True,
        smooth_zero_extension_verified=True,
    )
    data.update(kwargs)
    return Section9CorrectionExtensionWitness(**data)


def test_admits_complete_stage_and_checks_exact_schedule_lower_bound():
    cert = admit_section9_correction_stage(
        _hypothesis(), [_witness("p"), _witness("A"), _witness("B")]
    )

    # Independent closed-form oracle: a_3 >= a_1 * 2^(3-1) = 8*4 = 32.
    assert cert.stage == 3
    assert cert.scale == 40
    assert cert.schedule_lower_bound == 32
    assert cert.components == ("A", "B", "p")
    assert cert.formal_admission_ready
    assert cert.common_domain_verified
    assert cert.schedule_membership_witnessed
    assert cert.cutoff_support_verified
    assert cert.smooth_zero_extension_verified
    assert cert.component_completeness_verified

    # Formal admission is intentionally not a construction or convergence proof.
    assert not cert.source_theorems_machine_verified
    assert not cert.actual_correction_field_values_verified
    assert not cert.infinite_schedule_constructed_here
    assert not cert.eq_9_21_sum_constructed
    assert not cert.proposition_9_9_verified
    assert not cert.paper_exact_velocity_available


def test_rejects_scale_below_doubling_lower_bound_even_with_membership_assertion():
    rows = [_witness(name, scale=31) for name in ("A", "B", "p")]
    with pytest.raises(ValueError, match="doubling-schedule lower bound"):
        admit_section9_correction_stage(_hypothesis(), rows)


def test_rejects_missing_duplicate_or_incoherent_component_data():
    with pytest.raises(ValueError, match="cover exactly A, B and p"):
        admit_section9_correction_stage(_hypothesis(), [_witness("A"), _witness("B")])

    with pytest.raises(ValueError, match="duplicate Section 9 component"):
        admit_section9_correction_stage(
            _hypothesis(), [_witness("A"), _witness("A"), _witness("B"), _witness("p")]
        )

    with pytest.raises(ValueError, match="one stage"):
        admit_section9_correction_stage(
            _hypothesis(), [_witness("A"), _witness("B", stage=4), _witness("p")]
        )

    with pytest.raises(ValueError, match="same stage cutoff scale"):
        admit_section9_correction_stage(
            _hypothesis(), [_witness("A"), _witness("B", scale=48), _witness("p")]
        )

    with pytest.raises(ValueError, match="common q_big"):
        admit_section9_correction_stage(
            _hypothesis(),
            [_witness("A"), _witness("B", q_big=Fraction(1, 5)), _witness("p")],
        )


def test_rejects_unverified_extension_facts_and_nonrigorous_evidence():
    for field_name in (
        "common_domain_verified",
        "schedule_membership_verified",
        "cutoff_support_verified",
        "smooth_zero_extension_verified",
    ):
        with pytest.raises(ValueError, match=field_name):
            _witness("A", **{field_name: False})

    with pytest.raises(ValueError, match="kind must be"):
        _witness("A", kind="sampled")
    with pytest.raises(ValueError, match="provenance must be nonempty"):
        _witness("A", provenance="   ")


def test_rejects_witnesses_from_a_different_global_common_domain():
    hypothesis = _hypothesis(q_big=Fraction(1, 3))
    with pytest.raises(ValueError, match="local-sum hypothesis common q_big"):
        admit_section9_correction_stage(
            hypothesis, [_witness("A"), _witness("B"), _witness("p")]
        )


def test_component_and_stage_inputs_fail_closed():
    with pytest.raises(ValueError, match="component must be"):
        _witness("velocity")
    with pytest.raises(ValueError, match="stage must be a positive integer"):
        _witness("A", stage=0)
    with pytest.raises(ValueError, match="scale must be a positive integer"):
        _witness("A", scale=0)
