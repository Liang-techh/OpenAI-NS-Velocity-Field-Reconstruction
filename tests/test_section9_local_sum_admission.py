from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_local_sum_admission import (
    Section9LocalSumHypothesis,
    certify_eq_9_21_local_finiteness,
)


def _hypothesis(**kwargs):
    data = dict(
        q_big=Fraction(1, 4),
        first_scale=8,
        kind="paper-derived",
        provenance="Lemma 9.7 + Lemma 5.4 infinite doubling schedule",
    )
    data.update(kwargs)
    return Section9LocalSumHypothesis(**data)


def _oracle_active_stages(first_scale, q_lower, limit=20):
    # Independent direct enumeration of the worst-case doubling sequence.
    scales = [first_scale]
    for _ in range(1, limit):
        scales.append(2 * scales[-1])
    return [j for j, scale in enumerate(scales, start=1) if scale * q_lower < 1]


def test_exact_local_finiteness_bound_matches_direct_doubling_oracle():
    lower = Fraction(1, 40)
    cert = certify_eq_9_21_local_finiteness(
        _hypothesis(), q_lower=lower, q_upper=Fraction(1, 5)
    )
    assert _oracle_active_stages(8, lower) == [1, 2, 3]
    assert cert.max_potentially_active_stage == 3
    assert cert.first_guaranteed_inactive_stage == 4
    assert cert.common_domain_arithmetic_verified
    assert cert.local_finiteness_from_doubling_verified
    assert cert.paper_cutoff_support_rule_used
    assert not cert.source_theorems_machine_verified
    assert not cert.actual_correction_fields_verified
    assert not cert.eq_9_21_sum_constructed
    assert not cert.proposition_9_9_verified
    assert not cert.paper_exact_velocity_available


def test_cutoff_support_edge_is_strict_and_exact():
    # Here a_3*q_lower = 32/32 = 1 exactly, so chi(a_3 q)=0 by the
    # paper cutoff support rule. Only stages 1 and 2 can remain active.
    lower = Fraction(1, 32)
    cert = certify_eq_9_21_local_finiteness(
        _hypothesis(), q_lower=lower, q_upper=Fraction(1, 6)
    )
    assert _oracle_active_stages(8, lower) == [1, 2]
    assert cert.max_potentially_active_stage == 2
    assert cert.first_guaranteed_inactive_stage == 3


def test_strip_beyond_first_cutoff_has_no_positive_active_stage():
    hypothesis = Section9LocalSumHypothesis(
        q_big=Fraction(3, 4),
        first_scale=2,
        kind="lean-derived",
        provenance="machine theorem reference",
    )
    cert = certify_eq_9_21_local_finiteness(
        hypothesis, q_lower=Fraction(1, 2), q_upper=Fraction(2, 3)
    )
    assert cert.max_potentially_active_stage == 0
    assert cert.first_guaranteed_inactive_stage == 1


def test_common_domain_and_compact_strip_fail_closed():
    with pytest.raises(ValueError, match="q_lower < q_upper"):
        certify_eq_9_21_local_finiteness(
            _hypothesis(), q_lower=Fraction(1, 8), q_upper=Fraction(1, 8)
        )
    with pytest.raises(ValueError, match="q_upper < q_big"):
        certify_eq_9_21_local_finiteness(
            _hypothesis(), q_lower=Fraction(1, 8), q_upper=Fraction(1, 4)
        )
    with pytest.raises(ValueError, match="positive"):
        certify_eq_9_21_local_finiteness(
            _hypothesis(), q_lower=0, q_upper=Fraction(1, 8)
        )


def test_hypothesis_rejects_unrigorous_metadata_and_invalid_domain_data():
    with pytest.raises(ValueError, match="strict condition"):
        _hypothesis(q_big=Fraction(1, 8))
    with pytest.raises(ValueError, match="q_big<=1"):
        _hypothesis(q_big=Fraction(5, 4))
    with pytest.raises(ValueError, match="positive integer"):
        _hypothesis(first_scale=True)
    with pytest.raises(ValueError, match="kind must be one of"):
        _hypothesis(kind="sampled")
    with pytest.raises(ValueError, match="provenance"):
        _hypothesis(provenance="   ")
