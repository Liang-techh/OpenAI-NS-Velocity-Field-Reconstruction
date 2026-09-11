import pytest

from openai_ns_reconstruction.section9_flat_remainder import (
    Section9FlatRemainderPowerWitness,
    certify_flat_remainder_power_ladder,
)
from openai_ns_reconstruction.section9_stage_certificate import CertifiedBoundDatum


def _witness(power, *, stage=12, derivative_order=3, q_upper=0.125):
    return Section9FlatRemainderPowerWitness(
        stage=stage,
        derivative_order=derivative_order,
        power=power,
        q_upper=q_upper,
        constant_bound=CertifiedBoundDatum(
            upper_bound=2.0 + power,
            kind="paper-derived",
            provenance=f"Lemma 9.8 flat-remainder bound for N={power}",
        ),
    )


def test_finite_flat_remainder_ladder_preserves_exact_power_family_and_fail_closed_flags():
    # Supply the rows deliberately out of order.  The production gate may
    # normalize their order, but it must neither invent nor drop a power.
    supplied = (_witness(3), _witness(0), _witness(2), _witness(1))
    cert = certify_flat_remainder_power_ladder(supplied, max_power=3)

    assert cert.powers == (0, 1, 2, 3)
    assert tuple(w.constant_bound.upper_bound for w in cert.witnesses) == (2.0, 3.0, 4.0, 5.0)
    assert cert.stage == 12
    assert cert.derivative_order == 3
    assert cert.q_upper == 0.125
    assert cert.formal_ladder_ready
    assert cert.contiguous_finite_power_ladder_verified
    assert cert.common_q_domain_verified
    assert not cert.source_theorems_machine_verified
    assert not cert.flat_remainder_all_orders_verified
    assert not cert.actual_section9_sequence_verified
    assert not cert.paper_exact_velocity_available


def test_missing_duplicate_or_extra_power_fails_closed():
    with pytest.raises(ValueError, match="exactly powers"):
        certify_flat_remainder_power_ladder((_witness(0), _witness(2)), max_power=2)

    with pytest.raises(ValueError, match="duplicate"):
        certify_flat_remainder_power_ladder((_witness(0), _witness(1), _witness(1)), max_power=1)

    with pytest.raises(ValueError, match="exactly powers"):
        certify_flat_remainder_power_ladder((_witness(0), _witness(1), _witness(2)), max_power=1)


def test_mixing_stage_derivative_order_or_q_domain_fails_closed():
    with pytest.raises(ValueError, match="one Section 9 stage"):
        certify_flat_remainder_power_ladder((_witness(0), _witness(1, stage=13)), max_power=1)

    with pytest.raises(ValueError, match="one derivative order"):
        certify_flat_remainder_power_ladder(
            (_witness(0), _witness(1, derivative_order=4)), max_power=1
        )

    with pytest.raises(ValueError, match="one common q domain"):
        certify_flat_remainder_power_ladder(
            (_witness(0), _witness(1, q_upper=0.0625)), max_power=1
        )


def test_gate_rejects_unprovenanced_or_sampled_inputs_and_invalid_domains():
    with pytest.raises(ValueError, match="kind must be one of"):
        CertifiedBoundDatum(1.0, "sampled", "grid maximum")

    with pytest.raises(ValueError, match="provenance"):
        CertifiedBoundDatum(1.0, "paper-derived", "   ")

    with pytest.raises(ValueError, match="q_upper"):
        _witness(0, q_upper=1.0)

    with pytest.raises(ValueError, match="nonnegative integer"):
        Section9FlatRemainderPowerWitness(
            stage=0,
            derivative_order=0,
            power=True,
            q_upper=0.1,
            constant_bound=CertifiedBoundDatum(1.0, "lean-derived", "proof artifact"),
        )
