import math
from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_endpoint_majorant_adapter import (
    OFFICIAL_SECTION10_Q_MAX,
    Section9UniformResidualEnvelopeWitness,
    derive_section9_endpoint_majorant_from_uniform_envelope,
)
from openai_ns_reconstruction.section9_stage_certificate import CertifiedBoundDatum


def _datum(value: float, provenance: str, *, kind: str = "paper-derived"):
    return CertifiedBoundDatum(
        upper_bound=value,
        kind=kind,
        provenance=provenance,
    )


def _witness(**overrides):
    kwargs = dict(
        endpoint_derivative_degree=2,
        stage=100,
        spatial_window=3,
        h=Fraction(1, 4),
        derivative_loss=Fraction(1, 1),
        log_power=2,
        leading_constant=_datum(2.0, "Eq. (9.18) C_jm theorem"),
        flat_remainder_uniform_bound=_datum(
            0.125, "uniform flat-remainder theorem on source revision"
        ),
        valid_q_upper=Fraction(1, 2),
        source_id="section9-eq-9.21-residual",
        source_revision="fixture-r1",
        residual_envelope_uniform_certified=True,
        flat_remainder_uniform_certified=True,
        official_support_covered_certified=True,
    )
    kwargs.update(overrides)
    return Section9UniformResidualEnvelopeWitness(**kwargs)


def test_adapter_derives_uniform_alpha_zero_majorant_on_official_support():
    witness = _witness()
    record = derive_section9_endpoint_majorant_from_uniform_envelope(witness)

    # sigma_100 = 1/5 + 10 = 51/5 and h=1/4, so beta=51/20-1=31/20.
    assert record.residual_exponent == Fraction(31, 20)
    assert record.official_q_max == Fraction(5, 16)
    assert record.source_key == ("section9-eq-9.21-residual", "fixture-r1")
    assert record.section9_derivative_order == 3
    assert record.bridge.majorant.singularity_exponent == 0.0
    assert record.bridge.majorant.valid_from == 0.75
    assert record.bridge.localization_transfer.certified
    assert record.formal_adapter_ready

    expected = (
        witness.leading_constant.upper_bound * record.power_log_supremum
        + witness.flat_remainder_uniform_bound.upper_bound
    )
    assert record.derived_uniform_coefficient == pytest.approx(expected)
    assert record.bridge.majorant.coefficient == pytest.approx(expected)
    assert "q<=5/16" in record.bridge.evidence_provenance

    assert not record.actual_section9_sequence_verified
    assert not record.source_majorant_derived_from_actual_residual_verified
    assert not record.endpoint_limit_constructed
    assert not record.paper_exact_velocity_available


def test_power_log_supremum_dominates_dense_independent_grid():
    record = derive_section9_endpoint_majorant_from_uniform_envelope(_witness())
    beta = float(record.residual_exponent)
    p = 2
    q_max = float(OFFICIAL_SECTION10_Q_MAX)

    # Independent regression only: the production adapter uses the analytic
    # stationary point, not this sampled grid.
    for k in range(1, 4000):
        q = q_max * k / 4000.0
        value = q**beta * (1.0 + abs(math.log(q))) ** p
        assert value <= record.power_log_supremum * (1.0 + 2e-14)


def test_zero_exponent_is_admitted_only_without_log_growth():
    stage = 0
    h = Fraction(1, 4)
    # h*sigma_0 = (1/4)*(1/5) = 1/20.
    witness = _witness(
        stage=stage,
        h=h,
        derivative_loss=Fraction(1, 20),
        log_power=0,
    )
    record = derive_section9_endpoint_majorant_from_uniform_envelope(witness)
    assert record.residual_exponent == 0
    assert record.power_log_supremum == 1.0

    with pytest.raises(ValueError, match="beta=0"):
        derive_section9_endpoint_majorant_from_uniform_envelope(
            _witness(
                stage=stage,
                h=h,
                derivative_loss=Fraction(1, 20),
                log_power=1,
            )
        )


def test_adapter_rejects_negative_exponent_or_insufficient_q_domain():
    with pytest.raises(ValueError, match="leading exponent"):
        derive_section9_endpoint_majorant_from_uniform_envelope(
            _witness(stage=0, derivative_loss=1, log_power=0)
        )

    with pytest.raises(ValueError, match="5/16"):
        _witness(valid_q_upper=Fraction(5, 16))


def test_adapter_rejects_mixed_evidence_or_missing_uniform_theorem_fact():
    with pytest.raises(ValueError, match="same evidence kind"):
        _witness(
            flat_remainder_uniform_bound=_datum(
                0.125,
                "rigorous external remainder",
                kind="rigorous-external",
            )
        )

    with pytest.raises(ValueError, match="residual_envelope_uniform_certified"):
        _witness(residual_envelope_uniform_certified=False)

    with pytest.raises(ValueError, match="source_revision"):
        _witness(source_revision="  ")
