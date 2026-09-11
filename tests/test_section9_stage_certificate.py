from fractions import Fraction
import math

import pytest

from openai_ns_reconstruction.section9_stage_certificate import (
    CertifiedBoundDatum,
    check_eq_9_17_pointwise,
    check_eq_9_18_pointwise,
)


def _certified(
    value: float, provenance: str = "independent test enclosure"
) -> CertifiedBoundDatum:
    return CertifiedBoundDatum(value, "certified-numerical", provenance)


def test_eq_917_pointwise_check_uses_paper_exponent_and_log_factor():
    # Independent paper substitution, not a call back into the production
    # exponent helpers: h=1/200, m=2, j=10065 gives g_j-ell_m=1.
    h = Fraction(1, 200)
    q = 0.25
    constant = 3.0
    log_power = 2
    exponent = Fraction(1)
    oracle_rhs = constant * q ** float(exponent) * (1.0 + abs(math.log(q))) ** log_power

    passed = check_eq_9_17_pointwise(
        stage=10065,
        derivative_order=2,
        h=h,
        q=q,
        constant=constant,
        log_power=log_power,
        correction_bound=_certified(0.99 * oracle_rhs),
    )
    failed = check_eq_9_17_pointwise(
        stage=10065,
        derivative_order=2,
        h=h,
        q=q,
        constant=constant,
        log_power=log_power,
        correction_bound=_certified(1.01 * oracle_rhs),
    )

    assert passed.equation == "9.17"
    assert passed.exponent == exponent
    assert passed.arithmetic_check_passed is True
    assert failed.arithmetic_check_passed is False
    assert passed.status == "formal-structure"
    assert passed.uniform_in_q_verified is False
    assert passed.actual_section9_sequence_verified is False
    assert passed.paper_exact_velocity_available is False


def test_eq_918_pointwise_check_keeps_flat_remainder_separate():
    # Independent paper substitution: for h=1/200, K_m=3/2, j=6998,
    # h*(1/5+j/10)-K_m = 2.
    h = Fraction(1, 200)
    q = 0.25
    constant = 2.0
    log_power = 1
    flat = 1.0e-3
    exponent = Fraction(2)
    leading = constant * q**2 * (1.0 + abs(math.log(q)))
    oracle_rhs = leading + flat

    passed = check_eq_9_18_pointwise(
        stage=6998,
        derivative_order=4,
        h=h,
        q=q,
        derivative_loss=Fraction(3, 2),
        constant=constant,
        log_power=log_power,
        residual_bound=_certified(0.999 * oracle_rhs, "residual interval enclosure"),
        flat_remainder_bound=_certified(flat, "Eq. (9.5) remainder enclosure"),
    )
    failed = check_eq_9_18_pointwise(
        stage=6998,
        derivative_order=4,
        h=h,
        q=q,
        derivative_loss=Fraction(3, 2),
        constant=constant,
        log_power=log_power,
        residual_bound=_certified(1.001 * oracle_rhs, "residual interval enclosure"),
        flat_remainder_bound=_certified(flat, "Eq. (9.5) remainder enclosure"),
    )

    assert passed.equation == "9.18"
    assert passed.exponent == exponent
    assert passed.arithmetic_check_passed is True
    assert failed.arithmetic_check_passed is False
    assert passed.flat_remainder_upper_bound == flat
    assert passed.flat_remainder_provenance == "Eq. (9.5) remainder enclosure"
    assert passed.flat_remainder_all_orders_verified is False
    assert passed.paper_exact_velocity_available is False


def test_log_domain_check_does_not_underflow_for_tiny_q():
    # q**2000 underflows in binary64.  A zero certified upper bound is still a
    # valid arithmetic lower endpoint and must compare without evaluating that
    # power in ordinary scale.
    certificate = check_eq_9_18_pointwise(
        stage=4000198,
        derivative_order=0,
        h=Fraction(1, 200),
        q=1.0e-250,
        derivative_loss=0,
        constant=1.0,
        log_power=0,
        residual_bound=_certified(0.0),
        flat_remainder_bound=_certified(0.0),
    )
    assert certificate.exponent == 2000
    assert math.isfinite(certificate.leading_rhs_log)
    assert certificate.arithmetic_check_passed is True


def test_certificate_rejects_sampled_or_unprovenanced_data():
    with pytest.raises(ValueError, match="kind"):
        CertifiedBoundDatum(1.0, "sampled", "grid maximum")
    with pytest.raises(ValueError, match="provenance"):
        CertifiedBoundDatum(1.0, "paper-derived", "   ")
    with pytest.raises(ValueError):
        CertifiedBoundDatum(float("nan"), "paper-derived", "paper")


def test_certificate_rejects_invalid_paper_parameters():
    datum = _certified(1.0)
    with pytest.raises(ValueError):
        check_eq_9_17_pointwise(
            stage=0,
            derivative_order=0,
            h=Fraction(1, 200),
            q=1.0,
            constant=1.0,
            log_power=0,
            correction_bound=datum,
        )
    with pytest.raises(ValueError):
        check_eq_9_17_pointwise(
            stage=0,
            derivative_order=0,
            h=Fraction(1, 200),
            q=0.5,
            constant=0.0,
            log_power=0,
            correction_bound=datum,
        )
    for bad_loss in (-1, True, float("nan")):
        with pytest.raises(ValueError):
            check_eq_9_18_pointwise(
                stage=0,
                derivative_order=0,
                h=Fraction(1, 200),
                q=0.5,
                derivative_loss=bad_loss,
                constant=1.0,
                log_power=0,
                residual_bound=datum,
                flat_remainder_bound=_certified(0.0),
            )
