from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_cutoff_schedule import (
    build_slow_borel_cutoff_schedule,
)
from openai_ns_reconstruction.background_tail_order import (
    FixedPrefixTailOrderCertificate,
    chart_tail_exponent_exact,
    certify_schedule_physical_tail,
    minimal_prefix_for_chart_target,
    minimal_prefix_for_physical_target,
    physical_tail_exponent_exact,
)


def _unit_bound_schedule(h: float, max_order: int):
    return build_slow_borel_cutoff_schedule(
        h,
        [[1.0] * (j + 3) for j in range(1, max_order + 1)],
    )


def test_exact_fixed_prefix_exponents_match_hand_algebra():
    h = 0.125
    J = 7
    M = 2
    b = -0.25

    assert chart_tail_exponent_exact(h, J, M) == Fraction(-1, 1)
    assert physical_tail_exponent_exact(h, J, M, b) == Fraction(-13, 4)

    cert = FixedPrefixTailOrderCertificate(h, J, M, b)
    assert cert.dyadic_prefactor_exact == Fraction(1, 128)
    assert cert.worst_chart_exponent_exact == Fraction(-1, 1)
    assert cert.physical_exponent_exact == Fraction(-13, 4)
    assert cert.chart_exponent_exact(0) == Fraction(1, 1)


def test_minimal_chart_prefix_is_exact_and_minimal():
    # h=1/8, M=2, P=1 requires (J+1)/8 - 2 >= 1, hence J>=23.
    cert = minimal_prefix_for_chart_target(0.125, 2, 1.0)
    assert cert.prefix_order == 23
    assert cert.worst_chart_exponent_exact == Fraction(1, 1)
    assert cert.certifies_chart_target(1.0)

    previous = FixedPrefixTailOrderCertificate(0.125, 22, 2)
    assert not previous.certifies_chart_target(1.0)


def test_minimal_physical_prefix_is_exact_and_minimal():
    # (J+1)/8 + 1/4 - 4 >= 1/2 gives J>=33.
    cert = minimal_prefix_for_physical_target(0.125, 2, 0.25, 0.5)
    assert cert.prefix_order == 33
    assert cert.physical_exponent_exact == Fraction(1, 2)
    assert cert.certifies_physical_target(0.5)

    previous = FixedPrefixTailOrderCertificate(0.125, 32, 2, 0.25)
    assert not previous.certifies_physical_target(0.5)


def test_minimum_prefix_order_is_respected_without_weakening_target():
    cert = minimal_prefix_for_physical_target(
        0.125,
        0,
        0.0,
        0.25,
        minimum_prefix_order=9,
    )
    assert cert.prefix_order == 9
    assert cert.physical_exponent_exact == Fraction(5, 4)
    assert cert.certifies_physical_target(0.25)


def test_constructed_schedule_gate_fails_closed_when_prefix_is_too_short():
    short = _unit_bound_schedule(0.125, 2)
    with pytest.raises(ValueError, match="does not reach"):
        certify_schedule_physical_tail(short, 0, 0.0, 0.5)

    enough = _unit_bound_schedule(0.125, 3)
    cert = certify_schedule_physical_tail(enough, 0, 0.0, 0.5)
    assert cert.prefix_order == 3
    assert cert.physical_exponent_exact == Fraction(1, 2)
    assert cert.prefix_order <= enough.max_order


def test_derivative_budget_and_invalid_inputs_fail_closed():
    with pytest.raises(ValueError, match="M <= J\+3"):
        FixedPrefixTailOrderCertificate(0.125, 0, 4)
    with pytest.raises(ValueError, match="m <= J\+3"):
        chart_tail_exponent_exact(0.125, 0, 4)
    with pytest.raises(ValueError, match="0 < h < 1/2"):
        minimal_prefix_for_chart_target(0.5, 0, 1.0)
    with pytest.raises(ValueError, match="nonnegative integer"):
        minimal_prefix_for_chart_target(0.125, -1, 1.0)
    with pytest.raises(ValueError, match="finite"):
        minimal_prefix_for_physical_target(0.125, 0, 0.0, float("nan"))
    with pytest.raises(TypeError):
        certify_schedule_physical_tail(object(), 0, 0.0, 1.0)
