import math
from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_truncation_residual import (
    convolution_at,
    omitted_slow_tail_majorant,
    recurrence_truncation_breakdown,
    slow_order_exact,
    slow_weight,
)


def _fixture(order: int = 3):
    # Pure algebra fixture, not paper coefficient data.
    pair = [
        [0.20, -0.10, 0.05, 0.03],
        [0.07, 0.04, -0.02, 0.01],
        [-0.03, 0.06, 0.02, -0.04],
        [0.01, -0.05, 0.08, 0.02],
    ]
    shifted = [0.11, -0.09, 0.07, -0.05]
    pair = [row[: order + 1] for row in pair[: order + 1]]
    shifted = shifted[: order + 1]
    linear = []
    for n in range(order + 1):
        conv = sum(pair[i][n - i] for i in range(n + 1))
        previous = 0.0 if n == 0 else shifted[n - 1]
        linear.append(previous - conv)
    return linear, pair, shifted


def _weights(q: float, h: float, b: float, order: int):
    top = max(2 * order, order + 1)
    return [q ** (b + 2.0 * n * h) for n in range(top + 1)]


def _direct_lhs(weights, linear, pair, shifted, order):
    linear_sum = sum(weights[n] * linear[n] for n in range(order + 1))
    pair_sum = sum(
        weights[i + j] * pair[i][j]
        for i in range(order + 1)
        for j in range(order + 1)
    )
    shifted_sum = sum(weights[n + 1] * shifted[n] for n in range(order + 1))
    return linear_sum + pair_sum - shifted_sum


def _direct_tail(weights, pair, shifted, order):
    return (
        sum(
            weights[i + j] * pair[i][j]
            for i in range(order + 1)
            for j in range(order + 1)
            if i + j > order
        )
        - weights[order + 1] * shifted[order]
    )


def test_recurrence_truncation_matches_independent_direct_expansion():
    N = 3
    q, h, b = 0.5, 0.125, 0.25
    linear, pair, shifted = _fixture(N)
    weights = _weights(q, h, b, N)

    result = recurrence_truncation_breakdown(N, weights, linear, pair, shifted)
    direct = _direct_lhs(weights, linear, pair, shifted, N)
    direct_tail = _direct_tail(weights, pair, shifted, N)

    assert result.lhs == pytest.approx(direct, rel=0.0, abs=2e-16)
    assert result.rhs == pytest.approx(direct, rel=0.0, abs=2e-16)
    assert result.retained_recurrence == pytest.approx(0.0, abs=2e-16)
    assert result.max_recurrence_abs == pytest.approx(0.0, abs=2e-16)
    assert result.omitted_remainder == pytest.approx(direct_tail, rel=0.0, abs=2e-16)
    assert abs(result.algebraic_defect) <= 4e-16


def test_nonzero_retained_recurrence_is_not_silently_dropped():
    N = 2
    linear, pair, shifted = _fixture(N)
    linear[1] += 0.125
    weights = _weights(0.75, 0.1, -0.2, N)

    result = recurrence_truncation_breakdown(N, weights, linear, pair, shifted)
    direct = _direct_lhs(weights, linear, pair, shifted, N)

    assert result.recurrence_values[1] == pytest.approx(0.125)
    assert result.retained_recurrence == pytest.approx(weights[1] * 0.125)
    assert result.rhs == pytest.approx(direct, rel=2e-15, abs=2e-15)
    assert result.lhs != pytest.approx(result.omitted_remainder, rel=1e-12, abs=1e-12)


def test_zero_recurrence_tail_obeys_first_omitted_slow_order_majorant():
    N = 3
    h, b = 0.125, 0.25
    linear, pair, shifted = _fixture(N)

    for q in (0.75, 0.5, 0.25, 0.125):
        weights = _weights(q, h, b, N)
        result = recurrence_truncation_breakdown(N, weights, linear, pair, shifted)
        majorant = omitted_slow_tail_majorant(q, h, b, N, pair, shifted)

        # b + 2(N+1)h = 1/4 + 1 = 5/4 exactly for this fixture.
        assert majorant.exponent_exact == Fraction(5, 4)
        assert result.max_recurrence_abs <= 2e-16
        assert abs(result.omitted_remainder) <= majorant.binary64_bound()
        assert abs(result.lhs) <= majorant.binary64_bound() + 5e-16


def test_majorant_prefactor_is_independent_hand_l1_tail():
    N = 2
    _, pair, shifted = _fixture(N)
    majorant = omitted_slow_tail_majorant(0.5, 0.125, 0.0, N, pair, shifted)
    manual = sum(
        abs(pair[i][j])
        for i in range(N + 1)
        for j in range(N + 1)
        if i + j > N
    ) + abs(shifted[N])
    assert majorant.coefficient_l1 == pytest.approx(manual, rel=0.0, abs=1e-16)
    assert majorant.exponent_exact == Fraction(3, 4)


def test_order_zero_keeps_the_final_shifted_term():
    pair = [[0.3]]
    shifted = [0.4]
    linear = [-0.3]  # recurrence_0 = L_0 + K_00 = 0
    weights = [1.0, 0.25]
    result = recurrence_truncation_breakdown(0, weights, linear, pair, shifted)
    assert result.recurrence_values == pytest.approx((0.0,))
    assert result.pair_tail == 0.0
    assert result.shifted_tail == pytest.approx(-0.1)
    assert result.lhs == pytest.approx(-0.1)


def test_convolution_helper_and_exact_slow_order():
    _, pair, _ = _fixture(3)
    assert convolution_at(pair, 3, 2) == pytest.approx(
        pair[0][2] + pair[1][1] + pair[2][0]
    )
    assert slow_order_exact(0.125, 7) == Fraction(7, 4)
    assert slow_weight(0.5, 0.125, 0.25, 4) == pytest.approx(0.5 ** 1.25)


def test_invalid_shapes_domains_and_underflow_fail_closed():
    linear, pair, shifted = _fixture(2)
    weights = _weights(0.5, 0.125, 0.0, 2)

    with pytest.raises(ValueError, match="weights must have length"):
        recurrence_truncation_breakdown(2, weights[:-1], linear, pair, shifted)
    with pytest.raises(ValueError, match="pair coefficients"):
        recurrence_truncation_breakdown(2, weights, linear, pair[:-1], shifted)
    bad_pair = [row[:] for row in pair]
    bad_pair[0][0] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        recurrence_truncation_breakdown(2, weights, linear, bad_pair, shifted)
    with pytest.raises(ValueError, match="0 < q <= 1"):
        omitted_slow_tail_majorant(1.1, 0.125, 0.0, 2, pair, shifted)
    with pytest.raises(ValueError, match="index <= order"):
        convolution_at(pair, 2, 3)
    with pytest.raises(OverflowError, match="underflows"):
        slow_weight(1e-300, 0.125, 0.0, 10)


def test_zero_coefficient_tail_has_exact_zero_binary64_bound():
    majorant = omitted_slow_tail_majorant(0.5, 0.125, 0.0, 1, [[0.0, 0.0], [0.0, 0.0]], [0.0, 0.0])
    assert majorant.coefficient_l1 == 0.0
    assert majorant.log_bound == -math.inf
    assert majorant.binary64_bound() == 0.0
