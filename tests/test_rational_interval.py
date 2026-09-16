from fractions import Fraction

import pytest

from openai_ns_reconstruction.rational_interval import (
    RationalInterval,
    dyadic_step,
    interval_add,
    interval_divide_positive,
    interval_multiply,
    interval_subtract,
    map_monotone,
    outward_round_dyadic,
    positive_cap,
)


def test_arithmetic_handles_sign_changes_and_zero() -> None:
    left = RationalInterval(Fraction(-2), Fraction(3))
    right = RationalInterval(Fraction(-5), Fraction(7))

    assert interval_add(left, right) == RationalInterval(Fraction(-7), Fraction(10))
    assert interval_subtract(left, right) == RationalInterval(Fraction(-9), Fraction(8))
    assert interval_multiply(left, right) == RationalInterval(Fraction(-15), Fraction(21))
    assert interval_multiply(
        left,
        RationalInterval(Fraction(0), Fraction(0)),
    ) == RationalInterval(Fraction(0), Fraction(0))


def test_positive_division_near_zero_is_outward_and_fail_closed() -> None:
    numerator = RationalInterval(Fraction(-1), Fraction(2))
    denominator = RationalInterval(Fraction(1, 10**12), Fraction(1, 10**6))

    assert interval_divide_positive(numerator, denominator) == RationalInterval(
        Fraction(-10**12),
        Fraction(2 * 10**12),
    )

    with pytest.raises(ValueError, match="strictly positive"):
        interval_divide_positive(
            numerator,
            RationalInterval(Fraction(0), Fraction(1)),
        )


def test_monotone_endpoint_propagation_in_both_directions() -> None:
    interval = RationalInterval(Fraction(2), Fraction(5))

    assert map_monotone(interval, lambda x: x * x) == RationalInterval(
        Fraction(4), Fraction(25)
    )
    assert map_monotone(
        interval,
        lambda x: Fraction(1, x),
        increasing=False,
    ) == RationalInterval(Fraction(1, 5), Fraction(1, 2))


def test_dyadic_rounding_is_strictly_outward() -> None:
    interval = RationalInterval(Fraction(-3, 10), Fraction(7, 10))
    rounded = outward_round_dyadic(interval, Fraction(1, 4))

    assert rounded == RationalInterval(Fraction(-1, 2), Fraction(3, 4))
    assert rounded.lower <= interval.lower <= interval.upper <= rounded.upper
    assert dyadic_step(Fraction(3, 10)) == Fraction(1, 4)


def test_caps_reject_bool_zero_and_negative() -> None:
    for value in (True, False, 0, -1):
        with pytest.raises(ValueError, match="positive integer"):
            positive_cap(value, "cap")

    assert positive_cap(3, "cap") == 3
