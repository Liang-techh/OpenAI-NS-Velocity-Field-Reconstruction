from decimal import (
    Decimal,
    DivisionByZero,
    InvalidOperation,
    Overflow,
    ROUND_CEILING,
    ROUND_FLOOR,
    Underflow,
    localcontext,
)
from fractions import Fraction

import pytest

from openai_ns_reconstruction.axis_phase_log_enclosure import (
    RationalLogAmplitudeEnclosure,
)


def test_huge_lambda_preserves_exact_scaled_width_and_directed_bounds() -> None:
    phase_lower = Fraction(2)
    phase_upper = phase_lower + Fraction(1, 10**800)
    enclosure = RationalLogAmplitudeEnclosure(
        phase_lower,
        phase_upper,
        Decimal("2E784"),
        Decimal("4E784"),
    )

    assert enclosure.lower == 0
    assert enclosure.width == Fraction(2, 10**16)
    assert enclosure.width > 0
    lower, upper = enclosure.decimal_bounds(96)
    assert Fraction(lower) <= enclosure.lower
    assert Fraction(upper) >= enclosure.upper
    assert lower == 0
    assert upper > lower
    assert enclosure.paper_exact is False


def test_decimal_bounds_ignore_hostile_ambient_context() -> None:
    enclosure = RationalLogAmplitudeEnclosure(
        Fraction(-7, 13),
        Fraction(-7, 13) + Fraction(1, 10**37),
        Decimal("1E784"),
        Decimal("-3.25E783"),
    )

    with localcontext() as context:
        context.prec = 3
        context.rounding = ROUND_CEILING
        context.Emax = 2
        context.Emin = -2
        context.traps[DivisionByZero] = True
        context.traps[InvalidOperation] = True
        context.traps[Overflow] = True
        context.traps[Underflow] = True
        lower, upper = enclosure.decimal_bounds(48)

    assert Fraction(lower) <= enclosure.lower
    assert Fraction(upper) >= enclosure.upper
    assert lower <= upper


def test_zero_phase_interval_maps_to_exact_singleton() -> None:
    enclosure = RationalLogAmplitudeEnclosure(
        Fraction(1, 2),
        Fraction(1, 2),
        Decimal("1E4"),
        Decimal("1E3"),
    )

    assert enclosure.lower == enclosure.upper == Fraction(4000)
    assert enclosure.width == 0
    assert enclosure.decimal_bounds() == (Decimal("4000"), Decimal("4000"))


def test_constructor_and_precision_validation() -> None:
    with pytest.raises(TypeError, match="phase_lower must be a Fraction"):
        RationalLogAmplitudeEnclosure(1, Fraction(2), Decimal(1), Decimal(0))
    with pytest.raises(TypeError, match="phase_upper must be a Fraction"):
        RationalLogAmplitudeEnclosure(Fraction(1), 2, Decimal(1), Decimal(0))
    with pytest.raises(ValueError, match="phase bounds must be ordered"):
        RationalLogAmplitudeEnclosure(Fraction(2), Fraction(1), Decimal(1), Decimal(0))
    with pytest.raises(TypeError, match="Lambda must be a Decimal"):
        RationalLogAmplitudeEnclosure(Fraction(0), Fraction(1), 1, Decimal(0))
    with pytest.raises(ValueError, match="Lambda must be positive"):
        RationalLogAmplitudeEnclosure(Fraction(0), Fraction(1), Decimal(0), Decimal(0))
    with pytest.raises(ValueError, match="Lambda must be finite"):
        RationalLogAmplitudeEnclosure(Fraction(0), Fraction(1), Decimal("Infinity"), Decimal(0))
    with pytest.raises(TypeError, match="log_C must be a Decimal"):
        RationalLogAmplitudeEnclosure(Fraction(0), Fraction(1), Decimal(1), 0)
    with pytest.raises(ValueError, match="log_C must be finite"):
        RationalLogAmplitudeEnclosure(Fraction(0), Fraction(1), Decimal(1), Decimal("NaN"))

    enclosure = RationalLogAmplitudeEnclosure(
        Fraction(0), Fraction(1), Decimal(1), Decimal(0)
    )
    with pytest.raises(ValueError, match="precision must be a positive integer"):
        enclosure.decimal_bounds(0)
    with pytest.raises(ValueError, match="precision must be a positive integer"):
        enclosure.decimal_bounds(-1)
    with pytest.raises(ValueError, match="precision must be a positive integer"):
        enclosure.decimal_bounds(True)
