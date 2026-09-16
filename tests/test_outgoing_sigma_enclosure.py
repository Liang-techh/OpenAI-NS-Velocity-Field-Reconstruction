from decimal import Decimal, localcontext
from fractions import Fraction

import pytest

from openai_ns_reconstruction.outgoing_sigma_enclosure import (
    RationalInterval,
    validated_exp_negative,
    validated_outgoing_sigma,
    validated_sigma_primitive,
)
from openai_ns_reconstruction.outgoing_tail import _sigma_primitive


def test_exact_sigma_and_primitive_branches() -> None:
    for x in (Fraction(-1), Fraction(0)):
        assert validated_outgoing_sigma(x).width == 0
        assert validated_sigma_primitive(x).width == 0
    for x in (Fraction(1), Fraction(2)):
        assert validated_outgoing_sigma(x) == RationalInterval(Fraction(1), Fraction(1))
        assert validated_sigma_primitive(x) == RationalInterval(x - Fraction(1, 2), x - Fraction(1, 2))
    assert validated_outgoing_sigma(Fraction(1, 2)) == RationalInterval(
        Fraction(1, 2), Fraction(1, 2)
    )


def test_exp_negative_reference_width_and_caps() -> None:
    tolerance = Fraction(1, 256)
    with localcontext() as context:
        context.prec = 140
        for x in (Fraction(0), Fraction(1), Fraction(2), Fraction(10)):
            interval = validated_exp_negative(x, absolute_tolerance=tolerance)
            reference = (-Decimal(x.numerator) / Decimal(x.denominator)).exp()
            assert interval.width <= tolerance
            assert interval.lower <= Fraction(reference) <= interval.upper
    with pytest.raises(ValueError, match="absolute_tolerance must be positive"):
        validated_exp_negative(Fraction(0), absolute_tolerance=Fraction(0))
    with pytest.raises(ValueError, match="max_terms must be a positive integer"):
        validated_exp_negative(Fraction(0), absolute_tolerance=tolerance, max_terms=True)
    with pytest.raises(ValueError, match="max_squarings must be a positive integer"):
        validated_exp_negative(Fraction(0), absolute_tolerance=tolerance, max_squarings=0)
    with pytest.raises(ArithmeticError, match="max_terms"):
        validated_exp_negative(Fraction(1), absolute_tolerance=tolerance, max_terms=1)
    with pytest.raises(ArithmeticError, match="max_squarings"):
        validated_exp_negative(Fraction(4), absolute_tolerance=tolerance, max_squarings=1)


def test_sigma_reflection_and_decimal_diagnostic() -> None:
    tolerance = Fraction(1, 256)
    left = validated_outgoing_sigma(Fraction(1, 4), absolute_tolerance=tolerance)
    right = validated_outgoing_sigma(Fraction(3, 4), absolute_tolerance=tolerance)
    assert right == RationalInterval(1 - left.upper, 1 - left.lower)
    with localcontext() as context:
        context.prec = 140
        for x, interval in ((Fraction(1, 4), left), (Fraction(3, 4), right)):
            t = Decimal(x.numerator) / Decimal(x.denominator)
            r = (-1 / (t * t) + 1 / ((1 - t) * (1 - t))).exp()
            reference = r / (1 + r)
            assert interval.width <= tolerance
            assert interval.lower <= Fraction(reference) <= interval.upper


def test_primitive_reflection_refinement_and_quadrature_diagnostic() -> None:
    loose = validated_sigma_primitive(Fraction(1, 2), absolute_tolerance=Fraction(1, 64))
    tight = validated_sigma_primitive(Fraction(1, 2), absolute_tolerance=Fraction(1, 256))
    left = validated_sigma_primitive(Fraction(1, 4), absolute_tolerance=Fraction(1, 256))
    reflected = validated_sigma_primitive(Fraction(3, 4), absolute_tolerance=Fraction(1, 256))
    assert loose.width <= Fraction(1, 64)
    assert tight.width <= Fraction(1, 256)
    assert max(loose.lower, tight.lower) <= min(loose.upper, tight.upper)
    assert reflected == RationalInterval(
        Fraction(1, 4) + left.lower,
        Fraction(1, 4) + left.upper,
    )
    assert validated_sigma_primitive(Fraction(3, 2)) == RationalInterval(Fraction(1), Fraction(1))

    # The existing Gauss rule is a diagnostic comparison only; the rational
    # interval above is the certified result under test.
    quadrature_value = _sigma_primitive(0.5)
    assert tight.lower <= Fraction.from_float(quadrature_value) <= tight.upper


def test_primitive_invalid_controls_fail_before_exact_branches() -> None:
    with pytest.raises(TypeError, match="x must be a Fraction"):
        validated_sigma_primitive(0)
    with pytest.raises(ValueError, match="absolute_tolerance must be positive"):
        validated_sigma_primitive(Fraction(0), absolute_tolerance=Fraction(0))
    with pytest.raises(ValueError, match="max_cells must be a positive integer"):
        validated_sigma_primitive(Fraction(0), max_cells=True)
    with pytest.raises(ValueError, match="max_terms must be a positive integer"):
        validated_sigma_primitive(Fraction(0), max_terms=0)
    with pytest.raises(ValueError, match="max_squarings must be a positive integer"):
        validated_sigma_primitive(Fraction(0), max_squarings=False)
    with pytest.raises(ArithmeticError, match="max_cells"):
        validated_sigma_primitive(
            Fraction(1, 2), absolute_tolerance=Fraction(1, 4096), max_cells=1
        )
