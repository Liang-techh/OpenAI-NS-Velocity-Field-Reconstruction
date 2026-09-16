from dataclasses import replace
from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction

import pytest

from openai_ns_reconstruction.axis_amplitude_log_scale import (
    AmplitudeLogSource,
)
from openai_ns_reconstruction.axis_phase_log_enclosure import (
    RationalLogAmplitudeEnclosure,
)


def _source() -> AmplitudeLogSource:
    return AmplitudeLogSource(
        h=Fraction(1, 4),
        j=Fraction(2, 5),
        sigma=Fraction(3, 7),
        eta=Fraction(0),
        enclosure=RationalLogAmplitudeEnclosure(
            phase_lower=Fraction(1, 3),
            phase_upper=Fraction(1, 3),
            Lambda=Decimal("1"),
            log_C=Decimal("0"),
        ),
        midpoint=Decimal("0.3"),
    )


def test_exact_rational_delta_and_amplitude_power() -> None:
    source = _source()

    assert source.delta_lower == Fraction(1, 30)
    assert source.delta_upper == Fraction(1, 30)

    power = source.power(3, Decimal("0.9"))
    assert power.delta_lower == Fraction(1, 10)
    assert power.delta_upper == Fraction(1, 10)
    assert power.q == 3
    assert power.paper_exact is False
    assert power.coefficient_arithmetic_certified is False


def test_compatible_source_and_exact_identity_mismatch() -> None:
    source = _source()

    source.assert_compatible(replace(source))

    with pytest.raises(ValueError, match="eta mismatch"):
        source.assert_compatible(replace(source, eta=Fraction(1, 7)))


def test_power_rejects_invalid_q_and_nonzero_zero_power_scale() -> None:
    source = _source()

    with pytest.raises(ValueError, match="nonnegative integer"):
        source.power(True, Decimal("0"))
    with pytest.raises(ValueError, match="nonnegative integer"):
        source.power(-1, Decimal("0"))
    with pytest.raises(ValueError, match="q=0"):
        source.power(0, Decimal("1"))

    zero_power = source.power(0, Decimal("0"))
    assert zero_power.delta_lower == 0
    assert zero_power.delta_upper == 0


def test_exact_deltas_ignore_hostile_decimal_context() -> None:
    with localcontext() as context:
        context.prec = 2
        context.rounding = ROUND_FLOOR
        source = _source()
        power = source.power(3, Decimal("0.9"))

    assert source.delta_lower == Fraction(1, 30)
    assert source.delta_upper == Fraction(1, 30)
    assert power.delta_lower == Fraction(1, 10)
    assert power.delta_upper == Fraction(1, 10)


def test_power_delta_uses_returned_scale_rounding_not_scaled_base_offset() -> None:
    # Synthetic caller-supplied interval: this is a conditional algebra
    # fixture, not the actual eta=0 phase enclosure.
    source = replace(
        _source(),
        enclosure=RationalLogAmplitudeEnclosure(
            phase_lower=Fraction(10, 3),
            phase_upper=Fraction(10, 3),
            Lambda=Decimal("1"),
            log_C=Decimal("0"),
        ),
        midpoint=Decimal("3.33"),
    )

    with localcontext() as context:
        context.prec = 2
        returned_scale = Decimal(3) * source.midpoint
    assert returned_scale == Decimal("10")

    power = source.power(3, returned_scale)
    assert 3 * source.delta_lower == Fraction(1, 100)
    assert power.delta_lower == 0
    assert power.delta_upper == 0
