from fractions import Fraction

import pytest

from openai_ns_reconstruction.outgoing_release_lag_enclosure import (
    validated_release_lag_enclosure,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _actual_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


def test_actual_release_lag_interval_contains_legacy_value() -> None:
    data = _actual_data()
    tolerance = Fraction(1, 256)

    result = validated_release_lag_enclosure(data, absolute_tolerance=tolerance)

    assert result.h == Fraction.from_float(data.h)
    assert result.lam == Fraction.from_float(data.core.lam)
    assert result.width <= tolerance
    assert result.cells == 128
    assert result.lag.lower > 0
    assert result.paper_exact is False
    assert result.global_axis_norm_certified is False

    legacy_value = Fraction.from_float(data.release_lag_at_ramp_end)
    assert result.lag.lower <= legacy_value <= result.lag.upper


def test_release_lag_refinement_overlap_and_positivity() -> None:
    data = _actual_data()
    loose = validated_release_lag_enclosure(
        data,
        absolute_tolerance=Fraction(1, 64),
    )
    tight = validated_release_lag_enclosure(
        data,
        absolute_tolerance=Fraction(1, 256),
    )

    assert loose.width <= Fraction(1, 64)
    assert tight.width <= Fraction(1, 256)
    assert tight.width < loose.width
    assert max(loose.lag.lower, tight.lag.lower) <= min(loose.lag.upper, tight.lag.upper)
    assert loose.lag.lower > 0
    assert tight.lag.lower > 0


def test_release_lag_invalid_controls_and_cell_cap_fail_closed() -> None:
    data = _actual_data()

    with pytest.raises(TypeError, match="data must be TailData"):
        validated_release_lag_enclosure(object())
    with pytest.raises(ValueError, match="absolute_tolerance must be positive"):
        validated_release_lag_enclosure(data, absolute_tolerance=Fraction(0))
    with pytest.raises(ValueError, match="max_cells must be a positive integer"):
        validated_release_lag_enclosure(data, max_cells=True)
    with pytest.raises(ValueError, match="max_terms must be a positive integer"):
        validated_release_lag_enclosure(data, max_terms=0)
    with pytest.raises(ValueError, match="max_squarings must be a positive integer"):
        validated_release_lag_enclosure(data, max_squarings=False)
    with pytest.raises(ArithmeticError, match="max_cells"):
        validated_release_lag_enclosure(
            data,
            absolute_tolerance=Fraction(1, 256),
            max_cells=1,
        )
