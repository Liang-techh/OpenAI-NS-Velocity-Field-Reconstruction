from fractions import Fraction

import pytest

from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.outgoing_tail_debt_enclosure import (
    validated_tail_debt_enclosure,
)


def _actual_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


def test_actual_tail_debt_interval_contains_legacy_value() -> None:
    data = _actual_data()
    tolerance = Fraction(1, 10**8)

    result = validated_tail_debt_enclosure(data, absolute_tolerance=tolerance)

    assert result.h == Fraction.from_float(data.h)
    assert result.width <= tolerance
    assert result.cells == 512
    assert Fraction(0) < result.rho.lower <= result.rho.upper < Fraction(1)
    assert result.paper_exact is False
    assert result.global_axis_norm_certified is False

    legacy_value = Fraction.from_float(data.tail_debt)
    assert result.debt.lower <= legacy_value <= result.debt.upper


def test_refinement_overlap_and_positive_factor_bounds() -> None:
    data = _actual_data()
    loose = validated_tail_debt_enclosure(
        data,
        absolute_tolerance=Fraction(1, 10**6),
    )
    tight = validated_tail_debt_enclosure(
        data,
        absolute_tolerance=Fraction(1, 10**8),
    )

    assert loose.width <= Fraction(1, 10**6)
    assert tight.width <= Fraction(1, 10**8)
    assert max(loose.debt.lower, tight.debt.lower) <= min(loose.debt.upper, tight.debt.upper)

    ratio_lower = loose.rho.lower / (1 - loose.rho.lower)
    ratio_upper = loose.rho.upper / (1 - loose.rho.upper)
    assert loose.debt.lower >= ratio_lower
    assert loose.debt.upper <= 27 * ratio_upper


def test_invalid_controls_and_cell_cap_fail_closed() -> None:
    data = _actual_data()

    with pytest.raises(TypeError, match="data must be TailData"):
        validated_tail_debt_enclosure(object())
    with pytest.raises(ValueError, match="absolute_tolerance must be positive"):
        validated_tail_debt_enclosure(data, absolute_tolerance=Fraction(0))
    with pytest.raises(ValueError, match="max_cells must be a positive integer"):
        validated_tail_debt_enclosure(data, max_cells=True)
    with pytest.raises(ValueError, match="max_terms must be a positive integer"):
        validated_tail_debt_enclosure(data, max_terms=0)
    with pytest.raises(ValueError, match="max_squarings must be a positive integer"):
        validated_tail_debt_enclosure(data, max_squarings=False)
    with pytest.raises(ArithmeticError, match="max_cells"):
        validated_tail_debt_enclosure(
            data,
            absolute_tolerance=Fraction(1, 10**8),
            max_cells=1,
        )
