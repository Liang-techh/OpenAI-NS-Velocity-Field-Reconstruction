from decimal import Decimal
from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_actual_shear_replay import (
    replay_actual_shear_identity,
    require_corrected_actual_shear_split,
    require_leading_only_equivalence,
)


def fixture_kwargs():
    return dict(
        frozen_shear_magnitude=Fraction(7, 5),
        delta_normal=Fraction(-2, 9),
        delta_tangent=Fraction(3, 7),
        flux_normal=Fraction(11, 13),
        flux_tangent=Fraction(-5, 8),
    )


def test_corrected_actual_shear_split_closes_exactly():
    result = require_corrected_actual_shear_split(**fixture_kwargs())
    assert result.identity_residual == 0
    assert result.exact_identity_verified
    assert result.lhs_actual_production == result.rhs_corrected_split
    assert result.leading_only_defect == result.retained_remainder
    assert result.retained_remainder == Fraction(2987, 6552)


def test_nonzero_remainder_cannot_be_dropped_as_leading_only():
    result = replay_actual_shear_identity(**fixture_kwargs())
    assert result.retained_remainder != 0
    assert not result.leading_only_is_valid
    with pytest.raises(ValueError, match="nonzero actual-shear remainder"):
        require_leading_only_equivalence(**fixture_kwargs())


def test_leading_only_equivalence_is_allowed_only_for_exact_zero_remainder():
    result = require_leading_only_equivalence(
        frozen_shear_magnitude=Fraction(4, 3),
        delta_normal=0,
        delta_tangent=0,
        flux_normal=Fraction(-9, 10),
        flux_tangent=Fraction(17, 19),
    )
    assert result.retained_remainder == 0
    assert result.leading_only_defect == 0


def test_checker_rejects_nonexact_or_invalid_theorem_inputs():
    with pytest.raises(ValueError, match="nonnegative"):
        replay_actual_shear_identity(
            frozen_shear_magnitude=-1,
            delta_normal=0,
            delta_tangent=0,
            flux_normal=1,
            flux_tangent=0,
        )
    for bad in (0.1, Decimal("0.1"), True):
        kwargs = fixture_kwargs()
        kwargs["delta_normal"] = bad
        with pytest.raises(TypeError):
            replay_actual_shear_identity(**kwargs)
