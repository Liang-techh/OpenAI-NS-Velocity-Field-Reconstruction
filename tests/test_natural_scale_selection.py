import math

import pytest

from openai_ns_reconstruction.natural_scale_selection import (
    JET_RADIAL_DERIVATIVE_BOUND_AT_FIVE,
    JET_VALUE_BOUND_AT_FIVE,
    contraction_threshold,
    normalization_threshold,
    select_natural_scale,
    stability_scale,
)


def test_exact_pinned_threshold_algebra_and_closed_jet_sums() -> None:
    B = 2.0
    L = 3.0

    # Independent closed forms for sum q^n and sum (n+1)q^n at q=5/20.
    q = 0.25
    assert JET_VALUE_BOUND_AT_FIVE == 1.0 / (1.0 - q)
    assert JET_RADIAL_DERIVATIVE_BOUND_AT_FIVE == 1.0 / (1.0 - q) ** 2

    contraction = contraction_threshold(B, L)
    stability = stability_scale(B)

    assert contraction == 1.0 + B + L
    assert stability == 1.0 + (14000.0 / 9.0) * B


def test_selection_takes_positive_profile_threshold_then_exact_C() -> None:
    witness = select_natural_scale(
        remainder_bound=2.0,
        remainder_lipschitz=3.0,
        phase_real_part_sup=0.01,
    )

    expected_lambda = max(1.0 + 2.0 + 3.0, 1.0 + (14000.0 / 9.0) * 2.0)
    assert witness.Lambda == expected_lambda
    assert witness.C == math.exp(expected_lambda * 0.01)
    assert witness.admits(witness.Lambda, witness.C)


def test_admissibility_recomputes_C_threshold_at_candidate_Lambda() -> None:
    witness = select_natural_scale(
        remainder_bound=0.01,
        remainder_lipschitz=0.25,
        phase_real_part_sup=0.2,
    )
    larger_lambda = witness.Lambda + 1.0
    exact_larger_C = normalization_threshold(larger_lambda, witness.phase_real_part_sup)

    assert witness.admits(larger_lambda, exact_larger_C)
    assert not witness.admits(larger_lambda, math.nextafter(exact_larger_C, 0.0))
    below_scale = math.nextafter(witness.Lambda, 0.0)
    assert not witness.admits(below_scale, witness.C * 10.0)


def test_fail_closed_for_uncertified_numeric_domains() -> None:
    with pytest.raises(ValueError, match="remainder_bound"):
        contraction_threshold(-1.0, 0.0)
    with pytest.raises(ValueError, match="remainder_bound"):
        stability_scale(float("nan"))
    with pytest.raises(ValueError, match="phase_real_part_sup"):
        normalization_threshold(2.0, -0.1)
    with pytest.raises(ArithmeticError, match="binary64"):
        normalization_threshold(1000.0, 1.0)
