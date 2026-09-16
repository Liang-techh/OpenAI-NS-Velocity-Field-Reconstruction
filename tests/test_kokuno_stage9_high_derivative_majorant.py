from fractions import Fraction

import pytest

from openai_ns_reconstruction.kokuno_stage9_high_derivative_majorant import (
    a38_counting_receipt,
    a39_euclidean_squared_majorant,
    a39_scalar_component_majorant,
    binomial_partition_sum,
    multiindex_binomial,
    normalize_spacetime_multiindex,
)


def test_m3_a38_counts_are_exact():
    receipt = a38_counting_receipt((2, 1, 0, 0))
    assert receipt.total_derivative_order == 3
    assert receipt.enumerated_binomial_sum == 8
    assert receipt.expected_binomial_sum == 8
    assert receipt.exact_partition_identity
    assert receipt.linear_term_count == 5
    assert receipt.mixed_product_weight == 48
    assert receipt.quadratic_product_weight == 24


def test_direct_enumeration_matches_two_to_total_order():
    for index in ((0, 0, 0, 0), (1, 0, 0, 0), (2, 1, 0, 2), (1, 2, 3, 1)):
        assert binomial_partition_sum(index) == 2 ** sum(index)


def test_m0_agrees_with_existing_public_counting_surface():
    majorant = a39_scalar_component_majorant((0, 0, 0, 0), Fraction(5, 3), Fraction(7, 5))
    assert majorant == Fraction(672, 25)


def test_m3_majorant_and_squared_euclidean_bound_are_exact():
    index = (2, 1, 0, 0)
    majorant = a39_scalar_component_majorant(index, Fraction(5, 3), Fraction(7, 5))
    assert majorant == Fraction(4151, 25)
    assert a39_euclidean_squared_majorant(index, Fraction(5, 3), Fraction(7, 5)) == Fraction(51692403, 625)


def test_under_counting_leibniz_partition_fails_by_exact_gap():
    v_bound = Fraction(5, 3)
    e_bound = Fraction(7, 5)
    correct = a39_scalar_component_majorant((2, 1, 0, 0), v_bound, e_bound)
    incorrect_half_partition = (Fraction(5) + 24 * v_bound + 12 * e_bound) * e_bound
    assert correct - incorrect_half_partition == Fraction(1988, 25)
    assert incorrect_half_partition < correct


def test_multiindex_binomial_rejects_non_subindex():
    with pytest.raises(ValueError):
        multiindex_binomial((1, 0, 0, 0), (2, 0, 0, 0))


def test_fail_closed_on_approximate_negative_or_malformed_inputs():
    with pytest.raises(TypeError):
        a39_scalar_component_majorant((1, 0, 0, 0), 1.0, Fraction(1))
    with pytest.raises(ValueError):
        a39_scalar_component_majorant((1, 0, 0, 0), Fraction(-1), Fraction(1))
    with pytest.raises(TypeError):
        normalize_spacetime_multiindex((True, 0, 0, 0))
    with pytest.raises(ValueError):
        normalize_spacetime_multiindex((1, 0, 0))
