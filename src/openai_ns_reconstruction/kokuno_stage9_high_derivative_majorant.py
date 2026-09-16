"""Exact clean-room replay of the public Stage-9 A38/A39 counting seam.

This module intentionally implements only finite combinatorial bookkeeping from
public formulas.  It does not execute or copy the Kokuno source checker.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import comb
from typing import Iterable, Tuple

SpacetimeMultiIndex = Tuple[int, int, int, int]


def _exact_fraction(value: int | Fraction, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an exact int or Fraction")
    return Fraction(value)


def normalize_spacetime_multiindex(index: Iterable[int]) -> SpacetimeMultiIndex:
    values = tuple(index)
    if len(values) != 4:
        raise ValueError("A38 uses a four-entry spacetime multi-index (t,x1,x2,x3)")
    for entry in values:
        if isinstance(entry, bool) or not isinstance(entry, int):
            raise TypeError("multi-index entries must be nonnegative integers")
        if entry < 0:
            raise ValueError("multi-index entries must be nonnegative")
    return values  # type: ignore[return-value]


def total_order(index: Iterable[int]) -> int:
    return sum(normalize_spacetime_multiindex(index))


def sub_multiindices(index: Iterable[int]) -> tuple[SpacetimeMultiIndex, ...]:
    alpha = normalize_spacetime_multiindex(index)
    return tuple(product(*(range(a + 1) for a in alpha)))  # type: ignore[return-value]


def multiindex_binomial(alpha: Iterable[int], beta: Iterable[int]) -> int:
    a = normalize_spacetime_multiindex(alpha)
    b = normalize_spacetime_multiindex(beta)
    if any(bi > ai for ai, bi in zip(a, b)):
        raise ValueError("beta must satisfy beta <= alpha componentwise")
    out = 1
    for ai, bi in zip(a, b):
        out *= comb(ai, bi)
    return out


def binomial_partition_sum(index: Iterable[int]) -> int:
    """Compute sum_{J<=I} binom(I,J) by direct finite enumeration."""

    alpha = normalize_spacetime_multiindex(index)
    return sum(multiindex_binomial(alpha, beta) for beta in sub_multiindices(alpha))


@dataclass(frozen=True)
class A38CountingReceipt:
    index: SpacetimeMultiIndex
    total_derivative_order: int
    enumerated_binomial_sum: int
    expected_binomial_sum: int
    linear_term_count: int
    mixed_product_weight: int
    quadratic_product_weight: int

    @property
    def exact_partition_identity(self) -> bool:
        return self.enumerated_binomial_sum == self.expected_binomial_sum


def a38_counting_receipt(index: Iterable[int]) -> A38CountingReceipt:
    """Replay the finite coefficient count behind public formulas A38/A39.

    For one Cartesian output component, A38 has five linear derivative terms.
    Its two cross products for each of three spatial directions contribute
    ``6 * sum_{J<=I} binom(I,J)`` and the quadratic product contributes
    ``3 * sum_{J<=I} binom(I,J)``.
    """

    alpha = normalize_spacetime_multiindex(index)
    m = sum(alpha)
    enumerated = binomial_partition_sum(alpha)
    expected = 1 << m
    if enumerated != expected:
        raise ArithmeticError("multi-index binomial partition identity failed")
    return A38CountingReceipt(
        index=alpha,
        total_derivative_order=m,
        enumerated_binomial_sum=enumerated,
        expected_binomial_sum=expected,
        linear_term_count=5,
        mixed_product_weight=6 * enumerated,
        quadratic_product_weight=3 * enumerated,
    )


def a39_scalar_component_majorant(
    index: Iterable[int],
    v_mplus2_bound: int | Fraction,
    e_mplus2_bound: int | Fraction,
) -> Fraction:
    """Exact rational A39 scalar-component majorant before the sqrt(3) factor."""

    counts = a38_counting_receipt(index)
    v_bound = _exact_fraction(v_mplus2_bound, name="v_mplus2_bound")
    e_bound = _exact_fraction(e_mplus2_bound, name="e_mplus2_bound")
    if v_bound < 0 or e_bound < 0:
        raise ValueError("certified norm bounds must be nonnegative")
    coefficient = (
        Fraction(counts.linear_term_count)
        + Fraction(counts.mixed_product_weight) * v_bound
        + Fraction(counts.quadratic_product_weight) * e_bound
    )
    return coefficient * e_bound


def a39_euclidean_squared_majorant(
    index: Iterable[int],
    v_mplus2_bound: int | Fraction,
    e_mplus2_bound: int | Fraction,
) -> Fraction:
    """Square the public sqrt(3) vector-output bound without approximation."""

    scalar = a39_scalar_component_majorant(index, v_mplus2_bound, e_mplus2_bound)
    return Fraction(3) * scalar * scalar
