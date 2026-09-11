"""Executable coefficient-jet product from the pinned ``AxisOperators`` algebra.

The official natural-axis fixed point is built in a compatible coefficient
space, not from unrelated coefficient samples.  Current main already exposes
the actual SchedulePressure-derived reference ``phi/u`` eta-jets through
``AxisCoefficientJetState``.  This module lands one genuine ``AxisOperators``
primitive on that representation: the product jet family.

For radial degree ``n`` and parameter derivative order ``m`` the pinned Lean
sources define

``J_{A B}[n,m] = sum_{i+j=n} sum_{k+l=m} binom(m,k) J_A[i,k] J_B[j,l]``.

This is exactly ``AxisWeightEstimates.jetProduct`` / ``AxisOperators.productFamily``.
The returned object remains fail-closed with respect to complete ``AxisSpace``
membership: no finite collection of evaluated products is promoted to the
global all-index weighted norm certificate required by ``ofJetFamily``.
Consequently this module does not yet claim the full ``coefficientOperators``
record, ``naturalRemainder``, a Picard iterate, or a paper-exact leading field.
"""

from __future__ import annotations

import math

from .axis_coefficient_reference_state import AxisCoefficientJetState


def _same_epsilon(left: AxisCoefficientJetState, right: AxisCoefficientJetState) -> float:
    """Require one common pinned coefficient-space weight scale."""

    if left.epsilon != right.epsilon:
        raise ValueError("coefficient product requires exactly matching epsilon")
    return left.epsilon


def product_jet(
    left: AxisCoefficientJetState,
    right: AxisCoefficientJetState,
    n: int,
    m: int,
    eta: float,
) -> float:
    """Evaluate the pinned radial-convolution/parameter-Leibniz product jet.

    Input validation for ``n``, ``m`` and the official eta window is delegated
    to the underlying actual-jet states.  Every summand is checked before the
    stable finite sum is returned so floating overflow cannot silently become a
    coefficient value.
    """

    _same_epsilon(left, right)

    # Force the same index/window guards even in the n=m=0 case before any
    # arithmetic is attempted.  These calls are also the first actual factors.
    left.jet(n, m, eta)
    right.jet(n, m, eta)

    terms: list[float] = []
    for i in range(n + 1):
        j = n - i
        for k in range(m + 1):
            l = m - k
            term = float(math.comb(m, k)) * left.jet(i, k, eta) * right.jet(j, l, eta)
            if not math.isfinite(term):
                raise ArithmeticError("coefficient product term must remain finite")
            terms.append(term)

    value = math.fsum(terms)
    if not math.isfinite(value):
        raise ArithmeticError("coefficient product jet must remain finite")
    return value


def axis_coefficient_product(
    left: AxisCoefficientJetState,
    right: AxisCoefficientJetState,
) -> AxisCoefficientJetState:
    """Return a lazy actual-jet state for the pinned coefficient product.

    The constructor accepts no replacement coefficient table or independent
    epsilon.  The product therefore remains tied to the actual derivative
    families of its operands.  ``AxisCoefficientJetState`` intentionally keeps
    ``global_axis_norm_certified=False`` until the all-index weighted estimate
    is materialized separately.
    """

    epsilon = _same_epsilon(left, right)

    def provider(n: int, m: int, eta: float) -> float:
        return product_jet(left, right, n, m, eta)

    return AxisCoefficientJetState(
        epsilon=epsilon,
        origin=f"pinned AxisOperators product of ({left.origin}) and ({right.origin})",
        _jet_provider=provider,
    )
