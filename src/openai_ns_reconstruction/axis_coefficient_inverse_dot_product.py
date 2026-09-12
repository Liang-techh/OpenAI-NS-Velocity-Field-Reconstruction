r"""Executable pinned inverse radial-dot product on coefficient eta-jets.

This module materializes ``AxisOperators.inverseDotProduct`` from the pinned
official Lean source on the landed actual ``SchedulePressure``-derived
coefficient states.  In the evaluated natural-axis algebra this is the regular
radial inverse ``J_r`` applied to ``A * (Y * partial_Y B)``.  The exact
coefficient family is

``J_out[0,m](eta) = 0``

and, for ``n >= 1``,

``J_out[n,m](eta) = sum_{i+j=n-1} sum_{k+l=m}``
``    j * binom(m,k) J_A[i,k](eta) J_B[j,l](eta)``
``    / radialDivisor(r,n-1)``.

This is the pinned ``differentialFamily`` specialization ``p=0`` and
``d(j)=j``.  The radial dot is therefore the Euler derivative ``Y*partial_Y``
on the second argument; no standalone radial derivative table is accepted.
No caller-supplied divisor table, epsilon, integration constant, or replacement
coefficient data are accepted.

The returned lazy state remains fail-closed for global ``AxisSpace``
membership.  This increment does not instantiate the complete
``coefficientOperators`` record, evaluate ``naturalRemainder``, or claim a
paper-exact leading field.
"""

from __future__ import annotations

import math

from .axis_coefficient_reference_state import AxisCoefficientJetState
from .axis_coefficient_regular_inverse import radial_divisor


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _same_epsilon(left: AxisCoefficientJetState, right: AxisCoefficientJetState) -> float:
    if left.epsilon != right.epsilon:
        raise ValueError("inverse dot-product requires exactly matching epsilon")
    return left.epsilon


def inverse_dot_product_jet(
    left: AxisCoefficientJetState,
    right: AxisCoefficientJetState,
    r: int,
    n: int,
    m: int,
    eta: float,
) -> float:
    """Evaluate one pinned ``AxisOperators.inverseDotProduct`` coefficient jet.

    ``r`` is validated by the landed pinned ``radial_divisor`` helper.  The
    explicit index checks occur before arithmetic so booleans cannot be
    silently accepted as Python integers.  Both source states retain the common
    pinned eta-window guard.
    """

    _same_epsilon(left, right)
    n = _index(n, "n")
    m = _index(m, "m")
    # Validate the theorem-domain positive integer r even on the zero row.
    radial_divisor(r, 0)

    if n == 0:
        # Keep source-window and parameter-jet validation active although the
        # zero-datum regular inverse has an identically zero first row.
        left.jet(0, m, eta)
        right.jet(0, m, eta)
        return 0.0

    radial_degree = n - 1
    divisor = radial_divisor(r, radial_degree)
    terms: list[float] = []
    for i in range(radial_degree + 1):
        j = radial_degree - i
        radial_factor = float(j)
        for k in range(m + 1):
            l = m - k
            term = (
                radial_factor
                * float(math.comb(m, k))
                / divisor
                * left.jet(i, k, eta)
                * right.jet(j, l, eta)
            )
            if not math.isfinite(term):
                raise ArithmeticError("inverse dot-product term must remain finite")
            terms.append(term)

    value = math.fsum(terms)
    if not math.isfinite(value):
        raise ArithmeticError("inverse dot-product jet must remain finite")
    return value


def axis_coefficient_inverse_dot_product(
    left: AxisCoefficientJetState,
    right: AxisCoefficientJetState,
    r: int,
) -> AxisCoefficientJetState:
    """Return the lazy pinned inverse radial-dot-product coefficient family."""

    epsilon = _same_epsilon(left, right)
    # Fail before constructing a state if r is outside the theorem domain.
    radial_divisor(r, 0)

    def provider(n: int, m: int, eta: float) -> float:
        return inverse_dot_product_jet(left, right, r, n, m, eta)

    return AxisCoefficientJetState(
        epsilon=epsilon,
        origin=(
            f"pinned AxisOperators inverseDotProduct(r={r}) of "
            f"({left.origin}) and ({right.origin})"
        ),
        _jet_provider=provider,
    )
