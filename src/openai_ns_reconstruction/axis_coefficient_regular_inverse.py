r"""Executable pinned regular radial inverse on coefficient eta-jets.

This module materializes ``AxisOperators.regularInverse`` from the pinned
official Lean source on the already landed actual ``SchedulePressure``-derived
coefficient states.  The pinned coefficient map is

``J_out[0,m](eta) = 0``

and, for ``n >= 1``,

``J_out[n,m](eta) = J_in[n-1,m](eta) / radialDivisor(r, n-1)``,

where

``radialDivisor(r, k) = (k + 1) * (k + r)`` and ``r >= 1``.

At the evaluated-profile level this is the zero-datum regular inverse of
``Y * d_Y^2 + r * d_Y``.  No caller-supplied divisor table, integration
constant, epsilon, or replacement coefficients are accepted.

The returned lazy state remains fail-closed for global ``AxisSpace``
membership.  This increment does not certify the all-index weighted norm,
instantiate the complete ``coefficientOperators`` record, or evaluate
``naturalRemainder``.
"""

from __future__ import annotations

import math

from .axis_coefficient_reference_state import AxisCoefficientJetState


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def radial_divisor(r: int, n: int) -> float:
    """Return the pinned ``radialDivisor r n = (n+1)(n+r)``."""

    r = _positive_integer(r, "r")
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError("n must be a nonnegative integer")
    value = float((n + 1) * (n + r))
    if not math.isfinite(value) or value <= 0.0:
        raise ArithmeticError("regular-inverse radial divisor must be finite and positive")
    return value


def regular_inverse_jet(
    state: AxisCoefficientJetState,
    r: int,
    n: int,
    m: int,
    eta: float,
) -> float:
    """Evaluate one pinned ``AxisOperators.regularInverse`` coefficient jet."""

    r = _positive_integer(r, "r")
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError("n must be a nonnegative integer")

    if n == 0:
        # Route through the source state so m and eta retain exactly the same
        # pinned validation even though the zero-datum output row vanishes.
        state.jet(0, m, eta)
        return 0.0

    source = state.jet(n - 1, m, eta)
    value = source / radial_divisor(r, n - 1)
    if not math.isfinite(value):
        raise ArithmeticError("coefficient regular-inverse jet must remain finite")
    return value


def axis_coefficient_regular_inverse(
    state: AxisCoefficientJetState,
    r: int,
) -> AxisCoefficientJetState:
    """Return the lazy pinned ``AxisOperators.regularInverse`` jet family.

    ``r`` is restricted to the theorem's positive integer parameter domain.
    The radial divisor and zero datum are fixed by the pinned Lean definitions.
    """

    r = _positive_integer(r, "r")

    def provider(n: int, m: int, eta: float) -> float:
        return regular_inverse_jet(state, r, n, m, eta)

    return AxisCoefficientJetState(
        epsilon=state.epsilon,
        origin=f"pinned AxisOperators regularInverse(r={r}) of ({state.origin})",
        _jet_provider=provider,
    )
