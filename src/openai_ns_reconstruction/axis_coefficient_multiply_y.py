r"""Executable coefficient-jet multiplication by the radial variable ``Y``.

This module materializes the pinned ``AxisOperators.mulY`` row map on the
already landed actual ``SchedulePressure``-derived coefficient eta-jets.  At
the pinned official Lean commit, ``multiplyYScale`` is zero in row zero and one
in every positive row, while ``multiplyYData`` reads the predecessor radial
row with no parameter-jet shift.  Hence

``J_out[0,m](eta) = 0``
``J_out[n,m](eta) = J_in[n-1,m](eta)`` for ``n >= 1``.

This is exactly the coefficient action of ``F(Y,eta) -> Y * F(Y,eta)``.
The returned lazy state remains fail-closed for global ``AxisSpace``
membership: no all-index weighted norm certificate, full
``coefficientOperators`` record, or ``naturalRemainder`` evaluation is claimed
here.
"""

from __future__ import annotations

import math

from .axis_coefficient_reference_state import AxisCoefficientJetState


def multiply_y_jet(
    state: AxisCoefficientJetState,
    n: int,
    m: int,
    eta: float,
) -> float:
    """Evaluate one pinned ``AxisOperators.mulY`` coefficient jet."""

    # Match AxisCoefficientJetState's fail-closed index semantics before the
    # zero-row branch; in particular bool must not silently act as an integer.
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError("n must be a nonnegative integer")

    if n == 0:
        # Route through the source state so m and eta retain the exact same
        # pinned validation even though the output row is identically zero.
        state.jet(0, m, eta)
        return 0.0

    value = state.jet(n - 1, m, eta)
    if not math.isfinite(value):
        raise ArithmeticError("coefficient multiply-Y jet must remain finite")
    return value


def axis_coefficient_multiply_y(
    state: AxisCoefficientJetState,
) -> AxisCoefficientJetState:
    """Return the lazy pinned ``AxisOperators.mulY`` jet family.

    No caller-supplied radial scale, epsilon, shift, or replacement coefficient
    table is accepted.  The row-zero condition and predecessor-row shift are
    fixed by the pinned ``multiplyYData`` definition.
    """

    def provider(n: int, m: int, eta: float) -> float:
        return multiply_y_jet(state, n, m, eta)

    return AxisCoefficientJetState(
        epsilon=state.epsilon,
        origin=f"pinned AxisOperators mulY of ({state.origin})",
        _jet_provider=provider,
    )
