"""Executable coefficient-jet radial average from pinned ``AxisOperators``.

The official natural-axis fixed-point remainder uses the regular radial average
operator on the compatible coefficient space.  In the pinned Lean sources
``AxisOperators.averageData`` is the row operator

``J_avg[n,m](eta) = J[n,m](eta) / (n + 1)``.

This is the coefficient form of ``\int_0^1 F(t Y, eta) dt`` for a radial power
series ``F(Y,eta) = sum_n a_n(eta) Y^n``.  The operator leaves the parameter-jet
order unchanged and therefore preserves the adjacent-jet compatibility already
carried by ``AxisCoefficientJetState``.

The returned state intentionally remains fail-closed for global ``AxisSpace``
membership.  This module materializes the exact jet family used by the pinned
operator, but it does not certify the all-index weighted supremum, instantiate
the full ``coefficientOperators`` record, or evaluate ``naturalRemainder``.
"""

from __future__ import annotations

import math

from .axis_coefficient_reference_state import AxisCoefficientJetState


def average_jet(
    state: AxisCoefficientJetState,
    n: int,
    m: int,
    eta: float,
) -> float:
    """Evaluate one pinned radial-average coefficient jet.

    The underlying state performs the official nonnegative-index and eta-window
    validation.  For every radial degree ``n``, ``AxisOperators.averageData``
    uses the exact scalar ``1 / (n + 1)`` and the same source radial row.
    """

    source = state.jet(n, m, eta)
    scale = 1.0 / float(n + 1)
    value = scale * source
    if not math.isfinite(value):
        raise ArithmeticError("coefficient average jet must remain finite")
    return value


def axis_coefficient_average(state: AxisCoefficientJetState) -> AxisCoefficientJetState:
    """Return the lazy pinned ``AxisOperators.average`` jet family.

    No caller-supplied scale, epsilon, or replacement coefficient table is
    accepted.  ``AxisCoefficientJetState`` deliberately keeps
    ``global_axis_norm_certified=False`` until the official all-index weighted
    estimate is materialized separately.
    """

    def provider(n: int, m: int, eta: float) -> float:
        return average_jet(state, n, m, eta)

    return AxisCoefficientJetState(
        epsilon=state.epsilon,
        origin=f"pinned AxisOperators average of ({state.origin})",
        _jet_provider=provider,
    )
