r"""Executable coefficient-jet radial primitive from pinned ``AxisOperators``.

The official natural-axis fixed-point remainder uses the zero-at-axis radial
primitive on the compatible coefficient space.  At the pinned official Lean
commit, ``AxisOperators.primitiveData`` is the ``rowData`` operator with

``primitiveScale 0 = 0`` and ``primitiveScale (n + 1) = 1 / (n + 1)``.

Consequently, for actual coefficient jets ``J[n,m](eta)``, the output jets are

``J_prim[0,m](eta) = 0``
``J_prim[n,m](eta) = J[n-1,m](eta) / n`` for ``n >= 1``.

This is exactly the coefficient map for the radial primitive
``P[F](Y,eta) = integral_0^Y F(s,eta) ds``.  Parameter-jet order is unchanged,
so adjacent eta-jet compatibility is inherited from the input state.

The returned state remains fail-closed for global ``AxisSpace`` membership.
This module materializes the pinned jet family only; it does not certify the
all-index weighted norm, instantiate the complete ``coefficientOperators``
record, or evaluate ``naturalRemainder``.
"""

from __future__ import annotations

import math

from .axis_coefficient_reference_state import AxisCoefficientJetState


def primitive_jet(
    state: AxisCoefficientJetState,
    n: int,
    m: int,
    eta: float,
) -> float:
    """Evaluate one pinned zero-at-axis radial-primitive coefficient jet."""

    if n == 0:
        # Route through the underlying state so the official index/window guards
        # remain active even though the primitive's zeroth radial row vanishes.
        state.jet(0, m, eta)
        return 0.0

    source = state.jet(n - 1, m, eta)
    value = source / float(n)
    if not math.isfinite(value):
        raise ArithmeticError("coefficient primitive jet must remain finite")
    return value


def axis_coefficient_primitive(state: AxisCoefficientJetState) -> AxisCoefficientJetState:
    """Return the lazy pinned ``AxisOperators.primitive`` jet family.

    No caller-supplied integration constant, epsilon, radial scale, or
    replacement coefficient table is accepted.  The integration constant is
    fixed by the pinned ``primitiveScale 0 = 0`` row.
    """

    def provider(n: int, m: int, eta: float) -> float:
        return primitive_jet(state, n, m, eta)

    return AxisCoefficientJetState(
        epsilon=state.epsilon,
        origin=f"pinned AxisOperators primitive of ({state.origin})",
        _jet_provider=provider,
    )
