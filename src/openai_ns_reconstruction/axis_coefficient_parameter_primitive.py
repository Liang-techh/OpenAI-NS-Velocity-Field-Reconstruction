r"""Executable coefficient-jet parameter primitive from pinned ``AxisOperators``.

The pinned natural-axis contraction does not apply parameter differentiation as
an independently bounded operator.  Instead, ``AxisWeightEstimates`` packages
one eta derivative together with the zero-at-axis radial primitive:

``parameterPrimitiveJet(f)[0,m] = 0``
``parameterPrimitiveJet(f)[n,m] = f[n-1,m+1] / n`` for ``n >= 1``.

Thus, for actual compatible coefficient jets ``J[n,m](eta)``, this module
materializes the coefficient map of

``P_eta[F](Y,eta) = integral_0^Y partial_eta F(s,eta) ds``.

The ``m + 1`` offset is essential: it is the actual next parameter derivative,
not a fitted finite difference or a caller-supplied derivative table.  The
returned state remains fail-closed for global ``AxisSpace`` membership.  This
module does not certify the all-index weighted norm, instantiate the complete
``coefficientOperators`` record, or evaluate ``naturalRemainder``.
"""

from __future__ import annotations

import math

from .axis_coefficient_reference_state import AxisCoefficientJetState


def parameter_primitive_jet(
    state: AxisCoefficientJetState,
    n: int,
    m: int,
    eta: float,
) -> float:
    """Evaluate one pinned radial-primitive-plus-eta-derivative jet."""

    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError("n must be a nonnegative integer")

    if n == 0:
        # Preserve the underlying state's exact validation for m and eta even
        # though the zero-at-axis row vanishes identically.
        state.jet(0, m, eta)
        return 0.0

    # ``state.jet`` performs the pinned m/window validation.  The m+1 access is
    # the exact AxisWeightEstimates.parameterPrimitiveJet shift.
    source = state.jet(n - 1, m + 1, eta)
    value = source / float(n)
    if not math.isfinite(value):
        raise ArithmeticError("coefficient parameter primitive jet must remain finite")
    return value


def axis_coefficient_parameter_primitive(
    state: AxisCoefficientJetState,
) -> AxisCoefficientJetState:
    """Return the lazy pinned ``AxisOperators.parameterPrimitive`` jet family.

    No caller-supplied eta derivative, integration constant, epsilon, radial
    scale, or replacement coefficient table is accepted.
    """

    def provider(n: int, m: int, eta: float) -> float:
        return parameter_primitive_jet(state, n, m, eta)

    return AxisCoefficientJetState(
        epsilon=state.epsilon,
        origin=f"pinned AxisOperators parameterPrimitive of ({state.origin})",
        _jet_provider=provider,
    )
