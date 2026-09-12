"""Executable coefficient-jet image of the pinned natural-axis linear operator.

The official ``AxisResolvent.naturalOperator`` is

``Q(A) = (1/2) * J_2(chi * A)``

where ``J_2`` is the zero-datum regular radial inverse and
``chi = H^2/(H^2+sigma^2)`` is the actual theorem-selected natural-axis
multiplier.  The landed reference angular jets satisfy the pinned closed form

``phi0[n] = (-chi/2)^n / (n! (n+1)!)``.

In particular ``chi = -4 * phi0[1]`` together with all eta derivatives.  This
module uses that exact identity to expose ``chi`` in the same lazy
``AxisCoefficientJetState`` representation as the already-landed reference
pair, then composes the landed coefficient ``product`` and ``j2`` operators to
materialize ``Q``.

No caller-supplied sigma, chi table, epsilon, divisor, or replacement operator
is accepted.  This remains a jet-level executable object: the global all-index
weighted ``AxisSpace`` norm certificate and the infinite alternating resolvent
series are still separate blockers, so this module does not claim a genuine
``naturalResolvent`` or paper-exact leading profile.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .axis_coefficient_operators import (
    AxisCoefficientOperators,
    actual_schedule_coefficient_operators,
)
from .axis_coefficient_reference_state import (
    ActualScheduleReferenceAxisState,
    AxisCoefficientJetState,
)


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _scaled_state(
    state: AxisCoefficientJetState,
    scalar: float,
    *,
    origin: str,
) -> AxisCoefficientJetState:
    scalar = float(scalar)
    if not math.isfinite(scalar):
        raise ValueError("scalar must be finite")

    def provider(n: int, m: int, eta: float) -> float:
        value = scalar * state.jet(n, m, eta)
        if not math.isfinite(value):
            raise ArithmeticError("scaled coefficient jet must remain finite")
        return value

    return AxisCoefficientJetState(
        epsilon=state.epsilon,
        origin=origin,
        _jet_provider=provider,
    )


def actual_schedule_chi_coefficient_state(
    reference: ActualScheduleReferenceAxisState,
) -> AxisCoefficientJetState:
    """Return the actual radially constant ``chi`` coefficient state.

    ``AxisReference.reference_coefficient`` gives
    ``phi0[1] = -chi/4``.  The landed angular reference family evaluates
    arbitrary actual eta derivatives analytically, so the same identity gives
    every parameter jet of ``chi`` without numerical differentiation or a new
    pressure/sigma choice.
    """

    if not isinstance(reference, ActualScheduleReferenceAxisState):
        raise TypeError("reference must be ActualScheduleReferenceAxisState")

    def provider(n: int, m: int, eta: float) -> float:
        n = _index(n, "n")
        m = _index(m, "m")
        if n != 0:
            # Preserve the exact reference-state index/window validation even
            # though every positive radial row is mathematically zero.
            reference.phi.jet(0, m, eta)
            return 0.0
        value = -4.0 * reference.phi.jet(1, m, eta)
        if not math.isfinite(value):
            raise ArithmeticError("chi coefficient jet must remain finite")
        return value

    return AxisCoefficientJetState(
        epsilon=reference.epsilon,
        origin="actual SchedulePressure chi via pinned phi0[1] = -chi/4",
        _jet_provider=provider,
    )


@dataclass(frozen=True)
class AxisCoefficientNaturalOperator:
    """Jet-level executable image of ``AxisResolvent.naturalOperator``."""

    epsilon: float
    operators: AxisCoefficientOperators
    chi: AxisCoefficientJetState

    def __post_init__(self) -> None:
        epsilon = float(self.epsilon)
        if not math.isfinite(epsilon) or epsilon <= 0.0:
            raise ValueError("epsilon must be finite and positive")
        if self.operators.epsilon != epsilon or self.chi.epsilon != epsilon:
            raise ValueError("naturalOperator inputs must use exactly one epsilon")
        object.__setattr__(self, "epsilon", epsilon)

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def natural_resolvent_materialized(self) -> bool:
        return False

    def apply(self, state: AxisCoefficientJetState) -> AxisCoefficientJetState:
        """Apply the exact pinned map ``A -> (1/2) J_2(chi*A)``."""

        if not isinstance(state, AxisCoefficientJetState):
            raise TypeError("state must be AxisCoefficientJetState")
        if state.epsilon != self.epsilon:
            raise ValueError("naturalOperator requires exactly matching epsilon")
        product = self.operators.product(self.chi, state)
        inverse = self.operators.j2(product)
        return _scaled_state(
            inverse,
            0.5,
            origin=f"pinned AxisResolvent.naturalOperator applied to ({state.origin})",
        )

    def __call__(self, state: AxisCoefficientJetState) -> AxisCoefficientJetState:
        return self.apply(state)


def actual_schedule_natural_operator(
    reference: ActualScheduleReferenceAxisState,
) -> AxisCoefficientNaturalOperator:
    """Build the pinned ``naturalOperator`` from the actual schedule reference.

    The actual reference state is the sole input.  Its theorem-selected
    coefficient epsilon and analytic angular jet family determine both the
    operator scale and ``chi``.
    """

    if not isinstance(reference, ActualScheduleReferenceAxisState):
        raise TypeError("reference must be ActualScheduleReferenceAxisState")
    operators = actual_schedule_coefficient_operators(reference)
    chi = actual_schedule_chi_coefficient_state(reference)
    return AxisCoefficientNaturalOperator(
        epsilon=reference.epsilon,
        operators=operators,
        chi=chi,
    )
