"""Coefficientwise exact action of the pinned natural-axis resolvent.

The official ``AxisResolvent.naturalResolvent`` is the alternating series

``R(A) = sum_{k >= 0} (-Q)^k A``

for ``Q = AxisResolvent.naturalOperator``.  The Lean proof does not assume a
small operator norm: ``Q`` increases the radial filtration by one and its
powers satisfy a factorial majorant.

For an individual radial coefficient row ``n`` the filtration statement is
stronger than a numerical truncation estimate: ``Q^k A`` has zero row ``n``
for every ``k > n``.  Therefore the infinite series is *exactly finite* at each
coefficient coordinate,

``R(A)[n,m](eta) = sum_{k=0}^n (-1)^k (Q^k A)[n,m](eta)``.

This module materializes that exact lazy coefficient-jet action on the actual
SchedulePressure-derived coefficient scale.  It introduces no tolerance,
iteration cutoff, fitted coefficient, or caller-selected theorem parameter.

The scope is deliberately representation-level.  The repository still does
not claim the global all-index weighted ``AxisSpace`` norm/membership object
or the Lean continuous-linear-map norm certificate, so this module does not
make the leading profile paper-exact by itself.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .axis_coefficient_natural_operator import (
    AxisCoefficientNaturalOperator,
    actual_schedule_natural_operator,
)
from .axis_coefficient_reference_state import (
    ActualScheduleReferenceAxisState,
    AxisCoefficientJetState,
)


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


@dataclass(frozen=True)
class AxisCoefficientNaturalResolvent:
    """Lazy coefficientwise image of ``AxisResolvent.naturalResolvent``.

    The returned state evaluates each radial row with the exact finite prefix
    forced by the pinned radial filtration.  No approximate infinite-series
    stopping rule is used.
    """

    epsilon: float
    operator: AxisCoefficientNaturalOperator

    def __post_init__(self) -> None:
        epsilon = float(self.epsilon)
        if not math.isfinite(epsilon) or epsilon <= 0.0:
            raise ValueError("epsilon must be finite and positive")
        if self.operator.epsilon != epsilon:
            raise ValueError("naturalResolvent must use the naturalOperator epsilon")
        object.__setattr__(self, "epsilon", epsilon)

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def coefficientwise_series_exact(self) -> bool:
        """The radial filtration makes every requested coefficient sum finite."""

        return True

    @property
    def natural_remainder_materialized(self) -> bool:
        return False

    @staticmethod
    def exact_term_count(radial_index: int) -> int:
        """Return the exact number of alternating terms needed for one row."""

        return _index(radial_index, "radial_index") + 1

    def apply(self, state: AxisCoefficientJetState) -> AxisCoefficientJetState:
        """Apply the pinned alternating resolvent coefficientwise exactly.

        ``Q`` raises radial order by one.  Hence row ``n`` receives exactly the
        terms ``Q^0`` through ``Q^n`` and all later terms vanish identically at
        that row.  The operator iterates are cached as lazy states so repeated
        coordinate requests reuse the same exact composition graph.
        """

        if not isinstance(state, AxisCoefficientJetState):
            raise TypeError("state must be AxisCoefficientJetState")
        if state.epsilon != self.epsilon:
            raise ValueError("naturalResolvent requires exactly matching epsilon")

        iterates: list[AxisCoefficientJetState] = [state]

        def q_power(k: int) -> AxisCoefficientJetState:
            while len(iterates) <= k:
                iterates.append(self.operator(iterates[-1]))
            return iterates[k]

        def provider(n: int, m: int, eta: float) -> float:
            n = _index(n, "n")
            m = _index(m, "m")
            total = 0.0
            sign = 1.0
            for k in range(n + 1):
                total += sign * q_power(k).jet(n, m, eta)
                sign = -sign
            if not math.isfinite(total):
                raise ArithmeticError("naturalResolvent coefficient jet must remain finite")
            return total

        return AxisCoefficientJetState(
            epsilon=self.epsilon,
            origin=(
                "pinned AxisResolvent.naturalResolvent coefficientwise finite filtration "
                f"applied to ({state.origin})"
            ),
            _jet_provider=provider,
        )

    def __call__(self, state: AxisCoefficientJetState) -> AxisCoefficientJetState:
        return self.apply(state)


def actual_schedule_natural_resolvent(
    reference: ActualScheduleReferenceAxisState,
) -> AxisCoefficientNaturalResolvent:
    """Build the pinned natural resolvent from the actual schedule reference."""

    if not isinstance(reference, ActualScheduleReferenceAxisState):
        raise TypeError("reference must be ActualScheduleReferenceAxisState")
    operator = actual_schedule_natural_operator(reference)
    return AxisCoefficientNaturalResolvent(
        epsilon=reference.epsilon,
        operator=operator,
    )
