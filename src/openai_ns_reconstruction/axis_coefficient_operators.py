"""Executable jet-level assembly of the pinned ``coefficientOperators`` record.

All coefficient primitives required by ``AxisContraction.coefficientOperators``
are already landed individually on ``AxisCoefficientJetState``.  This module
assembles those exact primitives into one operator record with the same field
layout as the pinned Lean definition:

``product, average, primitive, parameterPrimitive, mulY, j1, j2,``
``param1, param2, dot1, dot2, mixed1, mixed2``.

The factory is anchored to the actual SchedulePressure-derived reference state,
so the common theorem-selected epsilon is inherited rather than supplied by a
caller.  Every application rejects coefficient states from a different
coefficient-space scale.

This is deliberately a jet-level executable record, not a claim of full Lean
``AxisSpace`` membership: the global all-index weighted norm certificate and
continuous-linear/bilinear operator norms are still not materialized.  Hence
``naturalRemainder`` and a genuine Picard iterate remain downstream blockers and
``paper_exact`` stays false.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .axis_coefficient_average import axis_coefficient_average
from .axis_coefficient_inverse_dot_product import axis_coefficient_inverse_dot_product
from .axis_coefficient_inverse_mixed import axis_coefficient_inverse_mixed
from .axis_coefficient_inverse_param_product import axis_coefficient_inverse_param_product
from .axis_coefficient_multiply_y import axis_coefficient_multiply_y
from .axis_coefficient_parameter_primitive import axis_coefficient_parameter_primitive
from .axis_coefficient_primitive import axis_coefficient_primitive
from .axis_coefficient_product import axis_coefficient_product
from .axis_coefficient_reference_state import (
    ActualScheduleReferenceAxisState,
    AxisCoefficientJetState,
)
from .axis_coefficient_regular_inverse import axis_coefficient_regular_inverse


LEAN_FIELD_NAMES = (
    "product",
    "average",
    "primitive",
    "parameterPrimitive",
    "mulY",
    "j1",
    "j2",
    "param1",
    "param2",
    "dot1",
    "dot2",
    "mixed1",
    "mixed2",
)


@dataclass(frozen=True)
class AxisCoefficientOperators:
    """One common-epsilon executable image of Lean ``coefficientOperators``."""

    epsilon: float
    reference_origin: str

    def __post_init__(self) -> None:
        epsilon = float(self.epsilon)
        if not math.isfinite(epsilon) or epsilon <= 0.0:
            raise ValueError("epsilon must be finite and positive")
        if not isinstance(self.reference_origin, str) or not self.reference_origin.strip():
            raise ValueError("reference_origin must be a nonempty provenance label")
        object.__setattr__(self, "epsilon", epsilon)

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def lean_field_names(self) -> tuple[str, ...]:
        return LEAN_FIELD_NAMES

    def _state(self, state: AxisCoefficientJetState) -> AxisCoefficientJetState:
        if not isinstance(state, AxisCoefficientJetState):
            raise TypeError("operator arguments must be AxisCoefficientJetState instances")
        if state.epsilon != self.epsilon:
            raise ValueError("coefficientOperators requires exactly matching epsilon")
        return state

    def _pair(
        self,
        left: AxisCoefficientJetState,
        right: AxisCoefficientJetState,
    ) -> tuple[AxisCoefficientJetState, AxisCoefficientJetState]:
        return self._state(left), self._state(right)

    def product(
        self,
        left: AxisCoefficientJetState,
        right: AxisCoefficientJetState,
    ) -> AxisCoefficientJetState:
        left, right = self._pair(left, right)
        return axis_coefficient_product(left, right)

    def average(self, state: AxisCoefficientJetState) -> AxisCoefficientJetState:
        return axis_coefficient_average(self._state(state))

    def primitive(self, state: AxisCoefficientJetState) -> AxisCoefficientJetState:
        return axis_coefficient_primitive(self._state(state))

    def parameter_primitive(self, state: AxisCoefficientJetState) -> AxisCoefficientJetState:
        return axis_coefficient_parameter_primitive(self._state(state))

    def mul_y(self, state: AxisCoefficientJetState) -> AxisCoefficientJetState:
        return axis_coefficient_multiply_y(self._state(state))

    def j1(self, state: AxisCoefficientJetState) -> AxisCoefficientJetState:
        return axis_coefficient_regular_inverse(self._state(state), 1)

    def j2(self, state: AxisCoefficientJetState) -> AxisCoefficientJetState:
        return axis_coefficient_regular_inverse(self._state(state), 2)

    def param1(
        self,
        left: AxisCoefficientJetState,
        right: AxisCoefficientJetState,
    ) -> AxisCoefficientJetState:
        left, right = self._pair(left, right)
        return axis_coefficient_inverse_param_product(left, right, 1)

    def param2(
        self,
        left: AxisCoefficientJetState,
        right: AxisCoefficientJetState,
    ) -> AxisCoefficientJetState:
        left, right = self._pair(left, right)
        return axis_coefficient_inverse_param_product(left, right, 2)

    def dot1(
        self,
        left: AxisCoefficientJetState,
        right: AxisCoefficientJetState,
    ) -> AxisCoefficientJetState:
        left, right = self._pair(left, right)
        return axis_coefficient_inverse_dot_product(left, right, 1)

    def dot2(
        self,
        left: AxisCoefficientJetState,
        right: AxisCoefficientJetState,
    ) -> AxisCoefficientJetState:
        left, right = self._pair(left, right)
        return axis_coefficient_inverse_dot_product(left, right, 2)

    def mixed1(
        self,
        left: AxisCoefficientJetState,
        right: AxisCoefficientJetState,
    ) -> AxisCoefficientJetState:
        left, right = self._pair(left, right)
        return axis_coefficient_inverse_mixed(left, right, 1)

    def mixed2(
        self,
        left: AxisCoefficientJetState,
        right: AxisCoefficientJetState,
    ) -> AxisCoefficientJetState:
        left, right = self._pair(left, right)
        return axis_coefficient_inverse_mixed(left, right, 2)


def actual_schedule_coefficient_operators(
    reference: ActualScheduleReferenceAxisState,
) -> AxisCoefficientOperators:
    """Assemble the pinned record at the actual theorem-selected epsilon.

    The actual reference pair is required so callers cannot invent an unrelated
    coefficient-space scale.  Both reference halves are rechecked before the
    record is returned.
    """

    if not isinstance(reference, ActualScheduleReferenceAxisState):
        raise TypeError("reference must be ActualScheduleReferenceAxisState")
    if reference.phi.epsilon != reference.u.epsilon or reference.phi.epsilon != reference.epsilon:
        raise ValueError("actual reference states must share exactly one epsilon")
    return AxisCoefficientOperators(
        epsilon=reference.epsilon,
        reference_origin="actual SchedulePressure referencePair",
    )
