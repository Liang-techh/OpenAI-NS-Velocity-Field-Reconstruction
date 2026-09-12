"""Executable coefficient-jet composition of pinned ``naturalRemainder``.

This module translates ``AxisContraction.naturalRemainder`` literally onto the
actual SchedulePressure-derived coefficient scale.  The operator record, fixed
AxisData fields, and angular resolvent are all constructed from the landed
actual reference state; callers may only provide the theorem variables ``t``,
``a``, and the coefficient-state pair ``x``.

The implementation is intentionally fail-closed.  It does not claim global
``AxisSpace`` membership/norm certification, does not choose the manuscript's
final amplitude coefficient ``a`` or ``Lambda``/``C``, and therefore does not
materialize ``naturalRemainder(x0)`` or a Picard iterate by itself.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .axis_coefficient_data import (
    ActualScheduleAxisCoefficientData,
    actual_schedule_axis_coefficient_data,
)
from .axis_coefficient_natural_resolvent import (
    AxisCoefficientNaturalResolvent,
    actual_schedule_natural_resolvent,
)
from .axis_coefficient_operators import (
    AxisCoefficientOperators,
    actual_schedule_coefficient_operators,
)
from .axis_coefficient_reference_state import (
    ActualScheduleReferenceAxisState,
    AxisCoefficientJetState,
)


def _checked_state(
    state: AxisCoefficientJetState,
    epsilon: float,
    name: str,
) -> AxisCoefficientJetState:
    if not isinstance(state, AxisCoefficientJetState):
        raise TypeError(f"{name} must be AxisCoefficientJetState")
    if state.epsilon != epsilon:
        raise ValueError(f"{name} must use exactly matching epsilon")
    return state


def _linear_combination(
    epsilon: float,
    origin: str,
    *terms: tuple[float, AxisCoefficientJetState],
) -> AxisCoefficientJetState:
    checked: list[tuple[float, AxisCoefficientJetState]] = []
    for coefficient, state in terms:
        coefficient = float(coefficient)
        if not math.isfinite(coefficient):
            raise ValueError("linear-combination coefficients must be finite")
        checked.append((coefficient, _checked_state(state, epsilon, "linear-combination state")))

    def provider(n: int, m: int, eta: float) -> float:
        return sum(coefficient * state.jet(n, m, eta) for coefficient, state in checked)

    return AxisCoefficientJetState(
        epsilon=epsilon,
        origin=origin,
        _jet_provider=provider,
    )


@dataclass(frozen=True)
class AxisCoefficientNaturalRemainder:
    """Actual-scale executable image of pinned ``AxisContraction.naturalRemainder``."""

    epsilon: float
    operators: AxisCoefficientOperators
    data: ActualScheduleAxisCoefficientData
    resolvent: AxisCoefficientNaturalResolvent

    def __post_init__(self) -> None:
        epsilon = float(self.epsilon)
        if not math.isfinite(epsilon) or epsilon <= 0.0:
            raise ValueError("epsilon must be finite and positive")
        if self.operators.epsilon != epsilon:
            raise ValueError("coefficientOperators epsilon mismatch")
        if self.data.epsilon != epsilon:
            raise ValueError("AxisData epsilon mismatch")
        if self.resolvent.epsilon != epsilon:
            raise ValueError("naturalResolvent epsilon mismatch")
        object.__setattr__(self, "epsilon", epsilon)

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def actual_amplitude_materialized(self) -> bool:
        return False

    def apply(
        self,
        t: float,
        a: AxisCoefficientJetState,
        x: tuple[AxisCoefficientJetState, AxisCoefficientJetState],
    ) -> tuple[AxisCoefficientJetState, AxisCoefficientJetState]:
        """Evaluate the exact pinned remainder composition on lazy coefficient jets.

        ``t`` and ``a`` are exactly the theorem arguments of ``naturalRemainder``.
        This method does not select them; the later manuscript instantiation must
        supply ``t = 1/Lambda`` and the coefficient state of the normalized
        amplitude associated with the theorem-selected ``Lambda`` and ``C``.
        """

        t = float(t)
        if not math.isfinite(t):
            raise ValueError("t must be finite")
        a = _checked_state(a, self.epsilon, "a")
        if not isinstance(x, tuple) or len(x) != 2:
            raise TypeError("x must be a pair of AxisCoefficientJetState instances")
        phi = _checked_state(x[0], self.epsilon, "x.phi")
        u = _checked_state(x[1], self.epsilon, "x.u")

        O = self.operators
        d = self.data
        e = self.epsilon
        bu = O.average(u)

        angular_linear = _linear_combination(
            e,
            "pinned angularLinearCoefficient",
            (1.0, d.wStar),
            (d.h, d.one),
            (-2.0 * d.h, O.product(d.eta, d.uStar)),
        )
        angular_quadratic = O.product(d.d, d.normalizedGradient)
        average_coefficient = _linear_combination(
            e,
            "pinned averageCoefficient",
            (2.0 * d.D, d.eta),
        )
        angular_slow = _linear_combination(
            e,
            "pinned angularSlowCoefficient",
            (2.0 * d.h, d.eta),
        )
        axial_linear = _linear_combination(
            e,
            "pinned axialLinearCoefficient",
            (d.A, d.one),
            (-4.0 * d.A, O.product(d.eta, d.uStar)),
            (1.0, O.product(d.d, d.uStarEta)),
        )
        axial_quadratic = _linear_combination(
            e,
            "pinned axialQuadraticCoefficient",
            (2.0 * d.A, d.eta),
        )

        lin1 = _linear_combination(
            e,
            "pinned naturalRemainder.lin1",
            (1.0, O.j2(O.product(angular_linear, phi))),
            (1.0, O.dot2(d.wStar, phi)),
            (1.0, O.param2(phi, d.hStar)),
        )
        quad1 = O.j2(O.product(O.product(angular_quadratic, u), phi))
        slow1 = _linear_combination(
            e,
            "pinned naturalRemainder.slow1",
            (
                1.0,
                O.j2(
                    O.product(
                        _linear_combination(
                            e,
                            "pinned naturalRemainder.slow1.transport",
                            (1.0, O.product(average_coefficient, bu)),
                            (1.0, O.product(angular_slow, u)),
                        ),
                        phi,
                    )
                ),
            ),
            (1.0, O.param2(bu, O.product(d.d, phi))),
            (1.0, O.dot2(O.product(average_coefficient, bu), phi)),
            (1.0, O.product(d.d, O.mixed2(bu, phi))),
            (-1.0, O.param2(phi, O.product(d.d, u))),
        )

        lin2 = _linear_combination(
            e,
            "pinned naturalRemainder.lin2",
            (1.0, O.j1(O.product(axial_linear, u))),
            (1.0, O.dot1(d.wStar, u)),
            (1.0, O.param1(u, d.hStar)),
        )
        slow2 = _linear_combination(
            e,
            "pinned naturalRemainder.slow2",
            (1.0, O.j1(O.product(axial_quadratic, O.product(u, u)))),
            (1.0, O.dot1(O.product(average_coefficient, bu), u)),
            (1.0, O.product(d.d, O.mixed1(bu, u))),
            (-1.0, O.param1(u, O.product(d.d, u))),
        )

        source = O.product(O.product(a, a), O.product(phi, phi))
        pressure = O.j1(
            _linear_combination(
                e,
                "pinned naturalRemainder.pressure",
                (-1.0, O.product(_linear_combination(e, "4A eta", (4.0 * d.A, d.eta)), O.primitive(source))),
                (1.0, O.product(d.d, O.parameterPrimitive(source))),
                (-1.0, O.product(_linear_combination(e, "2 eta", (2.0, d.eta)), O.mulY(source))),
            )
        )

        first_input = O.product(
            d.inverseL,
            _linear_combination(
                e,
                "pinned naturalRemainder.angular input",
                (1.0, lin1),
                (1.0, quad1),
                (-t, slow1),
            ),
        )
        second = O.product(
            d.inverseL,
            _linear_combination(
                e,
                "pinned naturalRemainder.axial input",
                (1.0, lin2),
                (-t, slow2),
                (1.0, pressure),
            ),
        )
        return self.resolvent(first_input), second

    def __call__(
        self,
        t: float,
        a: AxisCoefficientJetState,
        x: tuple[AxisCoefficientJetState, AxisCoefficientJetState],
    ) -> tuple[AxisCoefficientJetState, AxisCoefficientJetState]:
        return self.apply(t, a, x)


def actual_schedule_natural_remainder(
    reference: ActualScheduleReferenceAxisState,
) -> AxisCoefficientNaturalRemainder:
    """Assemble the pinned remainder composition from the actual reference state."""

    if not isinstance(reference, ActualScheduleReferenceAxisState):
        raise TypeError("reference must be ActualScheduleReferenceAxisState")
    return AxisCoefficientNaturalRemainder(
        epsilon=reference.epsilon,
        operators=actual_schedule_coefficient_operators(reference),
        data=actual_schedule_axis_coefficient_data(reference),
        resolvent=actual_schedule_natural_resolvent(reference),
    )
