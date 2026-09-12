"""Mixed-scale angular half of pinned ``naturalRemainder`` at ``x0``.

The theorem-selected scale satisfies ``t = 1/Lambda``.  On the current
SchedulePressure bounds ``Lambda`` is far outside binary64, so forming ``t`` as
``float`` would silently erase a mathematically nonzero term.  At the actual
contraction centre ``x0 = referencePair`` the angular branch is

    R_phi(x0) = naturalResolvent(
        inverseL * (lin1 + quad1 - slow1 / Lambda)
    ).

Both multiplication by ``inverseL`` and ``naturalResolvent`` are linear.
Therefore this module evaluates the two ordinary numerator states separately
and retains the theorem-selected denominator explicitly:

    naturalResolvent(inverseL * (lin1 + quad1))
      - naturalResolvent(inverseL * slow1) / Lambda.

This is an algebraic decomposition of the pinned formula, not a numerical
series approximation.  The landed natural resolvent is coefficientwise exact
because radial filtration makes each requested row a finite sum.  No surrogate
Lambda, fitted coefficient table, arbitrary cutoff, or float underflow of
``1/Lambda`` is permitted.

The result remains ``formal-structure``: this module closes only the angular
mixed-scale seam.  It does not combine the angular and axial branches into a
complete ``naturalRemainder(x0)``, apply the outer ``1/(2*Lambda)`` Picard
factor, or certify the global weighted ``AxisSpace`` fixed point.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math
from typing import Iterable

from .axis_coefficient_amplitude import actual_schedule_amplitude_log_state
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
    WINDOW_LEFT,
    WINDOW_RIGHT,
    actual_schedule_reference_axis_state,
)
from .outgoing_tail import TailData


_DECIMAL_PRECISION = 96


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _eta_in_window(value: float) -> float:
    eta = float(value)
    if not math.isfinite(eta) or not WINDOW_LEFT <= eta <= WINDOW_RIGHT:
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]")
    return eta


def _decimal_from_float(value: float, name: str) -> Decimal:
    value = float(value)
    if not math.isfinite(value):
        raise ArithmeticError(f"{name} must be finite")
    return Decimal.from_float(value)


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
    terms: Iterable[tuple[float, AxisCoefficientJetState]],
) -> AxisCoefficientJetState:
    checked: list[tuple[float, AxisCoefficientJetState]] = []
    for coefficient, state in terms:
        coefficient = float(coefficient)
        if not math.isfinite(coefficient):
            raise ValueError("linear-combination coefficients must be finite")
        checked.append((coefficient, _checked_state(state, epsilon, origin)))

    def provider(n: int, m: int, eta: float) -> float:
        return sum(coefficient * state.jet(n, m, eta) for coefficient, state in checked)

    return AxisCoefficientJetState(
        epsilon=epsilon,
        origin=origin,
        _jet_provider=provider,
    )


@dataclass(frozen=True)
class MixedScaleAngularCoefficientJet:
    """One angular remainder jet with the nonzero ``1/Lambda`` scale explicit."""

    ordinary_base: Decimal
    inverse_lambda_numerator: Decimal
    Lambda: Decimal

    def __post_init__(self) -> None:
        for name, value in (
            ("ordinary_base", self.ordinary_base),
            ("inverse_lambda_numerator", self.inverse_lambda_numerator),
            ("Lambda", self.Lambda),
        ):
            if not isinstance(value, Decimal) or not value.is_finite():
                raise ValueError(f"{name} must be a finite Decimal")
        if self.Lambda <= 0:
            raise ValueError("Lambda must be positive")

    def inverse_lambda_term_decimal(self) -> Decimal:
        """Return the slow contribution at the shared 96-digit Decimal precision."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(self.inverse_lambda_numerator / self.Lambda)


@dataclass(frozen=True)
class ActualScheduleReferenceWideAngularRemainderState:
    """The theorem-selected mixed-scale angular remainder at ``referencePair``."""

    reference: ActualScheduleReferenceAxisState
    operators: AxisCoefficientOperators
    axis_data: ActualScheduleAxisCoefficientData
    resolvent: AxisCoefficientNaturalResolvent
    ordinary_base: AxisCoefficientJetState
    ordinary_slow: AxisCoefficientJetState
    Lambda: Decimal

    def __post_init__(self) -> None:
        epsilon = self.reference.epsilon
        if self.operators.epsilon != epsilon:
            raise ValueError("coefficientOperators epsilon mismatch")
        if self.axis_data.epsilon != epsilon:
            raise ValueError("AxisData epsilon mismatch")
        if self.resolvent.epsilon != epsilon:
            raise ValueError("naturalResolvent epsilon mismatch")
        if self.ordinary_base.epsilon != epsilon or self.ordinary_slow.epsilon != epsilon:
            raise ValueError("ordinary angular components must use the reference epsilon")
        if not isinstance(self.Lambda, Decimal) or not self.Lambda.is_finite() or self.Lambda <= 0:
            raise ValueError("Lambda must be a finite positive Decimal")

    @property
    def epsilon(self) -> float:
        return self.reference.epsilon

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def mixed_scale_angular_remainder_complete(self) -> bool:
        return True

    def jet(self, n: int, m: int, eta: float) -> MixedScaleAngularCoefficientJet:
        """Return one exact two-scale decomposition of the angular remainder jet."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        base = _decimal_from_float(self.ordinary_base.jet(n, m, eta), "ordinary base jet")
        slow = _decimal_from_float(self.ordinary_slow.jet(n, m, eta), "ordinary slow jet")
        return MixedScaleAngularCoefficientJet(
            ordinary_base=base,
            inverse_lambda_numerator=-slow,
            Lambda=self.Lambda,
        )


def actual_schedule_reference_wide_angular_remainder_state(
    data: TailData,
    j: float,
    *,
    phase_samples: int = 4001,
) -> ActualScheduleReferenceWideAngularRemainderState:
    """Materialize the angular ``naturalRemainder(referencePair)`` scale split."""

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")

    reference = actual_schedule_reference_axis_state(data, j)
    operators = actual_schedule_coefficient_operators(reference)
    axis_data = actual_schedule_axis_coefficient_data(reference)
    resolvent = actual_schedule_natural_resolvent(reference)
    amplitude = actual_schedule_amplitude_log_state(
        data,
        j,
        phase_samples=phase_samples,
    )
    if amplitude.epsilon != reference.epsilon:
        raise RuntimeError("amplitude/reference epsilon mismatch")

    O = operators
    d = axis_data
    e = reference.epsilon
    phi = reference.phi
    u = reference.u
    bu = O.average(u)

    angular_linear = _linear_combination(
        e,
        "pinned angularLinearCoefficient",
        (
            (1.0, d.wStar),
            (d.h, d.one),
            (-2.0 * d.h, O.product(d.eta, d.uStar)),
        ),
    )
    angular_quadratic = O.product(d.d, d.normalizedGradient)
    average_coefficient = _linear_combination(
        e,
        "pinned averageCoefficient",
        ((2.0 * d.D, d.eta),),
    )
    angular_slow = _linear_combination(
        e,
        "pinned angularSlowCoefficient",
        ((2.0 * d.h, d.eta),),
    )

    lin1 = _linear_combination(
        e,
        "pinned naturalRemainder.lin1 at x0",
        (
            (1.0, O.j2(O.product(angular_linear, phi))),
            (1.0, O.dot2(d.wStar, phi)),
            (1.0, O.param2(phi, d.hStar)),
        ),
    )
    quad1 = O.j2(O.product(O.product(angular_quadratic, u), phi))
    slow1_transport = _linear_combination(
        e,
        "pinned naturalRemainder.slow1.transport at x0",
        (
            (1.0, O.product(average_coefficient, bu)),
            (1.0, O.product(angular_slow, u)),
        ),
    )
    slow1 = _linear_combination(
        e,
        "pinned naturalRemainder.slow1 at x0",
        (
            (1.0, O.j2(O.product(slow1_transport, phi))),
            (1.0, O.param2(bu, O.product(d.d, phi))),
            (1.0, O.dot2(O.product(average_coefficient, bu), phi)),
            (1.0, O.product(d.d, O.mixed2(bu, phi))),
            (-1.0, O.param2(phi, O.product(d.d, u))),
        ),
    )

    base_input = O.product(
        d.inverseL,
        _linear_combination(
            e,
            "pinned naturalRemainder.angular base input at x0",
            ((1.0, lin1), (1.0, quad1)),
        ),
    )
    slow_input = O.product(d.inverseL, slow1)
    ordinary_base = resolvent(base_input)
    ordinary_slow = resolvent(slow_input)

    return ActualScheduleReferenceWideAngularRemainderState(
        reference=reference,
        operators=operators,
        axis_data=axis_data,
        resolvent=resolvent,
        ordinary_base=ordinary_base,
        ordinary_slow=ordinary_slow,
        Lambda=amplitude.Lambda,
    )
