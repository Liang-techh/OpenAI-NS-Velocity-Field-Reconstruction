"""Mixed-scale axial half of pinned ``naturalRemainder`` at ``x0``.

The actual SchedulePressure theorem chain selects a finite mathematical
``Lambda`` and the amplitude

    a(eta) = exp(Lambda * realPhase(eta)) / C.

On the current conservative theorem bounds, neither ``1/Lambda`` nor the
``a^2`` pressure contribution should be forced through the legacy binary64
coefficient backend.  The landed wide-pressure module already propagates the
nonzero ``a^2 * phi0^2`` source through the pinned pressure chain.  This module
closes exactly the next seam at the contraction centre ``x0 = referencePair``:

    R_u(x0) = inverseL * (lin2 - (1/Lambda) * slow2 + pressure).

The three summands are kept as a *decomposition* rather than rounded into one
floating-point number:

* ``ordinary_base`` is ``inverseL * lin2`` from the existing coefficient
  operators;
* ``inverse_lambda_numerator / Lambda`` is ``-inverseL * slow2 / Lambda`` and
  keeps the theorem-selected Decimal ``Lambda`` as an explicit denominator;
* ``pressure`` is the signed-log ``inverseL * pressure`` coefficient jet with
  the genuine common ``a^2`` scale preserved.

No surrogate amplitude, caller-selected ``Lambda/C``, fitted derivative table,
or radial cutoff is accepted.  This is still a formal-structure bridge: it does
not yet assemble the angular half, apply the outer ``1/(2*Lambda)`` Picard
scale, or materialize ``x1`` / the fixed point.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math
from typing import Iterable

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_data import (
    ActualScheduleAxisCoefficientData,
    actual_schedule_axis_coefficient_data,
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
from .axis_coefficient_wide_natural_pressure import (
    ActualScheduleReferenceWidePressureLogState,
    actual_schedule_reference_wide_pressure_log_state,
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
class MixedScaleAxialCoefficientJet:
    """One coefficient jet of ``R_u(x0)`` without destructive scale collapse.

    The represented real number is exactly the three-term theorem expression

      ordinary_base + inverse_lambda_numerator / Lambda + pressure,

    where ``pressure`` is retained in signed-log form.  The class deliberately
    exposes no method that silently rounds all three scales into binary64.
    """

    ordinary_base: Decimal
    inverse_lambda_numerator: Decimal
    Lambda: Decimal
    pressure: SignedLogCoefficientJet

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
        if not isinstance(self.pressure, SignedLogCoefficientJet):
            raise TypeError("pressure must be SignedLogCoefficientJet")

    def inverse_lambda_term_decimal(self) -> Decimal:
        """Return the slow term at the shared 96-digit Decimal working precision."""

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(self.inverse_lambda_numerator / self.Lambda)


@dataclass(frozen=True)
class ActualScheduleReferenceWideAxialRemainderState:
    """The theorem-selected mixed-scale axial remainder at ``referencePair``."""

    reference: ActualScheduleReferenceAxisState
    operators: AxisCoefficientOperators
    axis_data: ActualScheduleAxisCoefficientData
    wide_pressure: ActualScheduleReferenceWidePressureLogState
    ordinary_base: AxisCoefficientJetState
    ordinary_slow: AxisCoefficientJetState
    Lambda: Decimal

    def __post_init__(self) -> None:
        epsilon = self.reference.epsilon
        if self.operators.epsilon != epsilon:
            raise ValueError("coefficientOperators epsilon mismatch")
        if self.axis_data.epsilon != epsilon:
            raise ValueError("AxisData epsilon mismatch")
        if self.wide_pressure.epsilon != epsilon:
            raise ValueError("wide pressure epsilon mismatch")
        if self.ordinary_base.epsilon != epsilon or self.ordinary_slow.epsilon != epsilon:
            raise ValueError("ordinary axial components must use the reference epsilon")
        if not isinstance(self.Lambda, Decimal) or not self.Lambda.is_finite() or self.Lambda <= 0:
            raise ValueError("Lambda must be a finite positive Decimal")
        if self.Lambda != self.wide_pressure.amplitude.Lambda:
            raise ValueError("Lambda must be the actual theorem-selected amplitude scale")

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
    def mixed_scale_axial_remainder_complete(self) -> bool:
        return True

    def pressure_normalized_factor(self, n: int, m: int, eta: float) -> Decimal:
        """Return ``inverseL * pressure`` after factoring out the common ``a^2``.

        The product uses the exact finite radial convolution and eta Leibniz
        rule.  Since ``wide_pressure.normalized_factor(i,k)`` represents
        ``d_eta^k pressure_i / a^2`` at the same eta, every term shares the
        same nonzero amplitude-square scale and can be accumulated in Decimal.
        """

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        inverse_l = self.axis_data.inverseL
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            total = Decimal(0)
            for i in range(n + 1):
                j = n - i
                for k in range(m + 1):
                    total += (
                        Decimal(math.comb(m, k))
                        * self.wide_pressure.normalized_factor(i, k, eta)
                        * _decimal_from_float(
                            inverse_l.jet(j, m - k, eta),
                            "inverseL coefficient jet",
                        )
                    )
            return +total

    def pressure_jet_log(self, n: int, m: int, eta: float) -> SignedLogCoefficientJet:
        """Return one full ``inverseL * pressure`` jet in signed-log form."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        factor = self.pressure_normalized_factor(n, m, eta)
        if factor == 0:
            return SignedLogCoefficientJet.zero()
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            log_scale = +(Decimal(2) * self.wide_pressure.amplitude.log_amplitude(eta))
            log_factor = +abs(factor).ln()
        return SignedLogCoefficientJet(
            sign=1 if factor > 0 else -1,
            log_scale=log_scale,
            log_factor=log_factor,
        )

    def jet(self, n: int, m: int, eta: float) -> MixedScaleAxialCoefficientJet:
        """Return the exact three-scale decomposition of one axial remainder jet."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        base = _decimal_from_float(self.ordinary_base.jet(n, m, eta), "ordinary base jet")
        slow = _decimal_from_float(self.ordinary_slow.jet(n, m, eta), "ordinary slow jet")
        return MixedScaleAxialCoefficientJet(
            ordinary_base=base,
            inverse_lambda_numerator=-slow,
            Lambda=self.Lambda,
            pressure=self.pressure_jet_log(n, m, eta),
        )


def actual_schedule_reference_wide_axial_remainder_state(
    data: TailData,
    j: float,
    *,
    phase_samples: int = 4001,
) -> ActualScheduleReferenceWideAxialRemainderState:
    """Materialize ``R_u(referencePair)`` as a theorem-selected scale decomposition."""

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")

    reference = actual_schedule_reference_axis_state(data, j)
    operators = actual_schedule_coefficient_operators(reference)
    axis_data = actual_schedule_axis_coefficient_data(reference)
    wide_pressure = actual_schedule_reference_wide_pressure_log_state(
        data,
        j,
        phase_samples=phase_samples,
    )
    if wide_pressure.epsilon != reference.epsilon:
        raise RuntimeError("wide pressure/reference epsilon mismatch")

    O = operators
    d = axis_data
    e = reference.epsilon
    u = reference.u
    bu = O.average(u)

    average_coefficient = _linear_combination(
        e,
        "pinned averageCoefficient",
        ((2.0 * d.D, d.eta),),
    )
    axial_linear = _linear_combination(
        e,
        "pinned axialLinearCoefficient",
        (
            (d.A, d.one),
            (-4.0 * d.A, O.product(d.eta, d.uStar)),
            (1.0, O.product(d.d, d.uStarEta)),
        ),
    )
    axial_quadratic = _linear_combination(
        e,
        "pinned axialQuadraticCoefficient",
        ((2.0 * d.A, d.eta),),
    )

    lin2 = _linear_combination(
        e,
        "pinned naturalRemainder.lin2 at x0",
        (
            (1.0, O.j1(O.product(axial_linear, u))),
            (1.0, O.dot1(d.wStar, u)),
            (1.0, O.param1(u, d.hStar)),
        ),
    )
    slow2 = _linear_combination(
        e,
        "pinned naturalRemainder.slow2 at x0",
        (
            (1.0, O.j1(O.product(axial_quadratic, O.product(u, u)))),
            (1.0, O.dot1(O.product(average_coefficient, bu), u)),
            (1.0, O.product(d.d, O.mixed1(bu, u))),
            (-1.0, O.param1(u, O.product(d.d, u))),
        ),
    )

    ordinary_base = O.product(d.inverseL, lin2)
    ordinary_slow = O.product(d.inverseL, slow2)
    Lambda = wide_pressure.amplitude.Lambda
    return ActualScheduleReferenceWideAxialRemainderState(
        reference=reference,
        operators=operators,
        axis_data=axis_data,
        wide_pressure=wide_pressure,
        ordinary_base=ordinary_base,
        ordinary_slow=ordinary_slow,
        Lambda=Lambda,
    )
