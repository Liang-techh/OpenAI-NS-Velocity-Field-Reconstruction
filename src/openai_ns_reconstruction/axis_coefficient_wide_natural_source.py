"""Wide/log coefficient state for the theorem-selected natural source at x0.

The pinned ``AxisContraction.naturalRemainder`` uses the pressure source

    source = a^2 * phi^2,

where ``a`` is ``NaturalAxisCoefficients.realAmplitude`` and ``phi`` is the
angular fixed-point field.  At the contraction center ``x0 = referencePair``
the angular field is already materialized as actual SchedulePressure coefficient
jets, while the theorem-selected ``a`` is far below binary64 on the current
conservative scale.

This module bridges exactly that first nonlinear amplitude seam without
inventing a surrogate amplitude.  It keeps the common ``a(eta)^2`` magnitude in
signed-log form and carries every eta derivative in a separate Decimal factor.
Because ``a`` has radial degree zero, radial convolution with ``a^2`` leaves the
radial row of ``phi0^2`` unchanged.  Parameter derivatives are evaluated by the
exact Leibniz rule together with the complete Bell recurrence for
``d^k(a^2) / a^2``.

The result is still formal-structure only.  It does not yet push this wide state
through ``primitive``/``parameterPrimitive``/``mulY``/``j1``, does not combine
the resulting pressure term with ordinary binary64 remainder terms, and does
not instantiate the theorem-selected ``t = 1/Lambda`` or the Picard step.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math

from .axis_coefficient_amplitude import (
    ActualScheduleAmplitudeLogState,
    SignedLogCoefficientJet,
    actual_schedule_amplitude_log_state,
)
from .axis_coefficient_operators import actual_schedule_coefficient_operators
from .axis_coefficient_reference_state import (
    AxisCoefficientJetState,
    WINDOW_LEFT,
    WINDOW_RIGHT,
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


@dataclass(frozen=True)
class ActualScheduleReferenceNaturalSourceLogState:
    """The exact coefficientwise ``a^2 * phi0^2`` source in split-log form.

    ``normalized_factor(n,m,eta)`` is defined by

      d_eta^m [a(eta)^2 * phi0^2_n(eta)]
        = a(eta)^2 * normalized_factor(n,m,eta).

    Thus all enormous/small common amplitude scale is carried once by
    ``2*log(a)`` while the derivative-specific factor remains independently
    representable in Decimal arithmetic.
    """

    amplitude: ActualScheduleAmplitudeLogState
    angular_square: AxisCoefficientJetState

    def __post_init__(self) -> None:
        if not isinstance(self.amplitude, ActualScheduleAmplitudeLogState):
            raise TypeError("amplitude must be ActualScheduleAmplitudeLogState")
        if not isinstance(self.angular_square, AxisCoefficientJetState):
            raise TypeError("angular_square must be AxisCoefficientJetState")
        if self.angular_square.epsilon != self.amplitude.epsilon:
            raise ValueError("source payload must use the actual amplitude epsilon")

    @property
    def epsilon(self) -> float:
        return self.amplitude.epsilon

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def x_is_reference_pair(self) -> bool:
        return True

    @property
    def radial_amplitude_degree_zero(self) -> bool:
        return True

    def _amplitude_square_bell_factor(self, order: int, eta: float) -> Decimal:
        """Return ``(d^order a^2) / a^2`` from the pinned gradient jets."""

        order = _index(order, "order")
        eta = _eta_in_window(eta)
        if order == 0:
            return Decimal(1)

        # If f = log(a), then a^2 = exp(2 f) and
        # (2 f)^(k+1) = 2 * Lambda * d_eta^k(realGradient).
        q: list[Decimal] = [Decimal(0)]
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for derivative_order in range(order):
                gradient_jet = self.amplitude.data.normalizedGradient.jet(
                    0, derivative_order, eta
                )
                q.append(
                    +(
                        Decimal(2)
                        * self.amplitude.Lambda
                        * _decimal_from_float(gradient_jet, "realGradient jet")
                    )
                )

            bell = [Decimal(1)]
            for n in range(order):
                total = Decimal(0)
                for k in range(n + 1):
                    total += (
                        Decimal(math.comb(n, k))
                        * q[k + 1]
                        * bell[n - k]
                    )
                bell.append(+total)
            return bell[order]

    def normalized_factor(self, n: int, m: int, eta: float) -> Decimal:
        """Return the factor after extracting the common ``a(eta)^2`` scale.

        The radial degree-zero property of ``a`` reduces the radial product to
        the corresponding row of ``phi0^2``.  Eta derivatives use the exact
        finite Leibniz sum; no finite-difference derivative table is accepted.
        """

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            total = Decimal(0)
            for k in range(m + 1):
                total += (
                    Decimal(math.comb(m, k))
                    * self._amplitude_square_bell_factor(k, eta)
                    * _decimal_from_float(
                        self.angular_square.jet(n, m - k, eta),
                        "phi0^2 jet",
                    )
                )
            return +total

    def jet_log(self, n: int, m: int, eta: float) -> SignedLogCoefficientJet:
        """Return one source jet without narrowing the theorem amplitude scale."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        factor = self.normalized_factor(n, m, eta)
        if factor == 0:
            return SignedLogCoefficientJet.zero()

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            log_scale = +(Decimal(2) * self.amplitude.log_amplitude(eta))
            log_factor = +abs(factor).ln()
        return SignedLogCoefficientJet(
            sign=1 if factor > 0 else -1,
            log_scale=log_scale,
            log_factor=log_factor,
        )

    def binary64_state(self) -> AxisCoefficientJetState:
        """Fail-closed projection into the legacy coefficient backend."""

        def provider(n: int, m: int, eta: float) -> float:
            return self.jet_log(n, m, eta).to_binary64()

        return AxisCoefficientJetState(
            epsilon=self.epsilon,
            origin=(
                "actual SchedulePressure theorem-selected natural source "
                "a^2*phi0^2; binary64 projection fails closed on range loss"
            ),
            _jet_provider=provider,
        )


def actual_schedule_reference_natural_source_log_state(
    data: TailData,
    j: float,
    *,
    phase_samples: int = 4001,
) -> ActualScheduleReferenceNaturalSourceLogState:
    """Materialize the pinned ``source`` at ``x0 = referencePair``.

    Callers provide only the already-required actual schedule input ``TailData``
    and ``j``.  Pressure/sigma/epsilon/Lambda/C, the angular reference field,
    and the coefficient product are all inherited from the landed theorem chain.
    """

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    amplitude = actual_schedule_amplitude_log_state(
        data,
        j,
        phase_samples=phase_samples,
    )
    operators = actual_schedule_coefficient_operators(amplitude.reference)
    angular_square = operators.product(
        amplitude.reference.phi,
        amplitude.reference.phi,
    )
    return ActualScheduleReferenceNaturalSourceLogState(
        amplitude=amplitude,
        angular_square=angular_square,
    )
