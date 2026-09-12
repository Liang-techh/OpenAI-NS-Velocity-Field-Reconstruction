"""Wide/log pressure component of pinned ``naturalRemainder`` at ``x0``.

The landed wide natural source carries the genuine theorem-selected
``a^2 * phi0^2`` source without narrowing its common ``a^2`` scale into
binary64.  This module advances exactly one downstream seam: the linear
pressure chain from ``AxisContraction.naturalRemainder``

    pressure = j1 (
      -(4 A eta) * primitive(source)
      + d * parameterPrimitive(source)
      -(2 eta) * mulY(source))

at the actual SchedulePressure reference pair.

Every operation is evaluated coefficientwise after factoring out the same
nonzero ``a(eta)^2`` scale.  The normalized factors therefore remain Decimal
quantities while the common amplitude magnitude stays in signed-log form.
No amplitude surrogate, caller-supplied Lambda/C, fitted derivative table, or
radial cutoff is accepted.

This is still formal-structure only.  It does not yet combine the wide pressure
with the ordinary axial remainder, instantiate the theorem-selected
``t = 1/Lambda`` inside that mixed-scale sum, or materialize a Picard iterate.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math
from typing import Callable

from .axis_coefficient_amplitude import SignedLogCoefficientJet
from .axis_coefficient_reference_state import (
    AxisCoefficientJetState,
    WINDOW_LEFT,
    WINDOW_RIGHT,
)
from .axis_coefficient_wide_natural_source import (
    ActualScheduleReferenceNaturalSourceLogState,
    actual_schedule_reference_natural_source_log_state,
)
from .outgoing_tail import TailData


_DECIMAL_PRECISION = 96
_WideProvider = Callable[[int, int, float], Decimal]


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
class ActualScheduleReferenceWidePressureLogState:
    """The exact pressure-chain contribution with common ``a^2`` factored out.

    For each coefficient jet, ``normalized_factor(n,m,eta)`` satisfies

      d_eta^m pressure_n(eta)
        = a(eta)^2 * normalized_factor(n,m,eta).

    This definition is valid because all operations in the pressure chain are
    linear in the wide source except multiplication by ordinary coefficient
    fields.  Their eta derivatives are handled by the exact finite Leibniz rule.
    """

    source: ActualScheduleReferenceNaturalSourceLogState

    def __post_init__(self) -> None:
        if not isinstance(self.source, ActualScheduleReferenceNaturalSourceLogState):
            raise TypeError("source must be ActualScheduleReferenceNaturalSourceLogState")

    @property
    def epsilon(self) -> float:
        return self.source.epsilon

    @property
    def amplitude(self):
        return self.source.amplitude

    @property
    def data(self):
        return self.source.amplitude.data

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def pressure_chain_complete(self) -> bool:
        return True

    def _ordinary_product_factor(
        self,
        ordinary: AxisCoefficientJetState,
        wide: _WideProvider,
        n: int,
        m: int,
        eta: float,
    ) -> Decimal:
        """Return the normalized product of one ordinary and one wide state."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        if not isinstance(ordinary, AxisCoefficientJetState):
            raise TypeError("ordinary factor must be AxisCoefficientJetState")
        if ordinary.epsilon != self.epsilon:
            raise ValueError("ordinary factor must use the actual source epsilon")

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            total = Decimal(0)
            for i in range(n + 1):
                j = n - i
                for k in range(m + 1):
                    total += (
                        Decimal(math.comb(m, k))
                        * wide(i, k, eta)
                        * _decimal_from_float(
                            ordinary.jet(j, m - k, eta),
                            "ordinary coefficient jet",
                        )
                    )
            return +total

    def primitive_factor(self, n: int, m: int, eta: float) -> Decimal:
        """Normalized pinned radial primitive of the wide source."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        if n == 0:
            self.source.normalized_factor(0, m, eta)
            return Decimal(0)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(self.source.normalized_factor(n - 1, m, eta) / Decimal(n))

    def parameter_primitive_factor(self, n: int, m: int, eta: float) -> Decimal:
        """Normalized pinned ``primitive(partial_eta source)`` coefficient jet."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        if n == 0:
            self.source.normalized_factor(0, m + 1, eta)
            return Decimal(0)
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(self.source.normalized_factor(n - 1, m + 1, eta) / Decimal(n))

    def multiply_y_factor(self, n: int, m: int, eta: float) -> Decimal:
        """Normalized pinned predecessor-row map for multiplication by ``Y``."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        if n == 0:
            self.source.normalized_factor(0, m, eta)
            return Decimal(0)
        return self.source.normalized_factor(n - 1, m, eta)

    def pressure_input_factor(self, n: int, m: int, eta: float) -> Decimal:
        """Factor before the final pinned ``j1`` regular inverse."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        d = self.data

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            eta_primitive = self._ordinary_product_factor(
                d.eta, self.primitive_factor, n, m, eta
            )
            d_parameter_primitive = self._ordinary_product_factor(
                d.d, self.parameter_primitive_factor, n, m, eta
            )
            eta_multiply_y = self._ordinary_product_factor(
                d.eta, self.multiply_y_factor, n, m, eta
            )
            return +(
                -Decimal(4) * _decimal_from_float(d.A, "A") * eta_primitive
                + d_parameter_primitive
                - Decimal(2) * eta_multiply_y
            )

    def normalized_factor(self, n: int, m: int, eta: float) -> Decimal:
        """Return the normalized coefficient jet after the final pinned ``j1``."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        if n == 0:
            self.pressure_input_factor(0, m, eta)
            return Decimal(0)

        # regularInverse with r=1 maps input row k=n-1 to output row n and
        # divides by radialDivisor(1,k)=(k+1)(k+1)=n^2.
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(
                self.pressure_input_factor(n - 1, m, eta)
                / Decimal(n * n)
            )

    def jet_log(self, n: int, m: int, eta: float) -> SignedLogCoefficientJet:
        """Return one wide pressure jet without narrowing the common scale."""

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
                "actual SchedulePressure theorem-selected natural pressure; "
                "binary64 projection fails closed on range loss"
            ),
            _jet_provider=provider,
        )


def actual_schedule_reference_wide_pressure_log_state(
    data: TailData,
    j: float,
    *,
    phase_samples: int = 4001,
) -> ActualScheduleReferenceWidePressureLogState:
    """Materialize the pinned wide pressure chain at ``x0 = referencePair``."""

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    source = actual_schedule_reference_natural_source_log_state(
        data,
        j,
        phase_samples=phase_samples,
    )
    return ActualScheduleReferenceWidePressureLogState(source=source)
