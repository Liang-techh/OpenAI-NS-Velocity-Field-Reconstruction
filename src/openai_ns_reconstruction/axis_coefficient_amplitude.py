"""Wide/log-domain theorem-selected natural-axis amplitude coefficient state.

The pinned natural-axis construction uses

    a(eta) = exp(Lambda * realPhase(eta)) / C

with Lambda and C chosen only after the actual SchedulePressure analytic and
contraction bounds are available.  On the current certified chain those scalar
choices are finite mathematical reals but can lie far outside binary64.  This
module therefore materializes the amplitude in signed-log coefficient-jet
coordinates instead of narrowing Lambda or C to float.

Only the scalar representation changes.  The parameter derivatives use
(realAmplitude)' = Lambda * realGradient * realAmplitude and the complete Bell
recurrence, while radial degree is exactly zero.  ``binary64_state`` is an
explicit projection for callers that need the existing coefficient backend; it
fails closed whenever a nonzero jet is outside binary64's representable range
rather than silently replacing it by zero or infinity.

The phase value still uses the landed numerical ``real_phase`` quadrature, so
this object remains ``formal-structure`` rather than a paper-exact certificate.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math
import sys

from .axis_coefficient_data import (
    ActualScheduleAxisCoefficientData,
    actual_schedule_axis_coefficient_data,
)
from .axis_coefficient_reference_state import (
    ActualScheduleReferenceAxisState,
    AxisCoefficientJetState,
    actual_schedule_reference_axis_state,
)
from .natural_axis import real_phase
from .natural_scale_selection_wide import WideNaturalScaleSelection
from .outgoing_tail import TailData
from .stage1_scale_chain_wide import diagnose_actual_schedule_scale_chain_wide


_DECIMAL_PRECISION = 96


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _finite_decimal(value: Decimal, name: str) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError(f"{name} must be a finite Decimal")
    return value


def _decimal_from_float(value: float, name: str) -> Decimal:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return Decimal.from_float(value)


def _sum_precision(a: Decimal, b: Decimal) -> int:
    """Precision sufficient to retain a small log correction beside a huge one."""

    if a == 0 or b == 0:
        gap = 0
    else:
        gap = abs(a.adjusted() - b.adjusted())
    digits = max(len(a.as_tuple().digits), len(b.as_tuple().digits))
    return max(_DECIMAL_PRECISION, gap + digits + 16)


with localcontext() as _ctx:
    _ctx.prec = _DECIMAL_PRECISION
    _BINARY64_MAX_LOG = Decimal.from_float(sys.float_info.max).ln()
    _BINARY64_MIN_SUBNORMAL_LOG = Decimal.from_float(math.ulp(0.0)).ln()


@dataclass(frozen=True)
class SignedLogCoefficientJet:
    """One real coefficient jet with a lossless split logarithmic magnitude.

    ``log_scale`` is the common amplitude logarithm and ``log_factor`` is the
    derivative-specific Bell-factor logarithm.  Keeping them separate matters
    on the current theorem scale: ``log_scale`` is about 10^784 in magnitude,
    while a derivative correction can be only about 10^3.  A fixed-precision
    Decimal sum would erase that correction even though both inputs are known.
    """

    sign: int
    log_scale: Decimal | None
    log_factor: Decimal | None

    def __post_init__(self) -> None:
        if self.sign not in (-1, 0, 1):
            raise ValueError("sign must be -1, 0, or 1")
        if self.sign == 0:
            if self.log_scale is not None or self.log_factor is not None:
                raise ValueError("zero jet must not carry logarithmic components")
            return
        if self.log_scale is None or self.log_factor is None:
            raise ValueError("nonzero jet requires both logarithmic components")
        _finite_decimal(self.log_scale, "log_scale")
        _finite_decimal(self.log_factor, "log_factor")

    @classmethod
    def zero(cls) -> "SignedLogCoefficientJet":
        return cls(sign=0, log_scale=None, log_factor=None)

    @property
    def log_abs(self) -> Decimal | None:
        """Return log(abs(value)) without losing the small factor correction."""

        if self.sign == 0:
            return None
        assert self.log_scale is not None
        assert self.log_factor is not None
        with localcontext() as ctx:
            ctx.prec = _sum_precision(self.log_scale, self.log_factor)
            return +(self.log_scale + self.log_factor)

    def to_binary64(self) -> float:
        """Project to a finite binary64 value, failing closed on range loss."""

        if self.sign == 0:
            return 0.0
        log_abs = self.log_abs
        assert log_abs is not None
        if log_abs > _BINARY64_MAX_LOG:
            raise ArithmeticError("nonzero amplitude jet overflows binary64")
        if log_abs < _BINARY64_MIN_SUBNORMAL_LOG:
            raise ArithmeticError("nonzero amplitude jet underflows binary64")
        magnitude = math.exp(float(log_abs))
        if magnitude == 0.0:
            raise ArithmeticError("nonzero amplitude jet underflows binary64")
        if not math.isfinite(magnitude):
            raise ArithmeticError("nonzero amplitude jet overflows binary64")
        return math.copysign(magnitude, float(self.sign))


@dataclass(frozen=True)
class ActualScheduleAmplitudeLogState:
    """The theorem-selected realAmplitude as a radial-degree-zero log state."""

    reference: ActualScheduleReferenceAxisState
    data: ActualScheduleAxisCoefficientData
    scale: WideNaturalScaleSelection
    phase_samples: int = 4001

    def __post_init__(self) -> None:
        samples = int(self.phase_samples)
        if samples < 3:
            raise ValueError("phase_samples must be at least 3")
        if samples % 2 == 0:
            samples += 1
        if self.data.epsilon != self.reference.epsilon:
            raise ValueError("amplitude data/reference epsilon mismatch")
        _finite_decimal(self.scale.Lambda, "scale.Lambda")
        _finite_decimal(self.scale.C.exponent_upper, "scale.C.exponent_upper")
        object.__setattr__(self, "phase_samples", samples)

    @property
    def epsilon(self) -> float:
        return self.reference.epsilon

    @property
    def Lambda(self) -> Decimal:
        return self.scale.Lambda

    @property
    def C_exponent(self) -> Decimal:
        """Return ``log C`` for the pinned symbolic C selection."""

        return self.scale.C.exponent_upper

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def phase_is_numerically_evaluated(self) -> bool:
        return True

    @property
    def radial_degree_zero(self) -> bool:
        return True

    def log_amplitude(self, eta: float) -> Decimal:
        """Evaluate ``Lambda*realPhase(eta) - log(C)`` without forming C."""

        eta = float(eta)
        if not math.isfinite(eta):
            raise ValueError("eta must be finite")
        # The landed real_phase implementation returns 0 exactly here. Preserve
        # the symbolic normalization exponent without passing a no-op through
        # the working Decimal precision.
        if eta == 0.0:
            return self.C_exponent.copy_negate()
        phase = real_phase(
            self.data.h,
            self.data.j,
            self.data.sigma,
            eta,
            samples=self.phase_samples,
        )
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            return +(self.Lambda * _decimal_from_float(phase, "realPhase") - self.C_exponent)

    def _bell_factor(self, order: int, eta: float) -> Decimal:
        """Return B_m for d^m exp(f)=exp(f) B_m(f',...,f^(m))."""

        order = _index(order, "order")
        if order == 0:
            return Decimal(1)

        # q[k] = f^(k) for k>=1, where f' = Lambda * realGradient.
        q: list[Decimal] = [Decimal(0)]
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            for derivative_order in range(order):
                gradient_jet = self.data.normalizedGradient.jet(
                    0, derivative_order, float(eta)
                )
                q.append(
                    +(self.Lambda * _decimal_from_float(gradient_jet, "realGradient jet"))
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

    def jet_log(self, n: int, m: int, eta: float) -> SignedLogCoefficientJet:
        """Return the (n,m) amplitude jet without narrowing its magnitude."""

        n = _index(n, "n")
        m = _index(m, "m")
        if n != 0:
            return SignedLogCoefficientJet.zero()

        factor = self._bell_factor(m, float(eta))
        if factor == 0:
            return SignedLogCoefficientJet.zero()

        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            log_factor = +abs(factor).ln()
        return SignedLogCoefficientJet(
            sign=1 if factor > 0 else -1,
            log_scale=self.log_amplitude(float(eta)),
            log_factor=log_factor,
        )

    def binary64_state(self) -> AxisCoefficientJetState:
        """Return a fail-closed projection into the existing binary64 backend.

        No underflowing nonzero amplitude coefficient is silently turned into
        zero.  Consequently this adapter may raise during jet evaluation on the
        current theorem-selected scale; that is the representation boundary the
        next coefficient-backend step must remove.
        """

        def provider(n: int, m: int, eta: float) -> float:
            return self.jet_log(n, m, eta).to_binary64()

        return AxisCoefficientJetState(
            epsilon=self.epsilon,
            origin=(
                "actual SchedulePressure theorem-selected realAmplitude; "
                "binary64 projection fails closed on range loss"
            ),
            _jet_provider=provider,
        )


def actual_schedule_amplitude_log_state(
    data: TailData,
    j: float,
    *,
    phase_samples: int = 4001,
) -> ActualScheduleAmplitudeLogState:
    """Bind realAmplitude to the actual schedule and pinned wide Lambda/C choice.

    Callers cannot replace sigma, epsilon, Lambda, C, pressure data, or the
    coefficient table.  If the actual-schedule wide chain does not reach the
    theorem's scale selection this factory fails closed.
    """

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    reference = actual_schedule_reference_axis_state(data, j)
    diagnostic = diagnose_actual_schedule_scale_chain_wide(data, j)
    if diagnostic.upstream_obstruction is not None:
        raise RuntimeError(diagnostic.upstream_obstruction)
    if diagnostic.scale is None:
        raise RuntimeError("actual-schedule chain did not produce the pinned wide scale")
    axis_data = actual_schedule_axis_coefficient_data(reference)
    return ActualScheduleAmplitudeLogState(
        reference=reference,
        data=axis_data,
        scale=diagnostic.scale,
        phase_samples=phase_samples,
    )
