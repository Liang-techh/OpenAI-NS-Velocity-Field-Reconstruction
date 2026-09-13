"""Wide/log-domain theorem-selected natural-axis amplitude coefficient state.

The pinned natural-axis construction uses

    a(eta) = exp(Lambda * realPhase(eta)) / C

with Lambda and C chosen only after the actual SchedulePressure analytic and
contraction bounds are available.  On the current certified chain those scalar
choices are finite mathematical reals but can lie far outside binary64.  This
module therefore materializes the amplitude in signed-log coefficient-jet
coordinates instead of narrowing Lambda or C to float.

Only the scalar representation changes.  The default phase value is the
nearest Decimal midpoint of a validated exact-rational kernel interval at the
selected finite tolerance.  Its returned point is still numerical, and a
large ``Lambda`` can amplify the phase error substantially; callers needing
that uncertainty should use the explicit log-amplitude enclosure.  The
parameter derivatives use (realAmplitude)' = Lambda * realGradient *
realAmplitude and the complete Bell recurrence, while radial degree is exactly
zero.  ``binary64_state`` is an explicit projection for callers that need the
existing coefficient backend; it fails closed whenever a nonzero jet is outside
binary64's representable range rather than silently replacing it by zero or
infinity.  No path here claims a paper-exact parameter selection, global norm,
or total reconstruction.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
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
from .axis_coefficient_rational_data import RationalAxisCoefficientData
from .axis_phase_integral import PhaseIntegralResult, validated_phase_integral
from .axis_phase_log_enclosure import RationalLogAmplitudeEnclosure
from .natural_axis import real_phase
from .natural_scale_selection_wide import WideNaturalScaleSelection
from .outgoing_tail import TailData
from .stage1_scale_chain_wide import diagnose_actual_schedule_scale_chain_wide


_DECIMAL_PRECISION = 96
_PHASE_INITIAL_ORDER = 16
_PHASE_MAX_ORDER = 4096
_PHASE_MAX_CELLS = 4096
_PHASE_MAX_DEPTH = 128


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


@lru_cache(maxsize=16)
def _cached_phase_integral(
    h: Fraction,
    j: Fraction,
    sigma: Fraction,
    eta: Fraction,
    absolute_tolerance: Fraction,
    initial_order: int,
    max_order: int,
    max_cells: int,
    max_depth: int,
) -> PhaseIntegralResult:
    """Cache only exact scalar inputs and finite integrator settings."""

    return validated_phase_integral(
        h,
        j,
        sigma,
        eta,
        absolute_tolerance=absolute_tolerance,
        initial_order=initial_order,
        max_order=max_order,
        max_cells=max_cells,
        max_depth=max_depth,
    )


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
    phase_absolute_tolerance: Fraction = Fraction(1, 10**12)
    rational_data: RationalAxisCoefficientData = field(init=False, repr=False)

    def __post_init__(self) -> None:
        samples = int(self.phase_samples)
        if samples < 3:
            raise ValueError("phase_samples must be at least 3")
        if samples % 2 == 0:
            samples += 1
        if self.data.epsilon != self.reference.epsilon:
            raise ValueError("amplitude data/reference epsilon mismatch")
        if not isinstance(self.phase_absolute_tolerance, Fraction):
            raise TypeError("phase_absolute_tolerance must be a Fraction")
        if self.phase_absolute_tolerance <= 0:
            raise ValueError("phase_absolute_tolerance must be positive")
        _finite_decimal(self.scale.Lambda, "scale.Lambda")
        _finite_decimal(self.scale.C.exponent_upper, "scale.C.exponent_upper")
        object.__setattr__(self, "phase_samples", samples)
        object.__setattr__(self, "rational_data", RationalAxisCoefficientData(self.data))

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
    def phase_has_error_enclosure(self) -> bool:
        return True

    @property
    def radial_degree_zero(self) -> bool:
        return True

    def log_amplitude(self, eta: float) -> Decimal:
        """Return the Decimal96 midpoint of the default validated log interval."""

        eta = float(eta)
        if not math.isfinite(eta):
            raise ValueError("eta must be finite")
        if eta == 0.0:
            return self.C_exponent.copy_negate()

        return self.default_log_amplitude_enclosure(eta).midpoint_decimal()

    def phase_enclosure(self, eta: float) -> PhaseIntegralResult:
        """Return the cached default exact-rational phase integral interval."""

        eta = float(eta)
        if not math.isfinite(eta):
            raise ValueError("eta must be finite")
        eta_fraction = Fraction.from_float(eta)
        return _cached_phase_integral(
            self.rational_data.h,
            self.rational_data.j,
            self.rational_data.sigma,
            eta_fraction,
            self.phase_absolute_tolerance,
            _PHASE_INITIAL_ORDER,
            _PHASE_MAX_ORDER,
            _PHASE_MAX_CELLS,
            _PHASE_MAX_DEPTH,
        )

    def default_log_amplitude_enclosure(
        self,
        eta: float,
    ) -> RationalLogAmplitudeEnclosure:
        """Return the default validated interval transported to log amplitude."""

        phase = self.phase_enclosure(eta)
        return RationalLogAmplitudeEnclosure(
            phase.lower,
            phase.upper,
            self.Lambda,
            self.C_exponent,
        )

    def legacy_log_amplitude_diagnostic(self, eta: float) -> Decimal:
        """Evaluate the former 4001-point phase path for diagnostics only."""

        eta = float(eta)
        if not math.isfinite(eta):
            raise ValueError("eta must be finite")
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

    def log_amplitude_enclosure(
        self,
        eta: Fraction,
        *,
        absolute_log_tolerance: Fraction,
        initial_order: int = 16,
        max_order: int = 4096,
        max_cells: int = 4096,
        max_depth: int = 128,
    ) -> RationalLogAmplitudeEnclosure:
        """Enclose ``log(a(eta))`` from a validated exact phase interval.

        The selected finite ``Lambda`` and ``C`` exponent are used exactly as
        the affine map ``Lambda * phase - log(C)``.  The supplied positive
        ``absolute_log_tolerance`` is the requested half-width in log space;
        the phase integrator receives the exact tolerance divided by
        ``Lambda``.  This remains conditional on the selected finite
        parameters and on the integrator's phase enclosure, rather than a
        proof of the scale selection or of the full reconstruction.
        """

        if not isinstance(eta, Fraction):
            raise TypeError("eta must be a Fraction")
        if not isinstance(absolute_log_tolerance, Fraction):
            raise TypeError("absolute_log_tolerance must be a Fraction")
        if absolute_log_tolerance <= 0:
            raise ValueError("absolute_log_tolerance must be positive")

        # Keep the phase tolerance rational all the way to the validated
        # integrator.  Fraction(Decimal) is exact for the selected Lambda.
        phase_tolerance = absolute_log_tolerance / Fraction(self.Lambda)
        from .axis_phase_integral import validated_phase_integral

        result = validated_phase_integral(
            self.rational_data.h,
            self.rational_data.j,
            self.rational_data.sigma,
            eta,
            absolute_tolerance=phase_tolerance,
            initial_order=initial_order,
            max_order=max_order,
            max_cells=max_cells,
            max_depth=max_depth,
        )
        return RationalLogAmplitudeEnclosure(
            result.lower,
            result.upper,
            self.Lambda,
            self.C_exponent,
        )

    def _bell_factor(self, order: int, eta: float) -> Decimal:
        """Return B_m for d^m exp(f)=exp(f) B_m(f',...,f^(m))."""

        order = _index(order, "order")
        return _finite_decimal(
            self.rational_data.amplitude_power_bell_decimal(
                1,
                order,
                float(eta),
                self.Lambda,
            ),
            "amplitude Bell factor",
        )

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
    phase_absolute_tolerance: Fraction = Fraction(1, 10**12),
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
        phase_absolute_tolerance=phase_absolute_tolerance,
    )
