"""Executable scalar core of the pinned ``OutgoingTail.finalAngular`` schedule.

This module follows the formulas in the pinned official Lean files
``OutgoingSchedule.lean``, ``OutgoingTail.lean`` and ``SchedulePressure.lean``.
It instantiates the otherwise opaque step-derivative choice with the already
proved repository witness ``S = 32``. Consequently the resulting tail is a
constructive theorem-admissible realization of the printed formulas; it is not
claimed to be definitionally equal to Lean's particular ``Classical.choose``.

Only the scalar outgoing angular field and its clock weight are materialized
here. The all-real-line pressure integral, the low-|Z| uniform margin and the
later Lambda/C/fixed-point choices remain separate certification tasks.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
import math

from .quadrature import integrate
from .schedule_pressure import (
    EXPLICIT_FLATTEN_LENGTH,
    EXPLICIT_STEP_BOUND,
    outgoing_sigma,
    outgoing_sigma_derivative,
    schedule_shape_exponent,
)

_SIGMA_QUADRATURE_ORDER = 64
_TRANSITION_QUADRATURE_ORDER = 96
EXPLICIT_TAIL_COEFFICIENT = math.exp(-5.0) / (16.0 * (EXPLICIT_STEP_BOUND + 1.0))


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _sigma_primitive(x: float) -> float:
    """Integral of the pinned smooth step from 0 to ``x``.

    Symmetry gives integral_0^1 sigma = 1/2, so only the unit transition needs
    numerical quadrature. Values outside it use the exact plateau formulas.
    """
    x = _finite(x, "x")
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return x - 0.5
    return integrate(outgoing_sigma, 0.0, x, n=_SIGMA_QUADRATURE_ORDER)


@dataclass(frozen=True)
class OutgoingCoreParameters:
    """Executable counterpart of ``OutgoingSchedule.Parameters``."""

    P: float
    m: float
    lam: float
    wait: float

    def __post_init__(self) -> None:
        for name in ("P", "m", "lam", "wait"):
            object.__setattr__(self, name, _finite(getattr(self, name), name))
        if self.P <= 0.0:
            raise ValueError("P must be positive")
        if self.m <= 0.0:
            raise ValueError("m must be positive")
        if not 0.0 < self.lam < 0.1:
            raise ValueError("lam must satisfy 0 < lam < 1/10")
        if self.wait <= 25.0:
            raise ValueError("wait must be greater than 25")

    @property
    def drop_length(self) -> float:
        return math.exp(self.m) + 10.0

    @property
    def hold_start(self) -> float:
        return self.drop_length + 2.0

    @property
    def pulse_start(self) -> float:
        return self.hold_start + self.wait

    @property
    def pulse_length(self) -> float:
        return 13.0 / self.lam

    @property
    def endpoint(self) -> float:
        return self.pulse_start + self.pulse_length


@dataclass(frozen=True)
class TailData:
    """Executable counterpart of ``OutgoingTail.TailData`` using ``S=32``."""

    core: OutgoingCoreParameters
    h: float

    def __post_init__(self) -> None:
        if not isinstance(self.core, OutgoingCoreParameters):
            raise TypeError("core must be OutgoingCoreParameters")
        object.__setattr__(self, "h", _finite(self.h, "h"))
        if self.h <= 0.0:
            raise ValueError("h must be positive")
        if 2.0 * self.h >= self.core.lam:
            raise ValueError("TailData requires 2*h < lam")

    @property
    def flatten_end(self) -> float:
        return self.core.endpoint + EXPLICIT_FLATTEN_LENGTH

    @property
    def uniform_wait(self) -> float:
        return 30.0 * math.log(1.0 / self.core.lam)

    @property
    def release_start(self) -> float:
        return self.flatten_end + self.uniform_wait

    @property
    def long_hold(self) -> float:
        return 4.0 * math.log(1.0 / self.h)

    @property
    def second_ramp_start(self) -> float:
        return 1.0 + self.long_hold

    @property
    def ramp_end(self) -> float:
        return self.second_ramp_start + 1.0

    @property
    def rho(self) -> float:
        return EXPLICIT_TAIL_COEFFICIENT * self.h

    @cached_property
    def tail_debt(self) -> float:
        value = integrate(
            lambda t: math.exp((1.0 - self.h) * t) * tail_shape_derivative(self, t),
            1.0,
            3.0,
            n=_TRANSITION_QUADRATURE_ORDER,
        ) / (1.0 - self.rho)
        if not math.isfinite(value) or value <= 0.0:
            raise RuntimeError("numerical tailDebt evaluation did not remain positive")
        return value

    @cached_property
    def release_lag_at_ramp_end(self) -> float:
        q = initial_lag(self)
        s = self.second_ramp_start
        a_mid = release_rate_primitive(self, 1.0)

        left = integrate(
            lambda t: math.exp(release_rate_primitive(self, t)) * release_source(self, t),
            0.0,
            1.0,
            n=_TRANSITION_QUADRATURE_ORDER,
        )
        middle = math.exp(a_mid) * (1.0 - self.h) * (s - 1.0)
        right = integrate(
            lambda t: math.exp(release_rate_primitive(self, t)) * release_source(self, t),
            s,
            self.ramp_end,
            n=_TRANSITION_QUADRATURE_ORDER,
        )
        value = math.exp(-release_rate_primitive(self, self.ramp_end)) * (q + left + middle + right)
        if not math.isfinite(value) or value <= self.tail_debt:
            raise RuntimeError("releaseLag(rampEnd) must exceed the positive tail debt")
        return value

    @cached_property
    def decay_hold(self) -> float:
        value = math.log(self.release_lag_at_ramp_end / self.tail_debt) / (1.0 - self.h)
        if not math.isfinite(value) or value <= 0.0:
            raise RuntimeError("decayHold must be positive")
        return value

    @cached_property
    def tail_start(self) -> float:
        return self.release_start + self.ramp_end + self.decay_hold

    @cached_property
    def tail_end(self) -> float:
        return self.tail_start + 3.0


def log_amplitude(core: OutgoingCoreParameters, y: float) -> float:
    """``OutgoingSchedule.logAmplitude`` using the exact plateau reduction."""

    y = _finite(y, "y")
    return (
        0.1 * y
        - 0.6 * _sigma_primitive(y)
        - core.lam * _sigma_primitive(y - (core.drop_length + 1.0))
    )


def radial_amplitude(core: OutgoingCoreParameters, y: float) -> float:
    return core.P * math.exp(log_amplitude(core, y))


def release_slope(data: TailData, t: float) -> float:
    t = _finite(t, "t")
    return (
        -data.core.lam
        - (1.0 - data.core.lam) * outgoing_sigma(t)
        + (1.0 - data.h) * outgoing_sigma(t - data.second_ramp_start)
    )


def release_source(data: TailData, t: float) -> float:
    return -release_slope(data, t) - data.h


def release_rate_primitive(data: TailData, t: float) -> float:
    """Integral from 0 to t of ``1 + releaseSlope``."""

    t = _finite(t, "t")
    return (
        (1.0 - data.core.lam) * t
        - (1.0 - data.core.lam) * _sigma_primitive(t)
        + (1.0 - data.h) * _sigma_primitive(t - data.second_ramp_start)
    )


def initial_lag(data: TailData) -> float:
    return (data.core.lam - data.h) / (1.0 - data.core.lam)


def release_adjustment(data: TailData, t: float) -> float:
    """``primitive (releaseSlope + lam) t`` in closed plateau form."""

    t = _finite(t, "t")
    return (
        -(1.0 - data.core.lam) * _sigma_primitive(t)
        + (1.0 - data.h) * _sigma_primitive(t - data.second_ramp_start)
    )


def tail_shape(data: TailData, t: float) -> float:
    t = _finite(t, "t")
    return 1.0 - data.rho + data.rho * outgoing_sigma((t - 1.0) / 2.0)


def tail_shape_derivative(data: TailData, t: float) -> float:
    t = _finite(t, "t")
    return (data.rho / 2.0) * outgoing_sigma_derivative((t - 1.0) / 2.0)


def final_angular(data: TailData, y: float, eta: float) -> float:
    """The complete scalar ``OutgoingTail.finalAngular`` formula.

    The evaluation is arranged in logarithmic form to avoid multiplying a tiny
    radial amplitude by a large release factor. This is numerically equivalent
    to ``flattened * exp(releaseAdjustment) * tailShape/(1-rho)``.
    """

    y = _finite(y, "y")
    eta = _finite(eta, "eta")
    s = outgoing_sigma((y - data.core.endpoint) / EXPLICIT_FLATTEN_LENGTH)
    log_shape = math.log1p(eta * eta)
    angular_shape_log = -(1.0 - s) * log_shape - s * math.log(2.0)
    log_value = (
        math.log(data.core.P)
        + log_amplitude(data.core, y)
        + angular_shape_log
        + release_adjustment(data, y - data.release_start)
        + math.log(tail_shape(data, y - data.tail_start))
        - math.log1p(-data.rho)
    )
    return math.exp(log_value)


def clock_weight(data: TailData, y: float) -> float:
    """``SchedulePressure.clockWeight = finalAngular(y,0)^2``."""

    value = final_angular(data, y, 0.0)
    return value * value


def power_constant(data: TailData) -> float:
    """Stable evaluation of ``OutgoingTail.powerConstant``."""

    y = data.tail_start
    log_carrier = (
        math.log(data.core.P)
        + log_amplitude(data.core, y)
        - math.log(2.0)
        + release_adjustment(data, y - data.release_start)
    )
    return math.exp(log_carrier + (0.5 + data.h) * y - math.log1p(-data.rho))


def shape_exponent(data: TailData, y: float) -> float:
    return schedule_shape_exponent(y, data.core.endpoint)
