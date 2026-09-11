"""Executable evaluator for the pinned ``SchedulePressure.axisPressure`` datum.

The official Lean definition is

    axisPressure(d, eta) = -1/2 * integral_R finalAngular(d, (y, eta))^2 dy.

``SchedulePressure.angular_square_factorization`` rewrites the same integrand as

    clockWeight(d, y) * PressureDatum.kernel(shapeExponent(d, y), eta).

This module evaluates that *actual constructed schedule* for any executable
:class:`~openai_ns_reconstruction.outgoing_tail.TailData`.  It does not replace
the schedule with the earlier ideal-prefix witness.

The improper integral is reduced using identities proved in the pinned Lean
sources.  The ideal prefix ``y <= 0`` is integrated exactly, the eventual power
tail is integrated exactly, and every long constant-slope plateau is integrated
as an exponential from its left endpoint.  Only the six smooth transition
intervals are sent through the repository's cached Gauss--Legendre quadrature.
The result is therefore a reproducible numerical value of the paper-defined
integral, not a formal error bound or the still-missing uniform low-|Z| margin.
"""

from __future__ import annotations

import math

from .outgoing_tail import TailData, clock_weight, shape_exponent
from .quadrature import integrate


_PRESSURE_QUADRATURE_ORDER = 192


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _require_data(data: TailData) -> TailData:
    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    return data


def pressure_kernel(exponent: float, eta: float) -> float:
    """Stable real ``PressureDatum.kernel`` = ``(1+eta^2)^(-2 exponent)``.

    The logarithm is evaluated without forming ``eta**2`` when ``|eta| > 1``
    so every finite real ``eta`` remains a valid numerical input.
    """

    exponent = _finite(exponent, "exponent")
    eta = _finite(eta, "eta")
    if exponent < 0.0 or exponent > 1.0:
        raise ValueError("schedule exponent must lie in [0, 1]")
    if exponent == 0.0:
        return 1.0

    x = abs(eta)
    if x <= 1.0:
        log_base = math.log1p(x * x)
    else:
        inv = 1.0 / x
        log_base = 2.0 * math.log(x) + math.log1p(inv * inv)
    return math.exp(-2.0 * exponent * log_base)


def _density(data: TailData, y: float, eta: float) -> float:
    """Factorized ``finalAngular^2`` from ``SchedulePressure``."""

    exponent = shape_exponent(data, y)
    value = clock_weight(data, y) * pressure_kernel(exponent, eta)
    if not math.isfinite(value) or value < 0.0:
        raise ArithmeticError("schedule pressure density must remain finite and nonnegative")
    return value


def _weighted_density(data: TailData, y: float, eta: float) -> float:
    """``shapeExponent * finalAngular^2`` used by the pressure derivative."""

    exponent = shape_exponent(data, y)
    value = exponent * clock_weight(data, y) * pressure_kernel(exponent, eta)
    if not math.isfinite(value) or value < 0.0:
        raise ArithmeticError("weighted schedule density must remain finite and nonnegative")
    return value


def _transition_integral(data: TailData, eta: float, a: float, b: float, *, weighted: bool = False) -> float:
    if b < a:
        raise ArithmeticError("outgoing transition endpoints are not ordered")
    fn = _weighted_density if weighted else _density
    return integrate(lambda y: fn(data, y, eta), a, b, n=_PRESSURE_QUADRATURE_ORDER)


def _exponential_interval_mass(data: TailData, eta: float, a: float, b: float, rate: float) -> float:
    """Integrate a constant-shape interval with ``finalAngular ~ exp(-rate*y)``."""

    if b < a or rate <= 0.0 or not math.isfinite(rate):
        raise ArithmeticError("invalid constant-slope schedule interval")
    if a == b:
        return 0.0
    left = _density(data, a, eta)
    length = b - a
    factor = -math.expm1(-2.0 * rate * length) / (2.0 * rate)
    return left * factor


def _geometry(data: TailData) -> tuple[float, ...]:
    """Return the exact transition/plateau endpoints used by ``finalAngular``."""

    try:
        points = (
            0.0,
            1.0,
            data.core.drop_length + 1.0,
            data.core.hold_start,
            data.core.endpoint,
            data.flatten_end,
            data.release_start,
            data.release_start + 1.0,
            data.release_start + data.second_ramp_start,
            data.release_start + data.ramp_end,
            data.tail_start + 1.0,
            data.tail_end,
        )
    except OverflowError as exc:
        raise ValueError("derived outgoing schedule geometry must remain finite") from exc
    if not all(math.isfinite(v) for v in points):
        raise ValueError("derived outgoing schedule geometry must remain finite")
    if any(b < a for a, b in zip(points, points[1:])):
        raise ArithmeticError("derived outgoing schedule geometry is not ordered")
    return points


def ideal_prefix_mass(data: TailData, eta: float) -> float:
    """Exact mass on ``(-infinity, 0]`` from ``clockWeight_ideal``."""

    data = _require_data(data)
    eta = _finite(eta, "eta")
    return 5.0 * data.core.P * data.core.P * pressure_kernel(1.0, eta)


def schedule_pressure_mass(data: TailData, eta: float) -> float:
    """Evaluate ``integral_R finalAngular(d,(y,eta))^2 dy``.

    Long plateaus are not sampled numerically.  Their exact constant logarithmic
    slopes from ``OutgoingSchedule``/``OutgoingTail`` are used to integrate the
    squared amplitude exponentially.  This both respects the printed schedule
    and avoids accuracy loss when ``m`` or the holding times are large.
    """

    data = _require_data(data)
    eta = _finite(eta, "eta")
    (
        zero,
        first_end,
        drop_transition_start,
        hold_start,
        endpoint,
        flatten_end,
        release_start,
        release_first_end,
        release_second_start,
        release_ramp_end,
        tail_transition_start,
        tail_end,
    ) = _geometry(data)

    mass = ideal_prefix_mass(data, eta)
    mass += _transition_integral(data, eta, zero, first_end)
    mass += _exponential_interval_mass(data, eta, first_end, drop_transition_start, 0.5)
    mass += _transition_integral(data, eta, drop_transition_start, hold_start)
    mass += _exponential_interval_mass(data, eta, hold_start, endpoint, 0.5 + data.core.lam)
    mass += _transition_integral(data, eta, endpoint, flatten_end)
    mass += _exponential_interval_mass(data, eta, flatten_end, release_start, 0.5 + data.core.lam)
    mass += _transition_integral(data, eta, release_start, release_first_end)
    mass += _exponential_interval_mass(data, eta, release_first_end, release_second_start, 1.5)
    mass += _transition_integral(data, eta, release_second_start, release_ramp_end)
    mass += _exponential_interval_mass(data, eta, release_ramp_end, tail_transition_start, 0.5 + data.h)
    mass += _transition_integral(data, eta, tail_transition_start, tail_end)

    # ``finalAngular_eventual_power`` gives logarithmic amplitude slope
    # ``-(1/2+h)`` for every y >= tailEnd, hence squared-density rate ``1+2h``.
    mass += _density(data, tail_end, eta) / (1.0 + 2.0 * data.h)
    if not math.isfinite(mass) or mass <= 0.0:
        raise ArithmeticError("schedule pressure mass must remain finite and positive")
    return mass


def axis_pressure(data: TailData, eta: float) -> float:
    """Numerical value of the actual pinned ``SchedulePressure.axisPressure``."""

    return -0.5 * schedule_pressure_mass(data, eta)


def exponent_weighted_mass(data: TailData, eta: float) -> float:
    """Evaluate ``integral shapeExponent * finalAngular^2``.

    The exponent equals one before ``endpoint`` and zero after ``flattenEnd``;
    only the flattening transition needs the exponent-weighted quadrature.
    This is the integral appearing in ``SchedulePressure.axisPressure_hasDerivAt``.
    """

    data = _require_data(data)
    eta = _finite(eta, "eta")
    (
        zero,
        first_end,
        drop_transition_start,
        hold_start,
        endpoint,
        flatten_end,
        *_rest,
    ) = _geometry(data)

    mass = ideal_prefix_mass(data, eta)
    mass += _transition_integral(data, eta, zero, first_end)
    mass += _exponential_interval_mass(data, eta, first_end, drop_transition_start, 0.5)
    mass += _transition_integral(data, eta, drop_transition_start, hold_start)
    mass += _exponential_interval_mass(data, eta, hold_start, endpoint, 0.5 + data.core.lam)
    mass += _transition_integral(data, eta, endpoint, flatten_end, weighted=True)
    if not math.isfinite(mass) or mass <= 0.0:
        raise ArithmeticError("positive exponent mass must remain finite and positive")
    return mass


def axis_pressure_derivative(data: TailData, eta: float) -> float:
    """Derivative formula proved by ``SchedulePressure.axisPressure_hasDerivAt``."""

    data = _require_data(data)
    eta = _finite(eta, "eta")
    if eta == 0.0:
        return 0.0
    x = abs(eta)
    if x <= 1.0:
        prefactor = 2.0 * eta / (1.0 + eta * eta)
    else:
        inv = 1.0 / eta
        prefactor = (2.0 * inv) / (1.0 + inv * inv)
    return prefactor * exponent_weighted_mass(data, eta)
