"""Analytic eta-jets of the actual ``SchedulePressure.axisPressure`` datum.

``SchedulePressure.axisPressure_contDiff`` in the pinned Lean source proves that
all real parameter derivatives exist by holomorphic extension.  The existing
:mod:`schedule_axis_pressure` module materializes only the value and first
derivative.  This module keeps the same actual outgoing schedule and computes
arbitrary *finite* normalized Taylor jets without finite differences or fitted
coefficients.

For a fixed shape exponent ``a`` the pressure kernel is

    k_a(eta) = (1 + eta**2)**(-2*a).

We propagate the Taylor series of ``log(1+eta**2)`` and then ``exp`` exactly at
the level of truncated power-series algebra.  The y-integral is split using the
same theorem geometry as ``schedule_axis_pressure``: exponent 1 before the
flattening interval, variable exponent only on that compact interval, and
exponent 0 afterwards.  Constant-exponent clock masses reuse the landed
analytic plateau/tail reductions; only the compact flattening interval uses the
repository Gauss--Legendre quadrature.

The returned numbers are reproducible floating-point evaluations of the
paper-defined derivative integrals, not rigorous interval enclosures.  They do
not promote Stage 1 to ``paper-exact``.
"""

from __future__ import annotations

from functools import lru_cache
import math

from .outgoing_tail import TailData, clock_weight, shape_exponent
from .quadrature import integrate
from .schedule_axis_pressure import (
    _PRESSURE_QUADRATURE_ORDER,
    _exponential_interval_mass,
    _geometry,
    _transition_integral,
    ideal_prefix_mass,
)


def _index(value: int, name: str = "order") -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _reciprocal(a: tuple[float, ...]) -> tuple[float, ...]:
    if not a or not math.isfinite(a[0]) or a[0] == 0.0:
        raise ArithmeticError("Taylor reciprocal requires a finite nonzero constant term")
    out = [1.0 / a[0]]
    for n in range(1, len(a)):
        out.append(-sum(a[k] * out[n - k] for k in range(1, n + 1)) / a[0])
    if not all(math.isfinite(x) for x in out):
        raise ArithmeticError("Taylor reciprocal produced a non-finite coefficient")
    return tuple(out)


def pressure_kernel_normalized_taylor(
    exponent: float,
    eta: float,
    order: int,
) -> tuple[float, ...]:
    """Return ``k_a^(m)(eta)/m!`` for ``m=0..order``.

    This is truncated power-series algebra for the exact pinned
    ``PressureDatum.kernel``; no numerical differentiation is used.
    """

    order = _index(order)
    exponent = _finite(exponent, "exponent")
    eta = _finite(eta, "eta")
    if not 0.0 <= exponent <= 1.0:
        raise ValueError("schedule exponent must lie in [0, 1]")
    if exponent == 0.0:
        return (1.0,) + (0.0,) * order

    # Normalized Taylor coefficients of b(t)=1+(eta+t)^2.
    base = [0.0] * (order + 1)
    base[0] = 1.0 + eta * eta
    if order >= 1:
        base[1] = 2.0 * eta
    if order >= 2:
        base[2] = 1.0
    inverse = _reciprocal(tuple(base))

    # log(b)' = b'/b.  If c_n stores f^(n)/n!, then the normalized
    # derivative series has coefficient (n+1)c_{n+1}.
    log_base = [0.0] * (order + 1)
    log_base[0] = math.log(base[0])
    for n in range(order):
        quotient_n = 0.0
        for k in range(n + 1):
            derivative_k = (k + 1) * base[k + 1] if k + 1 < len(base) else 0.0
            quotient_n += derivative_k * inverse[n - k]
        log_base[n + 1] = quotient_n / float(n + 1)

    exponent_series = tuple(-2.0 * exponent * x for x in log_base)

    # If f=exp(g), normalized coefficients satisfy
    # n f_n = sum_{k=1}^n k g_k f_{n-k}.
    out = [0.0] * (order + 1)
    out[0] = math.exp(exponent_series[0])
    for n in range(1, order + 1):
        out[n] = sum(
            k * exponent_series[k] * out[n - k] for k in range(1, n + 1)
        ) / float(n)
    if not all(math.isfinite(x) for x in out):
        raise ArithmeticError("pressure-kernel Taylor jet must remain finite")
    return tuple(out)


@lru_cache(maxsize=32)
def _constant_exponent_clock_masses(data: TailData) -> tuple[float, float]:
    """Clock mass before flattening and after flattening.

    The middle interval ``[endpoint, flatten_end]`` is intentionally excluded
    because its shape exponent varies from one to zero.
    """

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
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

    pre = ideal_prefix_mass(data, 0.0)
    pre += _transition_integral(data, 0.0, zero, first_end)
    pre += _exponential_interval_mass(data, 0.0, first_end, drop_transition_start, 0.5)
    pre += _transition_integral(data, 0.0, drop_transition_start, hold_start)
    pre += _exponential_interval_mass(data, 0.0, hold_start, endpoint, 0.5 + data.core.lam)

    post = _exponential_interval_mass(data, 0.0, flatten_end, release_start, 0.5 + data.core.lam)
    post += _transition_integral(data, 0.0, release_start, release_first_end)
    post += _exponential_interval_mass(data, 0.0, release_first_end, release_second_start, 1.5)
    post += _transition_integral(data, 0.0, release_second_start, release_ramp_end)
    post += _exponential_interval_mass(data, 0.0, release_ramp_end, tail_transition_start, 0.5 + data.h)
    post += _transition_integral(data, 0.0, tail_transition_start, tail_end)
    post += clock_weight(data, tail_end) / (1.0 + 2.0 * data.h)

    if not math.isfinite(pre) or not math.isfinite(post) or pre <= 0.0 or post <= 0.0:
        raise ArithmeticError("constant-exponent clock masses must remain finite and positive")
    return pre, post


def axis_pressure_normalized_taylor(
    data: TailData,
    eta: float,
    order: int,
) -> tuple[float, ...]:
    """Return ``axisPressure^(m)(eta)/m!`` for ``m=0..order``.

    The derivative chain is the actual schedule integral from the pinned
    ``SchedulePressure``/``PressureDatum`` construction.  It is not a sampled
    polynomial fit and accepts no caller-supplied pressure surrogate.
    """

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    eta = _finite(eta, "eta")
    order = _index(order)
    pre_mass, post_mass = _constant_exponent_clock_masses(data)
    pre_kernel = pressure_kernel_normalized_taylor(1.0, eta, order)

    endpoint = data.core.endpoint
    flatten_end = data.flatten_end
    flatten = []
    for m in range(order + 1):
        value = integrate(
            lambda y, m=m: clock_weight(data, y)
            * pressure_kernel_normalized_taylor(shape_exponent(data, y), eta, order)[m],
            endpoint,
            flatten_end,
            n=_PRESSURE_QUADRATURE_ORDER,
        )
        flatten.append(value)

    out = []
    for m in range(order + 1):
        mass_coefficient = pre_mass * pre_kernel[m] + flatten[m]
        if m == 0:
            mass_coefficient += post_mass
        out.append(-0.5 * mass_coefficient)
    if not all(math.isfinite(x) for x in out):
        raise ArithmeticError("schedule pressure Taylor jet must remain finite")
    return tuple(out)
