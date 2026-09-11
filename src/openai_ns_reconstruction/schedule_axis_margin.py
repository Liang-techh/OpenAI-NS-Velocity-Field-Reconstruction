"""Constructive low-|Z| H^2 margin for the actual schedule axis pressure.

This module turns the compactness-only step in pinned
``NaturalAxisRange.low_Z_has_H_margin`` into an explicit conservative bound for
an executable ``OutgoingTail.TailData``.  It is tied to the actual
``SchedulePressure.axisPressure`` datum, not the earlier ideal-prefix surrogate.

The argument deliberately does not sample the compact set.  In real arithmetic:

1. The schedule ideal prefix and ``P >= 2`` give the ``PressureData`` hypotheses.
2. A closed schedule envelope gives an upper bound ``M`` for total clock mass.
3. The pressure kernel gives ``|P| <= M/2``, ``|P'| <= M`` and
   ``|P''| <= 4 M`` on ``[-1,1]``.
4. Those bounds give a global Lipschitz constant ``L_Z`` for
   ``NaturalAxisData.Z``.
5. At the unique H-root, pinned ``Z_at_root_lower`` gives ``Z > j/5``.  Hence
   ``|Z| <= j/10`` forces distance at least ``j/(10 L_Z)`` from that root.
6. The exact ``cross_identity`` then gives a quantitative lower bound for
   ``|H|`` and therefore a positive uniform ``H^2`` margin.

The returned floating-point numbers are an executable realization of that
analytic inequality chain.  They are intentionally conservative and are not an
interval-arithmetic or Lean proof object; Stage 1 remains ``formal-structure``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .natural_axis import A, D, H, axis_U, d
from .natural_axis_range import (
    CutoffParametersFromMargin,
    NaturalAxisRangeParameters,
    cutoff_parameters_from_margin,
)
from .outgoing_tail import TailData
from .schedule_axis_pressure import axis_pressure, axis_pressure_derivative


# e^-5 < 1/32, and 16 * (S + 1) = 16 * 33 = 528 for the landed S=32
# realization.  Therefore rho < h / (32*528) = h/16896.
_RHO_DENOMINATOR = 16896.0


@dataclass(frozen=True)
class ScheduleLowZMarginWitness:
    """Explicit conservative witness for ``low_Z_has_H_margin``.

    ``margin`` is the value that may be passed to
    :func:`cutoff_parameters_from_margin`; ``cutoff`` records the resulting
    theorem-side choice ``delta=j/10`` and ``sigma=sqrt(margin)/20``.
    """

    parameters: NaturalAxisRangeParameters
    release_start_upper: float
    rho_upper: float
    clock_mass_upper: float
    z_lipschitz_upper: float
    root_distance_lower: float
    h_root_slope_lower: float
    margin: float
    cutoff: CutoffParametersFromMargin


def _finite_positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ArithmeticError(f"derived {name} must remain finite and positive")
    return value


def _exp_m_upper(m: float) -> float:
    """Algebraic upper bound for ``exp(m)`` without using it as the proof step.

    For integer ``n > m > 0``, ``-log(1-m/n) > m/n`` implies
    ``exp(m) < (1-m/n)^(-n)``.  We use ``n = 2*ceil(m)+2`` so ``m/n < 1/2``.
    """

    m = float(m)
    if not math.isfinite(m) or m <= 0.0:
        raise ValueError("m must be finite and positive")
    n = 2 * math.ceil(m) + 2
    base = 1.0 - m / n
    if not 0.0 < base < 1.0:
        raise ArithmeticError("failed to form exp(m) upper-bound base")
    try:
        value = base ** (-n)
    except OverflowError as exc:
        raise ValueError("m is too large for a finite runtime margin witness") from exc
    return _finite_positive(value, "exp(m) upper bound")


def release_start_upper(data: TailData) -> float:
    """Conservative algebraic upper bound for the printed ``releaseStart``.

    The exact scalar schedule has

    ``releaseStart = exp(m)+13+wait+13/lam+330 log(2)+30 log(1/lam)``.

    Using ``log(2)<1``, ``log(1/lam)<1/lam`` and ``_exp_m_upper`` gives the
    bound below.  It avoids depending on tailDebt/decayHold or any quadrature.
    """

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    c = data.core
    value = _exp_m_upper(c.m) + 343.0 + c.wait + 43.0 / c.lam
    return _finite_positive(value, "release-start upper bound")


def clock_mass_upper(data: TailData) -> tuple[float, float, float]:
    """Return ``(M, release_upper, rho_upper)`` with ``integral clockWeight <= M``.

    The ideal prefix contributes exactly ``5 P^2``.  On
    ``0 <= y <= releaseStart`` the release/tail factors are inactive and
    ``logAmplitude <= 1/10``; hence ``clockWeight <= exp(1/5) P^2 < 5 P^2/4``.
    After ``releaseStart`` the combined carrier logarithmic slope is at most
    ``-(1/2+h)``.  The terminal taper contributes at most ``1/(1-rho)`` and
    ``rho < h/16896``.  Integrating this exponential tail gives the last term.
    """

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    h = float(data.h)
    release = release_start_upper(data)
    rho_upper = h / _RHO_DENOMINATOR
    if not 0.0 <= rho_upper < 1.0:
        raise ArithmeticError("rho upper bound must lie in [0,1)")
    p2 = data.core.P * data.core.P
    middle = 1.25 * p2 * release
    tail = 1.25 * p2 / ((1.0 - rho_upper) ** 2 * (1.0 + 2.0 * h))
    total = 5.0 * p2 + middle + tail
    return _finite_positive(total, "clock-mass upper bound"), release, rho_upper


def schedule_natural_axis_Z(data: TailData, j: float, eta: float) -> float:
    """Evaluate the actual schedule pressure inside ``NaturalAxisData.Z``."""

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    parameters = NaturalAxisRangeParameters(h=data.h, j=j)
    eta = float(eta)
    if not math.isfinite(eta) or not -1.0 <= eta <= 1.0:
        raise ValueError("eta must be finite and lie in [-1,1]")
    u = axis_U(parameters.j, eta)
    p = axis_pressure(data, eta)
    dp = axis_pressure_derivative(data, eta)
    value = (
        -A(parameters.h) * (1.0 - 2.0 * eta * u) * u
        - 4.0 * H(parameters.h, parameters.j, eta)
        - d(eta) * dp
        + 4.0 * A(parameters.h) * eta * p
    )
    if not math.isfinite(value):
        raise ArithmeticError("schedule NaturalAxisData.Z must remain finite")
    return value


def certify_schedule_low_Z_margin(data: TailData, j: float) -> ScheduleLowZMarginWitness:
    """Construct an explicit uniform low-|Z| margin for the actual schedule datum.

    Preconditions match the pinned theorem chain: ``0<h<=1/100``,
    ``0<j<=1/20`` and the ideal-prefix amplitude ``P>=2``.  ``TailData`` already
    enforces the outgoing schedule constraints including ``2h < lam < 1/10``.

    No sampled minimum is used.  The proof constants are intentionally coarse:
    they trade sharpness for a transparent executable witness.
    """

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    parameters = NaturalAxisRangeParameters(h=data.h, j=j)
    if data.core.P < 2.0:
        raise ValueError("the NaturalAxisRange pressure theorem requires schedule amplitude P >= 2")

    mass_upper, release_upper, rho_upper = clock_mass_upper(data)
    h = parameters.h
    j = parameters.j
    a = A(h)
    dd = D(h)

    # Global derivative envelope for the polynomial part of Z on [-1,1].
    # U=4 eta+j, F=1-2 eta U, and K=F U.
    u_max = 4.0 + j
    f_max = 1.0 + 2.0 * u_max
    f_prime_max = 2.0 * (8.0 + j)
    k_prime_max = f_prime_max * u_max + 4.0 * f_max
    h_prime_max = dd + 2.0 * u_max + 4.0
    polynomial_z_prime = a * k_prime_max + 4.0 * h_prime_max

    # Kernel inequalities from the actual pressure integral give
    # |P|<=M/2, |P'|<=M, |P''|<=4M.  Differentiating
    # Z=-A*K-4H-d*P'+4A*eta*P then costs (6+6A)M.
    z_lipschitz = polynomial_z_prime + (6.0 + 6.0 * a) * mass_upper
    z_lipschitz = _finite_positive(z_lipschitz, "Z Lipschitz upper bound")

    root_distance = j / (10.0 * z_lipschitz)
    root_distance = _finite_positive(root_distance, "root-distance lower bound")

    # From cross_identity at the unique root r:
    # H(eta)d(r)=(eta-r)[D(1+eta*r)+4d(eta)d(r)].
    # Since |r|<j/4, d(r)<=1 and d(eta),d(r)>=0,
    # |H(eta)| > D(1-j/4)|eta-r|.
    h_root_slope = dd * (1.0 - j / 4.0)
    h_root_slope = _finite_positive(h_root_slope, "H root-separation slope")

    # The theorem-side inequalities are strict, so half the squared analytic
    # lower envelope is still a valid non-strict positive margin and leaves a
    # generous floating realization slack.
    raw = h_root_slope * root_distance
    margin = _finite_positive(0.5 * raw * raw, "uniform H^2 margin")
    cutoff = cutoff_parameters_from_margin(parameters, margin)

    return ScheduleLowZMarginWitness(
        parameters=parameters,
        release_start_upper=release_upper,
        rho_upper=rho_upper,
        clock_mass_upper=mass_upper,
        z_lipschitz_upper=z_lipschitz,
        root_distance_lower=root_distance,
        h_root_slope_lower=h_root_slope,
        margin=margin,
        cutoff=cutoff,
    )
