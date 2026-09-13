"""Analytic parameter jets for the actual-schedule angular reference coefficient.

The pinned ``AxisCoefficientSpace.CoefficientSpace`` stores *actual* adjacent
parameter derivatives, not unrelated arrays.  The previously landed
:mod:`axis_reference_pair` materializes only the zeroth parameter jet of the
reference pair.  This module removes that limitation for the angular component.

For the official reference,

    phi0[n](eta) = (-chi(eta)/2)^n / (n! (n+1)!),

where ``chi = H^2/(H^2+sigma^2)`` and ``H`` is a cubic polynomial in ``eta``.
Consequently every parameter derivative is obtained exactly from finite Taylor
algebra around the requested point: no finite differences, sampled derivative
tables, fitted coefficients, or default-zero higher jets are used.

This is deliberately *not* yet a complete ``AxisCoefficientSpace`` state.  A
full membership witness additionally needs a uniform weighted bound for the
infinite jet family, and the axial reference still needs all of its parameter
jets.  Those remain explicit blockers before ``naturalRemainder`` may consume a
complete reference state.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .axis_reference_pair import ActualScheduleReferencePair, actual_schedule_reference_pair
from .natural_axis import D
from .outgoing_tail import TailData
from .schedule_analytic_neighborhood import certify_actual_schedule_analytic_inputs

_WINDOW_LEFT = -11.0 / 10.0
_WINDOW_RIGHT = 11.0 / 10.0


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _eta(value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or not _WINDOW_LEFT <= value <= _WINDOW_RIGHT:
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]")
    return value


def _convolve(left: list[float], right: list[float], order: int) -> list[float]:
    out = [0.0] * (order + 1)
    for k in range(order + 1):
        total = 0.0
        for i in range(k + 1):
            if i < len(left) and k - i < len(right):
                total += left[i] * right[k - i]
        out[k] = total
    return out


def _divide_series(numerator: list[float], denominator: list[float], order: int) -> list[float]:
    d0 = denominator[0]
    if not math.isfinite(d0) or d0 <= 0.0:
        raise ArithmeticError("H^2+sigma^2 must stay positive")
    quotient = [0.0] * (order + 1)
    quotient[0] = numerator[0] / d0
    for k in range(1, order + 1):
        correction = sum(denominator[i] * quotient[k - i] for i in range(1, k + 1))
        quotient[k] = (numerator[k] - correction) / d0
    if not all(math.isfinite(value) for value in quotient):
        raise ArithmeticError("analytic quotient jet must remain finite")
    return quotient


def _power_series(base: list[float], exponent: int, order: int) -> list[float]:
    exponent = _index(exponent, "exponent")
    result = [0.0] * (order + 1)
    result[0] = 1.0
    factor = list(base)
    power = exponent
    # Truncated binary exponentiation keeps large radial indices practical.
    while power:
        if power & 1:
            result = _convolve(result, factor, order)
        power >>= 1
        if power:
            factor = _convolve(factor, factor, order)
    return result


def _h_taylor(h: float, j: float, eta: float, order: int) -> list[float]:
    """Return Taylor coefficients ``H^(m)(eta)/m!`` through ``order``.

    The exact polynomial identity is

        H = j + (D+4) eta - j eta^2 - 4 eta^3.
    """

    coefficients = (float(j), D(h) + 4.0, -float(j), -4.0)
    out: list[float] = []
    for k in range(order + 1):
        value = 0.0
        for degree in range(k, len(coefficients)):
            value += (
                coefficients[degree]
                * math.comb(degree, k)
                * eta ** (degree - k)
            )
        out.append(value)
    return out


def _chi_taylor(reference: ActualScheduleReferencePair, eta: float, order: int) -> list[float]:
    h_series = _h_taylor(reference.data.h, reference.j, eta, order)
    h_squared = _convolve(h_series, h_series, order)
    denominator = list(h_squared)
    denominator[0] += reference.sigma * reference.sigma
    return _divide_series(h_squared, denominator, order)


def coefficient_weight(epsilon: float, n: int, m: int) -> float:
    """Pinned ``AxisWeightEstimates.weight epsilon n m`` in binary64.

    This helper exposes the normalization used by ``AxisCoefficientSpace``; it
    is not by itself a uniform norm certificate.
    """

    epsilon = float(epsilon)
    n = _index(n, "n")
    m = _index(m, "m")
    if not math.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("epsilon must be finite and positive")
    value = (
        (1.0 / 20.0) ** n
        * epsilon ** (-m)
        * math.factorial(m)
        * math.comb(n + m, m)
        / ((n + 1.0) ** 2 * (m + 1.0) ** 2)
    )
    if not math.isfinite(value) or value <= 0.0:
        raise ArithmeticError("coefficient weight is not representable as a positive binary64 value")
    return value


@dataclass(frozen=True)
class ActualScheduleAngularReferenceJets:
    """All actual parameter jets of the angular reference radial coefficients.

    ``epsilon`` is not caller supplied: it is the canonical
    ``rho/2`` produced by the landed actual-schedule analytic-neighborhood
    certificate, matching the radius loss used by the Stage-1 norm chain.
    """

    reference: ActualScheduleReferencePair
    epsilon: float

    def __post_init__(self) -> None:
        if not isinstance(self.reference, ActualScheduleReferencePair):
            raise TypeError("reference must be ActualScheduleReferencePair")
        epsilon = float(self.epsilon)
        if not math.isfinite(epsilon) or epsilon <= 0.0:
            raise ValueError("epsilon must be finite and positive")
        object.__setattr__(self, "epsilon", epsilon)

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def coefficient_space_membership_certified(self) -> bool:
        """Remain false until the infinite weighted supremum bound is closed."""

        return False

    def parameter_jet(self, n: int, m: int, eta: float) -> float:
        """Return ``partial_eta^m phi0[n](eta)`` analytically.

        Taylor coefficients are propagated through the exact rational formula
        for ``chi`` and the exact factorial reference formula.  Multiplication
        by ``m!`` converts the local Taylor coefficient to the actual derivative.
        """

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta(eta)
        chi_series = _chi_taylor(self.reference, eta, m)
        base = [-0.5 * value for value in chi_series]
        powered = _power_series(base, n, m)
        denominator = math.factorial(n) * math.factorial(n + 1)
        value = math.factorial(m) * powered[m] / denominator
        if not math.isfinite(value):
            raise ArithmeticError("angular reference parameter jet must remain finite")
        return value

    def normalized_jet(self, n: int, m: int, eta: float) -> float:
        """Return the stored-coordinate normalization ``jet/weight``."""

        return self.parameter_jet(n, m, eta) / coefficient_weight(self.epsilon, n, m)


def actual_schedule_angular_reference_jets(
    data: TailData,
    j: float,
) -> ActualScheduleAngularReferenceJets:
    """Build the actual-schedule angular derivative family with no free radius.

    The same analytic-neighborhood certificate used by the landed Stage-1 norm
    chain supplies ``rho`` and hence ``epsilon=rho/2``.  The reference itself
    supplies the theorem-side ``sigma`` from the certified actual schedule.
    """

    reference = actual_schedule_reference_pair(data, j)
    analytic = certify_actual_schedule_analytic_inputs(data, j)
    epsilon = float(analytic.neighborhood.radius) / 2.0
    return ActualScheduleAngularReferenceJets(reference=reference, epsilon=epsilon)
