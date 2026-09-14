"""Analytic parameter jets for the actual-schedule axial reference coefficient.

The pinned reference pair has exactly one nonzero axial radial coefficient,

    u0[1](eta) = -1/2 * inverseL(eta) * zStar(eta),

with ``inverseL = 1 / (1 - 2 h eta^2)`` and ``zStar`` built from the actual
``SchedulePressure.axisPressure`` datum.  This module materializes arbitrary
finite eta-jets of that coefficient by truncated Taylor-series algebra around
the requested eta.  Pressure derivatives come only from
:mod:`schedule_axis_pressure_jets`; production code never uses finite
differences, sampled fits, caller-supplied pressure rows, or hand-filled higher
parameter derivatives.

Radial degrees other than one are exact structural zeros because the official
regular inverse ``J_1`` shifts the radially constant source into degree one.
That statement is part of the reference-pair formula, not a default-zero
coefficient convention.

The pressure jet evaluator still carries its documented finite Gauss--Legendre
quadrature caveat.  Therefore this file provides an executable actual-schedule
jet family, not an interval enclosure and not yet an
``AxisCoefficientSpace.ofJetFamily`` membership witness.  The latter still
requires one uniform all-order weighted bound for the angular and axial jet
families.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .axis_reference_pair import ActualScheduleReferencePair, actual_schedule_reference_pair
from .axis_reference_parameter_jets import coefficient_weight
from .outgoing_tail import TailData
from .schedule_analytic_neighborhood import certify_actual_schedule_analytic_inputs
from .schedule_axis_pressure_jets import axis_pressure_normalized_taylor

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


def _finite_series(values: list[float], name: str) -> list[float]:
    if not all(math.isfinite(value) for value in values):
        raise ArithmeticError(f"{name} Taylor series must remain finite")
    return values


def _constant(value: float, order: int) -> list[float]:
    out = [0.0] * (order + 1)
    out[0] = float(value)
    return out


def _coordinate(eta: float, order: int) -> list[float]:
    out = [0.0] * (order + 1)
    out[0] = eta
    if order >= 1:
        out[1] = 1.0
    return out


def _add(left: list[float], right: list[float]) -> list[float]:
    if len(left) != len(right):
        raise ValueError("Taylor series must have the same truncation order")
    return _finite_series([a + b for a, b in zip(left, right)], "sum")


def _scale(series: list[float], scalar: float) -> list[float]:
    scalar = float(scalar)
    if not math.isfinite(scalar):
        raise ValueError("Taylor scalar must be finite")
    return _finite_series([scalar * value for value in series], "scaled")


def _multiply(left: list[float], right: list[float]) -> list[float]:
    if len(left) != len(right):
        raise ValueError("Taylor series must have the same truncation order")
    order = len(left) - 1
    out = [0.0] * (order + 1)
    for n in range(order + 1):
        out[n] = math.fsum(left[k] * right[n - k] for k in range(n + 1))
    return _finite_series(out, "product")


def _divide(numerator: list[float], denominator: list[float]) -> list[float]:
    if len(numerator) != len(denominator):
        raise ValueError("Taylor series must have the same truncation order")
    if not denominator or not math.isfinite(denominator[0]) or denominator[0] == 0.0:
        raise ArithmeticError("Taylor quotient requires a finite nonzero constant denominator")
    order = len(numerator) - 1
    out = [0.0] * (order + 1)
    out[0] = numerator[0] / denominator[0]
    for n in range(1, order + 1):
        correction = math.fsum(denominator[k] * out[n - k] for k in range(1, n + 1))
        out[n] = (numerator[n] - correction) / denominator[0]
    return _finite_series(out, "quotient")


def _axial_normalized_taylor(
    reference: ActualScheduleReferencePair,
    eta: float,
    order: int,
) -> tuple[float, ...]:
    """Return ``u0[1]^(m)(eta)/m!`` through ``order``.

    The implementation is a direct normalized-Taylor replay of the pinned
    ``NaturalAxisData.Z`` and ``referencePair`` formulas.  ``P'`` is obtained by
    shifting the actual pressure Taylor jet, so order ``m`` legitimately asks
    the pressure backend for order ``m+1``.
    """

    eta = _eta(eta)
    order = _index(order, "order")
    data = reference.data

    pressure = list(axis_pressure_normalized_taylor(data, eta, order + 1))
    pressure_prime = [(m + 1.0) * pressure[m + 1] for m in range(order + 1)]
    pressure = pressure[: order + 1]

    x = _coordinate(eta, order)
    one = _constant(1.0, order)
    u_axis = _add(_constant(reference.j, order), _scale(x, 4.0))
    x_squared = _multiply(x, x)
    d_series = _add(one, _scale(x_squared, -1.0))

    D = 0.5 - data.h
    A = 0.5 + data.h
    H_series = _add(_scale(x, D), _multiply(d_series, u_axis))

    one_minus_2xu = _add(one, _scale(_multiply(x, u_axis), -2.0))
    z_series = _scale(_multiply(one_minus_2xu, u_axis), -A)
    z_series = _add(z_series, _scale(H_series, -4.0))
    z_series = _add(z_series, _scale(_multiply(d_series, pressure_prime), -1.0))
    z_series = _add(z_series, _scale(_multiply(x, pressure), 4.0 * A))

    L_series = _add(one, _scale(x_squared, -2.0 * data.h))
    u_series = _scale(_divide(z_series, L_series), -0.5)
    return tuple(u_series)


@dataclass(frozen=True)
class ActualScheduleAxialReferenceJets:
    """Actual eta-derivative family of the axial reference coefficients."""

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
        """Stay false until a uniform infinite weighted bound is proved."""

        return False

    def normalized_parameter_taylor_coefficient(self, n: int, m: int, eta: float) -> float:
        """Return ``u0[n]^(m)(eta)/m!`` from actual schedule data."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta(eta)
        if n != 1:
            return 0.0
        value = _axial_normalized_taylor(self.reference, eta, m)[m]
        if not math.isfinite(value):
            raise ArithmeticError("axial reference normalized parameter jet must remain finite")
        return value

    def parameter_jet(self, n: int, m: int, eta: float) -> float:
        """Return ``partial_eta^m u0[n](eta)`` in binary64 when representable."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta(eta)
        if n != 1:
            return 0.0
        if m == 0:
            # Keep the landed zeroth reference row authoritative.
            return self.reference.u_coefficient(n, eta)
        normalized = self.normalized_parameter_taylor_coefficient(n, m, eta)
        try:
            value = float(math.factorial(m)) * normalized
        except OverflowError as exc:
            raise ArithmeticError("axial reference parameter jet is not representable in binary64") from exc
        if not math.isfinite(value):
            raise ArithmeticError("axial reference parameter jet is not representable in binary64")
        return value

    def normalized_jet(self, n: int, m: int, eta: float) -> float:
        """Return the stored-coordinate normalization ``jet/weight``."""

        return self.parameter_jet(n, m, eta) / coefficient_weight(self.epsilon, n, m)


def actual_schedule_axial_reference_jets(
    data: TailData,
    j: float,
) -> ActualScheduleAxialReferenceJets:
    """Build the actual-schedule axial derivative family with canonical radius."""

    reference = actual_schedule_reference_pair(data, j)
    analytic = certify_actual_schedule_analytic_inputs(data, j)
    epsilon = float(analytic.neighborhood.radius) / 2.0
    return ActualScheduleAxialReferenceJets(reference=reference, epsilon=epsilon)
