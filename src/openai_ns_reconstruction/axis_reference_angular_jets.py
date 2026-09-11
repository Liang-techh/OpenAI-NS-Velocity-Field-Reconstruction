"""Actual parameter jets for the angular half of the natural-axis reference pair.

The pinned ``AxisCoefficientSpace`` stores *actual* parameter derivatives, not
independent arrays.  PR #82 materialized only the zeroth parameter jet of
``AxisContraction.referencePair``.  This module advances that boundary for the
angular component ``phi0 = naturalResolvent(one)``.

``AxisReference.reference_coefficient`` gives

    phi0[n](eta) = (-chi(eta)/2)^n / (n! (n+1)!).

Here ``chi = H^2/(H^2+sigma^2)`` and ``H`` is a cubic polynomial.  We evaluate
arbitrary eta derivatives by truncated Taylor algebra, using the exact radial
recurrence rather than finite differences or fitted coefficients.  The
factory obtains ``sigma`` and the coefficient-space ``epsilon=rho/2`` only
from the landed actual-SchedulePressure certificate chain.

This is still ``formal-structure``: the axial reference jets (which require
higher actual derivatives of the schedule pressure), the global AxisSpace norm
realization, coefficient operators, ``naturalRemainder``, and the Picard fixed
point remain unresolved.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .axis_reference_pair import ActualScheduleReferencePair, actual_schedule_reference_pair, radial_divisor
from .natural_axis import D, H
from .outgoing_tail import TailData
from .schedule_analytic_neighborhood import (
    ActualScheduleAnalyticInputCertificate,
    certify_actual_schedule_analytic_inputs,
)


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _mul(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
    """Multiply normalized Taylor coefficients: a[k]=f^(k)/k!."""

    if len(a) != len(b):
        raise ValueError("Taylor jets must have the same order")
    n = len(a)
    return tuple(sum(a[k] * b[r - k] for k in range(r + 1)) for r in range(n))


def _scale(a: tuple[float, ...], c: float) -> tuple[float, ...]:
    return tuple(float(c) * x for x in a)


def _reciprocal(a: tuple[float, ...]) -> tuple[float, ...]:
    """Reciprocal normalized Taylor jet, solved coefficient by coefficient."""

    if not a or a[0] == 0.0 or not math.isfinite(a[0]):
        raise ArithmeticError("Taylor reciprocal requires a finite nonzero constant term")
    out = [1.0 / a[0]]
    for r in range(1, len(a)):
        out.append(-sum(a[k] * out[r - k] for k in range(1, r + 1)) / a[0])
    if not all(math.isfinite(x) for x in out):
        raise ArithmeticError("Taylor reciprocal produced a non-finite coefficient")
    return tuple(out)


def _h_taylor(h: float, j: float, eta: float, order: int) -> tuple[float, ...]:
    """Normalized Taylor jet of the exact cubic ``NaturalAxisData.H``."""

    order = _index(order, "order")
    eta = float(eta)
    if not math.isfinite(eta):
        raise ValueError("eta must be finite")
    values = [0.0] * (order + 1)
    values[0] = H(h, j, eta)
    if order >= 1:
        values[1] = D(h) + 4.0 - 12.0 * eta * eta - 2.0 * j * eta
    if order >= 2:
        values[2] = -12.0 * eta - j  # H''/2!
    if order >= 3:
        values[3] = -4.0  # H'''/3!
    return tuple(values)


def _chi_taylor(h: float, j: float, sigma: float, eta: float, order: int) -> tuple[float, ...]:
    Hjet = _h_taylor(h, j, eta, order)
    square = _mul(Hjet, Hjet)
    denominator = list(square)
    denominator[0] += float(sigma) * float(sigma)
    return _mul(square, _reciprocal(tuple(denominator)))


def axis_weight(epsilon: float, n: int, m: int) -> float:
    """Pinned ``AxisWeightEstimates.weight epsilon n m`` in binary64."""

    n = _index(n, "n")
    m = _index(m, "m")
    epsilon = float(epsilon)
    if not math.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("epsilon must be finite and positive")
    try:
        value = (
            (1.0 / 20.0) ** n
            * epsilon ** (-m)
            * float(math.factorial(m))
            * float(math.comb(n + m, m))
            / ((n + 1.0) ** 2 * (m + 1.0) ** 2)
        )
    except (OverflowError, ValueError) as exc:
        raise ArithmeticError("Axis weight is not representable in binary64") from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ArithmeticError("Axis weight underflowed or overflowed binary64")
    return value


@dataclass(frozen=True)
class ActualScheduleAngularReferenceJets:
    """Lazy actual derivative family for the angular ``referencePair`` field."""

    reference: ActualScheduleReferencePair
    analytic_inputs: ActualScheduleAnalyticInputCertificate

    def __post_init__(self) -> None:
        if self.analytic_inputs.low_z.parameters != self.reference.margin_witness.parameters:
            raise ValueError("reference and analytic-input certificates must use the same parameters")
        if self.analytic_inputs.neighborhood.sigma != self.reference.sigma:
            raise ValueError("reference and analytic-input certificates must use the same sigma")

    @property
    def epsilon(self) -> float:
        return self.analytic_inputs.neighborhood.epsilon

    @property
    def paper_exact(self) -> bool:
        return False

    def normalized_taylor(self, n: int, eta: float, order: int) -> tuple[float, ...]:
        """Return ``phi0[n]^(m)(eta)/m!`` for ``m=0..order``."""

        n = _index(n, "n")
        order = _index(order, "order")
        # Reuse the landed window guard and actual theorem-selected sigma.
        self.reference.chi0(eta)
        c = _scale(
            _chi_taylor(self.reference.data.h, self.reference.j, self.reference.sigma, eta, order),
            -0.5,
        )
        phi = (1.0,) + (0.0,) * order
        for k in range(n):
            phi = _scale(_mul(c, phi), 1.0 / radial_divisor(2, k))
        if not all(math.isfinite(x) for x in phi):
            raise ArithmeticError("angular reference Taylor jet must remain finite")
        return phi

    def parameter_jet(self, n: int, m: int, eta: float) -> float:
        """Return the actual ``m``-th eta derivative of radial coefficient ``n``."""

        m = _index(m, "m")
        coefficient = self.normalized_taylor(n, eta, m)[m]
        try:
            value = coefficient * float(math.factorial(m))
        except OverflowError as exc:
            raise ArithmeticError("parameter jet is not representable in binary64") from exc
        if not math.isfinite(value):
            raise ArithmeticError("parameter jet must remain finite")
        return value

    def normalized_axis_coordinate(self, n: int, m: int, eta: float) -> float:
        """Return the stored RawJets coordinate ``jet/weight`` at this point.

        This is a pointwise coordinate only; this module does not claim a new
        global norm proof beyond the already-landed resolvent majorant.
        """

        return self.parameter_jet(n, m, eta) / axis_weight(self.epsilon, n, m)


def actual_schedule_angular_reference_jets(
    data: TailData,
    j: float,
) -> ActualScheduleAngularReferenceJets:
    """Build angular reference jets with no caller-supplied sigma/rho/epsilon."""

    reference = actual_schedule_reference_pair(data, j)
    analytic_inputs = certify_actual_schedule_analytic_inputs(data, j)
    return ActualScheduleAngularReferenceJets(reference=reference, analytic_inputs=analytic_inputs)
