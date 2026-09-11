"""Actual parameter jets for the axial half of the natural-axis reference pair.

The pinned ``AxisContraction.referencePair`` has

    u0 = -(1/2) J_1(inverseL * zStar),

so only radial degree one is nonzero and

    u0[1](eta) = -(1/2) inverseL(eta) * zStar(eta).

PR #82 materialized this coefficient only at parameter order zero.  This module
uses the actual ``SchedulePressure.axisPressure`` derivative chain from
:mod:`schedule_axis_pressure_jets` and exact truncated Taylor algebra for the
polynomial/rational ``NaturalAxisData`` fields to produce arbitrary finite
eta-derivative jets.  No pressure surrogate, fitted coefficient, or caller
chosen sigma/Lambda/C enters this construction.

This closes the reference-pair parameter-jet gap only.  A compatible global
``AxisCoefficientSpace`` backend, ``coefficientOperators``,
``naturalRemainder`` and the genuine Picard fixed point are still unresolved;
therefore ``paper_exact`` remains false.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .axis_reference_angular_jets import axis_weight
from .axis_reference_pair import ActualScheduleReferencePair, actual_schedule_reference_pair
from .natural_axis import A
from .outgoing_tail import TailData
from .schedule_analytic_neighborhood import (
    ActualScheduleAnalyticInputCertificate,
    certify_actual_schedule_analytic_inputs,
)
from .schedule_axis_pressure_jets import axis_pressure_normalized_taylor


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _add(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
    if len(a) != len(b):
        raise ValueError("Taylor jets must have the same order")
    return tuple(x + y for x, y in zip(a, b))


def _scale(a: tuple[float, ...], c: float) -> tuple[float, ...]:
    return tuple(float(c) * x for x in a)


def _mul(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
    if len(a) != len(b):
        raise ValueError("Taylor jets must have the same order")
    return tuple(sum(a[k] * b[n - k] for k in range(n + 1)) for n in range(len(a)))


def _reciprocal(a: tuple[float, ...]) -> tuple[float, ...]:
    if not a or not math.isfinite(a[0]) or a[0] == 0.0:
        raise ArithmeticError("Taylor reciprocal requires a finite nonzero constant term")
    out = [1.0 / a[0]]
    for n in range(1, len(a)):
        out.append(-sum(a[k] * out[n - k] for k in range(1, n + 1)) / a[0])
    if not all(math.isfinite(x) for x in out):
        raise ArithmeticError("Taylor reciprocal produced a non-finite coefficient")
    return tuple(out)


def _constant(value: float, order: int) -> tuple[float, ...]:
    return (float(value),) + (0.0,) * order


def _eta_taylor(eta: float, order: int) -> tuple[float, ...]:
    eta = float(eta)
    if not math.isfinite(eta):
        raise ValueError("eta must be finite")
    out = [0.0] * (order + 1)
    out[0] = eta
    if order >= 1:
        out[1] = 1.0
    return tuple(out)


def _z_star_normalized_taylor(
    data: TailData,
    j: float,
    eta: float,
    order: int,
) -> tuple[float, ...]:
    """Normalized Taylor jet of the pinned ``NaturalAxisData.Z`` field."""

    order = _index(order, "order")
    eta_jet = _eta_taylor(eta, order)
    one = _constant(1.0, order)
    u = _add(_scale(eta_jet, 4.0), _constant(j, order))
    eta2 = _mul(eta_jet, eta_jet)
    d_jet = _add(one, _scale(eta2, -1.0))
    h_jet = _add(
        _scale(eta_jet, 0.5 - data.h),
        _mul(d_jet, u),
    )

    # Z uses p and p'.  To obtain Z through order m we therefore need the
    # pressure Taylor series through order m+1.
    pressure_plus_one = axis_pressure_normalized_taylor(data, eta, order + 1)
    pressure = pressure_plus_one[: order + 1]
    pressure_derivative = tuple(
        (m + 1) * pressure_plus_one[m + 1] for m in range(order + 1)
    )

    a = A(data.h)
    one_minus_2_eta_u = _add(one, _scale(_mul(eta_jet, u), -2.0))
    nonlinear = _scale(_mul(one_minus_2_eta_u, u), -a)
    h_term = _scale(h_jet, -4.0)
    pressure_derivative_term = _scale(_mul(d_jet, pressure_derivative), -1.0)
    pressure_term = _scale(_mul(eta_jet, pressure), 4.0 * a)
    z = _add(_add(nonlinear, h_term), _add(pressure_derivative_term, pressure_term))
    if not all(math.isfinite(x) for x in z):
        raise ArithmeticError("schedule NaturalAxisData.Z Taylor jet must remain finite")
    return z


@dataclass(frozen=True)
class ActualScheduleAxialReferenceJets:
    """Lazy actual derivative family for the axial ``referencePair`` field."""

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
        """Return ``u0[n]^(m)(eta)/m!`` for ``m=0..order``."""

        n = _index(n, "n")
        order = _index(order, "order")
        # Reuse the landed coefficient-window guard even though the axial
        # reference itself does not depend on sigma.
        self.reference.inverse_l(eta)
        if n != 1:
            return (0.0,) * (order + 1)

        eta_jet = _eta_taylor(eta, order)
        eta2 = _mul(eta_jet, eta_jet)
        l_jet = _add(_constant(1.0, order), _scale(eta2, -2.0 * self.reference.data.h))
        inverse_l = _reciprocal(l_jet)
        z_star = _z_star_normalized_taylor(
            self.reference.data,
            self.reference.j,
            eta,
            order,
        )
        out = _scale(_mul(inverse_l, z_star), -0.5)
        if not all(math.isfinite(x) for x in out):
            raise ArithmeticError("axial reference Taylor jet must remain finite")
        return out

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
        """Return the pointwise stored RawJets coordinate ``jet/weight``."""

        return self.parameter_jet(n, m, eta) / axis_weight(self.epsilon, n, m)


def actual_schedule_axial_reference_jets(
    data: TailData,
    j: float,
) -> ActualScheduleAxialReferenceJets:
    """Build axial reference jets from the actual theorem schedule chain."""

    reference = actual_schedule_reference_pair(data, j)
    analytic_inputs = certify_actual_schedule_analytic_inputs(data, j)
    return ActualScheduleAxialReferenceJets(reference=reference, analytic_inputs=analytic_inputs)
