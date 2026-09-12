"""Actual SchedulePressure fixed coefficient data for ``naturalRemainder``.

The pinned ``NaturalAxisCoefficients.CoefficientFamily.axisData`` record contains
three scalars ``A, D, h`` and ten radially constant coefficient fields:
``one, eta, d, inverseL, uStar, uStarEta, wStar, hStar,
normalizedGradient, zStar``.  ``naturalRemainder`` cannot be instantiated from
the landed reference pair and ``coefficientOperators`` until these fixed fields
live in the same executable coefficient-jet representation.

This module materializes exactly that representation from the already landed
actual SchedulePressure construction.  Every field has radial degree zero, as
proved by ``CoefficientFamily.coefficient_eq``; parameter derivatives are
computed by finite truncated-Taylor algebra from the official real-field
formulas, with ``zStar`` using the actual schedule-pressure Taylor chain.

This is still a representation-level object.  It does not certify the global
all-index weighted ``AxisSpace`` norm, execute ``naturalRemainder``, or claim a
paper-exact leading profile.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

from .axis_coefficient_reference_state import ActualScheduleReferenceAxisState, AxisCoefficientJetState
from .natural_axis import A, D
from .schedule_axis_pressure_jets import axis_pressure_normalized_taylor


NormalizedTaylorProvider = Callable[[float, int], tuple[float, ...]]


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


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


def _add(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
    if len(a) != len(b):
        raise ValueError("Taylor jets must have the same order")
    return tuple(x + y for x, y in zip(a, b))


def _scale(a: tuple[float, ...], scalar: float) -> tuple[float, ...]:
    return tuple(float(scalar) * x for x in a)


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


def _base_taylors(h: float, j: float, eta: float, order: int) -> dict[str, tuple[float, ...]]:
    eta_jet = _eta_taylor(eta, order)
    one = _constant(1.0, order)
    eta2 = _mul(eta_jet, eta_jet)
    d_jet = _add(one, _scale(eta2, -1.0))
    l_jet = _add(one, _scale(eta2, -2.0 * h))
    u_jet = _add(_scale(eta_jet, 4.0), _constant(j, order))
    h_jet = _add(_scale(eta_jet, D(h)), _mul(d_jet, u_jet))
    w_jet = _add(
        _add(one, _scale(d_jet, -4.0)),
        _scale(_mul(_mul(_constant(D(h), order), eta_jet), u_jet), -2.0),
    )
    return {
        "one": one,
        "eta": eta_jet,
        "d": d_jet,
        "L": l_jet,
        "uStar": u_jet,
        "hStar": h_jet,
        "wStar": w_jet,
    }


def _z_star_taylor(data, j: float, eta: float, order: int) -> tuple[float, ...]:
    base = _base_taylors(data.h, j, eta, order)
    eta_jet = base["eta"]
    one = base["one"]
    u_jet = base["uStar"]
    d_jet = base["d"]
    h_jet = base["hStar"]

    pressure_plus_one = axis_pressure_normalized_taylor(data, eta, order + 1)
    pressure = pressure_plus_one[: order + 1]
    pressure_derivative = tuple(
        (m + 1) * pressure_plus_one[m + 1] for m in range(order + 1)
    )

    one_minus_2_eta_u = _add(one, _scale(_mul(eta_jet, u_jet), -2.0))
    nonlinear = _scale(_mul(one_minus_2_eta_u, u_jet), -A(data.h))
    h_term = _scale(h_jet, -4.0)
    pressure_derivative_term = _scale(_mul(d_jet, pressure_derivative), -1.0)
    pressure_term = _scale(_mul(eta_jet, pressure), 4.0 * A(data.h))
    result = _add(_add(nonlinear, h_term), _add(pressure_derivative_term, pressure_term))
    if not all(math.isfinite(x) for x in result):
        raise ArithmeticError("zStar Taylor jet must remain finite")
    return result


def _field_taylor(
    name: str,
    *,
    data,
    j: float,
    sigma: float,
    eta: float,
    order: int,
) -> tuple[float, ...]:
    order = _index(order, "order")
    base = _base_taylors(data.h, j, eta, order)
    if name in {"one", "eta", "d", "uStar", "wStar", "hStar"}:
        return base[name]
    if name == "inverseL":
        return _reciprocal(base["L"])
    if name == "uStarEta":
        return _constant(4.0, order)
    if name == "normalizedGradient":
        denominator = _mul(base["hStar"], base["hStar"])
        denominator = list(denominator)
        denominator[0] += sigma * sigma
        gradient = _scale(
            _mul(_mul(base["L"], base["hStar"]), _reciprocal(tuple(denominator))),
            -1.0,
        )
        if not all(math.isfinite(x) for x in gradient):
            raise ArithmeticError("normalizedGradient Taylor jet must remain finite")
        return gradient
    if name == "zStar":
        return _z_star_taylor(data, j, eta, order)
    raise KeyError(f"unknown AxisData field {name!r}")


def _row_zero_state(
    *,
    epsilon: float,
    origin: str,
    taylor: NormalizedTaylorProvider,
) -> AxisCoefficientJetState:
    def jet(n: int, m: int, eta: float) -> float:
        n = _index(n, "n")
        m = _index(m, "m")
        if n != 0:
            return 0.0
        normalized = taylor(float(eta), m)[m]
        try:
            value = normalized * float(math.factorial(m))
        except OverflowError as exc:
            raise ArithmeticError("parameter jet is not representable in binary64") from exc
        if not math.isfinite(value):
            raise ArithmeticError("parameter jet must remain finite")
        return value

    return AxisCoefficientJetState(epsilon=epsilon, origin=origin, _jet_provider=jet)


LEAN_AXIS_DATA_FIELD_NAMES = (
    "A",
    "D",
    "h",
    "one",
    "eta",
    "d",
    "inverseL",
    "uStar",
    "uStarEta",
    "wStar",
    "hStar",
    "normalizedGradient",
    "zStar",
)


@dataclass(frozen=True)
class ActualScheduleAxisCoefficientData:
    """Executable image of pinned ``CoefficientFamily.axisData``."""

    A: float
    D: float
    h: float
    j: float
    sigma: float
    epsilon: float
    one: AxisCoefficientJetState
    eta: AxisCoefficientJetState
    d: AxisCoefficientJetState
    inverseL: AxisCoefficientJetState
    uStar: AxisCoefficientJetState
    uStarEta: AxisCoefficientJetState
    wStar: AxisCoefficientJetState
    hStar: AxisCoefficientJetState
    normalizedGradient: AxisCoefficientJetState
    zStar: AxisCoefficientJetState

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def lean_field_names(self) -> tuple[str, ...]:
        return LEAN_AXIS_DATA_FIELD_NAMES


def actual_schedule_axis_coefficient_data(
    reference: ActualScheduleReferenceAxisState,
) -> ActualScheduleAxisCoefficientData:
    """Materialize the actual fixed coefficient fields at the reference epsilon.

    The only input is the landed actual SchedulePressure reference-state object;
    callers cannot supply replacement pressure data, sigma, epsilon, or field
    tables.
    """

    if not isinstance(reference, ActualScheduleReferenceAxisState):
        raise TypeError("reference must be ActualScheduleReferenceAxisState")
    data = reference.reference.data
    j = float(reference.reference.j)
    sigma = float(reference.reference.sigma)
    epsilon = float(reference.epsilon)

    if reference.phi.epsilon != epsilon or reference.u.epsilon != epsilon:
        raise ValueError("reference pair must use exactly one coefficient epsilon")

    states: dict[str, AxisCoefficientJetState] = {}
    for name in (
        "one",
        "eta",
        "d",
        "inverseL",
        "uStar",
        "uStarEta",
        "wStar",
        "hStar",
        "normalizedGradient",
        "zStar",
    ):
        states[name] = _row_zero_state(
            epsilon=epsilon,
            origin=f"actual SchedulePressure AxisData.{name}",
            taylor=lambda eta, order, field=name: _field_taylor(
                field,
                data=data,
                j=j,
                sigma=sigma,
                eta=eta,
                order=order,
            ),
        )

    return ActualScheduleAxisCoefficientData(
        A=A(data.h),
        D=D(data.h),
        h=float(data.h),
        j=j,
        sigma=sigma,
        epsilon=epsilon,
        one=states["one"],
        eta=states["eta"],
        d=states["d"],
        inverseL=states["inverseL"],
        uStar=states["uStar"],
        uStarEta=states["uStarEta"],
        wStar=states["wStar"],
        hStar=states["hStar"],
        normalizedGradient=states["normalizedGradient"],
        zStar=states["zStar"],
    )
