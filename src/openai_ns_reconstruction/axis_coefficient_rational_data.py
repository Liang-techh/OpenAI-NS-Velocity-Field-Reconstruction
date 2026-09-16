"""Exact rational scalar jets for the actual natural-axis coefficient data.

``ActualScheduleAxisCoefficientData`` is the sole anchor.  Its selected
``h``, ``j``, and ``sigma`` values, and requested eta values, are interpreted
as exact binary64 rationals.  The fixed coefficient fields are then evaluated
by normalized Taylor algebra over :class:`fractions.Fraction`; no pressure
``zStar`` field is fabricated because its schedule-pressure dependency is
outside this provider.

This is an exact representation for the selected finite binary inputs.  It is
not a theorem-parameter selection proof, a pressure certificate, a global
``AxisSpace`` result, or a complete coefficient-space reconstruction.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import (
    Decimal,
    ROUND_CEILING,
    ROUND_FLOOR,
    ROUND_HALF_EVEN,
    localcontext,
)
from fractions import Fraction
from functools import lru_cache
import math

from .axis_coefficient_data import ActualScheduleAxisCoefficientData


_DECIMAL_PRECISION = 96
_WINDOW_LEFT = Fraction(-11, 10)
_WINDOW_RIGHT = Fraction(11, 10)

_FIXED_FIELD_NAMES = frozenset(
    {
        "one",
        "eta",
        "d",
        "inverseL",
        "uStar",
        "uStarEta",
        "wStar",
        "hStar",
        "normalizedGradient",
        "chi",
    }
)


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _binary_fraction(value: float, name: str) -> Fraction:
    """Convert one finite float to its exact binary rational value."""

    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise TypeError(f"{name} must be a real float")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return Fraction.from_float(value)


def _eta_fraction(value: float) -> Fraction:
    eta = _binary_fraction(value, "eta")
    if not _WINDOW_LEFT <= eta <= _WINDOW_RIGHT:
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]")
    return eta


def _constant(value: Fraction, order: int) -> tuple[Fraction, ...]:
    return (value,) + (Fraction(0),) * order


def _eta_taylor(eta: Fraction, order: int) -> tuple[Fraction, ...]:
    values = [Fraction(0)] * (order + 1)
    values[0] = eta
    if order >= 1:
        values[1] = Fraction(1)
    return tuple(values)


def _add(
    left: tuple[Fraction, ...],
    right: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    if len(left) != len(right):
        raise ValueError("Taylor jets must have the same order")
    return tuple(a + b for a, b in zip(left, right))


def _scale(
    values: tuple[Fraction, ...], scalar: Fraction,
) -> tuple[Fraction, ...]:
    return tuple(scalar * value for value in values)


def _mul(
    left: tuple[Fraction, ...],
    right: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    if len(left) != len(right):
        raise ValueError("Taylor jets must have the same order")
    return tuple(
        sum((left[k] * right[n - k] for k in range(n + 1)), Fraction(0))
        for n in range(len(left))
    )


def _reciprocal(values: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    if not values or values[0] == 0:
        raise ArithmeticError("Taylor reciprocal requires a finite nonzero constant term")
    out = [Fraction(1, 1) / values[0]]
    for n in range(1, len(values)):
        numerator = sum(
            (values[k] * out[n - k] for k in range(1, n + 1)),
            Fraction(0),
        )
        out.append(-numerator / values[0])
    return tuple(out)


def _base_taylors(
    h: Fraction,
    j: Fraction,
    eta: Fraction,
    order: int,
) -> dict[str, tuple[Fraction, ...]]:
    """Return the normalized Taylor families from ``axis_coefficient_data``."""

    eta_jet = _eta_taylor(eta, order)
    one = _constant(Fraction(1), order)
    D = Fraction(1, 2) - h
    eta2 = _mul(eta_jet, eta_jet)
    d = _add(one, _scale(eta2, Fraction(-1)))
    L = _add(one, _scale(eta2, -2 * h))
    u = _add(_scale(eta_jet, Fraction(4)), _constant(j, order))
    h_star = _add(
        _scale(eta_jet, D),
        _mul(d, u),
    )
    w_star = _add(
        _add(one, _scale(d, Fraction(-4))),
        _scale(
            _mul(_mul(_constant(D, order), eta_jet), u),
            Fraction(-2),
        ),
    )
    return {
        "one": one,
        "eta": eta_jet,
        "d": d,
        "L": L,
        "uStar": u,
        "hStar": h_star,
        "wStar": w_star,
    }


def _chi_taylor(
    h_star: tuple[Fraction, ...],
    sigma: Fraction,
) -> tuple[Fraction, ...]:
    """Return normalized Taylor coefficients of ``H^2/(H^2+sigma^2)``."""

    square = _mul(h_star, h_star)
    denominator = list(square)
    denominator[0] += sigma * sigma
    reciprocal = _reciprocal(tuple(denominator))
    out = [square[0] / denominator[0]]
    # ``chi = 1 - sigma^2/(H^2+sigma^2)`` preserves the exact positive-order
    # coefficients without subtracting nearly equal rational terms.
    out.extend(-sigma * sigma * reciprocal[n] for n in range(1, len(square)))
    return tuple(out)


@lru_cache(maxsize=64)
def _normalized_field_taylor(
    name: str,
    h: Fraction,
    j: Fraction,
    sigma: Fraction,
    eta: Fraction,
    order: int,
) -> tuple[Fraction, ...]:
    """Cached normalized Taylor coefficients for one finite scalar request."""

    base = _base_taylors(h, j, eta, order)
    if name in {"one", "eta", "d", "uStar", "wStar", "hStar"}:
        return base[name]
    if name == "inverseL":
        return _reciprocal(base["L"])
    if name == "uStarEta":
        return _constant(Fraction(4), order)
    if name == "chi":
        return _chi_taylor(base["hStar"], sigma)
    if name == "normalizedGradient":
        denominator = list(_mul(base["hStar"], base["hStar"]))
        denominator[0] += sigma * sigma
        return _scale(
            _mul(
                base["L"],
                _mul(base["hStar"], _reciprocal(tuple(denominator))),
            ),
            Fraction(-1),
        )
    raise KeyError(f"unknown rational AxisData field {name!r}")


def _fraction_to_decimal(value: Fraction, rounding: str) -> Decimal:
    with localcontext() as context:
        context.prec = _DECIMAL_PRECISION
        context.rounding = rounding
        return +(Decimal(value.numerator) / Decimal(value.denominator))


@dataclass(frozen=True)
class RationalJetEnclosure:
    """Directed Decimal96 enclosure of one exact rational jet.

    Only ``exact_value`` is accepted by the constructor.  The directed bounds
    are computed from it, so callers cannot forge an enclosure by supplying
    unrelated lower or upper values.
    """

    exact_value: Fraction
    lower: Decimal = field(init=False)
    upper: Decimal = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.exact_value, Fraction):
            raise TypeError("exact_value must be a Fraction")
        lower = _fraction_to_decimal(self.exact_value, ROUND_FLOOR)
        upper = _fraction_to_decimal(self.exact_value, ROUND_CEILING)
        if Fraction(lower) > self.exact_value or Fraction(upper) < self.exact_value:
            raise ArithmeticError("directed Decimal enclosure failed to bracket the rational")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)


@dataclass(frozen=True)
class RationalAxisCoefficientData:
    """Exact rational scalar fields anchored by actual coefficient data."""

    actual_data: ActualScheduleAxisCoefficientData
    _h: Fraction = field(init=False, repr=False)
    _j: Fraction = field(init=False, repr=False)
    _sigma: Fraction = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.actual_data, ActualScheduleAxisCoefficientData):
            raise TypeError("actual_data must be ActualScheduleAxisCoefficientData")
        h = _binary_fraction(self.actual_data.h, "AxisData.h")
        j = _binary_fraction(self.actual_data.j, "AxisData.j")
        sigma = _binary_fraction(self.actual_data.sigma, "AxisData.sigma")
        if sigma <= 0:
            raise ValueError("AxisData.sigma must be positive")
        object.__setattr__(self, "_h", h)
        object.__setattr__(self, "_j", j)
        object.__setattr__(self, "_sigma", sigma)

    @property
    def h(self) -> Fraction:
        return self._h

    @property
    def j(self) -> Fraction:
        return self._j

    @property
    def sigma(self) -> Fraction:
        return self._sigma

    @property
    def A(self) -> Fraction:
        return Fraction(1, 2) + self._h

    @property
    def D(self) -> Fraction:
        return Fraction(1, 2) - self._h

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    def _field_name(self, name: str) -> str:
        if not isinstance(name, str):
            raise TypeError("field name must be a string")
        if name == "zStar":
            raise NotImplementedError(
                "zStar requires the schedule-pressure dependency and is unsupported"
            )
        if name not in _FIXED_FIELD_NAMES:
            raise KeyError(f"unknown rational AxisData field {name!r}")
        return name

    def jet_fraction(self, name: str, n: int, m: int, eta: float) -> Fraction:
        """Return the exact actual ``m``-th eta derivative of a fixed field."""

        name = self._field_name(name)
        n = _index(n, "n")
        m = _index(m, "m")
        eta_fraction = _eta_fraction(eta)
        if n != 0:
            return Fraction(0)
        normalized = _normalized_field_taylor(
            name,
            self._h,
            self._j,
            self._sigma,
            eta_fraction,
            m,
        )
        return normalized[m] * math.factorial(m)

    def angular_reference_jet_fraction(self, n: int, m: int, eta: float) -> Fraction:
        """Return the exact eta derivative of the angular reference coefficient."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta_fraction = _eta_fraction(eta)
        chi = _normalized_field_taylor(
            "chi",
            self._h,
            self._j,
            self._sigma,
            eta_fraction,
            m,
        )
        c = _scale(chi, Fraction(-1, 2))
        power = _constant(Fraction(1), m)
        for _ in range(n):
            power = _mul(c, power)
        denominator = math.factorial(n) * math.factorial(n + 1)
        return power[m] * math.factorial(m) / denominator

    def jet_decimal(self, name: str, n: int, m: int, eta: float) -> Decimal:
        """Return ``jet_fraction`` rounded to nearest Decimal96."""

        return _fraction_to_decimal(
            self.jet_fraction(name, n, m, eta),
            ROUND_HALF_EVEN,
        )

    def angular_reference_jet_decimal(self, n: int, m: int, eta: float) -> Decimal:
        """Return the angular reference jet rounded to nearest Decimal96."""

        return _fraction_to_decimal(
            self.angular_reference_jet_fraction(n, m, eta),
            ROUND_HALF_EVEN,
        )

    def angular_reference_jet_enclosure(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> RationalJetEnclosure:
        """Return a directed Decimal96 enclosure of an angular reference jet."""

        return RationalJetEnclosure(
            self.angular_reference_jet_fraction(n, m, eta)
        )

    def amplitude_power_bell_fraction(
        self,
        power: int,
        order: int,
        eta: float,
        Lambda: Decimal,
    ) -> Fraction:
        """Return ``(d_eta^order a^power) / a^power`` exactly.

        The logarithmic derivative is supplied by the exact rational
        ``normalizedGradient`` field.  The recurrence is the ordinary Bell
        recurrence for full eta derivatives divided by the common amplitude
        power, with no factorial normalization.  ``Lambda`` is retained as an
        exact Decimal rational; no amplitude, phase, or ``C`` is formed.
        """

        power = _index(power, "power")
        order = _index(order, "order")
        _eta_fraction(eta)
        if not isinstance(Lambda, Decimal) or not Lambda.is_finite():
            raise ValueError("Lambda must be a finite Decimal")
        if Lambda <= 0:
            raise ValueError("Lambda must be positive")
        lambda_fraction = Fraction(Lambda)
        q = [
            power
            * lambda_fraction
            * self.jet_fraction("normalizedGradient", 0, k, eta)
            for k in range(order)
        ]
        bell = [Fraction(1)]
        for n in range(order):
            bell.append(
                sum(
                    (
                        math.comb(n, k)
                        * q[k]
                        * bell[n - k]
                        for k in range(n + 1)
                    ),
                    Fraction(0),
                )
            )
        return bell[order]

    def amplitude_power_bell_decimal(
        self,
        power: int,
        order: int,
        eta: float,
        Lambda: Decimal,
    ) -> Decimal:
        """Return the Bell jet rounded to nearest Decimal96."""

        return _fraction_to_decimal(
            self.amplitude_power_bell_fraction(power, order, eta, Lambda),
            ROUND_HALF_EVEN,
        )

    def amplitude_power_bell_enclosure(
        self,
        power: int,
        order: int,
        eta: float,
        Lambda: Decimal,
    ) -> RationalJetEnclosure:
        """Return a directed Decimal96 enclosure of an amplitude Bell jet."""

        return RationalJetEnclosure(
            self.amplitude_power_bell_fraction(power, order, eta, Lambda)
        )

    def jet_enclosure(
        self,
        name: str,
        n: int,
        m: int,
        eta: float,
    ) -> RationalJetEnclosure:
        """Return a directed Decimal96 enclosure of one fixed-field jet."""

        return RationalJetEnclosure(self.jet_fraction(name, n, m, eta))


__all__ = ["RationalAxisCoefficientData", "RationalJetEnclosure"]
