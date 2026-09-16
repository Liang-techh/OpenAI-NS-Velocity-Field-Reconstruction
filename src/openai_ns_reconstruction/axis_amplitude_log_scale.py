"""Exact phase/common-log-scale bookkeeping for the actual amplitude.

This module binds one exact rational phase-kernel identity to the affine
``RationalLogAmplitudeEnclosure`` selected by the caller.  It records the
finite Decimal midpoint that a numerical adapter actually returned and keeps
the exact midpoint discrepancy visible as a Fraction interval.  The later
coefficient product and summation layers are deliberately outside this type.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal
from fractions import Fraction

from .axis_phase_integral import (
    _WINDOW_LEFT,
    _WINDOW_RIGHT,
)
from .axis_phase_log_enclosure import RationalLogAmplitudeEnclosure


def _require_fraction(value: object, name: str) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be a Fraction")
    return value


def _require_decimal(value: object, name: str) -> Decimal:
    if not isinstance(value, Decimal):
        raise TypeError(f"{name} must be a Decimal")
    if not value.is_finite():
        raise ValueError(f"{name} must be finite")
    return value


def _require_power(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("q must be a nonnegative integer")
    return value


def _validate_kernel_identity(
    h: object,
    j: object,
    sigma: object,
    eta: object,
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    h = _require_fraction(h, "h")
    j = _require_fraction(j, "j")
    sigma = _require_fraction(sigma, "sigma")
    eta = _require_fraction(eta, "eta")
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    if not _WINDOW_LEFT <= eta <= _WINDOW_RIGHT:
        raise ValueError("eta must lie in the pinned window [-11/10,11/10]")
    return h, j, sigma, eta


@dataclass(frozen=True)
class AmplitudeLogSource:
    """One actual finite amplitude log midpoint with exact source deltas.

    ``midpoint`` is a finite Decimal returned by a numerical adapter.  It is
    intentionally not required to lie inside ``enclosure``: Decimal rounding
    can place a narrow rounded midpoint just outside its exact Fraction
    interval.  The two delta fields retain that discrepancy exactly.
    """

    h: Fraction
    j: Fraction
    sigma: Fraction
    eta: Fraction
    enclosure: RationalLogAmplitudeEnclosure
    midpoint: Decimal
    delta_lower: Fraction = field(init=False)
    delta_upper: Fraction = field(init=False)

    def __post_init__(self) -> None:
        h, j, sigma, eta = _validate_kernel_identity(
            self.h,
            self.j,
            self.sigma,
            self.eta,
        )
        if not isinstance(self.enclosure, RationalLogAmplitudeEnclosure):
            raise TypeError(
                "enclosure must be RationalLogAmplitudeEnclosure"
            )
        midpoint = _require_decimal(self.midpoint, "midpoint")
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "j", j)
        object.__setattr__(self, "sigma", sigma)
        object.__setattr__(self, "eta", eta)
        object.__setattr__(self, "midpoint", midpoint)
        object.__setattr__(
            self,
            "delta_lower",
            self.enclosure.lower - Fraction(midpoint),
        )
        object.__setattr__(
            self,
            "delta_upper",
            self.enclosure.upper - Fraction(midpoint),
        )

    @property
    def coefficient_arithmetic_certified(self) -> bool:
        return False

    @property
    def paper_exact(self) -> bool:
        return False

    def assert_compatible(self, other: "AmplitudeLogSource") -> None:
        """Require exact source, interval, scalar, and midpoint identity."""

        if not isinstance(other, AmplitudeLogSource):
            raise TypeError("other must be AmplitudeLogSource")
        for name in ("h", "j", "sigma", "eta"):
            if getattr(self, name) != getattr(other, name):
                raise ValueError(f"amplitude log source {name} mismatch")
        for name in ("phase_lower", "phase_upper", "lower", "upper", "Lambda", "log_C"):
            if getattr(self.enclosure, name) != getattr(other.enclosure, name):
                raise ValueError(f"amplitude log source enclosure {name} mismatch")
        if self.midpoint != other.midpoint:
            raise ValueError("amplitude log source midpoint mismatch")

    def power(
        self,
        q: int,
        returned_log_scale: Decimal,
    ) -> "AmplitudePowerLogScale":
        """Return the exact common-log delta interval for an amplitude power."""

        q = _require_power(q)
        returned_log_scale = _require_decimal(
            returned_log_scale,
            "returned_log_scale",
        )
        if q == 0 and returned_log_scale != Decimal(0):
            raise ValueError("q=0 requires returned_log_scale exactly 0")
        return AmplitudePowerLogScale(
            source=self,
            q=q,
            returned_log_scale=returned_log_scale,
        )


@dataclass(frozen=True)
class AmplitudePowerLogScale:
    """Exact common-log delta interval for one nonnegative amplitude power."""

    source: AmplitudeLogSource
    q: int
    returned_log_scale: Decimal
    delta_lower: Fraction = field(init=False)
    delta_upper: Fraction = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.source, AmplitudeLogSource):
            raise TypeError("source must be AmplitudeLogSource")
        q = _require_power(self.q)
        returned_log_scale = _require_decimal(
            self.returned_log_scale,
            "returned_log_scale",
        )
        if q == 0 and returned_log_scale != Decimal(0):
            raise ValueError("q=0 requires returned_log_scale exactly 0")
        returned = Fraction(returned_log_scale)
        object.__setattr__(self, "q", q)
        object.__setattr__(self, "returned_log_scale", returned_log_scale)
        object.__setattr__(
            self,
            "delta_lower",
            q * self.source.enclosure.lower - returned,
        )
        object.__setattr__(
            self,
            "delta_upper",
            q * self.source.enclosure.upper - returned,
        )

    @property
    def coefficient_arithmetic_certified(self) -> bool:
        return False

    @property
    def paper_exact(self) -> bool:
        return False


__all__ = ["AmplitudeLogSource", "AmplitudePowerLogScale"]
