"""Exact rational reference source and pressure coefficient factors.

This module closes the chosen-kernel arithmetic seam at the reference point
``x0``.  It accepts the already-bound
:class:`ActualScheduleAmplitudeLogState`, derives the angular reference
coefficients from its exact rational data, and evaluates the pinned
``a(eta)^2 * phi0(eta)^2`` source and pressure chain entirely with
:class:`fractions.Fraction`.

The common ``a(eta)^2`` magnitude is kept in the source's exact phase/log
interval.  Each returned ``ReferenceCoefficientLogEnclosure`` separately
encloses the exact signed normalized source or pressure factor.  This is a
chosen finite-kernel arithmetic path; it does not certify downstream sums,
global coefficient-space bounds, or a paper-exact reconstruction.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
import math

from .axis_amplitude_derivative_enclosure import (
    RationalLogInterval,
    rational_log_enclosure,
)
from .axis_amplitude_log_scale import AmplitudeLogSource
from .axis_coefficient_amplitude import ActualScheduleAmplitudeLogState


_DEFAULT_LOG_TOLERANCE = Fraction(1, 10**96)
_WINDOW_LEFT = Fraction(-11, 10)
_WINDOW_RIGHT = Fraction(11, 10)


def _index(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _positive_fraction(value: object, name: str) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be a Fraction")
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _positive_cap(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError("max_terms must be a positive integer")
    return value


def _eta(value: object) -> float:
    try:
        eta = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]") from exc
    if not math.isfinite(eta):
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]")
    eta_fraction = Fraction.from_float(eta)
    if not _WINDOW_LEFT <= eta_fraction <= _WINDOW_RIGHT:
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]")
    return eta


def _validate_actual_identity(amplitude: ActualScheduleAmplitudeLogState) -> None:
    """Require that the rational scalar source is the amplitude's actual state."""

    if not isinstance(amplitude, ActualScheduleAmplitudeLogState):
        raise TypeError("amplitude must be ActualScheduleAmplitudeLogState")
    reference = amplitude.reference
    data = amplitude.data
    rational = amplitude.rational_data
    reference_pair = reference.reference

    expected_h = Fraction.from_float(reference_pair.data.h)
    expected_j = Fraction.from_float(reference_pair.j)
    expected_sigma = Fraction.from_float(reference_pair.sigma)
    if rational.h != expected_h or rational.j != expected_j or rational.sigma != expected_sigma:
        raise ValueError("amplitude rational data is incompatible with its actual reference identity")
    if rational.A != Fraction(1, 2) + rational.h:
        raise ValueError("rational A must be the exact half-plus-h scalar")
    if data.epsilon != reference.epsilon:
        raise ValueError("amplitude data must use the actual reference epsilon")


@dataclass(frozen=True)
class ReferenceCoefficientLogEnclosure:
    """Log enclosure for one exact normalized source or pressure coefficient."""

    kind: str
    n: int
    m: int
    exact_factor: Fraction
    source: AmplitudeLogSource
    log_factor: RationalLogInterval | None
    sign: int = field(init=False)
    total_log_lower: Fraction | None = field(init=False)
    total_log_upper: Fraction | None = field(init=False)

    def __post_init__(self) -> None:
        if self.kind not in {"source", "pressure"}:
            raise ValueError("kind must be 'source' or 'pressure'")
        n = _index(self.n, "n")
        m = _index(self.m, "m")
        if not isinstance(self.exact_factor, Fraction):
            raise TypeError("exact_factor must be a Fraction")
        if not isinstance(self.source, AmplitudeLogSource):
            raise TypeError("source must be an AmplitudeLogSource")

        if self.exact_factor == 0:
            if self.log_factor is not None:
                raise ValueError("zero coefficient factor must not carry a log interval")
            sign = 0
            lower = None
            upper = None
        else:
            if not isinstance(self.log_factor, RationalLogInterval):
                raise TypeError("nonzero coefficient factor requires a RationalLogInterval")
            sign = 1 if self.exact_factor > 0 else -1
            lower = 2 * self.source.enclosure.lower + self.log_factor.lower
            upper = 2 * self.source.enclosure.upper + self.log_factor.upper
            if lower > upper:
                raise ArithmeticError("coefficient log bounds are unordered")

        object.__setattr__(self, "n", n)
        object.__setattr__(self, "m", m)
        object.__setattr__(self, "sign", sign)
        object.__setattr__(self, "total_log_lower", lower)
        object.__setattr__(self, "total_log_upper", upper)

    @property
    def total_log_width(self) -> Fraction | None:
        if self.total_log_lower is None or self.total_log_upper is None:
            return None
        return self.total_log_upper - self.total_log_lower

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def downstream_arithmetic_certified(self) -> bool:
        return False


@dataclass
class _ExactReferenceCache:
    """Bounded per-call cache for exact rational convolution terms."""

    amplitude: ActualScheduleAmplitudeLogState
    eta: float
    bells: dict[int, Fraction] = field(default_factory=dict)
    phi: dict[tuple[int, int], Fraction] = field(default_factory=dict)
    fields: dict[tuple[str, int, int], Fraction] = field(default_factory=dict)
    source: dict[tuple[int, int], Fraction] = field(default_factory=dict)

    @property
    def rational(self):
        return self.amplitude.rational_data

    def bell_factor(self, order: int) -> Fraction:
        if order not in self.bells:
            self.bells[order] = self.rational.amplitude_power_bell_fraction(
                2,
                order,
                self.eta,
                self.amplitude.Lambda,
            )
        return self.bells[order]

    def phi_factor(self, n: int, m: int) -> Fraction:
        key = (n, m)
        if key not in self.phi:
            self.phi[key] = self.rational.angular_reference_jet_fraction(
                n,
                m,
                self.eta,
            )
        return self.phi[key]

    def field_factor(self, name: str, n: int, m: int) -> Fraction:
        key = (name, n, m)
        if key not in self.fields:
            self.fields[key] = self.rational.jet_fraction(name, n, m, self.eta)
        return self.fields[key]

    def source_factor(self, n: int, m: int) -> Fraction:
        key = (n, m)
        if key in self.source:
            return self.source[key]

        total = Fraction(0)
        for k in range(m + 1):
            phi_square = Fraction(0)
            remaining = m - k
            for i in range(n + 1):
                for ell in range(remaining + 1):
                    phi_square += (
                        math.comb(remaining, ell)
                        * self.phi_factor(i, ell)
                        * self.phi_factor(n - i, remaining - ell)
                    )
            total += math.comb(m, k) * self.bell_factor(k) * phi_square
        self.source[key] = total
        return total

    def primitive(self, n: int, m: int) -> Fraction:
        if n == 0:
            return Fraction(0)
        return self.source_factor(n - 1, m) / n

    def parameter_primitive(self, n: int, m: int) -> Fraction:
        if n == 0:
            return Fraction(0)
        return self.source_factor(n - 1, m + 1) / n

    def multiply_y(self, n: int, m: int) -> Fraction:
        if n == 0:
            return Fraction(0)
        return self.source_factor(n - 1, m)

    def ordinary_product(
        self,
        name: str,
        wide,
        n: int,
        m: int,
    ) -> Fraction:
        total = Fraction(0)
        for i in range(n + 1):
            for k in range(m + 1):
                total += (
                    math.comb(m, k)
                    * self.field_factor(name, i, k)
                    * wide(n - i, m - k)
                )
        return total

    def pressure_input(self, n: int, m: int) -> Fraction:
        eta_primitive = self.ordinary_product("eta", self.primitive, n, m)
        d_parameter_primitive = self.ordinary_product(
            "d", self.parameter_primitive, n, m
        )
        eta_multiply_y = self.ordinary_product("eta", self.multiply_y, n, m)
        return (
            -4 * self.rational.A * eta_primitive
            + d_parameter_primitive
            - 2 * eta_multiply_y
        )

    def pressure_factor(self, n: int, m: int) -> Fraction:
        if n == 0:
            return Fraction(0)
        return self.pressure_input(n - 1, m) / (n * n)


@dataclass(frozen=True)
class ActualScheduleReferenceRationalPressure:
    """Exact chosen-kernel source and pressure factors at the actual reference."""

    amplitude: ActualScheduleAmplitudeLogState

    def __post_init__(self) -> None:
        _validate_actual_identity(self.amplitude)

    @property
    def epsilon(self) -> float:
        return self.amplitude.epsilon

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    def _cache(self, eta: object) -> _ExactReferenceCache:
        return _ExactReferenceCache(self.amplitude, _eta(eta))

    def source_factor(self, n: int, m: int, eta: float) -> Fraction:
        n = _index(n, "n")
        m = _index(m, "m")
        cache = self._cache(eta)
        return cache.source_factor(n, m)

    def pressure_factor(self, n: int, m: int, eta: float) -> Fraction:
        n = _index(n, "n")
        m = _index(m, "m")
        cache = self._cache(eta)
        return cache.pressure_factor(n, m)

    def normalized_source_factor(self, n: int, m: int, eta: float) -> Fraction:
        return self.source_factor(n, m, eta)

    def normalized_pressure_factor(self, n: int, m: int, eta: float) -> Fraction:
        return self.pressure_factor(n, m, eta)

    def _log_enclosure(
        self,
        kind: str,
        n: int,
        m: int,
        eta: float,
        *,
        absolute_log_factor_tolerance: Fraction,
        max_terms: int,
    ) -> ReferenceCoefficientLogEnclosure:
        if kind not in {"source", "pressure"}:
            raise ValueError("kind must be 'source' or 'pressure'")
        n = _index(n, "n")
        m = _index(m, "m")
        tolerance = _positive_fraction(
            absolute_log_factor_tolerance,
            "absolute_log_factor_tolerance",
        )
        max_terms = _positive_cap(max_terms)
        eta = _eta(eta)

        source = self.amplitude.log_amplitude_source(eta)
        exact_factor = (
            self.source_factor(n, m, eta)
            if kind == "source"
            else self.pressure_factor(n, m, eta)
        )
        if exact_factor == 0:
            return ReferenceCoefficientLogEnclosure(
                kind=kind,
                n=n,
                m=m,
                exact_factor=exact_factor,
                source=source,
                log_factor=None,
            )
        log_factor = rational_log_enclosure(
            abs(exact_factor),
            absolute_tolerance=tolerance,
            max_terms=max_terms,
        )
        return ReferenceCoefficientLogEnclosure(
            kind=kind,
            n=n,
            m=m,
            exact_factor=exact_factor,
            source=source,
            log_factor=log_factor,
        )

    def source_log_enclosure(
        self,
        n: int,
        m: int,
        eta: float,
        *,
        absolute_log_factor_tolerance: Fraction = _DEFAULT_LOG_TOLERANCE,
        max_terms: int = 1024,
    ) -> ReferenceCoefficientLogEnclosure:
        return self._log_enclosure(
            "source",
            n,
            m,
            eta,
            absolute_log_factor_tolerance=absolute_log_factor_tolerance,
            max_terms=max_terms,
        )

    def pressure_log_enclosure(
        self,
        n: int,
        m: int,
        eta: float,
        *,
        absolute_log_factor_tolerance: Fraction = _DEFAULT_LOG_TOLERANCE,
        max_terms: int = 1024,
    ) -> ReferenceCoefficientLogEnclosure:
        return self._log_enclosure(
            "pressure",
            n,
            m,
            eta,
            absolute_log_factor_tolerance=absolute_log_factor_tolerance,
            max_terms=max_terms,
        )

    def log_enclosure(
        self,
        kind: str,
        n: int,
        m: int,
        eta: float,
        *,
        absolute_log_factor_tolerance: Fraction = _DEFAULT_LOG_TOLERANCE,
        max_terms: int = 1024,
    ) -> ReferenceCoefficientLogEnclosure:
        return self._log_enclosure(
            kind,
            n,
            m,
            eta,
            absolute_log_factor_tolerance=absolute_log_factor_tolerance,
            max_terms=max_terms,
        )


# Explicit aliases keep the public name discoverable alongside the existing
# ``ActualScheduleReferenceWidePressureLogState`` naming without accepting a
# second, potentially incompatible angular-square provider.
ActualScheduleRationalPressure = ActualScheduleReferenceRationalPressure


def actual_schedule_reference_rational_pressure(
    amplitude: ActualScheduleAmplitudeLogState,
) -> ActualScheduleReferenceRationalPressure:
    """Bind exact reference source/pressure arithmetic to one actual amplitude."""

    return ActualScheduleReferenceRationalPressure(amplitude=amplitude)


actual_schedule_rational_pressure = actual_schedule_reference_rational_pressure


__all__ = [
    "ActualScheduleReferenceRationalPressure",
    "ActualScheduleRationalPressure",
    "ReferenceCoefficientLogEnclosure",
    "actual_schedule_reference_rational_pressure",
    "actual_schedule_rational_pressure",
]
