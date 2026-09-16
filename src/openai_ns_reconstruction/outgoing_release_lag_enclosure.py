"""Certified rational enclosure of the actual outgoing release lag.

For ``TailData`` write ``h = Fraction.from_float(data.h)`` and
``lambda = Fraction.from_float(data.core.lam)``. The finite release schedule
reduces to a one-dimensional monotone Darboux enclosure, so this module
certifies the selected scalar ``releaseLag(rampEnd)`` without using the
rounded geometry or cached floating-point lag stored on ``TailData``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .axis_amplitude_derivative_enclosure import rational_log_enclosure
from .outgoing_sigma_enclosure import (
    validated_exp_negative,
    validated_exp_positive_on_0_3,
    validated_outgoing_sigma,
)
from .outgoing_tail import TailData
from .rational_interval import (
    RationalInterval,
    ceil_grid as _ceil_grid,
    dyadic_step as _dyadic_step,
    floor_grid as _floor_grid,
    positive_cap as _positive_cap,
    positive_fraction as _positive_fraction,
)


_DEFAULT_TOLERANCE = Fraction(1, 256)
_DEFAULT_MAX_CELLS = 4096
_DEFAULT_MAX_TERMS = 1024
_DEFAULT_MAX_SQUARINGS = 4096


def _round_interval(interval: RationalInterval, step: Fraction) -> RationalInterval:
    lower = _floor_grid(interval.lower, step)
    upper = _ceil_grid(interval.upper, step)
    if lower > upper:
        raise ArithmeticError("rounded release-lag interval became unordered")
    return RationalInterval(lower, upper)


def _round_unit_interval(interval: RationalInterval, step: Fraction) -> RationalInterval:
    rounded = _round_interval(interval, step)
    lower = max(Fraction(0), rounded.lower)
    upper = min(Fraction(1), rounded.upper)
    if lower > upper:
        raise ArithmeticError("rounded unit interval became unordered")
    return RationalInterval(lower, upper)


def _require_actual_data(data: object) -> tuple[TailData, Fraction, Fraction]:
    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    h = Fraction.from_float(data.h)
    lam = Fraction.from_float(data.core.lam)
    if not Fraction(0) < 2 * h < lam < Fraction(1, 10):
        raise ValueError("TailData must satisfy 0 < 2*h < lambda < 1/10")
    return data, h, lam


def _f_value_interval(
    value: Fraction,
    *,
    c: Fraction,
    k: Fraction,
    node_tolerance: Fraction,
    max_terms: int,
    max_squarings: int,
    common_step: Fraction,
) -> RationalInterval:
    """Enclose ``exp(-c*value) + exp(k*value)`` at one exact node value."""

    if not Fraction(0) <= value <= Fraction(1, 2):
        raise ValueError("release-lag primitive value must lie in [0, 1/2]")
    negative = validated_exp_negative(
        c * value,
        absolute_tolerance=node_tolerance,
        max_terms=max_terms,
        max_squarings=max_squarings,
    )
    positive = validated_exp_positive_on_0_3(
        k * value,
        absolute_tolerance=node_tolerance,
        max_terms=max_terms,
        max_squarings=max_squarings,
    )
    raw = RationalInterval(
        negative.lower + positive.lower,
        negative.upper + positive.upper,
    )
    return _round_interval(raw, common_step)


def _integral_f_levels(
    *,
    c: Fraction,
    k: Fraction,
    node_tolerance: Fraction,
    max_cells: int,
    max_terms: int,
    max_squarings: int,
):
    """Yield Darboux bounds for the release integral at powers of two cells."""

    common_step = _dyadic_step(node_tolerance / 8)
    sigma_cache: dict[Fraction, RationalInterval] = {}
    f_cache: dict[Fraction, RationalInterval] = {}

    def sigma(node: Fraction) -> RationalInterval:
        if node not in sigma_cache:
            raw = validated_outgoing_sigma(
                node,
                absolute_tolerance=node_tolerance,
                max_terms=max_terms,
                max_squarings=max_squarings,
            )
            sigma_cache[node] = _round_unit_interval(raw, common_step)
        return sigma_cache[node]

    def f_value(value: Fraction) -> RationalInterval:
        if value not in f_cache:
            f_cache[value] = _f_value_interval(
                value,
                c=c,
                k=k,
                node_tolerance=node_tolerance,
                max_terms=max_terms,
                max_squarings=max_squarings,
                common_step=common_step,
            )
        return f_cache[value]

    cells = 1
    while cells <= max_cells:
        width = Fraction(1, cells)
        sigma_nodes = [sigma(Fraction(index, cells)) for index in range(cells + 1)]

        prefixes = [RationalInterval(Fraction(0), Fraction(0))]
        for index in range(1, cells):
            previous = prefixes[-1]
            lower = previous.lower + width * sigma_nodes[index - 1].lower
            upper = previous.upper + width * sigma_nodes[index].upper
            prefixes.append(
                RationalInterval(
                    max(Fraction(0), lower),
                    min(Fraction(1, 2), upper),
                )
            )
        prefixes.append(RationalInterval(Fraction(1, 2), Fraction(1, 2)))

        f_nodes = []
        for prefix in prefixes:
            f_nodes.append(
                RationalInterval(
                    f_value(prefix.lower).lower,
                    f_value(prefix.upper).upper,
                )
            )

        lower = Fraction(0)
        upper = Fraction(0)
        for index in range(cells):
            lower += width * f_nodes[index].lower
            upper += width * f_nodes[index + 1].upper
        yield RationalInterval(lower, upper), cells
        cells *= 2


@dataclass(frozen=True)
class ReleaseLagEnclosure:
    """Exact rational bounds for ``releaseLag(rampEnd)``."""

    h: Fraction
    lam: Fraction
    lag: RationalInterval
    cells: int

    def __post_init__(self) -> None:
        if not isinstance(self.h, Fraction) or not isinstance(self.lam, Fraction):
            raise TypeError("h and lam must be Fractions")
        if not Fraction(0) < 2 * self.h < self.lam < Fraction(1, 10):
            raise ValueError("h and lam must satisfy 0 < 2*h < lam < 1/10")
        if not isinstance(self.lag, RationalInterval):
            raise TypeError("lag must be a RationalInterval")
        if self.lag.lower <= 0:
            raise ValueError("release lag interval must be positive")
        if isinstance(self.cells, bool) or not isinstance(self.cells, int) or self.cells <= 0:
            raise ValueError("cells must be a positive integer")

    @property
    def release_lag(self) -> RationalInterval:
        return self.lag

    @property
    def width(self) -> Fraction:
        return self.lag.width

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False


def validated_release_lag_enclosure(
    data: TailData,
    *,
    absolute_tolerance: Fraction = _DEFAULT_TOLERANCE,
    max_cells: int = _DEFAULT_MAX_CELLS,
    max_terms: int = _DEFAULT_MAX_TERMS,
    max_squarings: int = _DEFAULT_MAX_SQUARINGS,
) -> ReleaseLagEnclosure:
    """Return a certified rational interval for the actual release lag."""

    data, h, lam = _require_actual_data(data)
    tolerance = _positive_fraction(absolute_tolerance, "absolute_tolerance")
    max_cells = _positive_cap(max_cells, "max_cells")
    max_terms = _positive_cap(max_terms, "max_terms")
    max_squarings = _positive_cap(max_squarings, "max_squarings")

    epsilon = min(tolerance, Fraction(1))
    c = 1 - lam
    k = 1 - h
    a = c / 2
    B = (c + k) / 2
    q0 = (lam - h) / c

    log_interval = rational_log_enclosure(
        1 / h,
        absolute_tolerance=epsilon / 1024,
        max_terms=max_terms,
    )
    L = RationalInterval(4 * log_interval.lower, 4 * log_interval.upper)
    if L.lower <= 0:
        raise ArithmeticError("release-lag log interval lost positivity")

    node_tolerance = min(
        Fraction(1, 1024),
        epsilon / (1024 * (1 + L.upper)),
    )
    exp_a = validated_exp_positive_on_0_3(
        a,
        absolute_tolerance=node_tolerance,
        max_terms=max_terms,
        max_squarings=max_squarings,
    )
    exp_minus_B = validated_exp_negative(
        B,
        absolute_tolerance=node_tolerance,
        max_terms=max_terms,
        max_squarings=max_squarings,
    )

    for integral, cells in _integral_f_levels(
        c=c,
        k=k,
        node_tolerance=node_tolerance,
        max_cells=max_cells,
        max_terms=max_terms,
        max_squarings=max_squarings,
    ):
        inner = RationalInterval(L.lower + integral.lower, L.upper + integral.upper)
        release_integral = RationalInterval(
            exp_a.lower * inner.lower,
            exp_a.upper * inner.upper,
        )
        lag_argument = RationalInterval(
            q0 + 1 + k * release_integral.lower,
            q0 + 1 + k * release_integral.upper,
        )
        lag = RationalInterval(
            exp_minus_B.lower * lag_argument.lower - 1,
            exp_minus_B.upper * lag_argument.upper - 1,
        )
        if lag.lower > 0 and lag.width <= tolerance:
            return ReleaseLagEnclosure(h=h, lam=lam, lag=lag, cells=cells)

    raise ArithmeticError("release lag enclosure did not reach requested tolerance before max_cells")


release_lag_enclosure = validated_release_lag_enclosure
validated_release_lag = validated_release_lag_enclosure


__all__ = [
    "ReleaseLagEnclosure",
    "release_lag_enclosure",
    "validated_release_lag",
    "validated_release_lag_enclosure",
]
