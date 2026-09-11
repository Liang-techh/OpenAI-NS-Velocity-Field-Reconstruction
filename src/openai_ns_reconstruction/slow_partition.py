"""Executable squared slow partitions for Section 6.2.

The paper specifies two smooth squared partitions before Eq. (6.8): a dyadic
partition ``sum_ell chi_ell(q)^2 = 1`` with ``chi_ell`` supported where
``q / 2**(-ell) in [1/2, 2]``, and in every band a product partition
``sum_a chi_{ell,a}(R,Z,T)^2 = 1`` on a mesh of size ``S_*^-3`` whose support
extends by at most one mesh length from its grid point in each coordinate.

The manuscript deliberately leaves the seed bump unspecified.  This module
therefore makes one explicit paper-admissible C-infinity choice and records it
as an implementation choice, not as a uniquely paper-exact transition profile:

    b(s) = exp(-1/(1-s^2)) for |s|<1, and 0 otherwise.

Integer translates are normalized by the square root of their squared sum.
Because ``b`` has support in ``(-1,1)``, at most two integer translates are
nonzero at a point, so the normalization is local and the squared partition
identity can be checked without truncating an infinite sum.
"""
from __future__ import annotations

from dataclasses import dataclass
import itertools
import math

from .charts import DyadicChart
from .slow_labels import SlowBoxEnclosure, SlowLabel


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def seed_bump(s: float) -> float:
    """A fixed nonnegative C-infinity bump supported in ``[-1,1]``."""
    s = _finite(s, "s")
    if abs(s) >= 1.0:
        return 0.0
    return math.exp(-1.0 / (1.0 - s * s))


def local_integer_indices(x: float) -> tuple[int, int]:
    """The only integer translates that can be nonzero at ``x``."""
    x = _finite(x, "x")
    lo = math.floor(x)
    return lo, lo + 1


def normalized_translate(x: float, index: int) -> float:
    """Normalized translate ``b(x-index)/sqrt(sum_j b(x-j)^2)``.

    The denominator is evaluated from the two possible nonzero translates, so
    this is the exact local form of the infinite lattice normalization.
    """
    x = _finite(x, "x")
    if isinstance(index, bool) or not isinstance(index, int):
        raise ValueError("index must be an integer")
    i0, i1 = local_integer_indices(x)
    b0 = seed_bump(x - i0)
    b1 = seed_bump(x - i1)
    denom = math.hypot(b0, b1)
    if denom == 0.0 or not math.isfinite(denom):
        raise ArithmeticError("lattice normalization denominator vanished")
    return seed_bump(x - index) / denom


def dyadic_band_coordinate(q: float) -> float:
    """Return ``-log2(q)`` so band ``ell`` is centered at this coordinate."""
    q = _finite(q, "q")
    if q <= 0.0:
        raise ValueError("q must be positive")
    return -math.log2(q)


def dyadic_cutoff(ell: int, q: float) -> float:
    """Concrete ``chi_ell(q)`` with the paper's ``q/Q in [1/2,2]`` support."""
    if isinstance(ell, bool) or not isinstance(ell, int):
        raise ValueError("ell must be an integer")
    return normalized_translate(dyadic_band_coordinate(q), ell)


def active_dyadic_indices(q: float) -> tuple[int, int]:
    """The at-most-two bands that can carry nonzero dyadic cutoff at ``q``."""
    return local_integer_indices(dyadic_band_coordinate(q))


@dataclass(frozen=True)
class ProductSlowCutoff:
    """The product squared-partition factor ``chi_{ell,a}(R,Z,T)``.

    We choose the grid origin at zero.  Section 6.2 allows an arbitrary fixed
    translated mesh; the zero-origin choice is therefore implementation data,
    not a claim that the manuscript selected this particular translate.
    """

    chart: DyadicChart
    a: tuple[int, int, int]

    def __post_init__(self) -> None:
        if not isinstance(self.a, tuple) or len(self.a) != 3:
            raise ValueError("a must be a length-three integer tuple")
        if any(isinstance(v, bool) or not isinstance(v, int) for v in self.a):
            raise ValueError("a must be a length-three integer tuple")

    @property
    def mesh_size(self) -> float:
        return self.chart.S_star ** -3

    @property
    def center(self) -> tuple[float, float, float]:
        h = self.mesh_size
        return tuple(h * v for v in self.a)

    def value(self, R: float, Z: float, T: float) -> float:
        h = self.mesh_size
        point = tuple(_finite(v, "slow coordinate") for v in (R, Z, T))
        out = 1.0
        for x, idx in zip(point, self.a):
            out *= normalized_translate(x / h, idx)
        return out

    def support_contains(self, R: float, Z: float, T: float) -> bool:
        point = tuple(_finite(v, "slow coordinate") for v in (R, Z, T))
        return all(abs(x - c) <= self.mesh_size for x, c in zip(point, self.center))

    def enclosure(self, label: SlowLabel) -> SlowBoxEnclosure:
        """Bridge this actual product cutoff to the existing slow-label interface."""
        if label.ell != self.chart.ell or label.a != self.a:
            raise ValueError("label must identify this band and product-grid index")
        if self.center[0] - self.mesh_size <= 0.0:
            raise ValueError("paper active slow boxes must stay away from R=0")
        return SlowBoxEnclosure(
            chart=self.chart,
            label=label,
            center=self.center,
            half_width=(self.mesh_size,) * 3,
        )


def active_product_indices(chart: DyadicChart, R: float, Z: float, T: float):
    """Yield the at-most-eight product-grid indices nonzero at a slow point."""
    h = chart.S_star ** -3
    coords = tuple(_finite(v, "slow coordinate") / h for v in (R, Z, T))
    axes = [local_integer_indices(x) for x in coords]
    return tuple(itertools.product(*axes))


def product_partition_square_sum(chart: DyadicChart, R: float, Z: float, T: float) -> float:
    """Finite local evaluation of ``sum_a chi_{ell,a}^2``."""
    total = 0.0
    for a in active_product_indices(chart, R, Z, T):
        value = ProductSlowCutoff(chart, a).value(R, Z, T)
        total += value * value
    return total


def slow_label_cutoff(label: SlowLabel, q: float, R: float, Z: float, T: float) -> float:
    """Equation (6.9) cutoff ``eta_gamma=chi_ell(q) chi_{ell,a}``.

    The sign is intentionally absent from the value: Eq. (6.9) duplicates the
    same cutoff for the ``+`` and ``-`` wave families.
    """
    chart = DyadicChart(label.ell)
    return dyadic_cutoff(label.ell, q) * ProductSlowCutoff(chart, label.a).value(R, Z, T)
