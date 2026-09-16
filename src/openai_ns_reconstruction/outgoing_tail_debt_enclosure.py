"""Certified rational enclosure of the pinned outgoing ``tailDebt`` scalar.

For the actual :class:`~openai_ns_reconstruction.outgoing_tail.TailData`, write
``h = Fraction.from_float(data.h)`` and ``k = 1 - h``. Integration by parts
in the pinned tail transition gives

    tailDebt = rho / (1 - rho) * (exp(3 k) - 2 k J),
    J = integral_0^1 exp(k * (1 + 2 x)) * sigma(x) dx,

where ``rho = h * exp(-5) / 528``. The scalar is enclosed with exact
fractions, using the validated exponential and outgoing-sigma primitives.
This does not certify the complete pressure integral or any global norm.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from fractions import Fraction

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
    interval_divide_positive,
    interval_multiply,
    positive_cap as _positive_cap,
    positive_fraction as _positive_fraction,
)


_DEFAULT_TOLERANCE = Fraction(1, 10**8)
_DEFAULT_MAX_CELLS = 4096
_DEFAULT_MAX_TERMS = 1024
_DEFAULT_MAX_SQUARINGS = 4096
_RHO_DENOMINATOR = 528


def _require_tail_data(data: object) -> tuple[TailData, Fraction]:
    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    h = Fraction.from_float(data.h)
    if not Fraction(0) < h < Fraction(1, 20):
        raise ValueError("TailData h must satisfy 0 < h < 1/20")
    return data, h


@dataclass(frozen=True)
class TailDebtEnclosure:
    """Exact rational bounds for the selected finite ``tailDebt`` scalar."""

    h: Fraction
    rho: RationalInterval
    debt: RationalInterval
    cells: int

    def __post_init__(self) -> None:
        if not isinstance(self.h, Fraction):
            raise TypeError("h must be a Fraction")
        if not Fraction(0) < self.h < Fraction(1, 20):
            raise ValueError("h must satisfy 0 < h < 1/20")
        if not isinstance(self.rho, RationalInterval):
            raise TypeError("rho must be a RationalInterval")
        if not Fraction(0) <= self.rho.lower <= self.rho.upper < Fraction(1):
            raise ValueError("rho interval must lie in [0, 1)")
        if not isinstance(self.debt, RationalInterval):
            raise TypeError("debt must be a RationalInterval")
        if self.debt.lower <= 0:
            raise ValueError("tail debt interval must be positive")
        if isinstance(self.cells, bool) or not isinstance(self.cells, int) or self.cells <= 0:
            raise ValueError("cells must be a positive integer")

    @property
    def tail_debt(self) -> RationalInterval:
        return self.debt

    @property
    def width(self) -> Fraction:
        return self.debt.width

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False


def _rho_interval(
    h: Fraction,
    tolerance: Fraction,
    max_terms: int,
    max_squarings: int,
) -> RationalInterval:
    exp_minus_five = validated_exp_negative(
        Fraction(5),
        absolute_tolerance=tolerance,
        max_terms=max_terms,
        max_squarings=max_squarings,
    )
    rho = RationalInterval(
        h * exp_minus_five.lower / _RHO_DENOMINATOR,
        h * exp_minus_five.upper / _RHO_DENOMINATOR,
    )
    if rho.upper >= 1:
        raise ArithmeticError("rho interval is not contained in [0, 1)")
    return rho


def _tail_ratio(rho: RationalInterval) -> RationalInterval:
    denominator = RationalInterval(1 - rho.upper, 1 - rho.lower)
    return interval_divide_positive(rho, denominator)


def _integral_j_levels(
    k: Fraction,
    *,
    node_tolerance: Fraction,
    max_cells: int,
    max_terms: int,
    max_squarings: int,
) -> Iterator[tuple[RationalInterval, int]]:
    """Yield successively refined Darboux bounds for J.

    The integral itself is deliberately allowed to be wider than the final
    debt target: its coefficient is the small factor ``2 * k * rho / (1-rho)``.
    The caller checks the complete debt interval at every yielded level.
    """

    endpoint_tolerance = _positive_fraction(node_tolerance, "node_tolerance")
    product_step = _dyadic_step(endpoint_tolerance / 4)
    sigma_cache: dict[Fraction, RationalInterval] = {}
    exponential_cache: dict[Fraction, RationalInterval] = {}
    product_cache: dict[Fraction, RationalInterval] = {}

    def product(node: Fraction) -> RationalInterval:
        if node in product_cache:
            return product_cache[node]
        if node not in sigma_cache:
            sigma_cache[node] = validated_outgoing_sigma(
                node,
                absolute_tolerance=endpoint_tolerance,
                max_terms=max_terms,
                max_squarings=max_squarings,
            )
        if node not in exponential_cache:
            z = k * (1 + 2 * node)
            exponential_cache[node] = validated_exp_positive_on_0_3(
                z,
                absolute_tolerance=endpoint_tolerance,
                max_terms=max_terms,
                max_squarings=max_squarings,
            )
        raw = interval_multiply(sigma_cache[node], exponential_cache[node])
        rounded = RationalInterval(
            _floor_grid(raw.lower, product_step),
            _ceil_grid(raw.upper, product_step),
        )
        product_cache[node] = rounded
        return rounded

    cells = 1
    while cells <= max_cells:
        width = Fraction(1, cells)
        lower = Fraction(0)
        upper = Fraction(0)
        for index in range(cells):
            left = Fraction(index, cells)
            right = Fraction(index + 1, cells)
            lower += width * product(left).lower
            upper += width * product(right).upper
        result = RationalInterval(lower, upper)
        yield result, cells
        cells *= 2


def validated_tail_debt_enclosure(
    data: TailData,
    *,
    absolute_tolerance: Fraction = _DEFAULT_TOLERANCE,
    max_cells: int = _DEFAULT_MAX_CELLS,
    max_terms: int = _DEFAULT_MAX_TERMS,
    max_squarings: int = _DEFAULT_MAX_SQUARINGS,
) -> TailDebtEnclosure:
    """Return a certified rational interval for the actual schedule tail debt."""

    data, h = _require_tail_data(data)
    tolerance = _positive_fraction(absolute_tolerance, "absolute_tolerance")
    max_cells = _positive_cap(max_cells, "max_cells")
    max_terms = _positive_cap(max_terms, "max_terms")
    max_squarings = _positive_cap(max_squarings, "max_squarings")
    epsilon = min(tolerance, Fraction(1))
    k = 1 - h

    component_tolerance = min(Fraction(1, 1024), epsilon / 1024)
    rho = _rho_interval(
        h,
        component_tolerance,
        max_terms,
        max_squarings,
    )
    ratio = _tail_ratio(rho)
    exp_3k = validated_exp_positive_on_0_3(
        3 * k,
        absolute_tolerance=component_tolerance,
        max_terms=max_terms,
        max_squarings=max_squarings,
    )
    for j_interval, cells in _integral_j_levels(
        k,
        node_tolerance=component_tolerance,
        max_cells=max_cells,
        max_terms=max_terms,
        max_squarings=max_squarings,
    ):
        factor = RationalInterval(
            exp_3k.lower - 2 * k * j_interval.upper,
            exp_3k.upper - 2 * k * j_interval.lower,
        )
        factor = RationalInterval(
            max(Fraction(1), factor.lower),
            min(Fraction(27), factor.upper),
        )
        if factor.lower <= 0 or factor.lower > factor.upper:
            raise ArithmeticError("tail debt integration-by-parts factor lost positivity")
        debt = interval_multiply(ratio, factor)
        if debt.width <= tolerance:
            return TailDebtEnclosure(h=h, rho=rho, debt=debt, cells=cells)

    raise ArithmeticError("tail debt enclosure did not reach requested tolerance before max_cells")


tail_debt_enclosure = validated_tail_debt_enclosure
validated_tail_debt = validated_tail_debt_enclosure


__all__ = [
    "TailDebtEnclosure",
    "tail_debt_enclosure",
    "validated_tail_debt",
    "validated_tail_debt_enclosure",
]
