"""Certified rational enclosure of the pinned outgoing ``tailDebt`` scalar.

For the actual :class:`~openai_ns_reconstruction.outgoing_tail.TailData`, write
``h = Fraction.from_float(data.h)`` and ``k = 1 - h``.  Integration by parts
in the pinned tail transition gives

    tailDebt = rho / (1 - rho) * (exp(3 k) - 2 k J),
    J = integral_0^1 exp(k * (1 + 2 x)) * sigma(x) dx,

where ``rho = h * exp(-5) / 528``.  The scalar is enclosed with exact
fractions, using the validated exponential and outgoing-sigma primitives.
This does not certify the complete pressure integral or any global norm.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from fractions import Fraction

from .outgoing_sigma_enclosure import (
    RationalInterval,
    validated_exp_negative,
    validated_outgoing_sigma,
)
from .outgoing_tail import TailData


_DEFAULT_TOLERANCE = Fraction(1, 10**8)
_DEFAULT_MAX_CELLS = 4096
_DEFAULT_MAX_TERMS = 1024
_DEFAULT_MAX_SQUARINGS = 4096
_RHO_DENOMINATOR = 528
_EXP_SAFE_UPPER = Fraction(1)
_EXP_SAFE_LOWER = Fraction(1, 27)


def _positive_fraction(value: object, name: str) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be a Fraction")
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _positive_cap(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _dyadic_step(bound: Fraction) -> Fraction:
    if not isinstance(bound, Fraction) or bound <= 0:
        raise ValueError("dyadic step bound must be positive")
    exponent = max(0, bound.denominator.bit_length() - bound.numerator.bit_length())
    step = Fraction(1, 1 << exponent)
    while step > bound:
        exponent += 1
        step = Fraction(1, 1 << exponent)
    while exponent > 0 and Fraction(1, 1 << (exponent - 1)) <= bound:
        exponent -= 1
        step = Fraction(1, 1 << exponent)
    return step


def _floor_grid(value: Fraction, step: Fraction) -> Fraction:
    return (value // step) * step


def _ceil_grid(value: Fraction, step: Fraction) -> Fraction:
    quotient = value / step
    return (-((-quotient.numerator) // quotient.denominator)) * step


def _require_tail_data(data: object) -> tuple[TailData, Fraction]:
    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    h = Fraction.from_float(data.h)
    if not Fraction(0) < h < Fraction(1, 20):
        raise ValueError("TailData h must satisfy 0 < h < 1/20")
    return data, h


def _exp_positive_on_0_3(
    z: Fraction,
    absolute_tolerance: Fraction,
    max_terms: int,
    max_squarings: int,
) -> RationalInterval:
    """Enclose ``exp(z)`` for ``0 <= z <= 3`` using ``exp(-z)``."""

    if not isinstance(z, Fraction):
        raise TypeError("z must be a Fraction")
    if not Fraction(0) <= z <= Fraction(3):
        raise ValueError("z must lie in [0, 3]")
    tolerance = _positive_fraction(absolute_tolerance, "absolute_tolerance")
    max_terms = _positive_cap(max_terms, "max_terms")
    max_squarings = _positive_cap(max_squarings, "max_squarings")

    negative = validated_exp_negative(
        z,
        absolute_tolerance=tolerance / 729,
        max_terms=max_terms,
        max_squarings=max_squarings,
    )
    # The exact value exp(-z) lies in [1/27, 1].  Intersecting the input
    # interval with this independently proved range prevents a rounded lower
    # endpoint from causing a false loss of positivity on inversion.
    inverse_lower = max(_EXP_SAFE_LOWER, negative.lower)
    inverse_upper = min(_EXP_SAFE_UPPER, negative.upper)
    if inverse_lower <= 0 or inverse_lower > inverse_upper:
        raise ArithmeticError("positive exponential input interval is invalid")
    result = RationalInterval(
        Fraction(1, 1) / inverse_upper,
        Fraction(1, 1) / inverse_lower,
    )
    if result.width > tolerance:
        raise ArithmeticError("positive exponential enclosure exceeded requested tolerance")
    return result


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
    return RationalInterval(
        rho.lower / (1 - rho.lower),
        rho.upper / (1 - rho.upper),
    )


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
            exponential_cache[node] = _exp_positive_on_0_3(
                z,
                endpoint_tolerance,
                max_terms,
                max_squarings,
            )
        sigma = sigma_cache[node]
        exponential = exponential_cache[node]
        raw = RationalInterval(
            sigma.lower * exponential.lower,
            sigma.upper * exponential.upper,
        )
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
    return


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

    # Keep the local function enclosures small, then let the final debt width
    # determine when the spatial Darboux refinement can stop.
    component_tolerance = min(Fraction(1, 1024), epsilon / 1024)
    rho = _rho_interval(
        h,
        component_tolerance,
        max_terms,
        max_squarings,
    )
    ratio = _tail_ratio(rho)
    exp_3k = _exp_positive_on_0_3(
        3 * k,
        component_tolerance,
        max_terms,
        max_squarings,
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
        # Integration by parts identifies this factor with
        # integral exp(k(1+2x)) d sigma, whose total sigma mass is one.
        # Since k > 0, it lies in [exp(k), exp(3k)] subset [1,27].
        factor = RationalInterval(
            max(Fraction(1), factor.lower),
            min(Fraction(27), factor.upper),
        )
        if factor.lower <= 0 or factor.lower > factor.upper:
            raise ArithmeticError("tail debt integration-by-parts factor lost positivity")
        debt = RationalInterval(
            ratio.lower * factor.lower,
            ratio.upper * factor.upper,
        )
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
