"""Validated exact-rational integration of the natural-axis phase kernel.

For exact rational ``h``, ``j``, and positive ``sigma``, this module evaluates
the oriented integral of

    g(x) = -(1 - 2 h x^2) H(h, j, x) / (H(h, j, x)^2 + sigma^2),
    H(h, j, x) = (1/2 - h) x + (1 - x^2) (4 x + j),

using finite Taylor expansions on rational cells.  Each accepted cell carries
an exact rational estimate and Cauchy tail bound.  The adaptive routine only
bisects when the denominator gate fails and increases Taylor order when the
cell error exceeds its exact rational tolerance budget.  Resource limits fail
closed instead of returning an uncertified result.

This validates the chosen rational kernel integral only.  It does not evaluate
the downstream quadrature-based phase, select theorem parameters, or certify
the amplitude, pressure, global norm, or reconstructed velocity.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from fractions import Fraction
import math


_WINDOW_LEFT = Fraction(-11, 10)
_WINDOW_RIGHT = Fraction(11, 10)
_HALF = Fraction(1, 2)


class NeedSubdivision(ArithmeticError):
    """The exact denominator gate failed on a cell, so it must be bisected."""


class PhaseIntegrationLimit(ArithmeticError):
    """A finite adaptive integration resource limit was reached."""


def _fraction(value: Fraction, name: str) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be a Fraction")
    return value


def _nonnegative_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _validate_accumulation(value: str) -> str:
    if value not in ("exact", "dyadic"):
        raise ValueError("accumulation must be exactly 'exact' or 'dyadic'")
    return value


def _largest_dyadic_at_most(value: Fraction) -> Fraction:
    """Return the largest positive power of two no greater than ``value``."""

    if value <= 0:
        raise ValueError("dyadic step bound must be positive")
    exponent = value.numerator.bit_length() - value.denominator.bit_length()
    if exponent >= 0:
        step = Fraction(1 << exponent)
    else:
        step = Fraction(1, 1 << (-exponent))
    while step > value:
        step /= 2
    while step * 2 <= value:
        step *= 2
    return step


def _floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def _ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def _dyadic_cell_enclosure(
    result: "PhaseCellResult",
    budget: Fraction,
) -> tuple[Fraction, Fraction]:
    """Round one cell enclosure outward and return its dyadic midpoint/radius."""

    step = _largest_dyadic_at_most(budget / 4)
    lower = result.integral_estimate - result.error_bound
    upper = result.integral_estimate + result.error_bound
    lower_units = _floor_fraction(lower / step)
    upper_units = _ceil_fraction(upper / step)
    rounded_lower = step * lower_units
    rounded_upper = step * upper_units
    return (
        (rounded_lower + rounded_upper) / 2,
        (rounded_upper - rounded_lower) / 2,
    )


def _window_endpoint(value: Fraction, name: str) -> Fraction:
    value = _fraction(value, name)
    if not _WINDOW_LEFT <= value <= _WINDOW_RIGHT:
        raise ValueError(f"{name} must lie in the pinned window [-11/10,11/10]")
    return value


def _validate_cell_inputs(
    h: Fraction,
    j: Fraction,
    sigma: Fraction,
    left: Fraction,
    right: Fraction,
    order: int,
) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction, int]:
    h = _fraction(h, "h")
    j = _fraction(j, "j")
    sigma = _fraction(sigma, "sigma")
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    left = _window_endpoint(left, "left")
    right = _window_endpoint(right, "right")
    if left > right:
        raise ValueError("cell endpoints must be ordered")
    order = _nonnegative_integer(order, "order")
    return h, j, sigma, left, right, order


def _poly_add(
    left: tuple[Fraction, ...],
    right: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    length = max(len(left), len(right))
    values = [Fraction(0)] * length
    for index, value in enumerate(left):
        values[index] += value
    for index, value in enumerate(right):
        values[index] += value
    return tuple(values)


def _poly_mul(
    left: tuple[Fraction, ...],
    right: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    values = [Fraction(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            values[left_index + right_index] += left_value * right_value
    return tuple(values)


def _poly_scale(
    values: tuple[Fraction, ...], scalar: Fraction,
) -> tuple[Fraction, ...]:
    return tuple(scalar * value for value in values)


def _poly_shift(
    values: tuple[Fraction, ...], center: Fraction,
) -> tuple[Fraction, ...]:
    """Return coefficients of ``values(center + z)`` in powers of ``z``."""

    shifted = [Fraction(0)] * len(values)
    for source_degree, source_value in enumerate(values):
        for degree in range(source_degree + 1):
            shifted[degree] += (
                source_value
                * math.comb(source_degree, degree)
                * center ** (source_degree - degree)
            )
    return tuple(shifted)


def _centered_kernel_polynomials(
    h: Fraction,
    j: Fraction,
    sigma: Fraction,
    center: Fraction,
) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...]]:
    D = _HALF - h
    L = (Fraction(1), Fraction(0), -2 * h)
    H = (j, Fraction(4) + D, -j, Fraction(-4))
    P = _poly_scale(_poly_mul(L, H), Fraction(-1))
    Q = list(_poly_mul(H, H))
    Q[0] += sigma * sigma
    return _poly_shift(P, center), _poly_shift(tuple(Q), center)


def _cell_result_from_polynomials(
    P: tuple[Fraction, ...],
    Q: tuple[Fraction, ...],
    radius: Fraction,
    order: int,
) -> "PhaseCellResult":
    rho = 2 * radius
    q0 = Q[0]
    if q0 <= 0:
        raise ArithmeticError("centered denominator constant must be positive")

    theta = sum(
        (abs(Q[k]) * rho**k / q0 for k in range(1, len(Q))),
        Fraction(0),
    )
    if theta >= _HALF:
        raise NeedSubdivision("cell denominator gate theta < 1/2 failed")

    numerator_majorant = sum(
        (abs(P[k]) * rho**k for k in range(len(P))),
        Fraction(0),
    )
    M = 2 * numerator_majorant / q0

    coefficients: list[Fraction] = [P[0] / q0]
    for degree in range(1, order + 1):
        p_degree = P[degree] if degree < len(P) else Fraction(0)
        q_sum = sum(
            (
                Q[index] * coefficients[degree - index]
                for index in range(1, min(degree, len(Q) - 1) + 1)
            ),
            Fraction(0),
        )
        coefficients.append((p_degree - q_sum) / q0)

    estimate = Fraction(0)
    for degree, coefficient in enumerate(coefficients):
        if degree % 2 == 0:
            estimate += coefficient * (2 * radius ** (degree + 1) / (degree + 1))

    # E_N = 2*r*M*2^(-N-1)/(1/2) = 2*r*M*2^(-N).
    error = 2 * radius * M / (2**order)
    return PhaseCellResult(
        integral_estimate=estimate,
        error_bound=error,
        theta=theta,
        order=order,
    )


@dataclass(frozen=True)
class PhaseCellResult:
    """Exact rational estimate and enclosure for one accepted phase cell."""

    integral_estimate: Fraction
    error_bound: Fraction
    lower: Fraction = field(init=False)
    upper: Fraction = field(init=False)
    theta: Fraction
    order: int

    def __post_init__(self) -> None:
        for name in ("integral_estimate", "error_bound", "theta"):
            _fraction(getattr(self, name), name)
        if self.error_bound < 0:
            raise ValueError("error must be nonnegative")
        if self.theta < 0 or self.theta >= _HALF:
            raise ValueError("theta must lie in [0,1/2)")
        _nonnegative_integer(self.order, "order")
        object.__setattr__(self, "lower", self.integral_estimate - self.error_bound)
        object.__setattr__(self, "upper", self.integral_estimate + self.error_bound)

    @property
    def estimate(self) -> Fraction:
        """Compatibility alias for the exact cell integral estimate."""

        return self.integral_estimate

    @property
    def error(self) -> Fraction:
        """Compatibility alias for the exact cell error bound."""

        return self.error_bound


def validated_phase_cell(
    h: Fraction,
    j: Fraction,
    sigma: Fraction,
    left: Fraction,
    right: Fraction,
    order: int = 16,
) -> PhaseCellResult:
    """Validate and integrate one ordered rational cell exactly.

    The cell is accepted only when the centered denominator majorant satisfies
    ``theta < 1/2``.  ``NeedSubdivision`` requests bisection by the adaptive
    caller; it is not a numerical failure of the rational kernel itself.
    """

    h, j, sigma, left, right, order = _validate_cell_inputs(
        h, j, sigma, left, right, order
    )
    if left == right:
        return PhaseCellResult(
            integral_estimate=Fraction(0),
            error_bound=Fraction(0),
            theta=Fraction(0),
            order=order,
        )

    center = (left + right) / 2
    radius = (right - left) / 2
    P, Q = _centered_kernel_polynomials(h, j, sigma, center)
    return _cell_result_from_polynomials(P, Q, radius, order)


@dataclass(frozen=True)
class PhaseIntegralResult:
    """Exact rational adaptive phase integral and its conditional enclosure."""

    integral_estimate: Fraction
    error_bound: Fraction
    lower: Fraction = field(init=False)
    upper: Fraction = field(init=False)
    cell_count: int
    max_order_used: int

    def __post_init__(self) -> None:
        for name in ("integral_estimate", "error_bound"):
            _fraction(getattr(self, name), name)
        if self.error_bound < 0:
            raise ValueError("error_bound must be nonnegative")
        _nonnegative_integer(self.cell_count, "cell_count")
        _nonnegative_integer(self.max_order_used, "max_order_used")
        object.__setattr__(self, "lower", self.integral_estimate - self.error_bound)
        object.__setattr__(self, "upper", self.integral_estimate + self.error_bound)

    @property
    def paper_exact(self) -> bool:
        return False


def _order_for_budget(
    error: Fraction,
    budget: Fraction,
    initial_order: int,
    max_order: int,
) -> int:
    if error <= budget:
        return initial_order
    order = initial_order
    halved = error
    while halved > budget:
        order += 1
        if order > max_order:
            raise PhaseIntegrationLimit(
                "maximum Taylor order reached before the cell tolerance was met"
            )
        halved /= 2
    return order


def validated_phase_integral(
    h: Fraction,
    j: Fraction,
    sigma: Fraction,
    eta: Fraction,
    *,
    absolute_tolerance: Fraction,
    initial_order: int = 16,
    max_order: int = 4096,
    max_cells: int = 4096,
    max_depth: int = 128,
    accumulation: str = "exact",
) -> PhaseIntegralResult:
    """Integrate the validated rational phase kernel from zero to ``eta``.

    Internal cells are always ordered from the smaller endpoint to the larger
    endpoint; a negative ``eta`` negates the final interval.  The tolerance
    budget of a cell is ``absolute_tolerance * length / abs(eta)``, so the
    budgets of a complete partition sum exactly to the requested tolerance.
    ``max_cells`` counts every distinct geometric cell attempted, including a
    rejected parent whose split produced accepted children.  Order retries on
    the same geometry do not consume another cell count.  Any cap or failed
    order gate raises :class:`PhaseIntegrationLimit`.

    ``accumulation="exact"`` retains the exact Taylor midpoint and Cauchy
    radius sums.  ``accumulation="dyadic"`` first requires each cell Cauchy
    error to be at most half its budget, then rounds that cell's exact
    enclosure outward to a power-of-two grid with step at most one quarter of
    the budget.  The returned midpoint and radius therefore remain an exact
    rational enclosure while their accumulated denominators stay dyadic.
    """

    h = _fraction(h, "h")
    j = _fraction(j, "j")
    sigma = _fraction(sigma, "sigma")
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    eta = _window_endpoint(eta, "eta")
    absolute_tolerance = _fraction(absolute_tolerance, "absolute_tolerance")
    if absolute_tolerance <= 0:
        raise ValueError("absolute_tolerance must be positive")
    initial_order = _nonnegative_integer(initial_order, "initial_order")
    max_order = _positive_integer(max_order, "max_order")
    max_cells = _positive_integer(max_cells, "max_cells")
    max_depth = _positive_integer(max_depth, "max_depth")
    if initial_order > max_order:
        raise PhaseIntegrationLimit("initial_order exceeds max_order")
    accumulation = _validate_accumulation(accumulation)

    if eta == 0:
        return PhaseIntegralResult(
            integral_estimate=Fraction(0),
            error_bound=Fraction(0),
            cell_count=0,
            max_order_used=0,
        )

    orientation = 1 if eta > 0 else -1
    left, right = (Fraction(0), eta) if orientation > 0 else (eta, Fraction(0))
    total_length = right - left
    pending: deque[tuple[Fraction, Fraction, int]] = deque([(left, right, 0)])
    attempted_cells = 0
    estimate = Fraction(0)
    error_bound = Fraction(0)
    max_order_used = 0
    accepted_cells = 0

    while pending:
        cell_left, cell_right, depth = pending.popleft()
        attempted_cells += 1
        if attempted_cells > max_cells:
            raise PhaseIntegrationLimit("maximum attempted cell count reached")

        try:
            result = validated_phase_cell(
                h,
                j,
                sigma,
                cell_left,
                cell_right,
                initial_order,
            )
        except NeedSubdivision:
            if depth >= max_depth:
                raise PhaseIntegrationLimit(
                    "maximum subdivision depth reached before the denominator gate passed"
                ) from None
            if attempted_cells + len(pending) + 2 > max_cells:
                raise PhaseIntegrationLimit(
                    "maximum attempted cell count would be exceeded by subdivision"
                ) from None
            midpoint = (cell_left + cell_right) / 2
            pending.append((cell_left, midpoint, depth + 1))
            pending.append((midpoint, cell_right, depth + 1))
            continue

        budget = absolute_tolerance * (cell_right - cell_left) / total_length
        order_budget = budget if accumulation == "exact" else budget / 2
        target_order = _order_for_budget(
            result.error_bound,
            order_budget,
            initial_order,
            max_order,
        )
        if target_order != initial_order:
            result = validated_phase_cell(
                h,
                j,
                sigma,
                cell_left,
                cell_right,
                target_order,
            )
        if accumulation == "exact":
            estimate += result.integral_estimate
            error_bound += result.error_bound
        else:
            midpoint, radius = _dyadic_cell_enclosure(result, budget)
            estimate += midpoint
            error_bound += radius
        max_order_used = max(max_order_used, result.order)
        accepted_cells += 1

    if orientation < 0:
        estimate = -estimate
    return PhaseIntegralResult(
        integral_estimate=estimate,
        error_bound=error_bound,
        cell_count=accepted_cells,
        max_order_used=max_order_used,
    )


__all__ = [
    "NeedSubdivision",
    "PhaseCellResult",
    "PhaseIntegrationLimit",
    "PhaseIntegralResult",
    "validated_phase_cell",
    "validated_phase_integral",
]
