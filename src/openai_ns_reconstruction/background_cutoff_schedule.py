"""Paper-admissible recursive cutoff-scale witnesses for Section 5.

The pinned Lean construction ``SlowBorelBase.exists_admissibleScales`` first
bounds the normalized template jets by constants ``C[j,m]`` and then invokes
``DiagonalScale.exists_diagonal_scales`` with logarithmic power ``p=0`` and
positive gain ``g(j)=2*h*j``. Thus, at every positive order ``j``, it is
enough to choose a local integer scale ``b_j`` such that

    C[j,m] * (1 / b_j) ** (h*j) <= 2**(-j)    for m <= j+2,

and finally take the recursive doubling envelope

    a_0     = max(1, B),
    a_{j+1} = max(b_{j+1}, 2*a_j).

This module materializes that numerical selection for a requested finite
prefix once the *actual* normalized jet bounds C[j,m] are supplied. It does
not manufacture those bounds from sampled profiles, and a finite prefix is
not the paper's completed infinite all-order background.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Callable, Sequence
import math
import sys


_LOG2 = math.log(2.0)
_LOG_FLOAT_MAX = math.log(sys.float_info.max)
_LOG_FLOAT_MIN_SUBNORMAL = math.log(math.ulp(0.0))
# Guard against accidentally allocating an arbitrarily large Python integer
# when an upstream analytic bound or h is malformed.  This is an implementation
# resource bound only; it is not a theorem-side restriction on admissible scales.
_MAX_EXACT_SCALE_BITS = 1_000_000


def _positive_float(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _nonnegative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def local_scale_from_template_bounds(
    h: float,
    order: int,
    jet_bounds: Sequence[float],
) -> int:
    """Choose one explicit local scale for SlowBorelBase at positive order j.

    ``jet_bounds[m]`` is the normalized-template constant ``C[j,m]`` for
    ``0 <= m <= j+2``. Since the pinned SlowBorelBase specialization has
    logarithmic exponent p=0, the worst q in ``0 < q <= 1/b_j`` is the edge
    q=1/b_j. The returned integer therefore certifies the entire punctured
    interval by monotonicity of q**(h*j).

    When the real threshold fits in binary64, the routine keeps the previous
    near-minimal ``ceil(exp(log_required))`` selection.  If that threshold is
    larger than binary64, it switches to an arbitrary-precision Python integer
    power-of-two witness and verifies the same inequality in log space.  This
    avoids turning a representation limit into a false theorem obstruction.
    Non-finite logarithmic thresholds, or thresholds exceeding the explicit
    implementation bit budget, fail closed instead of clipping the scale.
    """
    h = _positive_float(h, "h")
    order = _nonnegative_int(order, "order")
    if order == 0:
        raise ValueError("order zero is exempt from the positive-gain scale condition")

    bounds = tuple(_positive_float(v, f"jet_bounds[{m}]") for m, v in enumerate(jet_bounds))
    if len(bounds) != order + 3:
        raise ValueError("positive order j requires exactly j+3 bounds for m=0,...,j+2")

    cmax = max(bounds)
    exponent = h * order
    log_required = (order * _LOG2 + math.log(cmax)) / exponent
    if not math.isfinite(log_required):
        raise OverflowError("required cutoff scale has a non-finite logarithmic threshold")

    target_log = -order * _LOG2
    if log_required <= 0.0:
        candidate = 1
    elif log_required < _LOG_FLOAT_MAX:
        candidate = max(1, math.ceil(math.exp(log_required)))
        # A one-integer correction keeps the representable-range behavior
        # near-minimal despite transcendental rounding.
        while math.log(cmax) - exponent * math.log(candidate) > target_log:
            candidate += 1
        return candidate
    else:
        bits = max(1, math.ceil(log_required / _LOG2))
        if bits > _MAX_EXACT_SCALE_BITS:
            raise OverflowError(
                "required cutoff scale exceeds the configured exact-integer bit budget"
            )
        candidate = 1 << bits
        # Rounding in log_required/log(2) can undershoot by one bit.  Doubling
        # is theorem-safe and avoids ever converting the wide integer to float.
        while math.log(cmax) - exponent * math.log(candidate) > target_log:
            bits += 1
            if bits > _MAX_EXACT_SCALE_BITS:
                raise OverflowError(
                    "required cutoff scale exceeds the configured exact-integer bit budget"
                )
            candidate <<= 1
        return candidate

    while math.log(cmax) - exponent * math.log(candidate) > target_log:
        candidate += 1
    return candidate


@dataclass(frozen=True)
class SlowBorelCutoffSchedule:
    """A finite-prefix witness for the recursive Section-5 cutoff schedule."""

    h: float
    initial_lower_bound: int
    jet_bounds: tuple[tuple[float, ...], ...]
    local_scales: tuple[int, ...]
    scales: tuple[int, ...]

    def __post_init__(self) -> None:
        h = _positive_float(self.h, "h")
        lower = _nonnegative_int(self.initial_lower_bound, "initial_lower_bound")
        rows = tuple(tuple(_positive_float(v, "jet bound") for v in row) for row in self.jet_bounds)
        local = tuple(_nonnegative_int(v, "local scale") for v in self.local_scales)
        scales = tuple(_nonnegative_int(v, "scale") for v in self.scales)
        if len(local) != len(rows) + 1 or len(scales) != len(rows) + 1:
            raise ValueError("schedule must contain stage zero plus one entry per positive-order bound row")
        if local[0] != max(1, lower) or scales[0] != local[0]:
            raise ValueError("stage-zero scale must be max(1, initial_lower_bound)")
        for j, row in enumerate(rows, start=1):
            if len(row) != j + 3:
                raise ValueError("row j must contain exactly j+3 normalized template bounds")
            if local[j] < 1 or scales[j] < local[j]:
                raise ValueError("positive-order scales must dominate their local witnesses")
            if scales[j] < 2 * scales[j - 1]:
                raise ValueError("recursive cutoff scales must at least double")
            target_log = -j * _LOG2
            exponent = h * j
            for c in row:
                if math.log(c) - exponent * math.log(scales[j]) > target_log + 32 * sys.float_info.epsilon:
                    raise ValueError("scale does not certify the SlowBorelBase dyadic edge bound")
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "initial_lower_bound", lower)
        object.__setattr__(self, "jet_bounds", rows)
        object.__setattr__(self, "local_scales", local)
        object.__setattr__(self, "scales", scales)

    @property
    def max_order(self) -> int:
        return len(self.scales) - 1

    def reciprocal_support_log_edge(self, order: int) -> float:
        """Return ``log(1/a_j)`` without narrowing the integer scale to binary64."""
        order = _nonnegative_int(order, "order")
        if order > self.max_order:
            raise IndexError("order lies outside the constructed finite prefix")
        return -math.log(self.scales[order])

    def reciprocal_support_edge(self, order: int) -> float:
        """Return 1/a_j when that positive edge is representable in binary64.

        For a mathematically valid wide integer scale, ``1/a_j`` can be smaller
        than the least positive binary64 subnormal.  Returning 0.0 would turn a
        representation artifact into an exact support claim, so this accessor
        fails closed and directs callers to ``reciprocal_support_log_edge``.
        """
        log_edge = self.reciprocal_support_log_edge(order)
        if log_edge < _LOG_FLOAT_MIN_SUBNORMAL:
            raise OverflowError(
                "reciprocal support edge underflows binary64; use reciprocal_support_log_edge"
            )
        value = math.exp(log_edge)
        if value == 0.0:
            raise OverflowError(
                "reciprocal support edge underflows binary64; use reciprocal_support_log_edge"
            )
        return value

    def edge_log_margin(self, order: int, derivative_order: int) -> float:
        """Logarithmic margin in C[j,m] a_j^(-h*j) <= 2^(-j)."""
        order = _nonnegative_int(order, "order")
        derivative_order = _nonnegative_int(derivative_order, "derivative_order")
        if order == 0 or order > self.max_order:
            raise ValueError("margin is defined only for constructed positive orders")
        row = self.jet_bounds[order - 1]
        if derivative_order >= len(row):
            raise ValueError("derivative_order must satisfy m <= j+2")
        c = row[derivative_order]
        log_bound = math.log(c) - (self.h * order) * math.log(self.scales[order])
        return (-order * _LOG2) - log_bound

    def majorant_at(self, order: int, derivative_order: int, q: float) -> tuple[float, float]:
        """Evaluate the certified p=0 scalar majorant and its dyadic target.

        This checks only the numerical majorant implied by supplied C[j,m]; it
        does not evaluate or certify an actual coefficient derivative.
        """
        q = _positive_float(q, "q")
        order = _nonnegative_int(order, "order")
        derivative_order = _nonnegative_int(derivative_order, "derivative_order")
        if order == 0 or order > self.max_order:
            raise ValueError("majorant is defined only for constructed positive orders")
        if q > self.reciprocal_support_edge(order):
            raise ValueError("q lies outside the scale-controlled active interval")
        row = self.jet_bounds[order - 1]
        if derivative_order >= len(row):
            raise ValueError("derivative_order must satisfy m <= j+2")
        bound = row[derivative_order] * q ** (self.h * order)
        target = 0.5 ** order
        return bound, target


def build_slow_borel_cutoff_schedule(
    h: float,
    positive_order_jet_bounds: Sequence[Sequence[float]],
    *,
    initial_lower_bound: int = 0,
) -> SlowBorelCutoffSchedule:
    """Construct the pinned Lean doubling-envelope schedule for a finite prefix.

    Row zero of ``positive_order_jet_bounds`` corresponds to paper order j=1.
    Each row j must contain the actual normalized template bounds C[j,m] for
    m=0,...,j+2. Those constants are analytic inputs from compactness of the
    true coefficient derivatives; this function never estimates them from
    point samples.
    """
    h = _positive_float(h, "h")
    lower = _nonnegative_int(initial_lower_bound, "initial_lower_bound")
    rows = tuple(tuple(float(v) for v in row) for row in positive_order_jet_bounds)

    local = [max(1, lower)]
    scales = [local[0]]
    for j, row in enumerate(rows, start=1):
        b = local_scale_from_template_bounds(h, j, row)
        local.append(b)
        scales.append(max(b, 2 * scales[-1]))

    return SlowBorelCutoffSchedule(h, lower, rows, tuple(local), tuple(scales))


def build_slow_borel_cutoff_schedule_from_provider(
    h: float,
    max_order: int,
    jet_bound_provider: Callable[[int, int], float],
    *,
    initial_lower_bound: int = 0,
) -> SlowBorelCutoffSchedule:
    """Provider adapter for a requested finite prefix of actual C[j,m] bounds."""
    max_order = _nonnegative_int(max_order, "max_order")
    if not callable(jet_bound_provider):
        raise TypeError("jet_bound_provider must be callable")
    rows = [
        [jet_bound_provider(j, m) for m in range(j + 3)]
        for j in range(1, max_order + 1)
    ]
    return build_slow_borel_cutoff_schedule(
        h,
        rows,
        initial_lower_bound=initial_lower_bound,
    )
