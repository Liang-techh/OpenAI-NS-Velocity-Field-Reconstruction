"""Exact finite-order residual bookkeeping for the Section 5 slow expansion.

The pinned official Lean file ``NavierStokes/SlowExpansionResidual.lean``
expands a finite slow series exactly.  Once the retained coefficient equations
are separated from the omitted terms, the finite residual has the universal
form

    sum_{n<=N} w_n R_n
      + sum_{i,j<=N, i+j>N} w_{i+j} K_{ij}
      - w_{N+1} A_N,

where ``R_n = L_n + sum_{i+j=n} K_{ij} - previous(A)_n``.  In particular, if
the *actual* retained coefficient recurrences vanish, only the omitted pair
interactions and final shifted-viscosity term remain.

This module makes that algebra executable and provides the corresponding
finite-order slow-weight decay majorant.  It does not assert that caller data
satisfy Eqs. (5.3)--(5.6), does not construct the missing profile-dependent
hierarchy, and does not prove Proposition 5.3 / all-jets-flat residual decay.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Sequence
import math
import sys

from .coordinates import validate_h


def _nonnegative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _checked_vector(values: Sequence[float], length: int, name: str) -> tuple[float, ...]:
    if len(values) != length:
        raise ValueError(f"{name} must have length {length}")
    return tuple(_finite(value, f"{name}[{i}]") for i, value in enumerate(values))


def _checked_pair_matrix(
    values: Sequence[Sequence[float]], order: int
) -> tuple[tuple[float, ...], ...]:
    size = order + 1
    if len(values) != size:
        raise ValueError(f"pair coefficients must have {size} rows")
    return tuple(_checked_vector(row, size, f"pair[{i}]") for i, row in enumerate(values))


def _product(a: float, b: float, name: str) -> float:
    value = a * b
    if not math.isfinite(value):
        raise OverflowError(f"{name} is outside binary64 range")
    return value


def _fsum(values: Sequence[float], name: str) -> float:
    try:
        value = math.fsum(values)
    except OverflowError as exc:
        raise OverflowError(f"{name} is outside binary64 range") from exc
    if not math.isfinite(value):
        raise OverflowError(f"{name} is outside binary64 range")
    return value


def slow_order_exact(h: float, order: int) -> Fraction:
    """Return the manuscript slow order ``lambda_n = 2*n*h`` exactly for h's binary64 value."""
    h = validate_h(h)
    order = _nonnegative_int(order, "order")
    return 2 * order * Fraction.from_float(h)


def slow_weight(q: float, h: float, base_power: float, order: int) -> float:
    """Evaluate ``q**(base_power + 2*order*h)`` for ``0 < q <= 1``.

    This is a numerical evaluator for the theorem's weights.  It deliberately
    rejects underflow instead of returning a spurious zero that could be
    mistaken for an exact residual certificate.
    """
    q = _finite(q, "q")
    if not 0.0 < q <= 1.0:
        raise ValueError("q must satisfy 0 < q <= 1")
    exponent = Fraction.from_float(_finite(base_power, "base_power")) + slow_order_exact(h, order)
    try:
        exponent_float = float(exponent)
    except OverflowError as exc:
        raise OverflowError("slow-weight exponent is outside binary64 range") from exc
    log_value = exponent_float * math.log(q)
    log_max = math.log(sys.float_info.max)
    log_min = math.log(math.ulp(0.0))
    if log_value > log_max:
        raise OverflowError("slow weight overflows binary64")
    if log_value < log_min:
        raise OverflowError("slow weight underflows binary64; keep the result in log form")
    value = math.exp(log_value)
    if not math.isfinite(value) or value == 0.0:
        raise OverflowError("slow weight is outside binary64 range")
    return value


def convolution_at(pair: Sequence[Sequence[float]], order: int, index: int) -> float:
    """Return ``sum_{i+j=index} K_ij`` inside an order-N square truncation."""
    order = _nonnegative_int(order, "order")
    index = _nonnegative_int(index, "index")
    if index > order:
        raise ValueError("retained convolution index must satisfy index <= order")
    matrix = _checked_pair_matrix(pair, order)
    terms = [matrix[i][index - i] for i in range(index + 1)]
    return _fsum(terms, "retained convolution")


@dataclass(frozen=True)
class FiniteTruncationBreakdown:
    """Executable counterpart of ``recurrence_truncation`` at one finite order."""

    order: int
    lhs: float
    retained_recurrence: float
    pair_tail: float
    shifted_tail: float
    rhs: float
    recurrence_values: tuple[float, ...]

    @property
    def algebraic_defect(self) -> float:
        """Floating-point defect of the exact Lean identity; expected only to roundoff."""
        return self.lhs - self.rhs

    @property
    def max_recurrence_abs(self) -> float:
        return max((abs(value) for value in self.recurrence_values), default=0.0)

    @property
    def omitted_remainder(self) -> float:
        """The residual left when every retained recurrence is independently zero."""
        return self.pair_tail + self.shifted_tail


def recurrence_truncation_breakdown(
    order: int,
    weights: Sequence[float],
    linear: Sequence[float],
    pair: Sequence[Sequence[float]],
    shifted: Sequence[float],
) -> FiniteTruncationBreakdown:
    """Evaluate the exact finite recurrence/truncation split from the pinned Lean theorem.

    ``linear[n]`` is ``L_n``, ``pair[i][j]`` is ``K_ij`` and ``shifted[n]``
    is ``A_n``.  The weight vector must reach every index appearing in the
    finite product, hence through ``max(2N, N+1)``.

    No tolerance is used to decide that a recurrence vanishes.  Callers that
    need the ``recurrence_truncation_of_zero`` specialization must establish
    those coefficient identities independently.
    """
    N = _nonnegative_int(order, "order")
    max_weight_index = max(2 * N, N + 1)
    w = _checked_vector(weights, max_weight_index + 1, "weights")
    L = _checked_vector(linear, N + 1, "linear")
    A = _checked_vector(shifted, N + 1, "shifted")
    K = _checked_pair_matrix(pair, N)

    convolutions: list[float] = []
    recurrences: list[float] = []
    for n in range(N + 1):
        conv = _fsum([K[i][n - i] for i in range(n + 1)], f"convolution[{n}]")
        previous = 0.0 if n == 0 else A[n - 1]
        recurrence = _fsum([L[n], conv, -previous], f"recurrence[{n}]")
        convolutions.append(conv)
        recurrences.append(recurrence)

    linear_sum = _fsum(
        [_product(w[n], L[n], f"w[{n}]*L[{n}]") for n in range(N + 1)],
        "linear weighted sum",
    )
    pair_sum = _fsum(
        [
            _product(w[i + j], K[i][j], f"w[{i+j}]*K[{i},{j}]")
            for i in range(N + 1)
            for j in range(N + 1)
        ],
        "pair weighted sum",
    )
    shifted_sum = _fsum(
        [_product(w[n + 1], A[n], f"w[{n+1}]*A[{n}]") for n in range(N + 1)],
        "shifted weighted sum",
    )
    lhs = _fsum([linear_sum, pair_sum, -shifted_sum], "finite expansion residual")

    retained = _fsum(
        [
            _product(w[n], recurrences[n], f"w[{n}]*recurrence[{n}]")
            for n in range(N + 1)
        ],
        "retained recurrence sum",
    )
    pair_tail = _fsum(
        [
            _product(w[i + j], K[i][j], f"tail w[{i+j}]*K[{i},{j}]")
            for i in range(N + 1)
            for j in range(N + 1)
            if i + j > N
        ],
        "omitted pair tail",
    )
    shifted_tail = -_product(w[N + 1], A[N], f"w[{N+1}]*A[{N}]")
    rhs = _fsum([retained, pair_tail, shifted_tail], "truncation right-hand side")

    return FiniteTruncationBreakdown(
        order=N,
        lhs=lhs,
        retained_recurrence=retained,
        pair_tail=pair_tail,
        shifted_tail=shifted_tail,
        rhs=rhs,
        recurrence_values=tuple(recurrences),
    )


@dataclass(frozen=True)
class SlowTailMajorant:
    """Conditional finite-order majorant for the omitted terms only.

    If the retained coefficient recurrences are independently known to vanish,
    ``coefficient_l1 * q**exponent`` also bounds the full finite residual in
    this scalar recurrence row.  This object itself does not prove that premise.
    """

    q: float
    h: float
    base_power: float
    order: int
    exponent_exact: Fraction
    coefficient_l1: float
    log_bound: float

    def binary64_bound(self) -> float:
        """Return a slightly outward-rounded binary64 value when representable."""
        if self.coefficient_l1 == 0.0:
            return 0.0
        log_max = math.log(sys.float_info.max)
        log_min = math.log(math.ulp(0.0))
        if self.log_bound > log_max:
            raise OverflowError("slow tail majorant overflows binary64")
        if self.log_bound < log_min:
            raise OverflowError("slow tail majorant underflows binary64; use log_bound")
        value = math.exp(self.log_bound)
        if not math.isfinite(value) or value == 0.0:
            raise OverflowError("slow tail majorant is outside binary64 range")
        return math.nextafter(value, math.inf)


def omitted_slow_tail_majorant(
    q: float,
    h: float,
    base_power: float,
    order: int,
    pair: Sequence[Sequence[float]],
    shifted: Sequence[float],
) -> SlowTailMajorant:
    """Bound the omitted pair tail and last shifted term by their first slow order.

    For ``w_n = q**(b + 2*n*h)`` and ``0 < q <= 1``, every omitted pair has
    ``i+j >= N+1``.  Therefore

      |pairTail_N - w_{N+1} A_N|
        <= q**(b + 2*(N+1)*h)
           * (sum_{i+j>N}|K_ij| + |A_N|).

    The returned bound is only for the explicit omitted terms.  It becomes a
    bound for the full finite expansion residual only after the retained
    recurrence identities are established from the actual coefficient hierarchy.
    """
    N = _nonnegative_int(order, "order")
    q = _finite(q, "q")
    if not 0.0 < q <= 1.0:
        raise ValueError("q must satisfy 0 < q <= 1")
    h = validate_h(h)
    base = _finite(base_power, "base_power")
    K = _checked_pair_matrix(pair, N)
    A = _checked_vector(shifted, N + 1, "shifted")

    coefficient_l1 = _fsum(
        [abs(K[i][j]) for i in range(N + 1) for j in range(N + 1) if i + j > N]
        + [abs(A[N])],
        "omitted coefficient l1 bound",
    )
    exponent = Fraction.from_float(base) + slow_order_exact(h, N + 1)
    if coefficient_l1 == 0.0:
        log_bound = -math.inf
    else:
        try:
            exponent_float = float(exponent)
        except OverflowError as exc:
            raise OverflowError("slow-tail exponent is outside binary64 range") from exc
        log_bound = math.log(coefficient_l1) + exponent_float * math.log(q)
        if math.isnan(log_bound):
            raise OverflowError("slow tail majorant is not numerically representable")

    return SlowTailMajorant(
        q=q,
        h=h,
        base_power=base,
        order=N,
        exponent_exact=exponent,
        coefficient_l1=coefficient_l1,
        log_bound=log_bound,
    )
