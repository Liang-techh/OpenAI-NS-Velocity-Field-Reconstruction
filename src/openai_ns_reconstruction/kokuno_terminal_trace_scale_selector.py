"""Clean-room exact audit of the public terminal-trace scale selection.

The public corrected reader bounds each differentiated cutoff-monomial term by
an exact finite constant and chooses the least increasing integer scale that
makes a finite family of such bounds at most ``2**(-j)``.  This module
reimplements only that finite arithmetic seam.  It does not import source
checker/proof code and it does not certify the terminal jets or full endpoint.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb, factorial
from typing import Mapping


class TerminalTraceScaleError(ValueError):
    """Raised when an exact finite terminal-trace scale cannot be selected."""


@dataclass(frozen=True)
class TerminalTraceConstraint:
    """One exact E.11-style terminal cutoff constraint."""

    spatial_order: int
    time_order: int
    cutoff_constant: Fraction
    jet_norm: Fraction
    bound: Fraction
    target: Fraction

    @property
    def satisfied(self) -> bool:
        return self.bound <= self.target


@dataclass(frozen=True)
class TerminalTraceScaleSelection:
    """Exact least-integer result for one finite terminal stage."""

    stage_index: int
    previous_scale: int
    scale: int
    constraints: tuple[TerminalTraceConstraint, ...]

    @property
    def max_target_ratio(self) -> Fraction:
        if not self.constraints:
            return Fraction(0, 1)
        return max(item.bound / item.target for item in self.constraints)


def _exact_nonnegative(value: Fraction, *, name: str) -> Fraction:
    if type(value) is not Fraction:
        raise TypeError(f"{name} must be fractions.Fraction")
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _exact_int(value: int, *, name: str, minimum: int) -> int:
    if type(value) is not int:
        raise TypeError(f"{name} must be int")
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return value


def initial_terminal_scale() -> int:
    """Return the public terminal-extension initialization ``b_0 = 1``."""

    return 1


def cutoff_derivative_constant(
    *,
    stage_index: int,
    time_order: int,
    cutoff_derivative_norms: Mapping[int, Fraction],
) -> Fraction:
    """Return the exact finite cutoff constant ``C_{j,m}``.

    The input table supplies exact certified upper bounds for
    ``||chi_0^(a)||_infinity``.  Only entries actually required by the finite
    sum are read.
    """

    j = _exact_int(stage_index, name="stage_index", minimum=0)
    m = _exact_int(time_order, name="time_order", minimum=0)
    if m > j:
        raise ValueError("time_order must not exceed stage_index")

    total = Fraction(0, 1)
    for a in range(m + 1):
        if m - a > j:
            continue
        try:
            norm = cutoff_derivative_norms[a]
        except KeyError as exc:
            raise TerminalTraceScaleError(
                f"missing cutoff derivative norm for order {a}"
            ) from exc
        norm = _exact_nonnegative(norm, name=f"cutoff_derivative_norms[{a}]")
        total += Fraction(comb(m, a), factorial(j - m + a)) * norm
    return total


def terminal_monomial_bound(
    *,
    stage_index: int,
    spatial_order: int,
    time_order: int,
    scale: int,
    cutoff_derivative_norms: Mapping[int, Fraction],
    jet_norms: Mapping[int, Fraction],
) -> TerminalTraceConstraint:
    """Evaluate one exact differentiated cutoff-monomial upper bound."""

    j = _exact_int(stage_index, name="stage_index", minimum=1)
    a0 = _exact_int(spatial_order, name="spatial_order", minimum=0)
    m = _exact_int(time_order, name="time_order", minimum=0)
    b = _exact_int(scale, name="scale", minimum=1)
    if a0 + m > j // 2:
        raise ValueError("constraint is outside the public finite E.11 index set")

    try:
        jet_norm = jet_norms[a0]
    except KeyError as exc:
        raise TerminalTraceScaleError(f"missing jet norm for order {a0}") from exc
    jet_norm = _exact_nonnegative(jet_norm, name=f"jet_norms[{a0}]")
    constant = cutoff_derivative_constant(
        stage_index=j,
        time_order=m,
        cutoff_derivative_norms=cutoff_derivative_norms,
    )
    exponent_gap = j - m
    if exponent_gap <= 0:
        raise AssertionError("E.11 finite index set must have m < j for j >= 1")
    bound = constant * jet_norm / (b**exponent_gap)
    target = Fraction(1, 2**j)
    return TerminalTraceConstraint(
        spatial_order=a0,
        time_order=m,
        cutoff_constant=constant,
        jet_norm=jet_norm,
        bound=bound,
        target=target,
    )


def terminal_constraints(
    *,
    stage_index: int,
    scale: int,
    cutoff_derivative_norms: Mapping[int, Fraction],
    jet_norms: Mapping[int, Fraction],
) -> tuple[TerminalTraceConstraint, ...]:
    """Return the complete finite E.11 constraint family for one stage."""

    j = _exact_int(stage_index, name="stage_index", minimum=1)
    b = _exact_int(scale, name="scale", minimum=1)
    half = j // 2
    constraints: list[TerminalTraceConstraint] = []
    for a0 in range(half + 1):
        for m in range(half - a0 + 1):
            constraints.append(
                terminal_monomial_bound(
                    stage_index=j,
                    spatial_order=a0,
                    time_order=m,
                    scale=b,
                    cutoff_derivative_norms=cutoff_derivative_norms,
                    jet_norms=jet_norms,
                )
            )
    return tuple(constraints)


def select_terminal_scale(
    *,
    stage_index: int,
    previous_scale: int,
    cutoff_derivative_norms: Mapping[int, Fraction],
    jet_norms: Mapping[int, Fraction],
    max_scale: int,
) -> TerminalTraceScaleSelection:
    """Select the least integer ``b_j > b_{j-1}`` satisfying all constraints.

    ``max_scale`` is an explicit fail-closed search cap.  The source existence
    argument says a sufficiently large integer exists for finite exact input
    bounds; this function does not silently replace a missing certified bound
    or an exhausted search by an approximate choice.
    """

    j = _exact_int(stage_index, name="stage_index", minimum=1)
    previous = _exact_int(previous_scale, name="previous_scale", minimum=1)
    cap = _exact_int(max_scale, name="max_scale", minimum=1)
    if cap <= previous:
        raise TerminalTraceScaleError("max_scale must exceed previous_scale")

    half = j // 2
    for order in range(half + 1):
        if order not in jet_norms:
            raise TerminalTraceScaleError(f"missing jet norm for order {order}")
        _exact_nonnegative(jet_norms[order], name=f"jet_norms[{order}]")
        if order not in cutoff_derivative_norms:
            raise TerminalTraceScaleError(
                f"missing cutoff derivative norm for order {order}"
            )
        _exact_nonnegative(
            cutoff_derivative_norms[order],
            name=f"cutoff_derivative_norms[{order}]",
        )

    for candidate in range(previous + 1, cap + 1):
        constraints = terminal_constraints(
            stage_index=j,
            scale=candidate,
            cutoff_derivative_norms=cutoff_derivative_norms,
            jet_norms=jet_norms,
        )
        if all(item.satisfied for item in constraints):
            return TerminalTraceScaleSelection(
                stage_index=j,
                previous_scale=previous,
                scale=candidate,
                constraints=constraints,
            )
    raise TerminalTraceScaleError(
        f"no exact terminal scale found in ({previous}, {cap}]"
    )
