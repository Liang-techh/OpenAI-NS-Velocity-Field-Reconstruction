"""Exact finite-prefix support/plateau certificates for Section 5 cutoffs.

The paper/official formalization chooses a strictly increasing (in fact at least
-doubling) integer scale sequence ``a_j`` and multiplies the positive-order
slow stages by a fixed smooth cutoff.  The executable cutoff on main satisfies

    chi(s) = 1  for s <= 1/2,
    chi(s) = 0  for s >= 1.

Consequently, at any fixed ``q > 0`` the scale schedule has three structural
regions: an exact uncut prefix, at most one transition order, and then an exact
zero tail.  The at-most-one transition statement uses only ``a_{j+1} >= 2 a_j``.

This module makes that finite-prefix statement executable using exact rational
comparisons against the *binary64 input value* of ``q``.  It does not claim an
infinite schedule or Proposition 5.3 residual flatness; those still require the
materialized coefficient hierarchy and analytic bounds upstream.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math

from .background_cutoff_schedule import SlowBorelCutoffSchedule


def _positive_finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _scaled_q_exact(scale: int, q: float) -> Fraction:
    """Return ``scale*q`` exactly for the supplied binary64 ``q`` value."""
    return scale * Fraction.from_float(q)


@dataclass(frozen=True)
class CutoffSupportCertificate:
    """Finite-prefix support classification at one positive similarity scale.

    ``plateau_through`` is the largest order whose positive-order cutoff is
    structurally equal to one.  Order zero is always uncut in ``background.py``,
    so ``plateau_through == 0`` means that no positive-order stage is on the
    cutoff plateau.

    ``transition_order`` is the unique order, if present in the constructed
    prefix, with ``1/2 < a_j q < 1``.  ``zero_from`` is the first order in the
    constructed prefix with ``a_j q >= 1``; every later constructed order is
    then also forced to have zero cutoff by monotonicity of the scales.

    A missing ``zero_from`` is deliberately *not* interpreted as an infinite
    tail statement: the finite prefix simply ended before the zero region was
    reached.
    """

    q: float
    max_order: int
    plateau_through: int
    transition_order: int | None
    zero_from: int | None

    def __post_init__(self) -> None:
        q = _positive_finite(self.q, "q")
        if isinstance(self.max_order, bool) or not isinstance(self.max_order, int) or self.max_order < 0:
            raise ValueError("max_order must be a nonnegative integer")
        if (
            isinstance(self.plateau_through, bool)
            or not isinstance(self.plateau_through, int)
            or not 0 <= self.plateau_through <= self.max_order
        ):
            raise ValueError("plateau_through must lie in the constructed prefix")
        for value, name in (
            (self.transition_order, "transition_order"),
            (self.zero_from, "zero_from"),
        ):
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, int)
                or not 1 <= value <= self.max_order
            ):
                raise ValueError(f"{name} must be a positive order inside the constructed prefix")
        if self.transition_order is not None and self.transition_order != self.plateau_through + 1:
            raise ValueError("transition order must immediately follow the exact plateau prefix")
        if self.zero_from is not None:
            expected = self.plateau_through + 1 if self.transition_order is None else self.transition_order + 1
            if self.zero_from != expected:
                raise ValueError("zero tail must immediately follow the plateau/transition region")
        object.__setattr__(self, "q", q)

    @property
    def tail_zero_certified_within_prefix(self) -> bool:
        return self.zero_from is not None

    @property
    def stable_truncation_order_within_prefix(self) -> int | None:
        """Smallest N whose sum through order N equals every longer constructed prefix.

        This is a finite-prefix statement only.  ``None`` means that the
        constructed schedule has not yet reached a forced-zero tail at this q.
        """
        return None if self.zero_from is None else self.zero_from - 1

    def status(self, order: int) -> str:
        """Return ``leading``, ``plateau``, ``transition`` or ``zero``."""
        if isinstance(order, bool) or not isinstance(order, int) or not 0 <= order <= self.max_order:
            raise ValueError("order must lie in the constructed prefix")
        if order == 0:
            return "leading"
        if order <= self.plateau_through:
            return "plateau"
        if self.transition_order == order:
            return "transition"
        if self.zero_from is not None and order >= self.zero_from:
            return "zero"
        raise ValueError("order lies beyond the certified part of this finite prefix")

    def certifies_prefix_stability(self, truncation_order: int) -> bool:
        """Whether truncating here is stable against longer *constructed* prefixes."""
        if isinstance(truncation_order, bool) or not isinstance(truncation_order, int):
            raise ValueError("truncation_order must be an integer")
        if not 0 <= truncation_order <= self.max_order:
            raise ValueError("truncation_order must lie in the constructed prefix")
        stable = self.stable_truncation_order_within_prefix
        return stable is not None and truncation_order >= stable


def classify_cutoff_support(
    schedule: SlowBorelCutoffSchedule,
    q: float,
) -> CutoffSupportCertificate:
    """Classify the exact cutoff geometry of a constructed finite schedule at ``q``.

    Comparisons use ``Fraction.from_float(q)`` and exact integer scales, so a
    value lying exactly on the executable boundaries ``a_j q = 1/2`` or ``1``
    is classified without a tolerance.  The result mirrors the support logic
    behind ``SlowBorelBase.slowStage_locally_zero`` and
    ``SlowBorelBase.slowSum_finite_at_scale`` but deliberately stops at the
    supplied finite prefix.
    """
    if not isinstance(schedule, SlowBorelCutoffSchedule):
        raise TypeError("schedule must be a SlowBorelCutoffSchedule")
    q = _positive_finite(q, "q")
    half = Fraction(1, 2)
    one = Fraction(1, 1)

    plateau_through = 0
    transition_order: int | None = None
    zero_from: int | None = None

    for order in range(1, schedule.max_order + 1):
        scaled = _scaled_q_exact(schedule.scales[order], q)
        if scaled <= half:
            if transition_order is not None:
                raise RuntimeError("doubling schedule cannot return to the cutoff plateau")
            plateau_through = order
            continue
        if scaled < one:
            if transition_order is not None:
                raise RuntimeError("doubling schedule cannot contain two transition orders")
            transition_order = order
            continue
        zero_from = order
        break

    if zero_from is not None:
        for order in range(zero_from, schedule.max_order + 1):
            if _scaled_q_exact(schedule.scales[order], q) < one:
                raise RuntimeError("schedule lost the forced-zero tail")

    return CutoffSupportCertificate(
        q=q,
        max_order=schedule.max_order,
        plateau_through=plateau_through,
        transition_order=transition_order,
        zero_from=zero_from,
    )


def exact_plateau_edge(schedule: SlowBorelCutoffSchedule, order: int) -> Fraction:
    """Return the exact rational edge ``q = 1/(2 a_order)`` for a positive stage."""
    if not isinstance(schedule, SlowBorelCutoffSchedule):
        raise TypeError("schedule must be a SlowBorelCutoffSchedule")
    if isinstance(order, bool) or not isinstance(order, int) or not 1 <= order <= schedule.max_order:
        raise ValueError("order must be a positive order inside the constructed prefix")
    return Fraction(1, 2 * schedule.scales[order])


def exact_zero_edge(schedule: SlowBorelCutoffSchedule, order: int) -> Fraction:
    """Return the exact rational edge ``q = 1/a_order`` for a positive stage."""
    if not isinstance(schedule, SlowBorelCutoffSchedule):
        raise TypeError("schedule must be a SlowBorelCutoffSchedule")
    if isinstance(order, bool) or not isinstance(order, int) or not 1 <= order <= schedule.max_order:
        raise ValueError("order must be a positive order inside the constructed prefix")
    return Fraction(1, schedule.scales[order])
