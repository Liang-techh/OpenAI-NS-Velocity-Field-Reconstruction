"""Target-driven finite Section 5 hierarchy closure from owned provider data.

This module composes three already-landed/stacked fail-closed ingredients:

* the pinned ``DiagonalScale.exists_uniform_tail_order`` order selector;
* hierarchy-owned ``C[j,m]`` derivative-bound rows feeding the recursive
  SlowBorel/DiagonalScale cutoff schedule; and
* exact retained-recurrence cancellation identities.

The important restriction is that callers do *not* pass ``C[j,m]`` arrays or a
list of formal identities.  One provider object owns the hierarchy identity,
source revision, coefficient state, every requested jet-bound row, and every
retained recurrence identity.  For a requested ordinary-jet derivative order
``m`` and target q-power ``N``, this module first selects the minimal truncation
order ``J`` satisfying the pinned scalar condition

    h * (J + 1) - m >= N,

with the additional theorem-side admissibility requirement ``m <= J + 3`` from
``SlowBorelBase.ordinary_tail_bound``.  It then queries exactly the finite
hierarchy prefix needed through order ``J`` and constructs the pinned recursive
cutoff schedule plus a contiguous exact recurrence-cancellation prefix.

This is intentionally not an all-order theorem.  If the real hierarchy API has
a finite derivative/order frontier, the provider query raises and no partial
certificate is returned.  Successful construction for one requested target
proves only that finite target-driven prefix; it does not prove that the
provider can answer every target, that an infinite diagonal schedule exists,
or that the actual PDE residual is all-jets-flat or super-algebraically small.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Protocol, runtime_checkable

from .background_hierarchy_cutoff_bounds import (
    HierarchyJetBoundProvider,
    HierarchySlowBorelPrefixWitness,
    build_slow_borel_schedule_from_hierarchy_certificate,
    certify_hierarchy_bound_prefix,
)
from .background_recurrence_cancellation import (
    ExactRetainedRecurrenceIdentity,
    FiniteRetainedCancellationCertificate,
)
from .background_uniform_tail_order import (
    PinnedUniformTailOrderCertificate,
    certify_pinned_uniform_tail_order,
)


def _nonnegative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _nonempty_text(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value.strip()


@runtime_checkable
class TargetDrivenHierarchyProvider(HierarchyJetBoundProvider, Protocol):
    """Owned hierarchy interface required for one target-driven finite prefix.

    Agent 2's eventual real hierarchy backend should implement this interface
    directly (or expose an object with the same attributes/methods).  The
    interface deliberately contains no setter or caller-supplied coefficient
    rows: all data must be produced by the hierarchy implementation itself.
    """

    coefficient_state_id: str

    def certified_retained_identity(self, order: int) -> ExactRetainedRecurrenceIdentity:
        ...


@dataclass(frozen=True)
class TargetDrivenHierarchyPrefixCertificate:
    """One finite hierarchy prefix selected by an exact requested tail power."""

    derivative_order: int
    coefficient_state_id: str
    target: PinnedUniformTailOrderCertificate
    cutoff: HierarchySlowBorelPrefixWitness
    cancellations: FiniteRetainedCancellationCertificate

    def __post_init__(self) -> None:
        derivative_order = _nonnegative_int(self.derivative_order, "derivative_order")
        coefficient_state_id = _nonempty_text(
            self.coefficient_state_id, "coefficient_state_id"
        )
        if not isinstance(self.target, PinnedUniformTailOrderCertificate):
            raise TypeError("target must be a PinnedUniformTailOrderCertificate")
        if not isinstance(self.cutoff, HierarchySlowBorelPrefixWitness):
            raise TypeError("cutoff must be a HierarchySlowBorelPrefixWitness")
        if not isinstance(self.cancellations, FiniteRetainedCancellationCertificate):
            raise TypeError("cancellations must be a FiniteRetainedCancellationCertificate")

        if self.target.derivative_loss != Fraction(derivative_order, 1):
            raise ValueError("target derivative loss must equal the ordinary-jet derivative order")
        if self.target.order < 1:
            raise ValueError("target-driven hierarchy prefix requires at least positive order one")
        if derivative_order > self.target.order + 3:
            raise ValueError("ordinary-tail theorem requires m <= J+3")
        if self.cutoff.bounds.max_order != self.target.order:
            raise ValueError("cutoff prefix order does not match the selected target order")
        if self.cancellations.max_order != self.target.order:
            raise ValueError("recurrence prefix order does not match the selected target order")

        bounds = self.cutoff.bounds
        cancellations = self.cancellations
        if bounds.hierarchy_id != cancellations.hierarchy_id:
            raise ValueError("cutoff bounds and recurrence identities use different hierarchies")
        if bounds.source_revision != cancellations.source_revision:
            raise ValueError("cutoff bounds and recurrence identities use different source revisions")
        if cancellations.coefficient_state_id != coefficient_state_id:
            raise ValueError("recurrence identities do not match the provider coefficient state")

        expected_power = self.target.h_exact * (self.target.order + 1) - derivative_order
        if self.target.tail_power_exact != expected_power:
            raise ValueError("target certificate does not match the ordinary-tail exponent")
        if expected_power < self.target.target_power:
            raise ValueError("selected finite prefix does not attain the requested target power")

        object.__setattr__(self, "derivative_order", derivative_order)
        object.__setattr__(self, "coefficient_state_id", coefficient_state_id)

    @property
    def hierarchy_id(self) -> str:
        return self.cutoff.bounds.hierarchy_id

    @property
    def source_revision(self) -> str:
        return self.cutoff.bounds.source_revision

    @property
    def selected_order(self) -> int:
        return self.target.order

    @property
    def ordinary_tail_power_exact(self) -> Fraction:
        return self.target.tail_power_exact

    @property
    def derivative_budget_margin(self) -> int:
        """Return the exact integer margin ``J+3-m`` from ordinary_tail_bound."""
        return self.selected_order + 3 - self.derivative_order

    @property
    def finite_target_prefix_only(self) -> bool:
        return True

    @property
    def all_order_hierarchy_verified(self) -> bool:
        return False

    @property
    def infinite_diagonal_schedule_verified(self) -> bool:
        return False

    @property
    def pde_residual_tail_verified(self) -> bool:
        return False

    @property
    def all_jets_flat(self) -> bool:
        return False

    @property
    def super_algebraic(self) -> bool:
        return False

    @property
    def paper_exact(self) -> bool:
        return False


def certify_target_driven_hierarchy_prefix(
    h: float,
    derivative_order: int,
    target_power: Fraction | int,
    provider: TargetDrivenHierarchyProvider,
    *,
    minimum_order: int = 0,
    initial_lower_bound: int = 0,
) -> TargetDrivenHierarchyPrefixCertificate:
    """Query exactly the hierarchy prefix required for one ordinary-jet target.

    The selected order is the minimal pinned target order subject to both the
    caller's finite ``minimum_order`` and the actual theorem-side condition
    ``m <= J+3``.  After selection, all ``C[j,k]`` rows and retained identities
    are obtained only through ``provider``.  A provider frontier therefore
    fails closed at the first missing hierarchy-owned artifact.
    """

    derivative_order = _nonnegative_int(derivative_order, "derivative_order")
    minimum_order = _nonnegative_int(minimum_order, "minimum_order")
    if not isinstance(provider, TargetDrivenHierarchyProvider):
        raise TypeError("provider must satisfy TargetDrivenHierarchyProvider")

    hierarchy_id = _nonempty_text(provider.hierarchy_id, "provider.hierarchy_id")
    source_revision = _nonempty_text(provider.source_revision, "provider.source_revision")
    coefficient_state_id = _nonempty_text(
        provider.coefficient_state_id, "provider.coefficient_state_id"
    )

    theorem_minimum = max(1, minimum_order, max(0, derivative_order - 3))
    target = certify_pinned_uniform_tail_order(
        h,
        derivative_order,
        target_power,
        minimum_order=theorem_minimum,
    )

    bounds = certify_hierarchy_bound_prefix(provider, target.order)
    cutoff = build_slow_borel_schedule_from_hierarchy_certificate(
        target.h,
        bounds,
        initial_lower_bound=initial_lower_bound,
    )

    identities: list[ExactRetainedRecurrenceIdentity] = []
    for order in range(target.order + 1):
        identity = provider.certified_retained_identity(order)
        if not isinstance(identity, ExactRetainedRecurrenceIdentity):
            raise TypeError("provider must return ExactRetainedRecurrenceIdentity values")
        if identity.order != order:
            raise ValueError("provider returned a retained recurrence witness for the wrong order")
        if identity.hierarchy_id != hierarchy_id:
            raise ValueError("provider returned a retained identity for a different hierarchy")
        if identity.source_revision != source_revision:
            raise ValueError("provider returned a retained identity for a different source revision")
        if identity.coefficient_state_id != coefficient_state_id:
            raise ValueError("provider returned a retained identity for a different coefficient state")
        identities.append(identity)

    cancellations = FiniteRetainedCancellationCertificate(
        hierarchy_id=hierarchy_id,
        source_revision=source_revision,
        coefficient_state_id=coefficient_state_id,
        max_order=target.order,
        identities=tuple(identities),
    )

    return TargetDrivenHierarchyPrefixCertificate(
        derivative_order=derivative_order,
        coefficient_state_id=coefficient_state_id,
        target=target,
        cutoff=cutoff,
        cancellations=cancellations,
    )
