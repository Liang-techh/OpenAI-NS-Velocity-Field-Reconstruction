"""Dependency-tracked finite-prefix C[j,m] bounds for Section 5 SlowBorel.

The existing ``background_cutoff_schedule`` layer correctly implements the
pinned SlowBorel/DiagonalScale scale inequality once numerical constants
``C[j,m]`` are supplied. This module closes one narrower provenance seam: it
requires every finite-prefix bound to carry an explicit hierarchy identity,
source revision, and coefficient dependency chain before those constants can
enter the cutoff scheduler.

This is intentionally *not* a derivation of the paper's uniform all-order
bounds. Until a genuine hierarchy bound provider proves every required
``C[j,m]`` row, the construction remains finite-prefix ``formal-structure``.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral
from typing import Protocol, runtime_checkable

from .background_cutoff_schedule import (
    SlowBorelCutoffSchedule,
    build_slow_borel_cutoff_schedule,
)


def _nonnegative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _positive_int(value: int, name: str) -> int:
    value = _nonnegative_int(value, name)
    if value == 0:
        raise ValueError(f"{name} must be positive")
    return value


def _nonempty_text(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value.strip()


@dataclass(frozen=True)
class HierarchyBoundDependency:
    """One coefficient artifact used to justify a normalized-template bound."""

    coefficient_order: int
    artifact: str
    provider: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "coefficient_order", _nonnegative_int(self.coefficient_order, "coefficient_order")
        )
        object.__setattr__(self, "artifact", _nonempty_text(self.artifact, "artifact"))
        object.__setattr__(self, "provider", _nonempty_text(self.provider, "provider"))


@dataclass(frozen=True)
class HierarchyJetBound:
    """A dependency-backed analytic upper bound for one ``C[j,m]`` entry."""

    order: int
    derivative_order: int
    bound: float
    hierarchy_id: str
    source_revision: str
    dependencies: tuple[HierarchyBoundDependency, ...]

    def __post_init__(self) -> None:
        order = _positive_int(self.order, "order")
        derivative_order = _nonnegative_int(self.derivative_order, "derivative_order")
        if derivative_order > order + 2:
            raise ValueError("derivative_order must satisfy m <= j+2")
        bound = float(self.bound)
        if not math.isfinite(bound) or bound <= 0.0:
            raise ValueError("bound must be finite and positive")
        hierarchy_id = _nonempty_text(self.hierarchy_id, "hierarchy_id")
        source_revision = _nonempty_text(self.source_revision, "source_revision")
        dependencies = tuple(self.dependencies)
        if not dependencies:
            raise ValueError("every C[j,m] bound must carry at least one hierarchy dependency")
        if not all(isinstance(dep, HierarchyBoundDependency) for dep in dependencies):
            raise TypeError("dependencies must contain HierarchyBoundDependency values")
        keys = [(dep.coefficient_order, dep.artifact, dep.provider) for dep in dependencies]
        if len(keys) != len(set(keys)):
            raise ValueError("hierarchy dependencies must be unique")
        if any(dep.coefficient_order > order for dep in dependencies):
            raise ValueError("a C[j,m] bound cannot depend on a future coefficient order")
        if not any(dep.coefficient_order == order for dep in dependencies):
            raise ValueError("a C[j,m] bound must bind the coefficient at its own order")
        object.__setattr__(self, "order", order)
        object.__setattr__(self, "derivative_order", derivative_order)
        object.__setattr__(self, "bound", bound)
        object.__setattr__(self, "hierarchy_id", hierarchy_id)
        object.__setattr__(self, "source_revision", source_revision)
        object.__setattr__(self, "dependencies", dependencies)


@runtime_checkable
class HierarchyJetBoundProvider(Protocol):
    """Protocol for an analytic hierarchy layer that owns ``C[j,m]`` proofs."""

    hierarchy_id: str
    source_revision: str

    def certified_jet_bound(self, order: int, derivative_order: int) -> HierarchyJetBound:
        ...


@dataclass(frozen=True)
class FiniteHierarchyBoundCertificate:
    """Exact triangular finite prefix of dependency-backed ``C[j,m]`` bounds."""

    hierarchy_id: str
    source_revision: str
    max_order: int
    entries: tuple[HierarchyJetBound, ...]

    def __post_init__(self) -> None:
        hierarchy_id = _nonempty_text(self.hierarchy_id, "hierarchy_id")
        source_revision = _nonempty_text(self.source_revision, "source_revision")
        max_order = _positive_int(self.max_order, "max_order")
        entries = tuple(self.entries)
        if not all(isinstance(entry, HierarchyJetBound) for entry in entries):
            raise TypeError("entries must contain HierarchyJetBound values")

        expected = {
            (j, m)
            for j in range(1, max_order + 1)
            for m in range(j + 3)
        }
        actual = [(entry.order, entry.derivative_order) for entry in entries]
        if len(actual) != len(set(actual)):
            raise ValueError("finite hierarchy bound certificate contains duplicate C[j,m] entries")
        if set(actual) != expected:
            missing = sorted(expected.difference(actual))
            extra = sorted(set(actual).difference(expected))
            raise ValueError(
                f"finite hierarchy bound certificate is not a complete triangular prefix; "
                f"missing={missing}, extra={extra}"
            )
        for entry in entries:
            if entry.hierarchy_id != hierarchy_id:
                raise ValueError("C[j,m] entry hierarchy_id does not match its certificate")
            if entry.source_revision != source_revision:
                raise ValueError("C[j,m] entry source_revision does not match its certificate")

        ordered = tuple(sorted(entries, key=lambda entry: (entry.order, entry.derivative_order)))
        object.__setattr__(self, "hierarchy_id", hierarchy_id)
        object.__setattr__(self, "source_revision", source_revision)
        object.__setattr__(self, "max_order", max_order)
        object.__setattr__(self, "entries", ordered)

    @property
    def paper_exact(self) -> bool:
        return False

    def entry(self, order: int, derivative_order: int) -> HierarchyJetBound:
        order = _positive_int(order, "order")
        derivative_order = _nonnegative_int(derivative_order, "derivative_order")
        if order > self.max_order or derivative_order > order + 2:
            raise ValueError("requested C[j,m] lies outside the certified finite prefix")
        index = sum(k + 3 for k in range(1, order)) + derivative_order
        entry = self.entries[index]
        if (entry.order, entry.derivative_order) != (order, derivative_order):
            raise RuntimeError("internal certificate ordering invariant was violated")
        return entry

    def rows(self) -> tuple[tuple[float, ...], ...]:
        return tuple(
            tuple(self.entry(j, m).bound for m in range(j + 3))
            for j in range(1, self.max_order + 1)
        )


def certify_hierarchy_bound_prefix(
    provider: HierarchyJetBoundProvider,
    max_order: int,
) -> FiniteHierarchyBoundCertificate:
    """Query exactly the SlowBorel triangular set from one hierarchy provider."""

    max_order = _positive_int(max_order, "max_order")
    if not isinstance(provider, HierarchyJetBoundProvider):
        raise TypeError("provider must satisfy HierarchyJetBoundProvider")
    hierarchy_id = _nonempty_text(provider.hierarchy_id, "provider.hierarchy_id")
    source_revision = _nonempty_text(provider.source_revision, "provider.source_revision")

    entries: list[HierarchyJetBound] = []
    for j in range(1, max_order + 1):
        for m in range(j + 3):
            entry = provider.certified_jet_bound(j, m)
            if not isinstance(entry, HierarchyJetBound):
                raise TypeError("provider must return HierarchyJetBound values")
            if (entry.order, entry.derivative_order) != (j, m):
                raise ValueError("provider returned a C[j,m] witness for the wrong index")
            if entry.hierarchy_id != hierarchy_id:
                raise ValueError("provider returned a C[j,m] witness for a different hierarchy")
            if entry.source_revision != source_revision:
                raise ValueError("provider returned a C[j,m] witness for a different source revision")
            entries.append(entry)

    return FiniteHierarchyBoundCertificate(
        hierarchy_id=hierarchy_id,
        source_revision=source_revision,
        max_order=max_order,
        entries=tuple(entries),
    )


@dataclass(frozen=True)
class HierarchySlowBorelPrefixWitness:
    """Cutoff schedule plus the exact hierarchy evidence that supplied its rows."""

    bounds: FiniteHierarchyBoundCertificate
    schedule: SlowBorelCutoffSchedule

    def __post_init__(self) -> None:
        if not isinstance(self.bounds, FiniteHierarchyBoundCertificate):
            raise TypeError("bounds must be a FiniteHierarchyBoundCertificate")
        if not isinstance(self.schedule, SlowBorelCutoffSchedule):
            raise TypeError("schedule must be a SlowBorelCutoffSchedule")
        if self.schedule.max_order != self.bounds.max_order:
            raise ValueError("schedule order does not match the hierarchy bound certificate")
        if self.schedule.jet_bounds != self.bounds.rows():
            raise ValueError("schedule C[j,m] rows do not match the hierarchy bound certificate")
        for j in range(1, self.bounds.max_order + 1):
            for m in range(j + 3):
                if self.schedule.edge_log_margin(j, m) < -64.0 * math.ulp(1.0):
                    raise ValueError("schedule does not certify one hierarchy-bound dyadic edge")

    @property
    def paper_exact(self) -> bool:
        return False

    def provenance(self, order: int, derivative_order: int) -> HierarchyJetBound:
        return self.bounds.entry(order, derivative_order)


def build_slow_borel_schedule_from_hierarchy_certificate(
    h: float,
    bounds: FiniteHierarchyBoundCertificate,
    *,
    initial_lower_bound: int = 0,
) -> HierarchySlowBorelPrefixWitness:
    """Feed only dependency-backed finite-prefix bounds into the landed scheduler."""

    if not isinstance(bounds, FiniteHierarchyBoundCertificate):
        raise TypeError("bounds must be a FiniteHierarchyBoundCertificate")
    schedule = build_slow_borel_cutoff_schedule(
        h,
        bounds.rows(),
        initial_lower_bound=initial_lower_bound,
    )
    return HierarchySlowBorelPrefixWitness(bounds=bounds, schedule=schedule)


def build_slow_borel_schedule_from_hierarchy_provider(
    h: float,
    max_order: int,
    provider: HierarchyJetBoundProvider,
    *,
    initial_lower_bound: int = 0,
) -> HierarchySlowBorelPrefixWitness:
    """One-shot fail-closed provider -> provenance certificate -> cutoff schedule."""

    bounds = certify_hierarchy_bound_prefix(provider, max_order)
    return build_slow_borel_schedule_from_hierarchy_certificate(
        h,
        bounds,
        initial_lower_bound=initial_lower_bound,
    )
