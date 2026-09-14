"""Outward-safe exact C[j,m] majorants for the target-driven Section 5 prefix.

The landed hierarchy-bound scheduler stores executable ``C[j,m]`` constants as
binary64 floats.  That is adequate only after a theorem-backed upper majorant
has been rounded *outward*: a nearest-float conversion may round a positive
real/rational upper bound downward before the recursive SlowBorel scale
inequality is evaluated.

This module closes that finite-prefix arithmetic seam without changing the
pinned theorem hypotheses.  The hierarchy provider must own an exact positive
integer/Fraction majorant for every requested ``C[j,m]`` row, together with the
same hierarchy/source/coefficient-state and coefficient dependency provenance
used by support and retained-recurrence evidence.  An internal adapter converts
each exact majorant to the least readily available finite binary64 value that is
verified to dominate it, then runs the existing target-driven dependency chain.

The exact rational is an implementation-level certified upper majorant; this
module does not claim that the source theorem's existential constant is itself
rational.  A provider with only a floating estimate cannot enter this path.
Successful construction remains one finite target-selected prefix and does not
prove total all-order hierarchy availability, an infinite diagonal schedule,
the actual PDE residual tail, all-jets flatness, or super-algebraic convergence.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from numbers import Integral
from typing import Protocol, runtime_checkable

from .background_hierarchy_cutoff_bounds import (
    HierarchyBoundDependency,
    HierarchyJetBound,
)
from .background_recurrence_cancellation import ExactRetainedRecurrenceIdentity
from .background_target_driven_dependency_chain import (
    TargetDrivenHierarchyDependencyChainCertificate,
    certify_target_driven_dependency_chain,
)
from .background_target_driven_support_prefix import HierarchyCoefficientSupport


def _positive_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _nonnegative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _nonempty_text(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value.strip()


def _positive_exact(value: Fraction | int, name: str) -> Fraction:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact integer/Fraction; booleans are forbidden")
    if isinstance(value, Fraction):
        result = value
    elif isinstance(value, Integral):
        result = Fraction(int(value), 1)
    else:
        raise TypeError(f"{name} must be an exact integer/Fraction; floats are forbidden")
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def outward_float_upper(value: Fraction | int) -> float:
    """Return a finite binary64 value proved to be >= the exact positive input."""

    exact = _positive_exact(value, "value")
    try:
        rounded = float(exact)
    except OverflowError as exc:
        raise OverflowError("exact majorant exceeds the finite binary64 scheduler range") from exc
    if not math.isfinite(rounded):
        raise OverflowError("exact majorant exceeds the finite binary64 scheduler range")

    if Fraction.from_float(rounded) < exact:
        rounded = math.nextafter(rounded, math.inf)
    if not math.isfinite(rounded) or Fraction.from_float(rounded) < exact:
        raise OverflowError("unable to construct a finite outward binary64 majorant")
    return rounded


@dataclass(frozen=True)
class ExactHierarchyJetMajorant:
    """Hierarchy-owned exact upper majorant for one finite-prefix ``C[j,m]``."""

    order: int
    derivative_order: int
    bound: Fraction | int
    hierarchy_id: str
    source_revision: str
    coefficient_state_id: str
    dependencies: tuple[HierarchyBoundDependency, ...]

    def __post_init__(self) -> None:
        order = _positive_int(self.order, "order")
        derivative_order = _nonnegative_int(self.derivative_order, "derivative_order")
        if derivative_order > order + 2:
            raise ValueError("derivative_order must satisfy m <= j+2")
        bound = _positive_exact(self.bound, "bound")
        hierarchy_id = _nonempty_text(self.hierarchy_id, "hierarchy_id")
        source_revision = _nonempty_text(self.source_revision, "source_revision")
        coefficient_state_id = _nonempty_text(
            self.coefficient_state_id, "coefficient_state_id"
        )
        dependencies = tuple(self.dependencies)
        if not dependencies:
            raise ValueError("every exact C[j,m] majorant must carry a hierarchy dependency")
        if not all(isinstance(dep, HierarchyBoundDependency) for dep in dependencies):
            raise TypeError("dependencies must contain HierarchyBoundDependency values")
        keys = [
            (dep.coefficient_order, dep.artifact, dep.provider)
            for dep in dependencies
        ]
        if len(keys) != len(set(keys)):
            raise ValueError("exact-majorant hierarchy dependencies must be unique")
        if any(dep.coefficient_order > order for dep in dependencies):
            raise ValueError("an exact C[j,m] majorant cannot depend on a future coefficient")
        if not any(dep.coefficient_order == order for dep in dependencies):
            raise ValueError("an exact C[j,m] majorant must bind its own coefficient order")

        object.__setattr__(self, "order", order)
        object.__setattr__(self, "derivative_order", derivative_order)
        object.__setattr__(self, "bound", bound)
        object.__setattr__(self, "hierarchy_id", hierarchy_id)
        object.__setattr__(self, "source_revision", source_revision)
        object.__setattr__(self, "coefficient_state_id", coefficient_state_id)
        object.__setattr__(self, "dependencies", dependencies)

    @property
    def outward_bound(self) -> float:
        return outward_float_upper(self.bound)

    def as_hierarchy_jet_bound(self) -> HierarchyJetBound:
        """Materialize the existing scheduler row with verified outward rounding."""

        return HierarchyJetBound(
            order=self.order,
            derivative_order=self.derivative_order,
            bound=self.outward_bound,
            hierarchy_id=self.hierarchy_id,
            source_revision=self.source_revision,
            dependencies=self.dependencies,
        )


@runtime_checkable
class TargetDrivenExactMajorantProvider(Protocol):
    """Hierarchy provider whose executable C rows originate in exact majorants."""

    hierarchy_id: str
    source_revision: str
    coefficient_state_id: str

    def certified_exact_jet_majorant(
        self, order: int, derivative_order: int
    ) -> ExactHierarchyJetMajorant:
        ...

    def certified_retained_identity(
        self, order: int
    ) -> ExactRetainedRecurrenceIdentity:
        ...

    def certified_coefficient_support(
        self, order: int
    ) -> HierarchyCoefficientSupport:
        ...


class _OutwardRoundedHierarchyProvider:
    """Private adapter from exact provider rows to the existing finite scheduler."""

    def __init__(self, provider: TargetDrivenExactMajorantProvider) -> None:
        self._provider = provider
        self.hierarchy_id = _nonempty_text(provider.hierarchy_id, "provider.hierarchy_id")
        self.source_revision = _nonempty_text(
            provider.source_revision, "provider.source_revision"
        )
        self.coefficient_state_id = _nonempty_text(
            provider.coefficient_state_id, "provider.coefficient_state_id"
        )
        self._exact: dict[tuple[int, int], ExactHierarchyJetMajorant] = {}

    def certified_jet_bound(self, order: int, derivative_order: int) -> HierarchyJetBound:
        exact = self._provider.certified_exact_jet_majorant(order, derivative_order)
        if not isinstance(exact, ExactHierarchyJetMajorant):
            raise TypeError("provider must return ExactHierarchyJetMajorant values")
        if (exact.order, exact.derivative_order) != (order, derivative_order):
            raise ValueError("provider returned an exact C[j,m] majorant for the wrong index")
        if exact.hierarchy_id != self.hierarchy_id:
            raise ValueError("exact C[j,m] majorant uses a different hierarchy")
        if exact.source_revision != self.source_revision:
            raise ValueError("exact C[j,m] majorant uses a different source revision")
        if exact.coefficient_state_id != self.coefficient_state_id:
            raise ValueError("exact C[j,m] majorant uses a different coefficient state")

        key = (order, derivative_order)
        previous = self._exact.get(key)
        if previous is not None and previous != exact:
            raise ValueError("provider returned inconsistent exact majorants for one C[j,m] row")
        self._exact[key] = exact
        return exact.as_hierarchy_jet_bound()

    def certified_retained_identity(
        self, order: int
    ) -> ExactRetainedRecurrenceIdentity:
        return self._provider.certified_retained_identity(order)

    def certified_coefficient_support(
        self, order: int
    ) -> HierarchyCoefficientSupport:
        return self._provider.certified_coefficient_support(order)

    def exact_entries(self, max_order: int) -> tuple[ExactHierarchyJetMajorant, ...]:
        expected = [
            (j, m)
            for j in range(1, max_order + 1)
            for m in range(j + 3)
        ]
        if set(self._exact) != set(expected):
            missing = sorted(set(expected).difference(self._exact))
            extra = sorted(set(self._exact).difference(expected))
            raise RuntimeError(
                "target-driven exact-majorant adapter did not consume the complete "
                f"triangular prefix; missing={missing}, extra={extra}"
            )
        return tuple(self._exact[key] for key in expected)


@dataclass(frozen=True)
class TargetDrivenExactMajorantChainCertificate:
    """Finite target chain whose scheduler inputs are outward-safe exact majorants."""

    chain: TargetDrivenHierarchyDependencyChainCertificate
    exact_majorants: tuple[ExactHierarchyJetMajorant, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.chain, TargetDrivenHierarchyDependencyChainCertificate):
            raise TypeError(
                "chain must be a TargetDrivenHierarchyDependencyChainCertificate"
            )
        exact_majorants = tuple(self.exact_majorants)
        if not all(
            isinstance(entry, ExactHierarchyJetMajorant) for entry in exact_majorants
        ):
            raise TypeError("exact_majorants must contain ExactHierarchyJetMajorant values")

        expected = [
            (j, m)
            for j in range(1, self.chain.selected_order + 1)
            for m in range(j + 3)
        ]
        actual = [(entry.order, entry.derivative_order) for entry in exact_majorants]
        if actual != expected:
            raise ValueError("exact majorants must be the complete ordered triangular prefix")

        prefix = self.chain.supported_prefix.prefix
        for entry in exact_majorants:
            if entry.hierarchy_id != self.chain.hierarchy_id:
                raise ValueError("exact majorant hierarchy does not match the dependency chain")
            if entry.source_revision != self.chain.source_revision:
                raise ValueError(
                    "exact majorant source revision does not match the dependency chain"
                )
            if entry.coefficient_state_id != self.chain.coefficient_state_id:
                raise ValueError(
                    "exact majorant coefficient state does not match the dependency chain"
                )
            runtime = prefix.cutoff.bounds.entry(
                entry.order, entry.derivative_order
            )
            if runtime.dependencies != entry.dependencies:
                raise ValueError(
                    "outward scheduler row dependencies differ from the exact majorant"
                )
            if Fraction.from_float(runtime.bound) < entry.bound:
                raise ValueError("binary64 scheduler row rounded the exact majorant downward")

        object.__setattr__(self, "exact_majorants", exact_majorants)

    @property
    def selected_order(self) -> int:
        return self.chain.selected_order

    @property
    def outward_binary64_majorants_verified(self) -> bool:
        return True

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


def certify_target_driven_exact_majorant_chain(
    h: float,
    derivative_order: int,
    target_power: Fraction | int,
    provider: TargetDrivenExactMajorantProvider,
    *,
    minimum_order: int = 0,
    initial_lower_bound: int = 0,
) -> TargetDrivenExactMajorantChainCertificate:
    """Run the existing finite target chain from provider-owned exact C majorants."""

    if not isinstance(provider, TargetDrivenExactMajorantProvider):
        raise TypeError("provider must satisfy TargetDrivenExactMajorantProvider")
    adapter = _OutwardRoundedHierarchyProvider(provider)
    chain = certify_target_driven_dependency_chain(
        h,
        derivative_order,
        target_power,
        adapter,
        minimum_order=minimum_order,
        initial_lower_bound=initial_lower_bound,
    )
    return TargetDrivenExactMajorantChainCertificate(
        chain=chain,
        exact_majorants=adapter.exact_entries(chain.selected_order),
    )
