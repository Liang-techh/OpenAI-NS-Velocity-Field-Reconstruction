"""Target-driven Section 5 prefix with hierarchy-owned common support evidence.

The pinned all-order construction does not use derivative bounds in isolation.
The actual positive-order coefficient fields are smooth *and* have a common
compact support box before the SlowBorel/DiagonalScale schedule is chosen.  In
the official construction this is provided by
``AssembledSlowBase.extendedCoefficient_support`` / ``...compactSupport``:
for every positive order and all four coefficient components,

    tsupport(coeff[n,i]) <= [-1, B^2/2] x [-outer, outer].

This module closes that support/provenance seam for the finite target-driven
prefix already constructed by ``background_target_driven_hierarchy_prefix``.
Callers cannot pass a support table.  One hierarchy provider must own each
support witness, and every support row must share an own-order coefficient
artifact with every C[j,m] row at that order.  The common box is represented
with exact ``Fraction`` arithmetic; binary floats are deliberately rejected.

This remains finite-prefix formal structure.  It neither proves that the
provider is total for arbitrary order nor promotes a finite supported prefix
to an infinite diagonal schedule, Proposition 5.3, all-jets flatness, or
super-algebraic residual convergence.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Protocol, runtime_checkable

from .background_target_driven_hierarchy_prefix import (
    TargetDrivenHierarchyPrefixCertificate,
    TargetDrivenHierarchyProvider,
    certify_target_driven_hierarchy_prefix,
)


def _positive_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _nonempty_text(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value.strip()


def _exact_fraction(value: Fraction | int, name: str) -> Fraction:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact integer/Fraction; booleans are forbidden")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value), 1)
    raise TypeError(f"{name} must be an exact integer/Fraction; floats are forbidden")


PINNED_COEFFICIENT_COMPONENTS = ("phi", "axial", "beta", "pressure")
PINNED_SUPPORT_THEOREM = "AssembledSlowBase.extendedCoefficient_support"
PINNED_COMPACT_SUPPORT_THEOREM = "AssembledSlowBase.extendedCoefficient_compactSupport"


@dataclass(frozen=True)
class HierarchySupportDependency:
    """One hierarchy-owned coefficient artifact used by a support proof."""

    coefficient_order: int
    artifact: str
    provider: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "coefficient_order", _positive_int(self.coefficient_order, "coefficient_order")
        )
        object.__setattr__(self, "artifact", _nonempty_text(self.artifact, "artifact"))
        object.__setattr__(self, "provider", _nonempty_text(self.provider, "provider"))

    @property
    def key(self) -> tuple[int, str, str]:
        return (self.coefficient_order, self.artifact, self.provider)


@dataclass(frozen=True)
class HierarchyCoefficientSupport:
    """Exact common-box support witness for one positive coefficient order.

    ``scheme_B`` and ``parameter_outer`` encode the two actual construction
    parameters appearing in the pinned support theorem.  The remaining box
    endpoints are checked rather than trusted: radial support must be exactly
    ``[-1, B^2/2]`` and the parameter support is exactly symmetric
    ``[-outer, outer]``.  ``outer > 1`` records that the common parameter
    window contains the physical ``[-1,1]`` strip used by the assembled base.
    """

    order: int
    hierarchy_id: str
    source_revision: str
    coefficient_state_id: str
    scheme_B: Fraction | int
    parameter_outer: Fraction | int
    radial_lower: Fraction | int
    radial_upper: Fraction | int
    components: tuple[str, ...]
    support_theorem: str
    compact_support_theorem: str
    dependencies: tuple[HierarchySupportDependency, ...]

    def __post_init__(self) -> None:
        order = _positive_int(self.order, "order")
        hierarchy_id = _nonempty_text(self.hierarchy_id, "hierarchy_id")
        source_revision = _nonempty_text(self.source_revision, "source_revision")
        coefficient_state_id = _nonempty_text(self.coefficient_state_id, "coefficient_state_id")
        scheme_B = _exact_fraction(self.scheme_B, "scheme_B")
        parameter_outer = _exact_fraction(self.parameter_outer, "parameter_outer")
        radial_lower = _exact_fraction(self.radial_lower, "radial_lower")
        radial_upper = _exact_fraction(self.radial_upper, "radial_upper")
        components = tuple(self.components)
        support_theorem = _nonempty_text(self.support_theorem, "support_theorem")
        compact_support_theorem = _nonempty_text(
            self.compact_support_theorem, "compact_support_theorem"
        )
        dependencies = tuple(self.dependencies)

        if scheme_B <= 0:
            raise ValueError("scheme_B must be positive")
        if parameter_outer <= 1:
            raise ValueError("parameter_outer must be strictly greater than one")
        if radial_lower != Fraction(-1, 1):
            raise ValueError("pinned radial support lower endpoint must equal -1 exactly")
        expected_upper = scheme_B * scheme_B / 2
        if radial_upper != expected_upper:
            raise ValueError("pinned radial support upper endpoint must equal B^2/2 exactly")
        if components != PINNED_COEFFICIENT_COMPONENTS:
            raise ValueError("support witness must cover exactly phi/axial/beta/pressure")
        if support_theorem != PINNED_SUPPORT_THEOREM:
            raise ValueError("support witness must cite the pinned extendedCoefficient_support theorem")
        if compact_support_theorem != PINNED_COMPACT_SUPPORT_THEOREM:
            raise ValueError(
                "support witness must cite the pinned extendedCoefficient_compactSupport theorem"
            )
        if not dependencies:
            raise ValueError("support witness must carry at least one hierarchy dependency")
        if not all(isinstance(dep, HierarchySupportDependency) for dep in dependencies):
            raise TypeError("dependencies must contain HierarchySupportDependency values")
        keys = [dep.key for dep in dependencies]
        if len(keys) != len(set(keys)):
            raise ValueError("support dependencies must be unique")
        if any(dep.coefficient_order > order for dep in dependencies):
            raise ValueError("support witness cannot depend on a future coefficient order")
        if not any(dep.coefficient_order == order for dep in dependencies):
            raise ValueError("support witness must bind the coefficient at its own order")

        object.__setattr__(self, "order", order)
        object.__setattr__(self, "hierarchy_id", hierarchy_id)
        object.__setattr__(self, "source_revision", source_revision)
        object.__setattr__(self, "coefficient_state_id", coefficient_state_id)
        object.__setattr__(self, "scheme_B", scheme_B)
        object.__setattr__(self, "parameter_outer", parameter_outer)
        object.__setattr__(self, "radial_lower", radial_lower)
        object.__setattr__(self, "radial_upper", radial_upper)
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "support_theorem", support_theorem)
        object.__setattr__(self, "compact_support_theorem", compact_support_theorem)
        object.__setattr__(self, "dependencies", dependencies)

    @property
    def exact_box(self) -> tuple[Fraction, Fraction, Fraction, Fraction]:
        """Return ``(X_lo, X_hi, eta_lo, eta_hi)`` exactly."""
        return (
            self.radial_lower,
            self.radial_upper,
            -self.parameter_outer,
            self.parameter_outer,
        )

    @property
    def own_order_dependency_keys(self) -> frozenset[tuple[int, str, str]]:
        return frozenset(dep.key for dep in self.dependencies if dep.coefficient_order == self.order)


@runtime_checkable
class TargetDrivenHierarchySupportProvider(TargetDrivenHierarchyProvider, Protocol):
    """Target-driven hierarchy provider that also owns coefficient support proofs."""

    def certified_coefficient_support(self, order: int) -> HierarchyCoefficientSupport:
        ...


@dataclass(frozen=True)
class FiniteHierarchySupportCertificate:
    """Common exact support box for every positive order in one finite prefix."""

    hierarchy_id: str
    source_revision: str
    coefficient_state_id: str
    max_order: int
    entries: tuple[HierarchyCoefficientSupport, ...]

    def __post_init__(self) -> None:
        hierarchy_id = _nonempty_text(self.hierarchy_id, "hierarchy_id")
        source_revision = _nonempty_text(self.source_revision, "source_revision")
        coefficient_state_id = _nonempty_text(self.coefficient_state_id, "coefficient_state_id")
        max_order = _positive_int(self.max_order, "max_order")
        entries = tuple(self.entries)
        if not all(isinstance(entry, HierarchyCoefficientSupport) for entry in entries):
            raise TypeError("entries must contain HierarchyCoefficientSupport values")
        if tuple(entry.order for entry in entries) != tuple(range(1, max_order + 1)):
            raise ValueError("support certificate must contain exactly the contiguous orders 1..J")

        reference = entries[0]
        for entry in entries:
            if entry.hierarchy_id != hierarchy_id:
                raise ValueError("support entry hierarchy_id does not match its certificate")
            if entry.source_revision != source_revision:
                raise ValueError("support entry source_revision does not match its certificate")
            if entry.coefficient_state_id != coefficient_state_id:
                raise ValueError("support entry coefficient_state_id does not match its certificate")
            if entry.exact_box != reference.exact_box:
                raise ValueError("positive-order coefficient supports do not share one common box")
            if entry.components != reference.components:
                raise ValueError("positive-order support witnesses cover different components")
            if entry.support_theorem != reference.support_theorem:
                raise ValueError("positive-order support witnesses cite different support theorems")
            if entry.compact_support_theorem != reference.compact_support_theorem:
                raise ValueError("positive-order support witnesses cite different compactness theorems")

        object.__setattr__(self, "hierarchy_id", hierarchy_id)
        object.__setattr__(self, "source_revision", source_revision)
        object.__setattr__(self, "coefficient_state_id", coefficient_state_id)
        object.__setattr__(self, "max_order", max_order)
        object.__setattr__(self, "entries", entries)

    @property
    def exact_box(self) -> tuple[Fraction, Fraction, Fraction, Fraction]:
        return self.entries[0].exact_box

    @property
    def all_positive_orders_share_common_support(self) -> bool:
        return True

    @property
    def paper_exact(self) -> bool:
        return False


@dataclass(frozen=True)
class TargetDrivenSupportedHierarchyPrefixCertificate:
    """Finite target prefix whose C rows and coefficients share owned support evidence."""

    prefix: TargetDrivenHierarchyPrefixCertificate
    supports: FiniteHierarchySupportCertificate

    def __post_init__(self) -> None:
        if not isinstance(self.prefix, TargetDrivenHierarchyPrefixCertificate):
            raise TypeError("prefix must be a TargetDrivenHierarchyPrefixCertificate")
        if not isinstance(self.supports, FiniteHierarchySupportCertificate):
            raise TypeError("supports must be a FiniteHierarchySupportCertificate")
        if self.supports.max_order != self.prefix.selected_order:
            raise ValueError("support prefix order does not match the target-driven prefix")
        if self.supports.hierarchy_id != self.prefix.hierarchy_id:
            raise ValueError("support prefix and derivative bounds use different hierarchies")
        if self.supports.source_revision != self.prefix.source_revision:
            raise ValueError("support prefix and derivative bounds use different source revisions")
        if self.supports.coefficient_state_id != self.prefix.coefficient_state_id:
            raise ValueError("support prefix and recurrence rows use different coefficient states")

        # Bind support and every C[j,m] row to at least one identical owned
        # coefficient artifact.  Sharing only a hierarchy label is not enough.
        for support in self.supports.entries:
            support_keys = support.own_order_dependency_keys
            for m in range(support.order + 3):
                bound = self.prefix.cutoff.bounds.entry(support.order, m)
                bound_keys = frozenset(
                    (dep.coefficient_order, dep.artifact, dep.provider)
                    for dep in bound.dependencies
                    if dep.coefficient_order == support.order
                )
                if support_keys.isdisjoint(bound_keys):
                    raise ValueError(
                        f"support order {support.order} does not share an own-order coefficient "
                        f"dependency with C[{support.order},{m}]"
                    )

    @property
    def selected_order(self) -> int:
        return self.prefix.selected_order

    @property
    def exact_support_box(self) -> tuple[Fraction, Fraction, Fraction, Fraction]:
        return self.supports.exact_box

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


def certify_target_driven_supported_hierarchy_prefix(
    h: float,
    derivative_order: int,
    target_power: Fraction | int,
    provider: TargetDrivenHierarchySupportProvider,
    *,
    minimum_order: int = 0,
    initial_lower_bound: int = 0,
) -> TargetDrivenSupportedHierarchyPrefixCertificate:
    """Build one finite target prefix and query its support only from ``provider``."""

    if not isinstance(provider, TargetDrivenHierarchySupportProvider):
        raise TypeError("provider must satisfy TargetDrivenHierarchySupportProvider")

    prefix = certify_target_driven_hierarchy_prefix(
        h,
        derivative_order,
        target_power,
        provider,
        minimum_order=minimum_order,
        initial_lower_bound=initial_lower_bound,
    )

    entries: list[HierarchyCoefficientSupport] = []
    for order in range(1, prefix.selected_order + 1):
        support = provider.certified_coefficient_support(order)
        if not isinstance(support, HierarchyCoefficientSupport):
            raise TypeError("provider must return HierarchyCoefficientSupport values")
        if support.order != order:
            raise ValueError("provider returned a support witness for the wrong coefficient order")
        if support.hierarchy_id != prefix.hierarchy_id:
            raise ValueError("provider returned a support witness for a different hierarchy")
        if support.source_revision != prefix.source_revision:
            raise ValueError("provider returned a support witness for a different source revision")
        if support.coefficient_state_id != prefix.coefficient_state_id:
            raise ValueError("provider returned a support witness for a different coefficient state")
        entries.append(support)

    supports = FiniteHierarchySupportCertificate(
        hierarchy_id=prefix.hierarchy_id,
        source_revision=prefix.source_revision,
        coefficient_state_id=prefix.coefficient_state_id,
        max_order=prefix.selected_order,
        entries=tuple(entries),
    )
    return TargetDrivenSupportedHierarchyPrefixCertificate(prefix=prefix, supports=supports)
