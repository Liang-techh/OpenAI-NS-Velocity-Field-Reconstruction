"""Bind target-driven Section 5 bounds, support, and recurrences to one coefficient artifact.

The target-driven hierarchy prefix already requires hierarchy-owned ``C[j,m]``
rows and exact retained recurrence identities, while the supported-prefix layer
requires each positive-order support witness to share an own-order coefficient
artifact with every derivative-bound row.  One provenance seam remained: an
exact retained recurrence identity could still name a different coefficient
artifact while reusing the same hierarchy/source/coefficient-state labels.

This module closes only that seam.  For every positive order in one finite,
target-selected prefix it requires a *single literal* ``(order, artifact,
provider)`` dependency to occur simultaneously in

* the compact-support witness,
* the exact retained recurrence identity, and
* every ``C[j,m]`` derivative-bound row for ``0 <= m <= j+2``.

The production constructor obtains the prefix through the existing
``TargetDrivenHierarchySupportProvider`` API; callers cannot supply a separate
``C[j,m]`` table, support table, or recurrence list.  Exact recurrence zero is
still enforced upstream with rational formal algebra, so no numerical
near-zero test is introduced here.

This remains a finite-prefix provenance certificate.  It does not prove that
the provider is total for arbitrary order, does not construct an infinite
diagonal schedule, and does not establish the actual PDE residual tail,
all-jets flatness, or super-algebraic convergence.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .background_target_driven_support_prefix import (
    TargetDrivenHierarchySupportProvider,
    TargetDrivenSupportedHierarchyPrefixCertificate,
    certify_target_driven_supported_hierarchy_prefix,
)


def _positive_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _nonempty_text(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value.strip()


def _dependency_key(dependency: object) -> tuple[int, str, str]:
    try:
        order = getattr(dependency, "coefficient_order")
        artifact = getattr(dependency, "artifact")
        provider = getattr(dependency, "provider")
    except AttributeError as exc:  # pragma: no cover - defensive contract guard
        raise TypeError("hierarchy dependencies must expose coefficient_order/artifact/provider") from exc
    return (
        _positive_int(order, "dependency.coefficient_order"),
        _nonempty_text(artifact, "dependency.artifact"),
        _nonempty_text(provider, "dependency.provider"),
    )


@dataclass(frozen=True)
class PositiveOrderCoefficientEvidenceLink:
    """One literal coefficient dependency shared by every evidence layer at an order."""

    order: int
    artifact: str
    provider: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "order", _positive_int(self.order, "order"))
        object.__setattr__(self, "artifact", _nonempty_text(self.artifact, "artifact"))
        object.__setattr__(self, "provider", _nonempty_text(self.provider, "provider"))

    @property
    def key(self) -> tuple[int, str, str]:
        return (self.order, self.artifact, self.provider)


@dataclass(frozen=True)
class TargetDrivenHierarchyDependencyChainCertificate:
    """Finite target prefix with one common owned coefficient artifact per order."""

    supported_prefix: TargetDrivenSupportedHierarchyPrefixCertificate
    links: tuple[PositiveOrderCoefficientEvidenceLink, ...]

    def __post_init__(self) -> None:
        if not isinstance(
            self.supported_prefix, TargetDrivenSupportedHierarchyPrefixCertificate
        ):
            raise TypeError(
                "supported_prefix must be a TargetDrivenSupportedHierarchyPrefixCertificate"
            )
        links = tuple(self.links)
        if not all(isinstance(link, PositiveOrderCoefficientEvidenceLink) for link in links):
            raise TypeError("links must contain PositiveOrderCoefficientEvidenceLink values")
        expected_orders = tuple(range(1, self.supported_prefix.selected_order + 1))
        if tuple(link.order for link in links) != expected_orders:
            raise ValueError("dependency links must cover exactly the contiguous positive orders 1..J")

        prefix = self.supported_prefix.prefix
        for link in links:
            support = self.supported_prefix.supports.entries[link.order - 1]
            identity = prefix.cancellations.identity(link.order)

            support_keys = support.own_order_dependency_keys
            recurrence_keys = frozenset(
                _dependency_key(dep)
                for dep in identity.dependencies
                if dep.coefficient_order == link.order
            )
            bound_key_sets = []
            for derivative_order in range(link.order + 3):
                bound = prefix.cutoff.bounds.entry(link.order, derivative_order)
                bound_key_sets.append(
                    frozenset(
                        _dependency_key(dep)
                        for dep in bound.dependencies
                        if dep.coefficient_order == link.order
                    )
                )

            if link.key not in support_keys:
                raise ValueError(
                    f"dependency link for order {link.order} is absent from the support witness"
                )
            if link.key not in recurrence_keys:
                raise ValueError(
                    f"dependency link for order {link.order} is absent from the retained recurrence"
                )
            if any(link.key not in keys for keys in bound_key_sets):
                raise ValueError(
                    f"dependency link for order {link.order} is absent from at least one C[j,m] row"
                )
            if not identity.exact_zero:
                raise ValueError("retained recurrence identity must be exact zero")

        object.__setattr__(self, "links", links)

    @property
    def selected_order(self) -> int:
        return self.supported_prefix.selected_order

    @property
    def hierarchy_id(self) -> str:
        return self.supported_prefix.prefix.hierarchy_id

    @property
    def source_revision(self) -> str:
        return self.supported_prefix.prefix.source_revision

    @property
    def coefficient_state_id(self) -> str:
        return self.supported_prefix.prefix.coefficient_state_id

    @property
    def retained_recurrences_exact(self) -> bool:
        return True

    @property
    def one_owned_artifact_per_positive_order(self) -> bool:
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


def bind_target_driven_dependency_chain(
    supported_prefix: TargetDrivenSupportedHierarchyPrefixCertificate,
) -> TargetDrivenHierarchyDependencyChainCertificate:
    """Bind an already provider-built supported prefix across all evidence layers."""

    if not isinstance(supported_prefix, TargetDrivenSupportedHierarchyPrefixCertificate):
        raise TypeError(
            "supported_prefix must be a TargetDrivenSupportedHierarchyPrefixCertificate"
        )

    prefix = supported_prefix.prefix
    links: list[PositiveOrderCoefficientEvidenceLink] = []
    for order in range(1, supported_prefix.selected_order + 1):
        support = supported_prefix.supports.entries[order - 1]
        identity = prefix.cancellations.identity(order)
        common_keys = set(support.own_order_dependency_keys)
        common_keys.intersection_update(
            _dependency_key(dep)
            for dep in identity.dependencies
            if dep.coefficient_order == order
        )
        for derivative_order in range(order + 3):
            bound = prefix.cutoff.bounds.entry(order, derivative_order)
            common_keys.intersection_update(
                _dependency_key(dep)
                for dep in bound.dependencies
                if dep.coefficient_order == order
            )
        if not common_keys:
            raise ValueError(
                f"order {order} has no single owned coefficient artifact shared by support, "
                "retained recurrence, and every C[j,m] row"
            )
        coefficient_order, artifact, provider = sorted(common_keys)[0]
        links.append(
            PositiveOrderCoefficientEvidenceLink(
                order=coefficient_order,
                artifact=artifact,
                provider=provider,
            )
        )

    return TargetDrivenHierarchyDependencyChainCertificate(
        supported_prefix=supported_prefix,
        links=tuple(links),
    )


def certify_target_driven_dependency_chain(
    h: float,
    derivative_order: int,
    target_power,
    provider: TargetDrivenHierarchySupportProvider,
    *,
    minimum_order: int = 0,
    initial_lower_bound: int = 0,
) -> TargetDrivenHierarchyDependencyChainCertificate:
    """Build one target-driven prefix and bind all positive-order evidence to it.

    ``target_power`` intentionally keeps the exact-input contract of the
    underlying target selector: integers/Fractions are accepted and binary
    floats are rejected there.  Missing hierarchy orders/support/recurrences
    propagate as failures; no partial certificate is returned.
    """

    supported = certify_target_driven_supported_hierarchy_prefix(
        h,
        derivative_order,
        target_power,
        provider,
        minimum_order=minimum_order,
        initial_lower_bound=initial_lower_bound,
    )
    return bind_target_driven_dependency_chain(supported)
