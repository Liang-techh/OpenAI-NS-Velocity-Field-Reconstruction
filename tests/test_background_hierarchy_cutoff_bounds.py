import math

import pytest

from openai_ns_reconstruction.background_hierarchy_cutoff_bounds import (
    FiniteHierarchyBoundCertificate,
    HierarchyBoundDependency,
    HierarchyJetBound,
    build_slow_borel_schedule_from_hierarchy_provider,
    certify_hierarchy_bound_prefix,
)


class _Provider:
    hierarchy_id = "section5-real-hierarchy-fixture"
    source_revision = "fixture-revision-1"

    def __init__(self) -> None:
        self.calls: list[tuple[int, int]] = []

    def certified_jet_bound(self, order: int, derivative_order: int) -> HierarchyJetBound:
        self.calls.append((order, derivative_order))
        dependencies = (
            HierarchyBoundDependency(
                coefficient_order=order,
                artifact=f"repaired-coefficient-{order}",
                provider=f"hierarchy.coefficient[{order}]",
            ),
        )
        if order > 1:
            dependencies += (
                HierarchyBoundDependency(
                    coefficient_order=order - 1,
                    artifact=f"repaired-coefficient-{order-1}",
                    provider=f"hierarchy.coefficient[{order-1}]",
                ),
            )
        return HierarchyJetBound(
            order=order,
            derivative_order=derivative_order,
            bound=1.0 + order + derivative_order / 10.0,
            hierarchy_id=self.hierarchy_id,
            source_revision=self.source_revision,
            dependencies=dependencies,
        )


def test_hierarchy_provider_binds_exact_triangular_prefix_into_cutoff_schedule() -> None:
    provider = _Provider()
    witness = build_slow_borel_schedule_from_hierarchy_provider(
        0.2,
        3,
        provider,
        initial_lower_bound=5,
    )

    assert provider.calls == [
        (1, 0), (1, 1), (1, 2), (1, 3),
        (2, 0), (2, 1), (2, 2), (2, 3), (2, 4),
        (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5),
    ]
    assert witness.bounds.rows() == witness.schedule.jet_bounds
    assert witness.bounds.hierarchy_id == provider.hierarchy_id
    assert witness.bounds.source_revision == provider.source_revision
    assert witness.paper_exact is False
    assert witness.bounds.paper_exact is False
    assert witness.schedule.scales[0] >= 5

    for j in range(1, 4):
        assert witness.schedule.scales[j] >= 2 * witness.schedule.scales[j - 1]
        for m in range(j + 3):
            entry = witness.provenance(j, m)
            assert entry.order == j
            assert entry.derivative_order == m
            assert any(dep.coefficient_order == j for dep in entry.dependencies)
            assert witness.schedule.edge_log_margin(j, m) >= -64.0 * math.ulp(1.0)


def test_certificate_rejects_missing_or_duplicate_triangular_entries() -> None:
    provider = _Provider()
    certificate = certify_hierarchy_bound_prefix(provider, 2)

    with pytest.raises(ValueError, match="complete triangular prefix"):
        FiniteHierarchyBoundCertificate(
            hierarchy_id=certificate.hierarchy_id,
            source_revision=certificate.source_revision,
            max_order=2,
            entries=certificate.entries[:-1],
        )

    with pytest.raises(ValueError, match="duplicate"):
        FiniteHierarchyBoundCertificate(
            hierarchy_id=certificate.hierarchy_id,
            source_revision=certificate.source_revision,
            max_order=2,
            entries=certificate.entries[:-1] + (certificate.entries[0],),
        )


def test_bound_rejects_future_or_unbound_own_order_dependencies() -> None:
    with pytest.raises(ValueError, match="future coefficient"):
        HierarchyJetBound(
            order=2,
            derivative_order=0,
            bound=3.0,
            hierarchy_id="h",
            source_revision="r",
            dependencies=(HierarchyBoundDependency(3, "future", "provider"),),
        )

    with pytest.raises(ValueError, match="own order"):
        HierarchyJetBound(
            order=2,
            derivative_order=0,
            bound=3.0,
            hierarchy_id="h",
            source_revision="r",
            dependencies=(HierarchyBoundDependency(1, "lower", "provider"),),
        )


@pytest.mark.parametrize("bad_bound", [0.0, -1.0, math.inf, math.nan])
def test_bound_rejects_nonpositive_or_nonfinite_values(bad_bound: float) -> None:
    with pytest.raises(ValueError, match="finite and positive"):
        HierarchyJetBound(
            order=1,
            derivative_order=0,
            bound=bad_bound,
            hierarchy_id="h",
            source_revision="r",
            dependencies=(HierarchyBoundDependency(1, "c1", "provider"),),
        )


def test_provider_identity_and_revision_mismatch_fail_closed() -> None:
    class WrongHierarchy(_Provider):
        def certified_jet_bound(self, order: int, derivative_order: int) -> HierarchyJetBound:
            entry = super().certified_jet_bound(order, derivative_order)
            if (order, derivative_order) == (1, 2):
                return HierarchyJetBound(
                    order=entry.order,
                    derivative_order=entry.derivative_order,
                    bound=entry.bound,
                    hierarchy_id="different-hierarchy",
                    source_revision=entry.source_revision,
                    dependencies=entry.dependencies,
                )
            return entry

    with pytest.raises(ValueError, match="different hierarchy"):
        certify_hierarchy_bound_prefix(WrongHierarchy(), 1)

    class WrongRevision(_Provider):
        def certified_jet_bound(self, order: int, derivative_order: int) -> HierarchyJetBound:
            entry = super().certified_jet_bound(order, derivative_order)
            if (order, derivative_order) == (1, 2):
                return HierarchyJetBound(
                    order=entry.order,
                    derivative_order=entry.derivative_order,
                    bound=entry.bound,
                    hierarchy_id=entry.hierarchy_id,
                    source_revision="different-revision",
                    dependencies=entry.dependencies,
                )
            return entry

    with pytest.raises(ValueError, match="different source revision"):
        certify_hierarchy_bound_prefix(WrongRevision(), 1)


def test_provider_wrong_index_fails_closed_before_schedule_construction() -> None:
    class WrongIndex(_Provider):
        def certified_jet_bound(self, order: int, derivative_order: int) -> HierarchyJetBound:
            entry = super().certified_jet_bound(order, derivative_order)
            if (order, derivative_order) == (1, 1):
                return HierarchyJetBound(
                    order=1,
                    derivative_order=0,
                    bound=entry.bound,
                    hierarchy_id=entry.hierarchy_id,
                    source_revision=entry.source_revision,
                    dependencies=entry.dependencies,
                )
            return entry

    with pytest.raises(ValueError, match="wrong index"):
        certify_hierarchy_bound_prefix(WrongIndex(), 1)
