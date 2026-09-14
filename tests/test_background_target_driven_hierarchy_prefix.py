from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_hierarchy_cutoff_bounds import (
    HierarchyBoundDependency,
    HierarchyJetBound,
)
from openai_ns_reconstruction.background_recurrence_cancellation import (
    ExactRetainedRecurrenceIdentity,
    FormalAtomTerm,
    RetainedRecurrenceDependency,
)
from openai_ns_reconstruction.background_target_driven_hierarchy_prefix import (
    TargetDrivenHierarchyPrefixCertificate,
    certify_target_driven_hierarchy_prefix,
)


class _OwnedHierarchyProvider:
    hierarchy_id = "section5-owned-hierarchy"
    source_revision = "hierarchy-revision-301"
    coefficient_state_id = "coefficient-state-301"

    def __init__(self, *, max_order: int | None = None) -> None:
        self.max_order = max_order
        self.bound_calls: list[tuple[int, int]] = []
        self.identity_calls: list[int] = []

    def certified_jet_bound(self, order: int, derivative_order: int) -> HierarchyJetBound:
        self.bound_calls.append((order, derivative_order))
        if self.max_order is not None and order > self.max_order:
            raise RuntimeError(f"hierarchy frontier reached at coefficient order {order}")
        dependencies = (
            HierarchyBoundDependency(
                coefficient_order=order,
                artifact=f"owned-repaired-coefficient-{order}",
                provider=f"hierarchy.coefficient[{order}]",
            ),
        )
        if order > 1:
            dependencies += (
                HierarchyBoundDependency(
                    coefficient_order=order - 1,
                    artifact=f"owned-repaired-coefficient-{order-1}",
                    provider=f"hierarchy.coefficient[{order-1}]",
                ),
            )
        return HierarchyJetBound(
            order=order,
            derivative_order=derivative_order,
            bound=1.0 + order + derivative_order / 16.0,
            hierarchy_id=self.hierarchy_id,
            source_revision=self.source_revision,
            dependencies=dependencies,
        )

    def certified_retained_identity(self, order: int) -> ExactRetainedRecurrenceIdentity:
        self.identity_calls.append(order)
        if self.max_order is not None and order > self.max_order:
            raise RuntimeError(f"hierarchy frontier reached at recurrence order {order}")
        atom = f"retained-row-{order}"
        return ExactRetainedRecurrenceIdentity(
            order=order,
            hierarchy_id=self.hierarchy_id,
            source_revision=self.source_revision,
            coefficient_state_id=self.coefficient_state_id,
            identity_source=f"hierarchy.retained_identity[{order}]",
            theorem_name="SlowExpansionResidual.retained_recurrence",
            dependencies=(
                RetainedRecurrenceDependency(
                    coefficient_order=order,
                    artifact=f"owned-repaired-coefficient-{order}",
                    provider=f"hierarchy.coefficient[{order}]",
                ),
            ),
            linear_terms=(FormalAtomTerm(atom, 1),),
            pair_terms=(),
            previous_shifted_terms=(FormalAtomTerm(atom, 1),),
        )


def test_target_selects_exact_needed_prefix_then_queries_only_owned_provider() -> None:
    provider = _OwnedHierarchyProvider()
    cert = certify_target_driven_hierarchy_prefix(
        0.01,
        derivative_order=0,
        target_power=Fraction(1, 50),
        provider=provider,
        initial_lower_bound=7,
    )

    # With the exact runtime binary64 h, J=1 is the minimal positive hierarchy
    # order attaining h*(J+1) >= 1/50.
    assert cert.selected_order == 1
    assert cert.ordinary_tail_power_exact >= Fraction(1, 50)
    assert cert.derivative_budget_margin == 4
    assert cert.hierarchy_id == provider.hierarchy_id
    assert cert.source_revision == provider.source_revision
    assert cert.coefficient_state_id == provider.coefficient_state_id

    assert provider.bound_calls == [(1, 0), (1, 1), (1, 2), (1, 3)]
    assert provider.identity_calls == [0, 1]
    assert cert.cutoff.bounds.rows() == cert.cutoff.schedule.jet_bounds
    assert cert.cutoff.schedule.scales[0] >= 7
    assert cert.cancellations.max_order == cert.selected_order
    assert all(identity.exact_zero for identity in cert.cancellations.identities)

    assert cert.finite_target_prefix_only
    assert not cert.all_order_hierarchy_verified
    assert not cert.infinite_diagonal_schedule_verified
    assert not cert.pde_residual_tail_verified
    assert not cert.all_jets_flat
    assert not cert.super_algebraic
    assert not cert.paper_exact


def test_provider_frontier_fails_before_any_partial_target_certificate_is_returned() -> None:
    provider = _OwnedHierarchyProvider(max_order=1)

    # This target requires J=2 for h=0.01.  The finite real-hierarchy provider
    # deliberately stops at order one, so the target-driven construction must
    # expose that exact upstream frontier rather than inventing C[2,m].
    with pytest.raises(RuntimeError, match="frontier reached at coefficient order 2"):
        certify_target_driven_hierarchy_prefix(
            0.01,
            derivative_order=0,
            target_power=Fraction(3, 100),
            provider=provider,
        )

    assert (2, 0) in provider.bound_calls
    assert provider.identity_calls == []


def test_cross_wired_recurrence_state_fails_closed() -> None:
    class WrongState(_OwnedHierarchyProvider):
        def certified_retained_identity(self, order: int) -> ExactRetainedRecurrenceIdentity:
            identity = super().certified_retained_identity(order)
            if order != 1:
                return identity
            atom = "wrong-state-row"
            return ExactRetainedRecurrenceIdentity(
                order=1,
                hierarchy_id=identity.hierarchy_id,
                source_revision=identity.source_revision,
                coefficient_state_id="different-coefficient-state",
                identity_source=identity.identity_source,
                theorem_name=identity.theorem_name,
                dependencies=identity.dependencies,
                linear_terms=(FormalAtomTerm(atom, 1),),
                pair_terms=(),
                previous_shifted_terms=(FormalAtomTerm(atom, 1),),
            )

    with pytest.raises(ValueError, match="different coefficient state"):
        certify_target_driven_hierarchy_prefix(
            0.01,
            derivative_order=0,
            target_power=Fraction(1, 50),
            provider=WrongState(),
        )


def test_certificate_rejects_derivative_loss_cross_wiring() -> None:
    provider = _OwnedHierarchyProvider()
    good = certify_target_driven_hierarchy_prefix(
        0.01,
        derivative_order=0,
        target_power=Fraction(1, 50),
        provider=provider,
    )

    with pytest.raises(ValueError, match="derivative loss"):
        TargetDrivenHierarchyPrefixCertificate(
            derivative_order=1,
            coefficient_state_id=good.coefficient_state_id,
            target=good.target,
            cutoff=good.cutoff,
            cancellations=good.cancellations,
        )


def test_float_target_is_rejected_and_cannot_be_rounded_into_theorem_gate() -> None:
    provider = _OwnedHierarchyProvider()
    with pytest.raises(TypeError, match="floats are forbidden"):
        certify_target_driven_hierarchy_prefix(
            0.01,
            derivative_order=0,
            target_power=0.02,
            provider=provider,
        )


def test_wrong_provider_shape_is_rejected_before_querying_rows() -> None:
    class MissingRecurrenceProvider:
        hierarchy_id = "h"
        source_revision = "r"
        coefficient_state_id = "s"

        def certified_jet_bound(self, order: int, derivative_order: int) -> HierarchyJetBound:
            raise AssertionError("must not be called")

    with pytest.raises(TypeError, match="TargetDrivenHierarchyProvider"):
        certify_target_driven_hierarchy_prefix(
            0.01,
            derivative_order=0,
            target_power=Fraction(1, 50),
            provider=MissingRecurrenceProvider(),
        )
