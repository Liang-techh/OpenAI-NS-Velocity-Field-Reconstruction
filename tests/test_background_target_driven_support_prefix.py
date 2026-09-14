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
from openai_ns_reconstruction.background_target_driven_support_prefix import (
    HierarchyCoefficientSupport,
    HierarchySupportDependency,
    PINNED_COEFFICIENT_COMPONENTS,
    PINNED_COMPACT_SUPPORT_THEOREM,
    PINNED_SUPPORT_THEOREM,
    TargetDrivenSupportedHierarchyPrefixCertificate,
    certify_target_driven_supported_hierarchy_prefix,
)


class _SupportedOwnedHierarchyProvider:
    hierarchy_id = "section5-owned-hierarchy"
    source_revision = "hierarchy-revision-301"
    coefficient_state_id = "coefficient-state-301"

    def __init__(
        self,
        *,
        max_order: int | None = None,
        support_max_order: int | None = None,
        support_B_by_order: dict[int, Fraction] | None = None,
    ) -> None:
        self.max_order = max_order
        self.support_max_order = support_max_order
        self.support_B_by_order = support_B_by_order or {}
        self.bound_calls: list[tuple[int, int]] = []
        self.identity_calls: list[int] = []
        self.support_calls: list[int] = []

    @staticmethod
    def _artifact(order: int) -> str:
        return f"owned-repaired-coefficient-{order}"

    @staticmethod
    def _provider(order: int) -> str:
        return f"hierarchy.coefficient[{order}]"

    def certified_jet_bound(self, order: int, derivative_order: int) -> HierarchyJetBound:
        self.bound_calls.append((order, derivative_order))
        if self.max_order is not None and order > self.max_order:
            raise RuntimeError(f"hierarchy frontier reached at coefficient order {order}")
        return HierarchyJetBound(
            order=order,
            derivative_order=derivative_order,
            bound=1.0 + order + derivative_order / 16.0,
            hierarchy_id=self.hierarchy_id,
            source_revision=self.source_revision,
            dependencies=(
                HierarchyBoundDependency(
                    coefficient_order=order,
                    artifact=self._artifact(order),
                    provider=self._provider(order),
                ),
            ),
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

    def certified_coefficient_support(self, order: int) -> HierarchyCoefficientSupport:
        self.support_calls.append(order)
        if self.support_max_order is not None and order > self.support_max_order:
            raise RuntimeError(f"support frontier reached at coefficient order {order}")
        B = self.support_B_by_order.get(order, Fraction(4, 1))
        return HierarchyCoefficientSupport(
            order=order,
            hierarchy_id=self.hierarchy_id,
            source_revision=self.source_revision,
            coefficient_state_id=self.coefficient_state_id,
            scheme_B=B,
            parameter_outer=Fraction(3, 1),
            radial_lower=Fraction(-1, 1),
            radial_upper=B * B / 2,
            components=PINNED_COEFFICIENT_COMPONENTS,
            support_theorem=PINNED_SUPPORT_THEOREM,
            compact_support_theorem=PINNED_COMPACT_SUPPORT_THEOREM,
            dependencies=(
                HierarchySupportDependency(
                    coefficient_order=order,
                    artifact=self._artifact(order),
                    provider=self._provider(order),
                ),
            ),
        )


def test_target_prefix_binds_every_C_row_to_one_exact_common_support_box() -> None:
    provider = _SupportedOwnedHierarchyProvider()
    cert = certify_target_driven_supported_hierarchy_prefix(
        0.01,
        derivative_order=0,
        target_power=Fraction(1, 50),
        provider=provider,
        initial_lower_bound=7,
    )

    assert cert.selected_order == 1
    assert cert.exact_support_box == (
        Fraction(-1, 1),
        Fraction(8, 1),
        Fraction(-3, 1),
        Fraction(3, 1),
    )
    assert cert.supports.all_positive_orders_share_common_support
    assert provider.support_calls == [1]
    assert provider.bound_calls == [(1, 0), (1, 1), (1, 2), (1, 3)]
    assert provider.identity_calls == [0, 1]

    assert cert.finite_target_prefix_only
    assert not cert.all_order_hierarchy_verified
    assert not cert.infinite_diagonal_schedule_verified
    assert not cert.pde_residual_tail_verified
    assert not cert.all_jets_flat
    assert not cert.super_algebraic
    assert not cert.paper_exact


def test_support_frontier_fails_closed_at_the_exact_requested_order() -> None:
    provider = _SupportedOwnedHierarchyProvider(support_max_order=1)

    # This target requires J=2.  Derivative/cancellation data exist through J,
    # but support ownership deliberately stops at one; no supported target
    # certificate may be returned.
    with pytest.raises(RuntimeError, match="support frontier reached at coefficient order 2"):
        certify_target_driven_supported_hierarchy_prefix(
            0.01,
            derivative_order=0,
            target_power=Fraction(3, 100),
            provider=provider,
        )

    assert provider.support_calls == [1, 2]


def test_support_and_C_rows_must_share_the_same_owned_coefficient_artifact() -> None:
    class CrossWiredSupport(_SupportedOwnedHierarchyProvider):
        def certified_coefficient_support(self, order: int) -> HierarchyCoefficientSupport:
            good = super().certified_coefficient_support(order)
            return HierarchyCoefficientSupport(
                order=good.order,
                hierarchy_id=good.hierarchy_id,
                source_revision=good.source_revision,
                coefficient_state_id=good.coefficient_state_id,
                scheme_B=good.scheme_B,
                parameter_outer=good.parameter_outer,
                radial_lower=good.radial_lower,
                radial_upper=good.radial_upper,
                components=good.components,
                support_theorem=good.support_theorem,
                compact_support_theorem=good.compact_support_theorem,
                dependencies=(
                    HierarchySupportDependency(
                        coefficient_order=order,
                        artifact=f"different-coefficient-{order}",
                        provider=f"different.provider[{order}]",
                    ),
                ),
            )

    with pytest.raises(ValueError, match="does not share an own-order coefficient dependency"):
        certify_target_driven_supported_hierarchy_prefix(
            0.01,
            derivative_order=0,
            target_power=Fraction(1, 50),
            provider=CrossWiredSupport(),
        )


def test_positive_orders_must_use_one_common_support_box() -> None:
    provider = _SupportedOwnedHierarchyProvider(
        support_B_by_order={1: Fraction(4, 1), 2: Fraction(5, 1)}
    )
    with pytest.raises(ValueError, match="do not share one common box"):
        certify_target_driven_supported_hierarchy_prefix(
            0.01,
            derivative_order=0,
            target_power=Fraction(3, 100),
            provider=provider,
        )


def test_binary_float_support_endpoints_are_rejected() -> None:
    with pytest.raises(TypeError, match="floats are forbidden"):
        HierarchyCoefficientSupport(
            order=1,
            hierarchy_id="h",
            source_revision="r",
            coefficient_state_id="s",
            scheme_B=4.0,
            parameter_outer=Fraction(3, 1),
            radial_lower=Fraction(-1, 1),
            radial_upper=Fraction(8, 1),
            components=PINNED_COEFFICIENT_COMPONENTS,
            support_theorem=PINNED_SUPPORT_THEOREM,
            compact_support_theorem=PINNED_COMPACT_SUPPORT_THEOREM,
            dependencies=(
                HierarchySupportDependency(1, "coefficient-1", "hierarchy.coefficient[1]"),
            ),
        )


def test_support_witness_rejects_nonpaper_box_or_theorem() -> None:
    base = dict(
        order=1,
        hierarchy_id="h",
        source_revision="r",
        coefficient_state_id="s",
        scheme_B=Fraction(4, 1),
        parameter_outer=Fraction(3, 1),
        radial_lower=Fraction(-1, 1),
        radial_upper=Fraction(8, 1),
        components=PINNED_COEFFICIENT_COMPONENTS,
        support_theorem=PINNED_SUPPORT_THEOREM,
        compact_support_theorem=PINNED_COMPACT_SUPPORT_THEOREM,
        dependencies=(HierarchySupportDependency(1, "coefficient-1", "provider-1"),),
    )

    with pytest.raises(ValueError, match=r"B\^2/2"):
        HierarchyCoefficientSupport(**{**base, "radial_upper": Fraction(9, 1)})
    with pytest.raises(ValueError, match="pinned extendedCoefficient_support"):
        HierarchyCoefficientSupport(**{**base, "support_theorem": "generic.support"})


def test_certificate_rejects_support_state_cross_wiring() -> None:
    provider = _SupportedOwnedHierarchyProvider()
    good = certify_target_driven_supported_hierarchy_prefix(
        0.01,
        derivative_order=0,
        target_power=Fraction(1, 50),
        provider=provider,
    )
    support = good.supports.entries[0]
    wrong = HierarchyCoefficientSupport(
        order=support.order,
        hierarchy_id=support.hierarchy_id,
        source_revision=support.source_revision,
        coefficient_state_id="different-state",
        scheme_B=support.scheme_B,
        parameter_outer=support.parameter_outer,
        radial_lower=support.radial_lower,
        radial_upper=support.radial_upper,
        components=support.components,
        support_theorem=support.support_theorem,
        compact_support_theorem=support.compact_support_theorem,
        dependencies=support.dependencies,
    )

    with pytest.raises(ValueError, match="different coefficient states"):
        TargetDrivenSupportedHierarchyPrefixCertificate(
            prefix=good.prefix,
            supports=type(good.supports)(
                hierarchy_id=good.supports.hierarchy_id,
                source_revision=good.supports.source_revision,
                coefficient_state_id="different-state",
                max_order=1,
                entries=(wrong,),
            ),
        )
