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
from openai_ns_reconstruction.background_target_driven_dependency_chain import (
    bind_target_driven_dependency_chain,
    certify_target_driven_dependency_chain,
)
from openai_ns_reconstruction.background_target_driven_support_prefix import (
    HierarchyCoefficientSupport,
    HierarchySupportDependency,
    PINNED_COEFFICIENT_COMPONENTS,
    PINNED_COMPACT_SUPPORT_THEOREM,
    PINNED_SUPPORT_THEOREM,
    certify_target_driven_supported_hierarchy_prefix,
)


class _Provider:
    hierarchy_id = "section5-owned-hierarchy"
    source_revision = "hierarchy-revision-agent2"
    coefficient_state_id = "coefficient-state-agent2"

    @staticmethod
    def _artifact(order: int) -> str:
        return f"assembled-coefficient-{order}"

    @staticmethod
    def _provider(order: int) -> str:
        return f"hierarchy.assembled[{order}]"

    def certified_jet_bound(self, order: int, derivative_order: int) -> HierarchyJetBound:
        return HierarchyJetBound(
            order=order,
            derivative_order=derivative_order,
            bound=2.0 + order + derivative_order / 16.0,
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
        atom = f"retained-row-{order}"
        return ExactRetainedRecurrenceIdentity(
            order=order,
            hierarchy_id=self.hierarchy_id,
            source_revision=self.source_revision,
            coefficient_state_id=self.coefficient_state_id,
            identity_source=f"hierarchy.retained[{order}]",
            theorem_name="SlowExpansionResidual.retained_recurrence",
            dependencies=(
                RetainedRecurrenceDependency(
                    coefficient_order=order,
                    artifact=self._artifact(order),
                    provider=self._provider(order),
                ),
            ),
            linear_terms=(FormalAtomTerm(atom, 1),),
            pair_terms=(),
            previous_shifted_terms=(FormalAtomTerm(atom, 1),),
        )

    def certified_coefficient_support(self, order: int) -> HierarchyCoefficientSupport:
        B = Fraction(4, 1)
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


def test_target_driven_chain_binds_one_literal_artifact_across_all_evidence() -> None:
    cert = certify_target_driven_dependency_chain(
        0.01,
        derivative_order=0,
        target_power=Fraction(1, 50),
        provider=_Provider(),
        initial_lower_bound=7,
    )

    assert cert.selected_order == 1
    assert [link.key for link in cert.links] == [
        (1, "assembled-coefficient-1", "hierarchy.assembled[1]")
    ]
    assert cert.retained_recurrences_exact
    assert cert.one_owned_artifact_per_positive_order
    assert cert.finite_target_prefix_only
    assert not cert.all_order_hierarchy_verified
    assert not cert.infinite_diagonal_schedule_verified
    assert not cert.pde_residual_tail_verified
    assert not cert.all_jets_flat
    assert not cert.super_algebraic
    assert not cert.paper_exact


def test_crosswired_recurrence_fails_even_when_supported_prefix_itself_is_valid() -> None:
    class CrossWiredRecurrence(_Provider):
        def certified_retained_identity(self, order: int) -> ExactRetainedRecurrenceIdentity:
            atom = f"retained-row-{order}"
            return ExactRetainedRecurrenceIdentity(
                order=order,
                hierarchy_id=self.hierarchy_id,
                source_revision=self.source_revision,
                coefficient_state_id=self.coefficient_state_id,
                identity_source=f"hierarchy.retained[{order}]",
                theorem_name="SlowExpansionResidual.retained_recurrence",
                dependencies=(
                    RetainedRecurrenceDependency(
                        coefficient_order=order,
                        artifact=f"different-coefficient-{order}",
                        provider=f"different.provider[{order}]",
                    ),
                ),
                linear_terms=(FormalAtomTerm(atom, 1),),
                pair_terms=(),
                previous_shifted_terms=(FormalAtomTerm(atom, 1),),
            )

    supported = certify_target_driven_supported_hierarchy_prefix(
        0.01,
        derivative_order=0,
        target_power=Fraction(1, 50),
        provider=CrossWiredRecurrence(),
    )
    assert supported.selected_order == 1

    with pytest.raises(ValueError, match="no single owned coefficient artifact"):
        bind_target_driven_dependency_chain(supported)


def test_pairwise_support_bound_matches_do_not_fake_one_global_dependency() -> None:
    class PairwiseOnly(_Provider):
        @staticmethod
        def _key_a(order: int) -> tuple[str, str]:
            return (f"coefficient-a-{order}", f"provider.a[{order}]")

        @staticmethod
        def _key_b(order: int) -> tuple[str, str]:
            return (f"coefficient-b-{order}", f"provider.b[{order}]")

        def certified_jet_bound(self, order: int, derivative_order: int) -> HierarchyJetBound:
            artifact, provider = (
                self._key_a(order) if derivative_order == 0 else self._key_b(order)
            )
            return HierarchyJetBound(
                order=order,
                derivative_order=derivative_order,
                bound=2.0 + order + derivative_order / 16.0,
                hierarchy_id=self.hierarchy_id,
                source_revision=self.source_revision,
                dependencies=(
                    HierarchyBoundDependency(order, artifact, provider),
                ),
            )

        def certified_coefficient_support(self, order: int) -> HierarchyCoefficientSupport:
            B = Fraction(4, 1)
            artifact_a, provider_a = self._key_a(order)
            artifact_b, provider_b = self._key_b(order)
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
                    HierarchySupportDependency(order, artifact_a, provider_a),
                    HierarchySupportDependency(order, artifact_b, provider_b),
                ),
            )

        def certified_retained_identity(self, order: int) -> ExactRetainedRecurrenceIdentity:
            atom = f"retained-row-{order}"
            if order == 0:
                artifact, provider = (self._artifact(0), self._provider(0))
            else:
                artifact, provider = self._key_a(order)
            return ExactRetainedRecurrenceIdentity(
                order=order,
                hierarchy_id=self.hierarchy_id,
                source_revision=self.source_revision,
                coefficient_state_id=self.coefficient_state_id,
                identity_source=f"hierarchy.retained[{order}]",
                theorem_name="SlowExpansionResidual.retained_recurrence",
                dependencies=(RetainedRecurrenceDependency(order, artifact, provider),),
                linear_terms=(FormalAtomTerm(atom, 1),),
                pair_terms=(),
                previous_shifted_terms=(FormalAtomTerm(atom, 1),),
            )

    # The support layer's pairwise condition is satisfied: every C row shares
    # either A or B with support.  The stronger chain correctly rejects it
    # because no one artifact is present in support, recurrence, and all C rows.
    supported = certify_target_driven_supported_hierarchy_prefix(
        0.01,
        derivative_order=0,
        target_power=Fraction(1, 50),
        provider=PairwiseOnly(),
    )
    with pytest.raises(ValueError, match="no single owned coefficient artifact"):
        bind_target_driven_dependency_chain(supported)


def test_nonzero_rational_recurrence_defect_cannot_enter_the_chain() -> None:
    class TinyDefect(_Provider):
        def certified_retained_identity(self, order: int) -> ExactRetainedRecurrenceIdentity:
            atom = f"retained-row-{order}"
            return ExactRetainedRecurrenceIdentity(
                order=order,
                hierarchy_id=self.hierarchy_id,
                source_revision=self.source_revision,
                coefficient_state_id=self.coefficient_state_id,
                identity_source=f"hierarchy.retained[{order}]",
                theorem_name="SlowExpansionResidual.retained_recurrence",
                dependencies=(
                    RetainedRecurrenceDependency(
                        coefficient_order=order,
                        artifact=self._artifact(order),
                        provider=self._provider(order),
                    ),
                ),
                linear_terms=(FormalAtomTerm(atom, 1),),
                pair_terms=(),
                previous_shifted_terms=(
                    FormalAtomTerm(atom, Fraction(10**30 - 1, 10**30)),
                ),
            )

    with pytest.raises(ValueError, match="not an exact zero identity"):
        certify_target_driven_dependency_chain(
            0.01,
            derivative_order=0,
            target_power=Fraction(1, 50),
            provider=TinyDefect(),
        )
