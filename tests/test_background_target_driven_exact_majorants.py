from fractions import Fraction
import math

import pytest

from openai_ns_reconstruction.background_hierarchy_cutoff_bounds import (
    HierarchyBoundDependency,
)
from openai_ns_reconstruction.background_recurrence_cancellation import (
    ExactRetainedRecurrenceIdentity,
    FormalAtomTerm,
    RetainedRecurrenceDependency,
)
from openai_ns_reconstruction.background_target_driven_exact_majorants import (
    ExactHierarchyJetMajorant,
    certify_target_driven_exact_majorant_chain,
    outward_float_upper,
)
from openai_ns_reconstruction.background_target_driven_support_prefix import (
    HierarchyCoefficientSupport,
    HierarchySupportDependency,
    PINNED_COEFFICIENT_COMPONENTS,
    PINNED_COMPACT_SUPPORT_THEOREM,
    PINNED_SUPPORT_THEOREM,
)


class _ExactProvider:
    hierarchy_id = "section5-owned-hierarchy"
    source_revision = "hierarchy-revision-agent2"
    coefficient_state_id = "coefficient-state-agent2"

    @staticmethod
    def _artifact(order: int) -> str:
        return f"assembled-coefficient-{order}"

    @staticmethod
    def _provider(order: int) -> str:
        return f"hierarchy.assembled[{order}]"

    def certified_exact_jet_majorant(
        self, order: int, derivative_order: int
    ) -> ExactHierarchyJetMajorant:
        if order == 1 and derivative_order == 0:
            # Nearest binary64 rounds this exact value downward to 1.0.
            bound = Fraction.from_float(1.0) + Fraction(1, 2**60)
        else:
            bound = Fraction(2 + order, 1) + Fraction(derivative_order, 16)
        return ExactHierarchyJetMajorant(
            order=order,
            derivative_order=derivative_order,
            bound=bound,
            hierarchy_id=self.hierarchy_id,
            source_revision=self.source_revision,
            coefficient_state_id=self.coefficient_state_id,
            dependencies=(
                HierarchyBoundDependency(
                    coefficient_order=order,
                    artifact=self._artifact(order),
                    provider=self._provider(order),
                ),
            ),
        )

    def certified_retained_identity(
        self, order: int
    ) -> ExactRetainedRecurrenceIdentity:
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


def test_outward_rounding_steps_above_a_downrounded_nearest_float() -> None:
    exact = Fraction.from_float(1.0) + Fraction(1, 2**60)
    assert float(exact) == 1.0

    rounded = outward_float_upper(exact)

    assert rounded == math.nextafter(1.0, math.inf)
    assert Fraction.from_float(rounded) >= exact


def test_target_chain_consumes_only_provider_owned_exact_majorants() -> None:
    cert = certify_target_driven_exact_majorant_chain(
        0.01,
        derivative_order=0,
        target_power=Fraction(1, 50),
        provider=_ExactProvider(),
        initial_lower_bound=7,
    )

    assert cert.selected_order == 1
    assert cert.outward_binary64_majorants_verified
    assert [
        (entry.order, entry.derivative_order)
        for entry in cert.exact_majorants
    ] == [(1, 0), (1, 1), (1, 2), (1, 3)]
    first = cert.exact_majorants[0]
    runtime = cert.chain.supported_prefix.prefix.cutoff.bounds.entry(1, 0)
    assert float(first.bound) == 1.0
    assert runtime.bound == math.nextafter(1.0, math.inf)
    assert Fraction.from_float(runtime.bound) >= first.bound
    assert cert.chain.links[0].key == (
        1,
        "assembled-coefficient-1",
        "hierarchy.assembled[1]",
    )
    assert cert.finite_target_prefix_only
    assert not cert.all_order_hierarchy_verified
    assert not cert.infinite_diagonal_schedule_verified
    assert not cert.pde_residual_tail_verified
    assert not cert.all_jets_flat
    assert not cert.super_algebraic
    assert not cert.paper_exact


def test_float_cannot_pose_as_an_exact_majorant() -> None:
    with pytest.raises(TypeError, match="exact integer/Fraction"):
        ExactHierarchyJetMajorant(
            order=1,
            derivative_order=0,
            bound=1.0,
            hierarchy_id="h",
            source_revision="r",
            coefficient_state_id="s",
            dependencies=(HierarchyBoundDependency(1, "a", "p"),),
        )


def test_exact_majorant_coefficient_state_must_match_provider() -> None:
    class WrongState(_ExactProvider):
        def certified_exact_jet_majorant(
            self, order: int, derivative_order: int
        ) -> ExactHierarchyJetMajorant:
            result = super().certified_exact_jet_majorant(order, derivative_order)
            return ExactHierarchyJetMajorant(
                order=result.order,
                derivative_order=result.derivative_order,
                bound=result.bound,
                hierarchy_id=result.hierarchy_id,
                source_revision=result.source_revision,
                coefficient_state_id="different-state",
                dependencies=result.dependencies,
            )

    with pytest.raises(ValueError, match="different coefficient state"):
        certify_target_driven_exact_majorant_chain(
            0.01,
            derivative_order=0,
            target_power=Fraction(1, 50),
            provider=WrongState(),
        )


def test_majorant_outside_finite_binary64_range_fails_closed() -> None:
    with pytest.raises(OverflowError, match="finite binary64 scheduler range"):
        outward_float_upper(Fraction(10**10000, 1))


def test_target_that_crosses_provider_frontier_fails_without_partial_certificate() -> None:
    class FirstOrderOnly(_ExactProvider):
        def certified_exact_jet_majorant(
            self, order: int, derivative_order: int
        ) -> ExactHierarchyJetMajorant:
            if order > 1:
                raise RuntimeError("real hierarchy frontier reached")
            return super().certified_exact_jet_majorant(order, derivative_order)

    # h=1/100 and N=3/100 select J=2 because h*(J+1) >= N.
    with pytest.raises(RuntimeError, match="real hierarchy frontier reached"):
        certify_target_driven_exact_majorant_chain(
            0.01,
            derivative_order=0,
            target_power=Fraction(3, 100),
            provider=FirstOrderOnly(),
        )
