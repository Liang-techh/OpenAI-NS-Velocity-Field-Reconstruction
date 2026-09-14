from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_hierarchy_cutoff_bounds import (
    HierarchyBoundDependency,
)
from openai_ns_reconstruction.background_joint_physical_jet_target_prefix import (
    certify_finite_joint_physical_jet_target_prefix,
)
from openai_ns_reconstruction.background_recurrence_cancellation import (
    ExactRetainedRecurrenceIdentity,
    FormalAtomTerm,
    RetainedRecurrenceDependency,
)
from openai_ns_reconstruction.background_target_driven_exact_majorants import (
    ExactHierarchyJetMajorant,
)
from openai_ns_reconstruction.background_target_driven_support_prefix import (
    HierarchyCoefficientSupport,
    HierarchySupportDependency,
    PINNED_COEFFICIENT_COMPONENTS,
    PINNED_COMPACT_SUPPORT_THEOREM,
    PINNED_SUPPORT_THEOREM,
)


class _PhysicalJointProvider:
    hierarchy_id = "section5-physical-joint-owned"
    source_revision = "agent2-physical-joint-source"
    coefficient_state_id = "agent2-physical-joint-state"

    def __init__(self, max_order: int | None = None) -> None:
        self.max_order = max_order

    @staticmethod
    def _artifact(order: int) -> str:
        return f"assembled-coefficient-{order}"

    @staticmethod
    def _provider(order: int) -> str:
        return f"hierarchy.assembled[{order}]"

    def _check_frontier(self, order: int) -> None:
        if self.max_order is not None and order > self.max_order:
            raise RuntimeError("real hierarchy frontier reached")

    def certified_exact_jet_majorant(
        self, order: int, derivative_order: int
    ) -> ExactHierarchyJetMajorant:
        self._check_frontier(order)
        return ExactHierarchyJetMajorant(
            order=order,
            derivative_order=derivative_order,
            bound=Fraction(6 + order, 1) + Fraction(derivative_order, 64),
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
        atom = f"retained-{order}"
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
        self._check_frontier(order)
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


def test_one_prefix_spends_physical_derivative_loss_for_every_requested_jet() -> None:
    cert = certify_finite_joint_physical_jet_target_prefix(
        0.125,
        max_derivative_order=1,
        target_physical_power=Fraction(1, 8),
        provider=_PhysicalJointProvider(),
        minimum_order=4,
        initial_lower_bound=3,
    )

    # Physical m=1 requires the ordinary request P+m=9/8:
    # (J+1)/8 - 1 >= 9/8 first holds at J=16.
    assert cert.selected_order == 16
    assert cert.ordinary_target_power_exact == Fraction(9, 8)
    assert cert.ordinary_tail_powers_exact == (Fraction(17, 8), Fraction(9, 8))
    assert cert.physical_tail_powers_exact == (Fraction(17, 8), Fraction(1, 8))
    assert cert.weakest_physical_tail_power_exact == Fraction(1, 8)
    assert cert.dyadic_prefactor_exact == Fraction(1, 65536)

    assert cert.joint_finite_physical_jet_target_arithmetic_verified
    assert cert.one_common_recursive_schedule_for_requested_physical_jets
    assert cert.retained_recurrences_exact

    prefix = cert.ordinary.chain.chain.supported_prefix.prefix
    assert prefix.derivative_order == 1
    assert prefix.target.target_power == Fraction(9, 8)
    assert prefix.cutoff.bounds.max_order == 16
    assert prefix.cancellations.max_order == 16

    assert cert.finite_physical_jet_prefix_only
    assert not cert.physical_chart_finite_bound_verified
    assert not cert.uncut_tail_germ_verified
    assert not cert.all_order_hierarchy_verified
    assert not cert.infinite_diagonal_schedule_verified
    assert not cert.pde_residual_tail_verified
    assert not cert.all_jets_flat
    assert not cert.super_algebraic
    assert not cert.paper_exact


def test_zero_derivative_budget_reduces_to_the_ordinary_target() -> None:
    cert = certify_finite_joint_physical_jet_target_prefix(
        0.125,
        max_derivative_order=0,
        target_physical_power=Fraction(3, 8),
        provider=_PhysicalJointProvider(),
    )

    assert cert.ordinary_target_power_exact == Fraction(3, 8)
    assert cert.physical_tail_powers_exact == cert.ordinary_tail_powers_exact
    assert cert.weakest_physical_tail_power_exact >= Fraction(3, 8)


def test_joint_physical_target_fails_closed_at_the_real_hierarchy_frontier() -> None:
    with pytest.raises(RuntimeError, match="real hierarchy frontier reached"):
        certify_finite_joint_physical_jet_target_prefix(
            0.125,
            max_derivative_order=1,
            target_physical_power=Fraction(1, 8),
            provider=_PhysicalJointProvider(max_order=15),
            minimum_order=4,
        )


def test_pinned_physical_tail_requires_h_strictly_below_one_half() -> None:
    with pytest.raises(ValueError, match="0<h<1/2"):
        certify_finite_joint_physical_jet_target_prefix(
            0.5,
            max_derivative_order=0,
            target_physical_power=Fraction(1, 2),
            provider=_PhysicalJointProvider(),
        )


def test_float_physical_target_cannot_cross_the_exact_gate() -> None:
    with pytest.raises(TypeError, match="floats are forbidden"):
        certify_finite_joint_physical_jet_target_prefix(
            0.125,
            max_derivative_order=1,
            target_physical_power=0.125,
            provider=_PhysicalJointProvider(),
        )
