from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_hierarchy_cutoff_bounds import HierarchyBoundDependency
from openai_ns_reconstruction.background_joint_physical_uncut_radius import (
    certify_finite_joint_physical_uncut_radius,
)
from openai_ns_reconstruction.background_recurrence_cancellation import (
    ExactRetainedRecurrenceIdentity,
    FormalAtomTerm,
    RetainedRecurrenceDependency,
)
from openai_ns_reconstruction.background_target_driven_exact_majorants import ExactHierarchyJetMajorant
from openai_ns_reconstruction.background_target_driven_support_prefix import (
    HierarchyCoefficientSupport,
    HierarchySupportDependency,
    PINNED_COEFFICIENT_COMPONENTS,
    PINNED_COMPACT_SUPPORT_THEOREM,
    PINNED_SUPPORT_THEOREM,
)


class _UncutProvider:
    hierarchy_id = "section5-uncut-owned"
    source_revision = "agent2-uncut-source"
    coefficient_state_id = "agent2-uncut-state"

    @staticmethod
    def _artifact(order: int) -> str:
        return f"assembled-coefficient-{order}"

    @staticmethod
    def _provider(order: int) -> str:
        return f"hierarchy.assembled[{order}]"

    def certified_exact_jet_majorant(self, order: int, derivative_order: int) -> ExactHierarchyJetMajorant:
        return ExactHierarchyJetMajorant(
            order=order,
            derivative_order=derivative_order,
            bound=Fraction(6 + order, 1) + Fraction(derivative_order, 64),
            hierarchy_id=self.hierarchy_id,
            source_revision=self.source_revision,
            coefficient_state_id=self.coefficient_state_id,
            dependencies=(HierarchyBoundDependency(
                coefficient_order=order,
                artifact=self._artifact(order),
                provider=self._provider(order),
            ),),
        )

    def certified_retained_identity(self, order: int) -> ExactRetainedRecurrenceIdentity:
        atom = f"retained-{order}"
        return ExactRetainedRecurrenceIdentity(
            order=order,
            hierarchy_id=self.hierarchy_id,
            source_revision=self.source_revision,
            coefficient_state_id=self.coefficient_state_id,
            identity_source=f"hierarchy.retained[{order}]",
            theorem_name="SlowExpansionResidual.retained_recurrence",
            dependencies=(RetainedRecurrenceDependency(
                coefficient_order=order,
                artifact=self._artifact(order),
                provider=self._provider(order),
            ),),
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
            dependencies=(HierarchySupportDependency(
                coefficient_order=order,
                artifact=self._artifact(order),
                provider=self._provider(order),
            ),),
        )


def _certificate():
    return certify_finite_joint_physical_uncut_radius(
        0.125,
        max_derivative_order=1,
        target_physical_power=Fraction(1, 8),
        provider=_UncutProvider(),
        minimum_order=4,
        initial_lower_bound=3,
    )


def test_actual_recursive_schedule_gives_exact_common_uncut_radius() -> None:
    cert = _certificate()

    assert cert.selected_order == 16
    assert cert.selected_scale == 51010836299776
    assert cert.uncut_radius_exact == Fraction(1, 102021672599552)
    assert cert.retained_plateau_edges_exact[-1] == Fraction(1, 2)
    assert all(edge <= Fraction(1, 2) for edge in cert.retained_plateau_edges_exact)

    q = cert.uncut_radius_exact / 2
    arguments = cert.certify_q(q)
    assert arguments[-1] == Fraction(1, 4)
    assert all(argument < Fraction(1, 2) for argument in arguments)

    assert cert.finite_uncut_plateau_arithmetic_verified
    assert cert.one_provider_owned_schedule_used
    assert not cert.function_level_uncut_germ_verified
    assert not cert.physical_chart_finite_bound_verified
    assert not cert.actual_slow_sum_materialized
    assert not cert.pde_residual_tail_verified
    assert not cert.all_order_hierarchy_verified
    assert not cert.infinite_diagonal_schedule_verified
    assert not cert.all_jets_flat
    assert not cert.super_algebraic
    assert not cert.paper_exact


def test_uncut_window_is_strict_and_exact_only() -> None:
    cert = _certificate()

    with pytest.raises(ValueError, match="strictly inside"):
        cert.certify_q(cert.uncut_radius_exact)
    with pytest.raises(ValueError, match="strictly inside"):
        cert.certify_q(cert.uncut_radius_exact * 2)
    with pytest.raises(TypeError, match="floats are forbidden"):
        cert.certify_q(float(cert.uncut_radius_exact / 2))


def test_uncut_radius_is_derived_from_same_selected_schedule() -> None:
    cert = _certificate()
    schedule = cert.schedule

    assert schedule.max_order == cert.physical.selected_order
    assert cert.selected_scale == schedule.scales[-1]
    assert cert.uncut_radius_exact == Fraction(1, 2 * schedule.scales[-1])
    assert all(b >= a for a, b in zip(schedule.scales, schedule.scales[1:]))
    assert all(b >= 2 * a for a, b in zip(schedule.scales, schedule.scales[1:]))
