from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_hierarchy_cutoff_bounds import HierarchyBoundDependency
from openai_ns_reconstruction.background_physical_tail_prefix_coherence import (
    FiniteIncreasingPhysicalTailCertificate,
    certify_increasing_physical_tail_pair,
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


class _Provider:
    hierarchy_id = "section5-prefix-coherence"
    source_revision = "agent2-prefix-coherence-source"
    coefficient_state_id = "agent2-prefix-coherence-state"

    @staticmethod
    def _artifact(order: int) -> str:
        return f"assembled-coefficient-{order}"

    @staticmethod
    def _provider(order: int) -> str:
        return f"hierarchy.assembled[{order}]"

    def certified_exact_jet_majorant(
        self, order: int, derivative_order: int
    ) -> ExactHierarchyJetMajorant:
        return ExactHierarchyJetMajorant(
            order=order,
            derivative_order=derivative_order,
            bound=Fraction(8 + order, 1) + Fraction(derivative_order, 128),
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

    def certified_retained_identity(self, order: int) -> ExactRetainedRecurrenceIdentity:
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


def test_stronger_power_request_extends_one_literal_recursive_schedule() -> None:
    cert = certify_increasing_physical_tail_pair(
        0.125,
        earlier_max_derivative_order=1,
        earlier_target_physical_power=Fraction(1, 8),
        later_max_derivative_order=1,
        later_target_physical_power=Fraction(3, 8),
        provider=_Provider(),
        minimum_order=4,
        initial_lower_bound=3,
    )

    assert cert.earlier.truncation_order == 16
    assert cert.later.truncation_order == 18
    assert cert.truncation_increment == 2
    assert cert.physical_power_gain_exact == Fraction(1, 4)
    assert cert.raw_first_omitted_power_gain_exact == Fraction(1, 2)
    assert cert.dyadic_prefactor_ratio_exact == Fraction(1, 4)

    schedule0 = cert.earlier.first_omitted.uncut.schedule
    schedule1 = cert.later.first_omitted.uncut.schedule
    assert schedule1.scales[: len(schedule0.scales)] == schedule0.scales
    assert schedule1.local_scales[: len(schedule0.local_scales)] == schedule0.local_scales
    assert cert.later.first_omitted.uncut_radius_exact <= cert.earlier.first_omitted.uncut_radius_exact


def test_larger_derivative_budget_preserves_existing_rows_and_gains_exact_power() -> None:
    cert = certify_increasing_physical_tail_pair(
        0.125,
        earlier_max_derivative_order=0,
        earlier_target_physical_power=Fraction(1, 8),
        later_max_derivative_order=1,
        later_target_physical_power=Fraction(1, 8),
        provider=_Provider(),
        minimum_order=4,
        initial_lower_bound=3,
    )

    assert cert.earlier.truncation_order == 4
    assert cert.later.truncation_order == 16
    assert cert.truncation_increment == 12
    assert cert.physical_power_gain_exact == Fraction(3, 2)

    row0_before = cert.earlier.first_omitted.jet_powers[0]
    row0_after = cert.later.first_omitted.jet_powers[0]
    assert row0_after.physical_tail_power - row0_before.physical_tail_power == Fraction(3, 2)
    assert row0_after.ordinary_tail_power - row0_before.ordinary_tail_power == Fraction(3, 2)


class _DriftingMajorantProvider(_Provider):
    def __init__(self) -> None:
        self._calls: dict[tuple[int, int], int] = {}

    def certified_exact_jet_majorant(
        self, order: int, derivative_order: int
    ) -> ExactHierarchyJetMajorant:
        key = (order, derivative_order)
        self._calls[key] = self._calls.get(key, 0) + 1
        base = super().certified_exact_jet_majorant(order, derivative_order)
        if key == (1, 0) and self._calls[key] >= 2:
            return ExactHierarchyJetMajorant(
                order=base.order,
                derivative_order=base.derivative_order,
                bound=base.bound + 1,
                hierarchy_id=base.hierarchy_id,
                source_revision=base.source_revision,
                coefficient_state_id=base.coefficient_state_id,
                dependencies=base.dependencies,
            )
        return base


def test_provider_that_changes_shared_majorant_across_requests_fails_closed() -> None:
    with pytest.raises(ValueError, match=r"shared exact C\[j,m\] majorants changed"):
        certify_increasing_physical_tail_pair(
            0.125,
            earlier_max_derivative_order=0,
            earlier_target_physical_power=Fraction(1, 8),
            later_max_derivative_order=1,
            later_target_physical_power=Fraction(1, 8),
            provider=_DriftingMajorantProvider(),
            minimum_order=4,
            initial_lower_bound=3,
        )


def test_dominating_request_may_not_reduce_derivative_budget() -> None:
    with pytest.raises(ValueError, match="later_max_derivative_order"):
        certify_increasing_physical_tail_pair(
            0.125,
            earlier_max_derivative_order=1,
            earlier_target_physical_power=Fraction(1, 8),
            later_max_derivative_order=0,
            later_target_physical_power=Fraction(1, 8),
            provider=_Provider(),
            minimum_order=4,
            initial_lower_bound=3,
        )


def test_dominating_request_may_not_reduce_target_power() -> None:
    with pytest.raises(ValueError, match="must not reduce the target physical power"):
        certify_increasing_physical_tail_pair(
            0.125,
            earlier_max_derivative_order=0,
            earlier_target_physical_power=Fraction(3, 8),
            later_max_derivative_order=0,
            later_target_physical_power=Fraction(1, 8),
            provider=_Provider(),
            minimum_order=4,
            initial_lower_bound=3,
        )


def test_truth_boundary_stays_finite_even_after_prefix_coherence() -> None:
    cert = certify_increasing_physical_tail_pair(
        0.125,
        earlier_max_derivative_order=0,
        earlier_target_physical_power=Fraction(1, 8),
        later_max_derivative_order=1,
        later_target_physical_power=Fraction(1, 8),
        provider=_Provider(),
        minimum_order=4,
        initial_lower_bound=3,
    )

    assert isinstance(cert, FiniteIncreasingPhysicalTailCertificate)
    assert cert.shared_schedule_prefix_verified
    assert cert.shared_coefficient_evidence_prefix_verified
    assert cert.increasing_truncation_exponent_gain_verified
    assert cert.finite_request_pair_only
    assert not cert.provider_totality_verified
    assert not cert.infinite_diagonal_schedule_verified
    assert not cert.actual_pde_residual_verified
    assert not cert.all_order_hierarchy_verified
    assert not cert.all_jets_flat
    assert not cert.super_algebraic
    assert not cert.paper_exact
