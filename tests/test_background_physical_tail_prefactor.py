from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_hierarchy_cutoff_bounds import HierarchyBoundDependency
from openai_ns_reconstruction.background_physical_tail_prefactor import (
    FinitePhysicalTailPrefactorCertificate,
    PhysicalChartExistentialJetConstant,
    PhysicalChartFiniteBoundShape,
    certify_finite_physical_tail_prefactor,
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
    hierarchy_id = "section5-physical-prefactor"
    source_revision = "agent2-physical-prefactor-source"
    coefficient_state_id = "agent2-physical-prefactor-state"

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


def _certificate() -> FinitePhysicalTailPrefactorCertificate:
    return certify_finite_physical_tail_prefactor(
        0.125,
        max_derivative_order=1,
        target_physical_power=Fraction(1, 8),
        provider=_Provider(),
        minimum_order=4,
        initial_lower_bound=3,
    )


def test_physical_chart_bound_uses_hierarchy_owned_support_interval() -> None:
    cert = _certificate()
    assert cert.chart_bound.x_lower == Fraction(-1, 1)
    assert cert.chart_bound.x_upper == Fraction(8, 1)
    assert cert.chart_bound.D_expression == (1, ("C_phys_0", "C_phys_1"))
    assert cert.chart_bound.existential_constants_positive
    assert not cert.chart_bound.numeric_D_materialized


def test_chain_rule_prefactor_shape_is_exact_for_every_requested_jet() -> None:
    cert = _certificate()
    assert cert.truncation_order == 16
    assert cert.first_omitted_order == 17

    row0, row1 = cert.rows
    assert row0.factorial_exact == 1
    assert row0.dyadic_prefactor_exact == Fraction(1, 65536)
    assert row0.D_power == 0
    assert row0.physical_power == Fraction(17, 8)
    assert row0.target_slack == Fraction(2, 1)

    assert row1.factorial_exact == 1
    assert row1.dyadic_prefactor_exact == Fraction(1, 65536)
    assert row1.D_power == 1
    assert row1.physical_power == Fraction(1, 8)
    assert row1.target_slack == Fraction(0, 1)


def test_second_derivative_budget_records_factorial_and_D_square_without_materializing_D() -> None:
    cert = certify_finite_physical_tail_prefactor(
        0.125,
        max_derivative_order=2,
        target_physical_power=Fraction(1, 8),
        provider=_Provider(),
        minimum_order=4,
        initial_lower_bound=3,
    )
    row = cert.rows[2]
    assert row.factorial_exact == 2
    assert row.D_power == 2
    assert row.physical_power >= Fraction(1, 8)
    assert cert.chart_bound.D_expression[0] == 1
    assert cert.chart_bound.D_expression[1] == ("C_phys_0", "C_phys_1", "C_phys_2")


def test_cross_wired_support_interval_is_rejected() -> None:
    cert = _certificate()
    bad_chart = PhysicalChartFiniteBoundShape(
        max_derivative_order=1,
        x_lower=Fraction(-2, 1),
        x_upper=Fraction(8, 1),
        jet_constants=(
            PhysicalChartExistentialJetConstant(0, "C_phys_0"),
            PhysicalChartExistentialJetConstant(1, "C_phys_1"),
        ),
    )
    with pytest.raises(ValueError, match="support box"):
        FinitePhysicalTailPrefactorCertificate(
            first_omitted=cert.first_omitted,
            chart_bound=bad_chart,
        )


def test_missing_or_relabelled_chart_constant_is_rejected() -> None:
    with pytest.raises(ValueError, match="exactly 0..M"):
        PhysicalChartFiniteBoundShape(
            max_derivative_order=1,
            x_lower=Fraction(-1, 1),
            x_upper=Fraction(8, 1),
            jet_constants=(PhysicalChartExistentialJetConstant(0, "C_phys_0"),),
        )
    with pytest.raises(ValueError, match="must be 'C_phys_1'"):
        PhysicalChartExistentialJetConstant(1, "C1")


def test_theorem_drift_is_rejected() -> None:
    cert = _certificate()
    with pytest.raises(ValueError, match="composition_theorem"):
        FinitePhysicalTailPrefactorCertificate(
            first_omitted=cert.first_omitted,
            chart_bound=cert.chart_bound,
            composition_theorem="wrong.theorem",
        )


def test_truth_boundary_does_not_promote_symbolic_finiteness_to_residual_convergence() -> None:
    cert = _certificate()
    assert cert.physical_chart_finite_bound_shape_verified
    assert cert.physical_composition_prefactor_shape_verified
    assert not cert.numeric_physical_chart_constant_materialized
    assert not cert.actual_pde_residual_verified
    assert not cert.all_order_hierarchy_verified
    assert not cert.infinite_diagonal_schedule_verified
    assert not cert.all_jets_flat
    assert not cert.super_algebraic
    assert not cert.paper_exact
