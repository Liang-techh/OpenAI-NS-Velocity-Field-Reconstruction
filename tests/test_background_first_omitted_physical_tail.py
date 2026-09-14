from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_first_omitted_physical_tail import (
    FiniteFirstOmittedPhysicalTailCertificate,
    certify_finite_first_omitted_physical_tail,
)
from openai_ns_reconstruction.background_hierarchy_cutoff_bounds import HierarchyBoundDependency
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


class _FirstOmittedProvider:
    hierarchy_id = "section5-first-omitted-owned"
    source_revision = "agent2-first-omitted-source"
    coefficient_state_id = "agent2-first-omitted-state"

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


def _certificate() -> FiniteFirstOmittedPhysicalTailCertificate:
    return certify_finite_first_omitted_physical_tail(
        0.125,
        max_derivative_order=1,
        target_physical_power=Fraction(1, 8),
        provider=_FirstOmittedProvider(),
        minimum_order=4,
        initial_lower_bound=3,
    )


def test_first_omitted_power_is_linked_to_same_physical_prefix() -> None:
    cert = _certificate()

    assert cert.truncation_order == 16
    assert cert.first_omitted_order == 17
    assert cert.h_exact == Fraction(1, 8)
    assert cert.raw_first_omitted_power_exact == Fraction(17, 4)
    assert cert.diagonal_tail_power_exact == Fraction(17, 8)
    assert cert.dyadic_prefactor_exact == Fraction(1, 65536)
    assert cert.uncut_radius_exact == Fraction(1, 102021672599552)

    rows = cert.jet_powers
    assert tuple(row.derivative_order for row in rows) == (0, 1)
    assert tuple(row.ordinary_tail_power for row in rows) == (
        Fraction(17, 8),
        Fraction(9, 8),
    )
    assert tuple(row.physical_tail_power for row in rows) == (
        Fraction(17, 8),
        Fraction(1, 8),
    )
    assert tuple(row.target_slack for row in rows) == (Fraction(2, 1), Fraction(0, 1))


def test_retained_exact_cancellation_stops_before_first_omitted_order() -> None:
    cert = _certificate()

    assert cert.retained_orders_exactly_cancelled == tuple(range(17))
    assert cert.first_omitted_order == cert.retained_orders_exactly_cancelled[-1] + 1
    assert cert.retained_recurrences_exact
    assert cert.first_omitted_order_exposed

    # The omitted order is deliberately not promoted to a cancelled recurrence
    # or a coefficient-majorant theorem merely because its q-power is known.
    assert not cert.first_omitted_recurrence_cancelled
    assert not cert.first_omitted_coefficient_majorant_verified


def test_truth_boundary_remains_finite_only() -> None:
    cert = _certificate()

    assert not cert.physical_chart_finite_bound_verified
    assert not cert.actual_pde_residual_verified
    assert not cert.all_order_hierarchy_verified
    assert not cert.infinite_diagonal_schedule_verified
    assert not cert.all_jets_flat
    assert not cert.super_algebraic
    assert not cert.paper_exact


def test_zero_derivative_budget_reduces_to_half_of_full_stage_power() -> None:
    cert = certify_finite_first_omitted_physical_tail(
        0.125,
        max_derivative_order=0,
        target_physical_power=Fraction(1, 2),
        provider=_FirstOmittedProvider(),
        minimum_order=4,
        initial_lower_bound=3,
    )

    row = cert.jet_powers[0]
    assert row.ordinary_tail_power == cert.diagonal_tail_power_exact
    assert row.physical_tail_power == cert.diagonal_tail_power_exact
    assert row.raw_first_omitted_power == 2 * row.physical_tail_power
    assert row.physical_tail_power >= Fraction(1, 2)


def test_theorem_identity_drift_is_rejected() -> None:
    base = _certificate()

    with pytest.raises(ValueError, match="stage_power_theorem"):
        FiniteFirstOmittedPhysicalTailCertificate(
            uncut=base.uncut,
            stage_power_theorem="wrong.theorem",
        )
    with pytest.raises(ValueError, match="ordinary_tail_theorem"):
        FiniteFirstOmittedPhysicalTailCertificate(
            uncut=base.uncut,
            ordinary_tail_theorem="wrong.theorem",
        )
    with pytest.raises(ValueError, match="physical_tail_theorem"):
        FiniteFirstOmittedPhysicalTailCertificate(
            uncut=base.uncut,
            physical_tail_theorem="wrong.theorem",
        )


def test_approximate_target_remains_forbidden_upstream() -> None:
    with pytest.raises(TypeError, match="floats are forbidden"):
        certify_finite_first_omitted_physical_tail(
            0.125,
            max_derivative_order=1,
            target_physical_power=0.125,
            provider=_FirstOmittedProvider(),
            minimum_order=4,
            initial_lower_bound=3,
        )
