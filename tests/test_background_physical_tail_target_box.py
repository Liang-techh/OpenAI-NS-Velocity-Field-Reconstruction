from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_hierarchy_cutoff_bounds import HierarchyBoundDependency
from openai_ns_reconstruction.background_physical_tail_cofinal_ladder import (
    certify_finite_cofinal_physical_tail_ladder,
)
from openai_ns_reconstruction.background_physical_tail_target_box import (
    FinitePhysicalTailTargetBoxCertificate,
    certify_finite_physical_tail_target_box,
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
    hierarchy_id = "section5-finite-target-box"
    source_revision = "agent2-finite-target-box-source"
    coefficient_state_id = "agent2-finite-target-box-state"

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
            bound=Fraction(20 + order, 1) + Fraction(derivative_order, 512),
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


def _ladder(max_level: int):
    return certify_finite_cofinal_physical_tail_ladder(
        0.25,
        max_level=max_level,
        provider=_Provider(),
        minimum_order=4,
        initial_lower_bound=3,
    )


def test_one_frontier_prefix_covers_complete_finite_jet_power_rectangle() -> None:
    ladder = _ladder(3)
    cert = certify_finite_physical_tail_target_box(
        ladder,
        max_derivative_order=2,
        max_integer_power=3,
    )

    assert cert.frontier_level == 3
    assert cert.frontier_certificate is ladder.levels[2]
    assert cert.target_count == 9
    assert tuple(
        (cell.derivative_order, cell.requested_integer_power)
        for cell in cert.cells
    ) == tuple((m, N) for m in range(3) for N in range(1, 4))

    expected_minimal_levels = (
        (1, 2, 3),
        (1, 2, 3),
        (2, 2, 3),
    )
    for m in range(3):
        for N in range(1, 4):
            cell = cert.cell(m, N)
            assert cell.minimal_canonical_level == expected_minimal_levels[m][N - 1]
            assert cell.frontier_certificate is cert.frontier_certificate
            assert cell.prefactor_row == cert.frontier_certificate.rows[m]
            assert cell.prefactor_row.derivative_order == m
            assert cell.prefactor_row.D_power == m
            assert cell.prefactor_row.physical_power >= N
            assert cell.physical_power_slack >= 0


def test_target_box_uses_one_dyadic_frontier_prefactor_for_every_cell() -> None:
    cert = certify_finite_physical_tail_target_box(
        _ladder(3),
        max_derivative_order=2,
        max_integer_power=3,
    )
    expected = cert.frontier_certificate.first_omitted.dyadic_prefactor_exact
    assert all(cell.prefactor_row.dyadic_prefactor_exact == expected for cell in cert.cells)


def test_box_beyond_materialized_cofinal_frontier_fails_closed() -> None:
    ladder = _ladder(2)
    with pytest.raises(ValueError, match="exceeds the materialized cofinal ladder"):
        certify_finite_physical_tail_target_box(
            ladder,
            max_derivative_order=2,
            max_integer_power=3,
        )


def test_invalid_finite_box_requests_fail_closed() -> None:
    ladder = _ladder(1)
    with pytest.raises(ValueError, match="max_integer_power must be positive"):
        certify_finite_physical_tail_target_box(
            ladder,
            max_derivative_order=0,
            max_integer_power=0,
        )
    with pytest.raises(ValueError, match="max_derivative_order must be a nonnegative integer"):
        certify_finite_physical_tail_target_box(
            ladder,
            max_derivative_order=True,  # type: ignore[arg-type]
            max_integer_power=1,
        )
    with pytest.raises(TypeError, match="ladder must be"):
        certify_finite_physical_tail_target_box(  # type: ignore[arg-type]
            object(),
            max_derivative_order=0,
            max_integer_power=1,
        )


def test_target_lookup_outside_box_fails_closed() -> None:
    cert = certify_finite_physical_tail_target_box(
        _ladder(2),
        max_derivative_order=1,
        max_integer_power=2,
    )
    with pytest.raises(ValueError, match="outside the certified finite target box"):
        cert.cell(2, 1)
    with pytest.raises(ValueError, match="outside the certified finite target box"):
        cert.cell(1, 3)


def test_truth_boundary_stays_finite_after_rectangular_quantifier_check() -> None:
    cert = certify_finite_physical_tail_target_box(
        _ladder(2),
        max_derivative_order=1,
        max_integer_power=2,
    )

    assert isinstance(cert, FinitePhysicalTailTargetBoxCertificate)
    assert cert.one_common_frontier_prefix_verified
    assert cert.finite_all_jet_integer_power_box_verified
    assert cert.finite_target_box_only
    assert not cert.provider_totality_verified
    assert not cert.infinite_diagonal_schedule_verified
    assert not cert.actual_pde_residual_verified
    assert not cert.all_order_hierarchy_verified
    assert not cert.all_jets_flat
    assert not cert.super_algebraic
    assert not cert.paper_exact
