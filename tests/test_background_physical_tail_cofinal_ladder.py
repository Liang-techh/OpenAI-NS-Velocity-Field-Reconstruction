from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_hierarchy_cutoff_bounds import HierarchyBoundDependency
from openai_ns_reconstruction.background_physical_tail_cofinal_ladder import (
    FiniteCofinalPhysicalTailLadderCertificate,
    canonical_level_for_physical_target,
    certify_finite_cofinal_physical_tail_ladder,
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
    hierarchy_id = "section5-cofinal-tail-ladder"
    source_revision = "agent2-cofinal-tail-ladder-source"
    coefficient_state_id = "agent2-cofinal-tail-ladder-state"

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
            bound=Fraction(12 + order, 1) + Fraction(derivative_order, 256),
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


def test_canonical_level_is_exact_least_positive_level_dominating_target() -> None:
    assert canonical_level_for_physical_target(0, Fraction(0, 1)) == 1
    assert canonical_level_for_physical_target(1, Fraction(9, 8)) == 2
    assert canonical_level_for_physical_target(4, Fraction(3, 2)) == 4
    assert canonical_level_for_physical_target(0, Fraction(-9, 8)) == 1

    with pytest.raises(TypeError, match="exact integer/Fraction"):
        canonical_level_for_physical_target(1, 1.25)  # type: ignore[arg-type]


def test_finite_canonical_ladder_uses_one_literal_recursive_schedule_family() -> None:
    cert = certify_finite_cofinal_physical_tail_ladder(
        0.25,
        max_level=2,
        provider=_Provider(),
        minimum_order=4,
        initial_lower_bound=3,
    )

    assert cert.truncation_orders == (11, 23)
    assert cert.levels[0].first_omitted.jet_powers[1].physical_tail_power == 1
    assert cert.levels[1].first_omitted.jet_powers[2].physical_tail_power == 2

    for previous, later in zip(cert.levels[:-1], cert.levels[1:], strict=True):
        schedule0 = previous.first_omitted.uncut.schedule
        schedule1 = later.first_omitted.uncut.schedule
        assert schedule1.local_scales[: len(schedule0.local_scales)] == schedule0.local_scales
        assert schedule1.scales[: len(schedule0.scales)] == schedule0.scales
        assert later.first_omitted.uncut_radius_exact <= previous.first_omitted.uncut_radius_exact


def test_materialized_level_dominates_an_arbitrary_exact_target_within_frontier() -> None:
    cert = certify_finite_cofinal_physical_tail_ladder(
        0.25,
        max_level=2,
        provider=_Provider(),
        minimum_order=4,
        initial_lower_bound=3,
    )

    selected = cert.certificate_for_target(1, Fraction(9, 8))
    assert selected is cert.levels[1]
    assert cert.required_level(1, Fraction(9, 8)) == 2
    assert all(
        row.physical_tail_power >= Fraction(9, 8)
        for row in selected.first_omitted.jet_powers[:2]
    )


def test_target_beyond_materialized_ladder_fails_closed() -> None:
    cert = certify_finite_cofinal_physical_tail_ladder(
        0.25,
        max_level=1,
        provider=_Provider(),
        minimum_order=4,
        initial_lower_bound=3,
    )

    with pytest.raises(ValueError, match="beyond the materialized finite canonical ladder"):
        cert.certificate_for_target(1, Fraction(9, 8))

    with pytest.raises(ValueError, match="max_level must be positive"):
        certify_finite_cofinal_physical_tail_ladder(
            0.25,
            max_level=0,
            provider=_Provider(),
            minimum_order=4,
            initial_lower_bound=3,
        )


class _LateDriftingProvider(_Provider):
    def __init__(self) -> None:
        self._calls: dict[tuple[int, int], int] = {}

    def certified_exact_jet_majorant(
        self, order: int, derivative_order: int
    ) -> ExactHierarchyJetMajorant:
        key = (order, derivative_order)
        self._calls[key] = self._calls.get(key, 0) + 1
        base = super().certified_exact_jet_majorant(order, derivative_order)
        # Levels 1 and 2 see the same row.  The drift appears only when the
        # deeper level 3 request re-asks for the retained prefix.
        if key == (1, 0) and self._calls[key] >= 3:
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


def test_late_provider_drift_is_rejected_at_the_first_incoherent_level() -> None:
    with pytest.raises(ValueError, match=r"shared exact C\[j,m\] majorants changed"):
        certify_finite_cofinal_physical_tail_ladder(
            0.25,
            max_level=3,
            provider=_LateDriftingProvider(),
            minimum_order=4,
            initial_lower_bound=3,
        )


def test_truth_boundary_remains_finite_after_canonical_cofinal_arithmetic() -> None:
    cert = certify_finite_cofinal_physical_tail_ladder(
        0.25,
        max_level=1,
        provider=_Provider(),
        minimum_order=4,
        initial_lower_bound=3,
    )

    assert isinstance(cert, FiniteCofinalPhysicalTailLadderCertificate)
    assert cert.canonical_target_cofinality_arithmetic_verified
    assert cert.finite_shared_schedule_ladder_verified
    assert cert.finite_ladder_only
    assert not cert.provider_totality_verified
    assert not cert.infinite_diagonal_schedule_verified
    assert not cert.actual_pde_residual_verified
    assert not cert.all_order_hierarchy_verified
    assert not cert.all_jets_flat
    assert not cert.super_algebraic
    assert not cert.paper_exact
