from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_hierarchy_cutoff_bounds import (
    HierarchyBoundDependency,
)
from openai_ns_reconstruction.background_joint_jet_target_prefix import (
    certify_finite_joint_jet_target_prefix,
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


class _JointProvider:
    hierarchy_id = "section5-joint-owned"
    source_revision = "agent2-joint-source"
    coefficient_state_id = "agent2-joint-state"

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
            bound=Fraction(4 + order, 1) + Fraction(derivative_order, 32),
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


def test_one_prefix_hits_one_target_for_every_requested_ordinary_jet() -> None:
    cert = certify_finite_joint_jet_target_prefix(
        0.125,
        max_derivative_order=1,
        target_power=Fraction(1, 8),
        provider=_JointProvider(),
        minimum_order=4,
        initial_lower_bound=3,
    )

    # (J+1)/8 - 1 >= 1/8 first holds at J=8 above Jmin=4.
    assert cert.selected_order == 8
    assert cert.ordinary_tail_powers_exact == (Fraction(9, 8), Fraction(1, 8))
    assert cert.weakest_tail_power_exact == Fraction(1, 8)
    assert cert.dyadic_prefactor_exact == Fraction(1, 256)
    assert cert.joint_finite_jet_target_arithmetic_verified
    assert cert.one_common_recursive_schedule_for_requested_jets

    prefix = cert.chain.chain.supported_prefix.prefix
    assert prefix.derivative_order == 1
    assert prefix.target.minimum_order == 4
    assert prefix.target.target_power == Fraction(1, 8)
    assert prefix.cutoff.bounds.max_order == 8
    assert prefix.cancellations.max_order == 8

    assert cert.finite_joint_jet_prefix_only
    assert not cert.uncut_tail_germ_verified
    assert not cert.all_order_hierarchy_verified
    assert not cert.infinite_diagonal_schedule_verified
    assert not cert.pde_residual_tail_verified
    assert not cert.all_jets_flat
    assert not cert.super_algebraic
    assert not cert.paper_exact


def test_minimum_order_can_dominate_without_weakening_the_target() -> None:
    cert = certify_finite_joint_jet_target_prefix(
        0.125,
        max_derivative_order=1,
        target_power=Fraction(1, 8),
        provider=_JointProvider(),
        minimum_order=12,
    )

    assert cert.selected_order == 12
    assert cert.ordinary_tail_powers_exact == (Fraction(13, 8), Fraction(5, 8))
    assert all(power >= cert.target_power for power in cert.ordinary_tail_powers_exact)


def test_joint_target_fails_closed_when_one_common_prefix_crosses_real_frontier() -> None:
    # The exact weakest-jet target selects J=8. A provider that stops at 7
    # must not return separate lower-order certificates for the easier jets.
    with pytest.raises(RuntimeError, match="real hierarchy frontier reached"):
        certify_finite_joint_jet_target_prefix(
            0.125,
            max_derivative_order=1,
            target_power=Fraction(1, 8),
            provider=_JointProvider(max_order=7),
            minimum_order=4,
        )


def test_float_target_cannot_cross_the_exact_joint_tail_gate() -> None:
    with pytest.raises(TypeError, match="floats are forbidden"):
        certify_finite_joint_jet_target_prefix(
            0.125,
            max_derivative_order=1,
            target_power=0.125,
            provider=_JointProvider(),
        )


def test_negative_derivative_budget_and_zero_target_fail_closed() -> None:
    with pytest.raises(ValueError, match="nonnegative integer"):
        certify_finite_joint_jet_target_prefix(
            0.125,
            max_derivative_order=-1,
            target_power=Fraction(1, 8),
            provider=_JointProvider(),
        )

    with pytest.raises(ValueError, match="positive"):
        certify_finite_joint_jet_target_prefix(
            0.125,
            max_derivative_order=0,
            target_power=0,
            provider=_JointProvider(),
        )
