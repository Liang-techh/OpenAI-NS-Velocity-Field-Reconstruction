"""Exact first-omitted-order power ledger for finite Section 5 physical tails.

This module closes one arithmetic seam between the finite retained-recurrence
certificate and the physical SlowBorel tail theorem.  The pinned formal source
keeps two distinct gains visible:

* a fixed positive stage ``j`` has its full slow power ``q^(2*h*j)``;
* summing the diagonal tail through prefix ``J`` spends half that slow gain to
  absorb coefficient/cutoff growth, giving ``q^(h*(J+1)-m)`` for an ordinary
  derivative of order ``m``;
* composition with the physical chart spends one additional q-power per
  physical derivative, giving ``q^(h*(J+1)-2*m)`` before weakening to a chosen
  target power ``P``.

The input is the already fail-closed finite physical-prefix/uncut-radius
certificate.  Consequently the truncation order, exact ``C[j,m]`` majorants,
common support, retained recurrence identities, recursive scale schedule, and
uncut window all come from one provider-owned chain.  There is no API for a
caller-supplied first-omitted exponent or coefficient table.

This is still finite theorem bookkeeping.  In particular, it does not bound
the coefficient multiplying the first omitted order, construct the existential
physical-chart constant ``D``, evaluate a PDE residual, materialize ``slowSum``,
or quantify over an infinite hierarchy.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .background_joint_physical_uncut_radius import (
    FiniteJointPhysicalUncutRadiusCertificate,
    certify_finite_joint_physical_uncut_radius,
)
from .background_target_driven_exact_majorants import TargetDrivenExactMajorantProvider


PINNED_FORMAL_REVISION = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_STAGE_POWER_THEOREM = "NavierStokes.SlowBorelBase.exists_stage_blown_power_bound"
PINNED_ORDINARY_TAIL_THEOREM = "NavierStokes.SlowBorelBase.ordinary_tail_bound"
PINNED_PHYSICAL_TAIL_THEOREM = "NavierStokes.SlowBorelBase.exists_physical_uncut_tail"


@dataclass(frozen=True)
class FirstOmittedPhysicalJetPower:
    """Exact power ledger for one requested physical derivative order."""

    derivative_order: int
    raw_first_omitted_power: Fraction
    diagonal_tail_power: Fraction
    ordinary_tail_power: Fraction
    physical_tail_power: Fraction
    target_power: Fraction

    def __post_init__(self) -> None:
        if isinstance(self.derivative_order, bool) or not isinstance(self.derivative_order, int):
            raise TypeError("derivative_order must be an integer")
        if self.derivative_order < 0:
            raise ValueError("derivative_order must be nonnegative")
        for name in (
            "raw_first_omitted_power",
            "diagonal_tail_power",
            "ordinary_tail_power",
            "physical_tail_power",
            "target_power",
        ):
            if not isinstance(getattr(self, name), Fraction):
                raise TypeError(f"{name} must be an exact Fraction")

        m = self.derivative_order
        if self.diagonal_tail_power * 2 != self.raw_first_omitted_power:
            raise ValueError("diagonal tail must spend exactly half the fixed-stage slow power")
        if self.ordinary_tail_power != self.diagonal_tail_power - m:
            raise ValueError("ordinary tail derivative loss must be exactly m q-powers")
        if self.physical_tail_power != self.ordinary_tail_power - m:
            raise ValueError("physical chart must spend exactly one further q-power per derivative")
        if self.physical_tail_power < self.target_power:
            raise ValueError("physical first-omitted power misses the requested target")

    @property
    def target_slack(self) -> Fraction:
        return self.physical_tail_power - self.target_power


@dataclass(frozen=True)
class FiniteFirstOmittedPhysicalTailCertificate:
    """Bind retained exact cancellation to the first omitted physical q-power."""

    uncut: FiniteJointPhysicalUncutRadiusCertificate
    formal_revision: str = PINNED_FORMAL_REVISION
    stage_power_theorem: str = PINNED_STAGE_POWER_THEOREM
    ordinary_tail_theorem: str = PINNED_ORDINARY_TAIL_THEOREM
    physical_tail_theorem: str = PINNED_PHYSICAL_TAIL_THEOREM

    def __post_init__(self) -> None:
        if not isinstance(self.uncut, FiniteJointPhysicalUncutRadiusCertificate):
            raise TypeError("uncut must be a FiniteJointPhysicalUncutRadiusCertificate")
        if self.formal_revision != PINNED_FORMAL_REVISION:
            raise ValueError("formal_revision does not match the pinned source")
        if self.stage_power_theorem != PINNED_STAGE_POWER_THEOREM:
            raise ValueError("stage_power_theorem does not match the pinned theorem")
        if self.ordinary_tail_theorem != PINNED_ORDINARY_TAIL_THEOREM:
            raise ValueError("ordinary_tail_theorem does not match the pinned theorem")
        if self.physical_tail_theorem != PINNED_PHYSICAL_TAIL_THEOREM:
            raise ValueError("physical_tail_theorem does not match the pinned theorem")

        physical = self.uncut.physical
        dependency_chain = physical.ordinary.chain.chain
        cancellations = dependency_chain.supported_prefix.prefix.cancellations
        J = physical.selected_order

        if dependency_chain.selected_order != J:
            raise ValueError("dependency chain ends at a different truncation order")
        if cancellations.max_order != J:
            raise ValueError("retained cancellation prefix ends at a different truncation order")
        if not dependency_chain.retained_recurrences_exact:
            raise ValueError("all retained recurrence identities must be exact zero")
        if len(dependency_chain.links) != J:
            raise ValueError("positive-order evidence links must cover every retained order 1..J")
        if physical.dyadic_prefactor_exact != Fraction(1, 2) ** J:
            raise ValueError("physical prefix lost the pinned dyadic tail prefactor")

        expected_raw = 2 * physical.h_exact * (J + 1)
        expected_diagonal = expected_raw / 2
        if expected_diagonal != physical.h_exact * (J + 1):
            raise RuntimeError("first-omitted half-gain arithmetic drifted")

        rows = self.jet_powers
        if len(rows) != physical.max_derivative_order + 1:
            raise RuntimeError("first-omitted physical ledger has the wrong jet budget")
        if tuple(row.ordinary_tail_power for row in rows) != physical.ordinary_tail_powers_exact:
            raise ValueError("ordinary tail powers drifted from the provider-owned physical prefix")
        if tuple(row.physical_tail_power for row in rows) != physical.physical_tail_powers_exact:
            raise ValueError("physical tail powers drifted from the provider-owned physical prefix")
        if any(row.target_power != physical.target_physical_power for row in rows):
            raise RuntimeError("physical target power drifted inside the first-omitted ledger")

        # The finite retained certificate must stop exactly before the omitted order.
        # We deliberately do not synthesize a recurrence identity for J+1.
        if self.first_omitted_order != cancellations.max_order + 1:
            raise RuntimeError("first omitted order is not immediately after the retained prefix")

    @property
    def truncation_order(self) -> int:
        return self.uncut.physical.selected_order

    @property
    def first_omitted_order(self) -> int:
        return self.truncation_order + 1

    @property
    def h_exact(self) -> Fraction:
        return self.uncut.physical.h_exact

    @property
    def raw_first_omitted_power_exact(self) -> Fraction:
        """Full fixed-stage slow power ``2*h*(J+1)`` before tail losses."""
        return 2 * self.h_exact * self.first_omitted_order

    @property
    def diagonal_tail_power_exact(self) -> Fraction:
        """Half-gain ``h*(J+1)`` retained by the diagonal tail estimate."""
        return self.raw_first_omitted_power_exact / 2

    @property
    def dyadic_prefactor_exact(self) -> Fraction:
        return self.uncut.physical.dyadic_prefactor_exact

    @property
    def uncut_radius_exact(self) -> Fraction:
        return self.uncut.uncut_radius_exact

    @property
    def jet_powers(self) -> tuple[FirstOmittedPhysicalJetPower, ...]:
        physical = self.uncut.physical
        raw = self.raw_first_omitted_power_exact
        diagonal = raw / 2
        return tuple(
            FirstOmittedPhysicalJetPower(
                derivative_order=m,
                raw_first_omitted_power=raw,
                diagonal_tail_power=diagonal,
                ordinary_tail_power=diagonal - m,
                physical_tail_power=diagonal - 2 * m,
                target_power=physical.target_physical_power,
            )
            for m in range(physical.max_derivative_order + 1)
        )

    @property
    def retained_orders_exactly_cancelled(self) -> tuple[int, ...]:
        return tuple(range(self.truncation_order + 1))

    @property
    def retained_recurrences_exact(self) -> bool:
        return True

    @property
    def first_omitted_order_exposed(self) -> bool:
        return True

    @property
    def first_omitted_coefficient_majorant_verified(self) -> bool:
        return False

    @property
    def first_omitted_recurrence_cancelled(self) -> bool:
        return False

    @property
    def physical_chart_finite_bound_verified(self) -> bool:
        return False

    @property
    def actual_pde_residual_verified(self) -> bool:
        return False

    @property
    def all_order_hierarchy_verified(self) -> bool:
        return False

    @property
    def infinite_diagonal_schedule_verified(self) -> bool:
        return False

    @property
    def all_jets_flat(self) -> bool:
        return False

    @property
    def super_algebraic(self) -> bool:
        return False

    @property
    def paper_exact(self) -> bool:
        return False


def certify_finite_first_omitted_physical_tail(
    h: float,
    max_derivative_order: int,
    target_physical_power: Fraction | int,
    provider: TargetDrivenExactMajorantProvider,
    *,
    minimum_order: int = 0,
    initial_lower_bound: int = 0,
) -> FiniteFirstOmittedPhysicalTailCertificate:
    """Build one provider-owned finite first-omitted physical-power ledger."""
    uncut = certify_finite_joint_physical_uncut_radius(
        h,
        max_derivative_order=max_derivative_order,
        target_physical_power=target_physical_power,
        provider=provider,
        minimum_order=minimum_order,
        initial_lower_bound=initial_lower_bound,
    )
    return FiniteFirstOmittedPhysicalTailCertificate(uncut=uncut)
