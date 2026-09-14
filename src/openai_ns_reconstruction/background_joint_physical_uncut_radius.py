"""Constructive uncut windows for finite joint physical SlowBorel tails.

The recursive SlowBorel/DiagonalScale schedule already carried by the Agent-7
finite physical-jet certificate is a positive doubling sequence ``a_j``.  The
pinned formal theorem ``GenericRealizationBounds.prefix_eventually_uncut``
prints an explicit plateau condition for a monotone positive schedule:

    0 < q < 1 / (2 * a_J).

On that interval every retained cutoff through stage ``J`` has
``|a_j q| < 1/2`` and is therefore exactly one.  This module computes that
radius as an exact ``Fraction`` from the provider-owned recursive schedule; it
does not accept a caller cutoff radius or a second schedule.

This closes only the finite theorem-hypothesis arithmetic.  Python does not
replay the Lean eventual-equality proof, materialize ``slowSum``, construct the
physical-chart finite constant, or prove an actual PDE residual estimate.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .background_joint_physical_jet_target_prefix import (
    FiniteJointPhysicalJetTargetPrefixCertificate,
    certify_finite_joint_physical_jet_target_prefix,
)
from .background_target_driven_exact_majorants import TargetDrivenExactMajorantProvider


PINNED_FORMAL_REVISION = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_UNCUT_THEOREM = "NavierStokes.GenericRealizationBounds.prefix_eventually_uncut"
PINNED_TAIL_THEOREM = "NavierStokes.GenericRealizationBounds.realized_tail_bound"
PINNED_CUTOFF_PLATEAU = "NavierStokes.SmoothCutoffs.scaledCutoff_one_of_abs_le"
PINNED_COORDINATE = "Prod.fst"


def _positive_exact(value: Fraction | int, name: str) -> Fraction:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact positive integer/Fraction")
    if isinstance(value, Fraction):
        result = value
    elif isinstance(value, Integral):
        result = Fraction(int(value), 1)
    else:
        raise TypeError(
            f"{name} must be an exact positive integer/Fraction; floats are forbidden"
        )
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


@dataclass(frozen=True)
class FiniteJointPhysicalUncutRadiusCertificate:
    """One exact uncut plateau window for a provider-owned finite prefix."""

    physical: FiniteJointPhysicalJetTargetPrefixCertificate
    formal_revision: str = PINNED_FORMAL_REVISION
    uncut_theorem: str = PINNED_UNCUT_THEOREM
    tail_theorem: str = PINNED_TAIL_THEOREM
    cutoff_plateau_theorem: str = PINNED_CUTOFF_PLATEAU
    coordinate_map: str = PINNED_COORDINATE

    def __post_init__(self) -> None:
        if not isinstance(self.physical, FiniteJointPhysicalJetTargetPrefixCertificate):
            raise TypeError(
                "physical must be a FiniteJointPhysicalJetTargetPrefixCertificate"
            )
        if self.formal_revision != PINNED_FORMAL_REVISION:
            raise ValueError("formal_revision does not match the pinned source")
        if self.uncut_theorem != PINNED_UNCUT_THEOREM:
            raise ValueError("uncut_theorem does not match the pinned formal theorem")
        if self.tail_theorem != PINNED_TAIL_THEOREM:
            raise ValueError("tail_theorem does not match the pinned formal theorem")
        if self.cutoff_plateau_theorem != PINNED_CUTOFF_PLATEAU:
            raise ValueError("cutoff plateau theorem does not match the pinned source")
        if self.coordinate_map != PINNED_COORDINATE:
            raise ValueError("the SlowBorel chart coordinate must remain Prod.fst")

        schedule = self.schedule
        if schedule.max_order != self.physical.selected_order:
            raise ValueError("recursive schedule does not end at the selected prefix")
        if len(schedule.scales) != self.physical.selected_order + 1:
            raise ValueError("recursive schedule has the wrong finite-prefix length")
        if any(scale <= 0 for scale in schedule.scales):
            raise ValueError("pinned uncut theorem requires positive cutoff scales")
        if any(
            schedule.scales[j + 1] < schedule.scales[j]
            for j in range(len(schedule.scales) - 1)
        ):
            raise ValueError("pinned uncut theorem requires a monotone cutoff schedule")
        if any(
            schedule.scales[j + 1] < 2 * schedule.scales[j]
            for j in range(len(schedule.scales) - 1)
        ):
            raise ValueError("recursive SlowBorel schedule lost its doubling envelope")

        radius = self.uncut_radius_exact
        if radius <= 0 or radius > Fraction(1, 2):
            raise RuntimeError("exact uncut radius is inconsistent with a_J >= 1")
        half = Fraction(1, 2)
        final_scale = self.selected_scale
        if Fraction(final_scale, 1) * radius != half:
            raise RuntimeError("exact plateau radius arithmetic drifted")
        if any(Fraction(scale, 1) * radius > half for scale in schedule.scales):
            raise RuntimeError("an earlier retained scale exceeds the plateau edge")

    @property
    def schedule(self):
        prefix = self.physical.ordinary.chain.chain.supported_prefix.prefix
        return prefix.cutoff.schedule

    @property
    def selected_order(self) -> int:
        return self.physical.selected_order

    @property
    def selected_scale(self) -> int:
        return self.schedule.scales[self.selected_order]

    @property
    def uncut_radius_exact(self) -> Fraction:
        """Return the exact open radius ``1/(2*a_J)`` from the real schedule."""
        return Fraction(1, 2 * self.selected_scale)

    @property
    def retained_plateau_edges_exact(self) -> tuple[Fraction, ...]:
        """Return exact ``a_j/(2*a_J) <= 1/2`` edge products for ``j<=J``."""
        radius = self.uncut_radius_exact
        return tuple(Fraction(scale, 1) * radius for scale in self.schedule.scales)

    def certify_q(self, q: Fraction | int) -> tuple[Fraction, ...]:
        """Check one exact ``0<q<1/(2*a_J)`` point against every retained stage.

        The returned products are exact theorem-side cutoff arguments ``a_j*q``.
        This is not a numerical evaluation of the cutoff function.
        """
        q_exact = _positive_exact(q, "q")
        if q_exact >= self.uncut_radius_exact:
            raise ValueError("q must lie strictly inside the exact uncut plateau")
        products = tuple(
            Fraction(scale, 1) * q_exact for scale in self.schedule.scales
        )
        if any(product >= Fraction(1, 2) for product in products):
            raise RuntimeError("recursive prefix left the pinned cutoff plateau")
        return products

    @property
    def finite_uncut_plateau_arithmetic_verified(self) -> bool:
        return True

    @property
    def one_provider_owned_schedule_used(self) -> bool:
        return True

    @property
    def function_level_uncut_germ_verified(self) -> bool:
        return False

    @property
    def physical_chart_finite_bound_verified(self) -> bool:
        return False

    @property
    def actual_slow_sum_materialized(self) -> bool:
        return False

    @property
    def pde_residual_tail_verified(self) -> bool:
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


def certify_finite_joint_physical_uncut_radius(
    h: float,
    max_derivative_order: int,
    target_physical_power: Fraction | int,
    provider: TargetDrivenExactMajorantProvider,
    *,
    minimum_order: int = 0,
    initial_lower_bound: int = 0,
) -> FiniteJointPhysicalUncutRadiusCertificate:
    """Select the physical prefix and derive its exact theorem-ready uncut radius."""
    physical = certify_finite_joint_physical_jet_target_prefix(
        h,
        max_derivative_order=max_derivative_order,
        target_physical_power=target_physical_power,
        provider=provider,
        minimum_order=minimum_order,
        initial_lower_bound=initial_lower_bound,
    )
    return FiniteJointPhysicalUncutRadiusCertificate(physical=physical)
