"""Finite joint physical-jet target prefixes for Section 5 SlowBorel tails.

The pinned theorem ``SlowBorelBase.exists_physical_uncut_tail`` reaches a
physical q-power ``P`` for one derivative order ``m`` by first invoking the
ordinary uncut-tail theorem at the stronger target ``P + m``.  The physical
chart chain rule then spends one further factor ``q^{-m}``.

For a finite physical derivative budget ``0 <= m <= M`` it is therefore enough
to choose one provider-owned ordinary prefix for the strongest request
``m=M`` and ordinary target ``P+M``.  The existing joint ordinary selector then
gives, for every requested ``m``,

    h * (J + 1) - m >= P + M,

hence after the physical chart loss

    h * (J + 1) - 2*m >= P.

This module records only that exact finite order arithmetic on top of the
hierarchy-owned exact-majorant/support/retained-recurrence chain.  It does not
construct the theorem's physical-chart derivative constant, the local uncut
germ/delta, an infinite hierarchy, or an actual PDE residual estimate.  A
finite successful request must not be promoted to all-jets flatness or
super-algebraic residual convergence.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .background_joint_jet_target_prefix import (
    FiniteJointJetTargetPrefixCertificate,
    certify_finite_joint_jet_target_prefix,
)
from .background_target_driven_exact_majorants import (
    TargetDrivenExactMajorantProvider,
)


PINNED_FORMAL_REVISION = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_PHYSICAL_TAIL_THEOREM = (
    "NavierStokes.SlowBorelBase.exists_physical_uncut_tail"
)


def _nonnegative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


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
class FiniteJointPhysicalJetTargetPrefixCertificate:
    """One finite provider-owned prefix shared by all physical jets ``m <= M``."""

    max_derivative_order: int
    minimum_order: int
    target_physical_power: Fraction | int
    ordinary: FiniteJointJetTargetPrefixCertificate
    formal_revision: str = PINNED_FORMAL_REVISION
    theorem_name: str = PINNED_PHYSICAL_TAIL_THEOREM

    def __post_init__(self) -> None:
        max_derivative_order = _nonnegative_int(
            self.max_derivative_order, "max_derivative_order"
        )
        minimum_order = _nonnegative_int(self.minimum_order, "minimum_order")
        target_physical_power = _positive_exact(
            self.target_physical_power, "target_physical_power"
        )
        if not isinstance(self.ordinary, FiniteJointJetTargetPrefixCertificate):
            raise TypeError(
                "ordinary must be a FiniteJointJetTargetPrefixCertificate"
            )
        if self.formal_revision != PINNED_FORMAL_REVISION:
            raise ValueError("formal_revision does not match the pinned source")
        if self.theorem_name != PINNED_PHYSICAL_TAIL_THEOREM:
            raise ValueError(
                "theorem_name does not match the pinned physical-tail theorem"
            )
        if self.ordinary.max_derivative_order != max_derivative_order:
            raise ValueError("ordinary prefix uses a different derivative budget")
        if self.ordinary.minimum_order != minimum_order:
            raise ValueError("ordinary prefix uses a different minimum order")

        expected_ordinary_target = target_physical_power + max_derivative_order
        if self.ordinary.target_power != expected_ordinary_target:
            raise ValueError(
                "physical tail must use the pinned ordinary target P+M"
            )

        h_exact = self.ordinary.h_exact
        if not Fraction(0, 1) < h_exact < Fraction(1, 2):
            raise ValueError(
                "pinned physical-tail theorem requires the exact range 0 < h < 1/2"
            )

        ordinary_powers = self.ordinary.ordinary_tail_powers_exact
        physical_powers = tuple(
            power - derivative_order
            for derivative_order, power in enumerate(ordinary_powers)
        )
        expected_physical_powers = tuple(
            h_exact * (self.ordinary.selected_order + 1) - 2 * derivative_order
            for derivative_order in range(max_derivative_order + 1)
        )
        if physical_powers != expected_physical_powers:
            raise RuntimeError("physical-chart derivative-loss arithmetic drifted")
        if any(power < target_physical_power for power in physical_powers):
            raise ValueError(
                "at least one requested physical jet misses target power P"
            )

        object.__setattr__(self, "max_derivative_order", max_derivative_order)
        object.__setattr__(self, "minimum_order", minimum_order)
        object.__setattr__(self, "target_physical_power", target_physical_power)

    @property
    def selected_order(self) -> int:
        return self.ordinary.selected_order

    @property
    def h_exact(self) -> Fraction:
        return self.ordinary.h_exact

    @property
    def ordinary_target_power_exact(self) -> Fraction:
        return self.ordinary.target_power

    @property
    def ordinary_tail_powers_exact(self) -> tuple[Fraction, ...]:
        return self.ordinary.ordinary_tail_powers_exact

    @property
    def physical_tail_powers_exact(self) -> tuple[Fraction, ...]:
        """Return ``h*(J+1)-2m`` for every requested physical derivative."""
        return tuple(
            power - derivative_order
            for derivative_order, power in enumerate(self.ordinary_tail_powers_exact)
        )

    @property
    def weakest_physical_tail_power_exact(self) -> Fraction:
        return self.physical_tail_powers_exact[-1]

    @property
    def dyadic_prefactor_exact(self) -> Fraction:
        return self.ordinary.dyadic_prefactor_exact

    @property
    def joint_finite_physical_jet_target_arithmetic_verified(self) -> bool:
        return True

    @property
    def one_common_recursive_schedule_for_requested_physical_jets(self) -> bool:
        return True

    @property
    def retained_recurrences_exact(self) -> bool:
        return self.ordinary.chain.chain.retained_recurrences_exact

    @property
    def finite_physical_jet_prefix_only(self) -> bool:
        return True

    @property
    def physical_chart_finite_bound_verified(self) -> bool:
        return False

    @property
    def uncut_tail_germ_verified(self) -> bool:
        return False

    @property
    def all_order_hierarchy_verified(self) -> bool:
        return False

    @property
    def infinite_diagonal_schedule_verified(self) -> bool:
        return False

    @property
    def pde_residual_tail_verified(self) -> bool:
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


def certify_finite_joint_physical_jet_target_prefix(
    h: float,
    max_derivative_order: int,
    target_physical_power: Fraction | int,
    provider: TargetDrivenExactMajorantProvider,
    *,
    minimum_order: int = 0,
    initial_lower_bound: int = 0,
) -> FiniteJointPhysicalJetTargetPrefixCertificate:
    """Select one exact-majorant prefix for all physical jets through ``M``.

    The pinned physical theorem asks the ordinary tail for target ``P+m`` for
    derivative order ``m``.  One common finite schedule for every ``m<=M`` is
    obtained by the stronger ordinary request ``P+M`` at the same finite
    derivative budget ``M``.  Missing hierarchy depth remains a hard failure.
    """

    max_derivative_order = _nonnegative_int(
        max_derivative_order, "max_derivative_order"
    )
    minimum_order = _nonnegative_int(minimum_order, "minimum_order")
    target_physical_power = _positive_exact(
        target_physical_power, "target_physical_power"
    )
    ordinary_target = target_physical_power + max_derivative_order
    ordinary = certify_finite_joint_jet_target_prefix(
        h,
        max_derivative_order=max_derivative_order,
        target_power=ordinary_target,
        provider=provider,
        minimum_order=minimum_order,
        initial_lower_bound=initial_lower_bound,
    )
    return FiniteJointPhysicalJetTargetPrefixCertificate(
        max_derivative_order=max_derivative_order,
        minimum_order=minimum_order,
        target_physical_power=target_physical_power,
        ordinary=ordinary,
    )
