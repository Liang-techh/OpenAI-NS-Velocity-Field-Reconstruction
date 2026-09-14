"""Finite joint-jet target prefixes for the Section 5 SlowBorel diagonal.

The pinned theorem ``SlowBorelBase.exists_ordinary_uncut_tail`` has a stronger
finite quantifier shape than the single-derivative target selector already
exposed by this reconstruction: for one finite derivative budget ``M`` and one
target power ``P``, it chooses one sufficiently late prefix ``J`` that works for
every ordinary derivative order ``m <= M``.

This module replays only the exact order arithmetic of that common-prefix
choice on top of the existing hierarchy-owned exact-majorant chain.  It asks
the provider for one prefix using the weakest derivative exponent, ``m=M``,
and requires ``J >= M`` and ``J >= Jmin``.  Consequently

    h * (J + 1) - m >= h * (J + 1) - M >= P

for every ``m <= M``.  The same provider-owned ``C[j,m]`` rows, exact common
support, exact retained recurrences, and recursive SlowBorel/DiagonalScale
schedule are therefore shared by the whole requested finite jet budget.

This is deliberately *not* the function-level Lean theorem.  In particular it
does not materialize the infinite ``slowSum``, prove the
``cutPrefix_eventually_uncut`` germ/``delta`` conclusion, or quantify over all
``M`` and ``P``.  A finite successful request must not be promoted to
all-jets-flatness or super-algebraic residual convergence.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .background_target_driven_exact_majorants import (
    TargetDrivenExactMajorantChainCertificate,
    TargetDrivenExactMajorantProvider,
    certify_target_driven_exact_majorant_chain,
)


PINNED_FORMAL_REVISION = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_JOINT_TAIL_THEOREM = (
    "NavierStokes.SlowBorelBase.exists_ordinary_uncut_tail"
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
class FiniteJointJetTargetPrefixCertificate:
    """One provider-owned finite prefix shared by all ordinary jets ``m <= M``."""

    max_derivative_order: int
    minimum_order: int
    target_power: Fraction | int
    chain: TargetDrivenExactMajorantChainCertificate
    formal_revision: str = PINNED_FORMAL_REVISION
    theorem_name: str = PINNED_JOINT_TAIL_THEOREM

    def __post_init__(self) -> None:
        max_derivative_order = _nonnegative_int(
            self.max_derivative_order, "max_derivative_order"
        )
        minimum_order = _nonnegative_int(self.minimum_order, "minimum_order")
        target_power = _positive_exact(self.target_power, "target_power")
        if not isinstance(self.chain, TargetDrivenExactMajorantChainCertificate):
            raise TypeError(
                "chain must be a TargetDrivenExactMajorantChainCertificate"
            )
        if self.formal_revision != PINNED_FORMAL_REVISION:
            raise ValueError("formal_revision does not match the pinned source")
        if self.theorem_name != PINNED_JOINT_TAIL_THEOREM:
            raise ValueError("theorem_name does not match the pinned joint-tail theorem")

        prefix = self.chain.chain.supported_prefix.prefix
        target = prefix.target
        expected_minimum = max(1, minimum_order, max_derivative_order)
        if prefix.derivative_order != max_derivative_order:
            raise ValueError(
                "underlying target prefix must be selected at the weakest jet m=M"
            )
        if target.target_power != target_power:
            raise ValueError("underlying target prefix uses a different target power")
        if target.minimum_order != expected_minimum:
            raise ValueError(
                "underlying target prefix does not enforce J >= max(Jmin,M,1)"
            )
        if self.chain.selected_order < minimum_order:
            raise ValueError("selected prefix is earlier than Jmin")
        if self.chain.selected_order < max_derivative_order:
            raise ValueError("selected prefix must satisfy J >= M")
        if prefix.ordinary_tail_power_exact < target_power:
            raise ValueError("weakest requested ordinary jet does not reach target P")

        exact_h = target.h_exact
        joint_powers = tuple(
            exact_h * (self.chain.selected_order + 1) - m
            for m in range(max_derivative_order + 1)
        )
        if joint_powers[-1] != prefix.ordinary_tail_power_exact:
            raise RuntimeError("joint-jet weakest exponent drifted from the target prefix")
        if any(power < target_power for power in joint_powers):
            raise ValueError("at least one requested ordinary jet misses target P")

        object.__setattr__(self, "max_derivative_order", max_derivative_order)
        object.__setattr__(self, "minimum_order", minimum_order)
        object.__setattr__(self, "target_power", target_power)

    @property
    def selected_order(self) -> int:
        return self.chain.selected_order

    @property
    def h_exact(self) -> Fraction:
        return self.chain.chain.supported_prefix.prefix.target.h_exact

    @property
    def ordinary_tail_powers_exact(self) -> tuple[Fraction, ...]:
        """Return ``h*(J+1)-m`` for every requested ``0 <= m <= M``."""
        return tuple(
            self.h_exact * (self.selected_order + 1) - m
            for m in range(self.max_derivative_order + 1)
        )

    @property
    def weakest_tail_power_exact(self) -> Fraction:
        return self.ordinary_tail_powers_exact[-1]

    @property
    def dyadic_prefactor_exact(self) -> Fraction:
        """Return the pinned ordinary-tail prefactor ``2^-J`` exactly."""
        return Fraction(1, 2) ** self.selected_order

    @property
    def joint_finite_jet_target_arithmetic_verified(self) -> bool:
        return True

    @property
    def one_common_recursive_schedule_for_requested_jets(self) -> bool:
        return True

    @property
    def finite_joint_jet_prefix_only(self) -> bool:
        return True

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


def certify_finite_joint_jet_target_prefix(
    h: float,
    max_derivative_order: int,
    target_power: Fraction | int,
    provider: TargetDrivenExactMajorantProvider,
    *,
    minimum_order: int = 0,
    initial_lower_bound: int = 0,
) -> FiniteJointJetTargetPrefixCertificate:
    """Select one exact-majorant prefix for all ordinary jets through ``M``.

    The finite request is reduced to the weakest exponent ``m=M``.  All
    coefficient bounds/support/recurrences are still obtained only through the
    hierarchy provider by ``certify_target_driven_exact_majorant_chain``.
    Missing hierarchy depth therefore propagates as a hard failure; no partial
    joint-jet certificate is returned.
    """

    max_derivative_order = _nonnegative_int(
        max_derivative_order, "max_derivative_order"
    )
    minimum_order = _nonnegative_int(minimum_order, "minimum_order")
    target_power = _positive_exact(target_power, "target_power")

    common_minimum = max(minimum_order, max_derivative_order)
    chain = certify_target_driven_exact_majorant_chain(
        h,
        derivative_order=max_derivative_order,
        target_power=target_power,
        provider=provider,
        minimum_order=common_minimum,
        initial_lower_bound=initial_lower_bound,
    )
    return FiniteJointJetTargetPrefixCertificate(
        max_derivative_order=max_derivative_order,
        minimum_order=minimum_order,
        target_power=target_power,
        chain=chain,
    )
