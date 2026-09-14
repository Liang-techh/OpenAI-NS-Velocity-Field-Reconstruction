"""Finite canonical target ladders for Section 5 physical SlowBorel tails.

The pinned all-order argument quantifies over every finite physical derivative
order and every requested decay power.  The preceding prefix-coherence layer
checks only one weaker/stronger request pair.  This module takes the next
quantifier step without pretending that a finite hierarchy is infinite: it
uses the canonical exact positive-target family

    level k >= 1 := derivatives 0..k, physical q-power at least k,

builds the first finitely many levels from one hierarchy provider, and requires
every adjacent level to extend the same coefficient evidence and the same
recursive SlowBorel/DiagonalScale schedule literally.

The family is cofinal for exact rational targets: a request (M,P) is dominated
by k = max(M, ceil(P), 1).  A finite ladder can discharge only requests whose
canonical level is already materialized.  Requests beyond that frontier fail
closed; no provider totality, infinite schedule, actual PDE residual, all-jets
flatness, or super-algebraic convergence is inferred from the finite ladder.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .background_physical_tail_prefactor import (
    FinitePhysicalTailPrefactorCertificate,
    certify_finite_physical_tail_prefactor,
)
from .background_physical_tail_prefix_coherence import (
    FiniteIncreasingPhysicalTailCertificate,
    PINNED_ADMISSIBLE_SCALES,
    PINNED_ENVELOPE,
    PINNED_FORMAL_REVISION,
    PINNED_PHYSICAL_TAIL,
    PINNED_SCALE_THEOREM,
)
from .background_target_driven_exact_majorants import TargetDrivenExactMajorantProvider


def _nonnegative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _positive_int(value: int, name: str) -> int:
    result = _nonnegative_int(value, name)
    if result == 0:
        raise ValueError(f"{name} must be positive")
    return result


def _exact_fraction(value: Fraction | int, name: str) -> Fraction:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact integer/Fraction")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value), 1)
    raise TypeError(f"{name} must be an exact integer/Fraction")


def _ceil_fraction(value: Fraction) -> int:
    """Exact integer ceiling, including negative rationals."""
    return -((-value.numerator) // value.denominator)


def canonical_level_for_physical_target(
    max_derivative_order: int,
    target_physical_power: Fraction | int,
) -> int:
    """Return the least positive canonical level dominating one exact target.

    Canonical level ``k>=1`` asks for all physical derivatives through ``k``
    and positive power ``k``.  Therefore ``k=max(M,ceil(P),1)`` is exactly the
    least level in this family dominating ``(M,P)``.  Using a positive level is
    deliberate: the landed executable physical-tail theorem correctly rejects
    nonpositive target powers.
    """

    M = _nonnegative_int(max_derivative_order, "max_derivative_order")
    P = _exact_fraction(target_physical_power, "target_physical_power")
    return max(M, _ceil_fraction(P), 1)


@dataclass(frozen=True)
class FiniteCofinalPhysicalTailLadderCertificate:
    """A finite prefix of the canonical positive cofinal physical-target family."""

    levels: tuple[FinitePhysicalTailPrefactorCertificate, ...]
    max_level: int
    formal_revision: str = PINNED_FORMAL_REVISION
    scale_theorem: str = PINNED_SCALE_THEOREM
    doubling_envelope: str = PINNED_ENVELOPE
    admissible_scales_theorem: str = PINNED_ADMISSIBLE_SCALES
    physical_tail_theorem: str = PINNED_PHYSICAL_TAIL

    def __post_init__(self) -> None:
        K = _positive_int(self.max_level, "max_level")
        if len(self.levels) != K:
            raise ValueError("finite canonical ladder must contain exactly levels 1..max_level")
        if self.formal_revision != PINNED_FORMAL_REVISION:
            raise ValueError("formal_revision does not match the pinned source")
        if self.scale_theorem != PINNED_SCALE_THEOREM:
            raise ValueError("scale_theorem does not match exists_diagonal_scales")
        if self.doubling_envelope != PINNED_ENVELOPE:
            raise ValueError("doubling_envelope does not match the pinned recursive envelope")
        if self.admissible_scales_theorem != PINNED_ADMISSIBLE_SCALES:
            raise ValueError("admissible_scales_theorem does not match the pinned source")
        if self.physical_tail_theorem != PINNED_PHYSICAL_TAIL:
            raise ValueError("physical_tail_theorem does not match the pinned source")

        previous: FinitePhysicalTailPrefactorCertificate | None = None
        common_minimum_order: int | None = None
        common_initial_lower_bound: int | None = None
        common_h: Fraction | None = None

        for k, certificate in enumerate(self.levels, start=1):
            if not isinstance(certificate, FinitePhysicalTailPrefactorCertificate):
                raise TypeError("every ladder level must be a FinitePhysicalTailPrefactorCertificate")

            physical = certificate.first_omitted.uncut.physical
            if physical.max_derivative_order != k:
                raise ValueError("ladder derivative budget is not the canonical level index")
            if physical.target_physical_power != Fraction(k, 1):
                raise ValueError("ladder physical target power is not the canonical level index")

            h_exact = certificate.first_omitted.h_exact
            if common_h is None:
                common_h = h_exact
            elif h_exact != common_h:
                raise ValueError("finite canonical ladder changed h across levels")

            if common_minimum_order is None:
                common_minimum_order = physical.minimum_order
            elif physical.minimum_order != common_minimum_order:
                raise ValueError("finite canonical ladder changed the minimum truncation order")

            schedule = certificate.first_omitted.uncut.schedule
            if common_initial_lower_bound is None:
                common_initial_lower_bound = schedule.initial_lower_bound
            elif schedule.initial_lower_bound != common_initial_lower_bound:
                raise ValueError("finite canonical ladder changed the initial scale lower bound")

            if len(certificate.first_omitted.jet_powers) != k + 1:
                raise RuntimeError("canonical level lost a requested physical derivative row")
            for m, row in enumerate(certificate.first_omitted.jet_powers):
                if row.derivative_order != m:
                    raise RuntimeError("canonical physical derivative rows are misaligned")
                if row.physical_tail_power < k:
                    raise ValueError("canonical level does not reach its advertised physical power")

            if previous is not None:
                # Reuse the strict two-request checker instead of weakening any
                # premise for the multi-level family.
                FiniteIncreasingPhysicalTailCertificate(earlier=previous, later=certificate)
            previous = certificate

    @property
    def truncation_orders(self) -> tuple[int, ...]:
        return tuple(level.truncation_order for level in self.levels)

    def required_level(
        self,
        max_derivative_order: int,
        target_physical_power: Fraction | int,
    ) -> int:
        return canonical_level_for_physical_target(max_derivative_order, target_physical_power)

    def certificate_for_target(
        self,
        max_derivative_order: int,
        target_physical_power: Fraction | int,
    ) -> FinitePhysicalTailPrefactorCertificate:
        """Return a materialized dominating level or fail at the finite frontier."""

        M = _nonnegative_int(max_derivative_order, "max_derivative_order")
        P = _exact_fraction(target_physical_power, "target_physical_power")
        k = canonical_level_for_physical_target(M, P)
        if k > self.max_level:
            raise ValueError(
                "requested physical target lies beyond the materialized finite canonical ladder"
            )

        certificate = self.levels[k - 1]
        rows = certificate.first_omitted.jet_powers
        if len(rows) <= M:
            raise RuntimeError("selected canonical level does not contain the requested derivative budget")
        for row in rows[: M + 1]:
            if row.physical_tail_power < P:
                raise RuntimeError("selected canonical level does not dominate the requested physical power")
        return certificate

    @property
    def canonical_target_cofinality_arithmetic_verified(self) -> bool:
        """The exact selector k=max(M,ceil(P),1) dominates every rational target."""
        return True

    @property
    def finite_shared_schedule_ladder_verified(self) -> bool:
        return True

    @property
    def finite_ladder_only(self) -> bool:
        return True

    @property
    def provider_totality_verified(self) -> bool:
        return False

    @property
    def infinite_diagonal_schedule_verified(self) -> bool:
        return False

    @property
    def actual_pde_residual_verified(self) -> bool:
        return False

    @property
    def all_order_hierarchy_verified(self) -> bool:
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


def certify_finite_cofinal_physical_tail_ladder(
    h: float,
    *,
    max_level: int,
    provider: TargetDrivenExactMajorantProvider,
    minimum_order: int = 0,
    initial_lower_bound: int = 0,
) -> FiniteCofinalPhysicalTailLadderCertificate:
    """Materialize levels ``1..max_level`` from one provider and one schedule family.

    Each level is constructed once.  Adjacent levels are then checked with the
    existing strict prefix-coherence certificate, so a stateful provider that
    changes any previously emitted coefficient/bound/support/recurrence row on
    a later request fails closed.
    """

    K = _positive_int(max_level, "max_level")
    levels: list[FinitePhysicalTailPrefactorCertificate] = []
    for k in range(1, K + 1):
        certificate = certify_finite_physical_tail_prefactor(
            h,
            max_derivative_order=k,
            target_physical_power=Fraction(k, 1),
            provider=provider,
            minimum_order=minimum_order,
            initial_lower_bound=initial_lower_bound,
        )
        levels.append(certificate)

    return FiniteCofinalPhysicalTailLadderCertificate(
        levels=tuple(levels),
        max_level=K,
    )
