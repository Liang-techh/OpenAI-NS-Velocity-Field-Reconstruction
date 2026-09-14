"""One exact Section 9 prefix for a finite ordinary-jet residual budget.

The pinned theorem ``ActualIterationLedger.eventually_residualRate_ge`` is
stated for one fixed derivative order at a time.  A finite ``C^M`` residual
check needs one *common* truncation index that works simultaneously for every
ordinary derivative order ``0 <= m <= M``.  This module constructs that common
index from the exact single-order selector without assuming that the loss is
monotone in ``m``.

Only finite exponent arithmetic is certified here.  The constructor does not
materialize correction fields, evaluate a Navier--Stokes residual, quantify
over all derivative orders, or prove convergence of the Eq. (9.21) summed
local field.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from numbers import Integral

from .section9_residual_rate_target import (
    PINNED_EVENTUALLY_GE_SYMBOL,
    PINNED_FORMAL_COMMIT,
    PINNED_FORMAL_REPOSITORY,
    PINNED_RESIDUAL_RATE_SYMBOL,
    Section9ResidualRateTargetCertificate,
    residual_rate,
    select_residual_rate_stage,
)


PINNED_LEDGER_FILE = "NavierStokes/ActualIterationLedger.lean"


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be a nonnegative integer")
    result = int(value)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


@dataclass(frozen=True)
class Section9FiniteJetResidualTargetCertificate:
    """Select one minimal prefix for all residual jets ``m <= M``.

    Each individual derivative target is first solved by the exact pinned
    affine ``residualRate`` formula.  The common stage is the maximum of those
    individually minimal stages.  We then check every derivative again at the
    common stage and check global minimality directly at ``J-1``.  No
    monotonicity in derivative order is assumed; this matters because the
    theorem-side parameter ``beta`` is kept as exact input rather than being
    silently given an extra sign hypothesis.
    """

    h: Fraction | int
    beta: Fraction | int
    max_derivative: int
    target: Fraction | int

    stage: int = field(init=False)
    per_order: tuple[Section9ResidualRateTargetCertificate, ...] = field(init=False)

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT

    def __post_init__(self) -> None:
        M = _natural(self.max_derivative, "max_derivative")
        if self.formal_repository != PINNED_FORMAL_REPOSITORY:
            raise ValueError("formal_repository must match the pinned source exactly")
        if self.formal_commit != PINNED_FORMAL_COMMIT:
            raise ValueError("formal_commit must match the pinned source exactly")

        rows = tuple(
            select_residual_rate_stage(
                h=self.h,
                beta=self.beta,
                derivative_order=m,
                target=self.target,
            )
            for m in range(M + 1)
        )
        # The single-order certificates normalize every scalar to Fraction and
        # already reject float/Decimal theorem inputs.
        canonical = rows[0]
        if any(
            row.h != canonical.h or row.beta != canonical.beta or row.target != canonical.target
            for row in rows
        ):
            raise AssertionError("finite derivative rows do not share one exact scalar budget")

        selected = max(row.stage for row in rows)
        object.__setattr__(self, "h", canonical.h)
        object.__setattr__(self, "beta", canonical.beta)
        object.__setattr__(self, "target", canonical.target)
        object.__setattr__(self, "max_derivative", M)
        object.__setattr__(self, "per_order", rows)
        object.__setattr__(self, "stage", selected)

        if any(value < self.target for value in self.selected_rates):
            raise AssertionError("common stage does not reach the target for every derivative order")
        if selected > 0 and not any(value < self.target for value in self.previous_rates):
            raise AssertionError("common stage is not the minimal simultaneous target-reaching stage")
        if selected > 0 and any(
            after - before != self.h / 10
            for before, after in zip(self.previous_rates, self.selected_rates)
        ):
            raise AssertionError("pinned residual-rate successor gain is not exactly h/10")

    @property
    def derivative_orders(self) -> tuple[int, ...]:
        return tuple(range(self.max_derivative + 1))

    @property
    def per_order_minimal_stages(self) -> tuple[int, ...]:
        return tuple(row.stage for row in self.per_order)

    @property
    def selected_rates(self) -> tuple[Fraction, ...]:
        return tuple(
            residual_rate(self.h, self.beta, self.stage, m)
            for m in self.derivative_orders
        )

    @property
    def previous_rates(self) -> tuple[Fraction, ...]:
        if self.stage == 0:
            return tuple()
        return tuple(
            residual_rate(self.h, self.beta, self.stage - 1, m)
            for m in self.derivative_orders
        )

    @property
    def active_orders(self) -> tuple[int, ...]:
        """Orders whose individual minimal stage equals the common stage."""
        return tuple(m for m, row in enumerate(self.per_order) if row.stage == self.stage)

    @property
    def joint_target_margin(self) -> Fraction:
        return min(rate - self.target for rate in self.selected_rates)

    @property
    def one_common_prefix_machine_checked(self) -> bool:
        return True

    @property
    def finite_derivative_budget_only(self) -> bool:
        return True

    @property
    def actual_stage_field_consumed(self) -> bool:
        return False

    @property
    def actual_residual_evaluated(self) -> bool:
        return False

    @property
    def all_derivative_orders_certified(self) -> bool:
        return False

    @property
    def infinite_correction_sequence_certified(self) -> bool:
        return False

    @property
    def eq_9_21_summed_field_certified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False


def select_finite_jet_residual_stage(
    *,
    h: Fraction | int,
    beta: Fraction | int,
    max_derivative: int,
    target: Fraction | int,
) -> Section9FiniteJetResidualTargetCertificate:
    """Return the minimal common stage reaching ``target`` for all ``m <= M``."""
    return Section9FiniteJetResidualTargetCertificate(
        h=h,
        beta=beta,
        max_derivative=max_derivative,
        target=target,
    )


__all__ = [
    "PINNED_EVENTUALLY_GE_SYMBOL",
    "PINNED_RESIDUAL_RATE_SYMBOL",
    "PINNED_LEDGER_FILE",
    "Section9FiniteJetResidualTargetCertificate",
    "select_finite_jet_residual_stage",
]
