"""Constructive exact stage selection for pinned Section 9 residual rates.

``ActualIterationLedger.lean`` proves that, for every fixed derivative order,
``residualRate h beta J m`` tends to ``+infinity`` when ``h > 0``.  The formal
proof is existential/eventual.  For exact rational inputs this module makes the
same arithmetic constructive: it returns the *smallest* finite prefix index
``J`` whose pinned residual-rate exponent reaches a requested target.

This is deliberately only exponent/rate arithmetic.  It does not manufacture a
correction stage, evaluate a Navier--Stokes residual, or turn a finite-prefix
bound into convergence of Eq. (9.21).  Float/Decimal inputs are rejected so a
rounded near-threshold value cannot change the selected stage.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from numbers import Integral

from .section9_iteration_ledger import (
    PINNED_FORMAL_COMMIT,
    PINNED_FORMAL_REPOSITORY,
    residual_wave,
)


PINNED_LEDGER_FILE = "NavierStokes/ActualIterationLedger.lean"
PINNED_GRAPH_FILE = "NavierStokes/PhysicalGraphBounds.lean"
PINNED_COORDINATE_FILE = "NavierStokes/CoordinateAlgebra.lean"
PINNED_RESIDUAL_LOSS_SYMBOL = "NavierStokes.ActualIterationLedger.residualLoss"
PINNED_RESIDUAL_RATE_SYMBOL = "NavierStokes.ActualIterationLedger.residualRate"
PINNED_TENDSTO_SYMBOL = "NavierStokes.ActualIterationLedger.residualRate_tendsto_atTop"
PINNED_EVENTUALLY_GE_SYMBOL = "NavierStokes.ActualIterationLedger.eventually_residualRate_ge"
PINNED_GRAPH_LOSS_SYMBOL = "NavierStokes.PhysicalGraphBounds.graphLoss"
PINNED_A_SYMBOL = "NavierStokes.CoordinateAlgebra.A"


def _exact_fraction(value: object, name: str) -> Fraction:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact integer/Fraction, not bool")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value), 1)
    raise TypeError(
        f"{name} must be an exact integer/Fraction; float/Decimal approximations are rejected"
    )


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be a nonnegative integer")
    result = int(value)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _ceil_fraction(value: Fraction) -> int:
    """Exact integer ceiling, without converting through binary floating point."""
    return -((-value.numerator) // value.denominator)


def coordinate_A(h: Fraction | int) -> Fraction:
    """Pinned ``CoordinateAlgebra.A h = 1/2 + h``."""
    hq = _exact_fraction(h, "h")
    return Fraction(1, 2) + hq


def graph_loss(derivative_order: int) -> Fraction:
    """Pinned graph-composition loss ``m * (m + 2)``."""
    m = _natural(derivative_order, "derivative_order")
    return Fraction(m * (m + 2), 1)


def residual_loss(
    h: Fraction | int,
    beta: Fraction | int,
    derivative_order: int,
) -> Fraction:
    """Replay ``ActualIterationLedger.residualLoss`` exactly."""
    hq = _exact_fraction(h, "h")
    betaq = _exact_fraction(beta, "beta")
    m = _natural(derivative_order, "derivative_order")
    return (
        graph_loss(m)
        + 1
        + (2 * coordinate_A(hq) + Fraction(1, 2))
        + betaq * m
    )


def residual_rate(
    h: Fraction | int,
    beta: Fraction | int,
    stage: int,
    derivative_order: int,
) -> Fraction:
    """Replay ``h * residualWave(stage) - residualLoss(h,beta,m)`` exactly."""
    hq = _exact_fraction(h, "h")
    betaq = _exact_fraction(beta, "beta")
    j = _natural(stage, "stage")
    m = _natural(derivative_order, "derivative_order")
    return hq * residual_wave(j) - residual_loss(hq, betaq, m)


@dataclass(frozen=True)
class Section9ResidualRateTargetCertificate:
    """Minimal exact finite-prefix index reaching one fixed derivative target.

    The pinned theorem ``eventually_residualRate_ge`` quantifies over every
    target real number.  For rational ``h``, ``beta`` and target, the affine
    formula can be solved exactly.  ``stage`` is computed internally and is
    then checked both for sufficiency and (unless it is zero) minimality.
    """

    h: Fraction | int
    beta: Fraction | int
    derivative_order: int
    target: Fraction | int
    stage: int = field(init=False)

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT

    def __post_init__(self) -> None:
        hq = _exact_fraction(self.h, "h")
        betaq = _exact_fraction(self.beta, "beta")
        targetq = _exact_fraction(self.target, "target")
        m = _natural(self.derivative_order, "derivative_order")
        if hq <= 0:
            raise ValueError("h must be strictly positive, matching residualRate_tendsto_atTop")
        if self.formal_repository != PINNED_FORMAL_REPOSITORY:
            raise ValueError("formal_repository must match the pinned source exactly")
        if self.formal_commit != PINNED_FORMAL_COMMIT:
            raise ValueError("formal_commit must match the pinned source exactly")

        object.__setattr__(self, "h", hq)
        object.__setattr__(self, "beta", betaq)
        object.__setattr__(self, "target", targetq)
        object.__setattr__(self, "derivative_order", m)

        # h * ((J + 7)/10) - residualLoss >= target
        # iff J >= 10 * (target + residualLoss) / h - 7.
        lower = 10 * (targetq + self.loss) / hq - 7
        selected = max(0, _ceil_fraction(lower))
        object.__setattr__(self, "stage", selected)

        if self.selected_rate < targetq:
            raise AssertionError("selected stage does not reach the requested residual-rate target")
        if selected > 0 and self.previous_rate is not None and self.previous_rate >= targetq:
            raise AssertionError("selected stage is not the minimal target-reaching stage")
        if selected > 0 and self.selected_rate - self.previous_rate != hq / 10:
            raise AssertionError("pinned residual-rate successor increment is not exactly h/10")

    @property
    def loss(self) -> Fraction:
        return residual_loss(self.h, self.beta, self.derivative_order)

    @property
    def selected_rate(self) -> Fraction:
        return residual_rate(self.h, self.beta, self.stage, self.derivative_order)

    @property
    def previous_rate(self) -> Fraction | None:
        if self.stage == 0:
            return None
        return residual_rate(self.h, self.beta, self.stage - 1, self.derivative_order)

    @property
    def target_margin(self) -> Fraction:
        return self.selected_rate - self.target

    @property
    def exact_minimal_stage_machine_checked(self) -> bool:
        return True

    @property
    def fixed_derivative_target_only(self) -> bool:
        return True

    @property
    def actual_stage_field_consumed(self) -> bool:
        return False

    @property
    def actual_residual_evaluated(self) -> bool:
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


def select_residual_rate_stage(
    *,
    h: Fraction | int,
    beta: Fraction | int,
    derivative_order: int,
    target: Fraction | int,
) -> Section9ResidualRateTargetCertificate:
    """Return the minimal exact stage reaching ``target`` for fixed ``m``."""
    return Section9ResidualRateTargetCertificate(
        h=h,
        beta=beta,
        derivative_order=derivative_order,
        target=target,
    )
