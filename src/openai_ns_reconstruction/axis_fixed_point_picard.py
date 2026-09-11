"""Constructive Picard bridge for the pinned natural-axis fixed point.

The pinned OpenAI theorem ``AxisContraction.exists_unique_natural_fixedPoint``
uses the map

``F(x) = x0 + (1 / (2 * Lambda)) * naturalRemainder(..., x)``

on the closed unit ball around ``x0``.  Its proof requires only

``(1 / (2*Lambda)) * remainderBound <= 1`` and
``(1 / (2*Lambda)) * remainderLip <= 1/2``.

Current main already carries the actual SchedulePressure analytic-input chain
through conservative wide ``remainderBound/remainderLip`` values and the pinned
wide ``Lambda`` choice.  This module closes the next *scalar theorem gate*: it
turns that actual-schedule chain into a fail-closed contraction certificate and
provides the exact Picard-map shape without forcing ``Lambda`` back through
binary64.

This is deliberately not the missing coefficient-space implementation.  A
caller must still provide the genuine ``AxisCoefficientSpace`` state algebra
and the actual ``naturalRemainder`` evaluator before any returned iterates can
be interpreted as the paper's ``phi/u`` fields.  No finite-dimensional or toy
state is promoted to paper-exact status here.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, ROUND_CEILING, localcontext
from typing import Callable, Generic, TypeVar, Union
import math

from .natural_scale_selection_wide import WideNaturalScaleSelection
from .outgoing_tail import TailData
from .stage1_scale_chain_wide import diagnose_actual_schedule_scale_chain_wide

DecimalLike = Union[Decimal, float, int]
StateT = TypeVar("StateT")
_DECIMAL_PRECISION = 96
_HALF = Decimal(1) / Decimal(2)
_ONE = Decimal(1)


def _decimal(value: DecimalLike, name: str) -> Decimal:
    if isinstance(value, Decimal):
        result = value
    elif isinstance(value, int):
        result = Decimal(value)
    else:
        numeric = float(value)
        if not math.isfinite(numeric):
            raise ValueError(f"{name} must be finite")
        result = Decimal.from_float(numeric)
    if not result.is_finite():
        raise ValueError(f"{name} must be finite")
    return result


def _nonnegative(value: DecimalLike, name: str) -> Decimal:
    result = _decimal(value, name)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _positive(value: DecimalLike, name: str) -> Decimal:
    result = _decimal(value, name)
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _ratio_up(numerator: Decimal, denominator: Decimal) -> Decimal:
    """Positive division rounded upward, used only for certificate upper bounds."""

    if numerator < 0 or denominator <= 0:
        raise ValueError("ratio requires nonnegative numerator and positive denominator")
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.rounding = ROUND_CEILING
        return +(numerator / denominator)


def _mul_up(*values: Decimal) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.rounding = ROUND_CEILING
        out = Decimal(1)
        for value in values:
            if value < 0:
                raise ValueError("certificate multiplication requires nonnegative factors")
            out *= value
        return +out


@dataclass(frozen=True)
class NaturalPicardContractionCertificate:
    """Wide certificate for the two hypotheses used by the Lean fixed-point step.

    ``inverse_two_lambda_upper`` is an upward-rounded upper bound for
    ``1/(2*Lambda)``.  Consequently the stored ``one_step_radius_upper`` and
    ``contraction_factor_upper`` are conservative upper bounds for the two
    products appearing in ``exists_fixedPoint_of_controlled``.
    """

    remainder_bound_upper: Decimal
    remainder_lipschitz_upper: Decimal
    Lambda: Decimal
    inverse_two_lambda_upper: Decimal
    one_step_radius_upper: Decimal
    contraction_factor_upper: Decimal

    def __post_init__(self) -> None:
        _nonnegative(self.remainder_bound_upper, "remainder_bound_upper")
        _nonnegative(self.remainder_lipschitz_upper, "remainder_lipschitz_upper")
        _positive(self.Lambda, "Lambda")
        _positive(self.inverse_two_lambda_upper, "inverse_two_lambda_upper")
        _nonnegative(self.one_step_radius_upper, "one_step_radius_upper")
        _nonnegative(self.contraction_factor_upper, "contraction_factor_upper")
        if self.one_step_radius_upper > _ONE:
            raise ValueError("fixed-point map is not certified to preserve the unit ball")
        if self.contraction_factor_upper > _HALF:
            raise ValueError("fixed-point map is not certified with contraction factor <= 1/2")

    @classmethod
    def from_scale(
        cls,
        scale: WideNaturalScaleSelection,
    ) -> "NaturalPicardContractionCertificate":
        """Derive the theorem gate from the pinned wide Lambda selection.

        The method refuses a scale object that no longer dominates its own
        pinned contraction threshold.  It does not accept an independent
        caller-supplied ``B``, ``L`` or ``Lambda``.
        """

        bound = _nonnegative(scale.remainder_bound, "scale.remainder_bound")
        lip = _nonnegative(scale.remainder_lipschitz, "scale.remainder_lipschitz")
        Lambda = _positive(scale.Lambda, "scale.Lambda")
        contraction = _positive(scale.contraction, "scale.contraction")
        if Lambda < contraction:
            raise ValueError("Lambda must dominate the pinned contraction threshold")

        denominator = _mul_up(Decimal(2), Lambda)
        inverse_upper = _ratio_up(Decimal(1), denominator)
        step_upper = _mul_up(inverse_upper, bound)
        q_upper = _mul_up(inverse_upper, lip)
        return cls(
            remainder_bound_upper=bound,
            remainder_lipschitz_upper=lip,
            Lambda=Lambda,
            inverse_two_lambda_upper=inverse_upper,
            one_step_radius_upper=step_upper,
            contraction_factor_upper=q_upper,
        )

    def tail_multiplier_upper(self, completed_iterations: int) -> Decimal:
        """Geometric multiplier for ``||x* - x_n||`` from the first Picard step.

        For ``x_0 = referencePair`` and contraction factor ``q``, Banach's
        estimate is ``q^n/(1-q) * ||x_1-x_0||``.  This method returns a
        conservative Decimal upper bound for the scalar multiplier only; the
        coefficient-space norm of the first step must come from the eventual
        genuine coefficient backend.
        """

        n = int(completed_iterations)
        if n < 0 or n != completed_iterations:
            raise ValueError("completed_iterations must be a nonnegative integer")
        q = self.contraction_factor_upper
        with localcontext() as ctx:
            ctx.prec = _DECIMAL_PRECISION
            ctx.rounding = ROUND_CEILING
            numerator = +(q**n)
            denominator = _ONE - q
            if denominator <= 0:
                raise ValueError("contraction factor must be < 1")
            return +(numerator / denominator)

    def tail_bound_from_first_step(
        self,
        first_step_norm_upper: DecimalLike,
        completed_iterations: int,
    ) -> Decimal:
        first = _nonnegative(first_step_norm_upper, "first_step_norm_upper")
        return _mul_up(self.tail_multiplier_upper(completed_iterations), first)


def actual_schedule_picard_certificate(
    data: TailData,
    j: float,
) -> NaturalPicardContractionCertificate:
    """Build the Picard theorem gate from the actual SchedulePressure chain.

    This path intentionally exposes no override for the remainder bounds or
    ``Lambda``.  If the upstream schedule/analytic/resolvent/wide-remainder
    chain has not landed, it fails closed instead of manufacturing a fixed-point
    budget.
    """

    diagnostic = diagnose_actual_schedule_scale_chain_wide(data, j)
    if diagnostic.upstream_obstruction is not None:
        raise RuntimeError(diagnostic.upstream_obstruction)
    if diagnostic.remainder is None or diagnostic.scale is None:
        raise RuntimeError("actual-schedule chain did not produce a wide remainder and scale")
    if diagnostic.scale.remainder_bound != diagnostic.remainder.remainder_bound_upper:
        raise RuntimeError("wide scale/remainder bound mismatch")
    if diagnostic.scale.remainder_lipschitz != diagnostic.remainder.remainder_lipschitz_upper:
        raise RuntimeError("wide scale/remainder Lipschitz mismatch")
    return NaturalPicardContractionCertificate.from_scale(diagnostic.scale)


@dataclass(frozen=True)
class PicardAlgebra(Generic[StateT]):
    """Operations needed to execute the exact fixed-point-map shape.

    ``scale_inverse_two_lambda(value, Lambda)`` is intentionally delegated to
    the state backend.  The future genuine coefficient-space implementation can
    therefore apply the scalar ``1/(2*Lambda)`` in its own exact/wide
    representation instead of narrowing the theorem scale to ``float``.
    """

    add: Callable[[StateT, StateT], StateT]
    scale_inverse_two_lambda: Callable[[StateT, Decimal], StateT]


def natural_picard_step(
    reference: StateT,
    current: StateT,
    *,
    certificate: NaturalPicardContractionCertificate,
    remainder: Callable[[StateT], StateT],
    algebra: PicardAlgebra[StateT],
) -> StateT:
    """Execute ``x0 + (1/(2*Lambda))*R(x)`` without choosing a surrogate state."""

    return algebra.add(
        reference,
        algebra.scale_inverse_two_lambda(remainder(current), certificate.Lambda),
    )


def iterate_natural_picard(
    reference: StateT,
    *,
    certificate: NaturalPicardContractionCertificate,
    remainder: Callable[[StateT], StateT],
    algebra: PicardAlgebra[StateT],
    iterations: int,
) -> tuple[StateT, ...]:
    """Return ``x_0,...,x_N`` for the pinned Picard map.

    This function is an execution adapter, not a certificate that ``remainder``
    is the official ``naturalRemainder``.  That identity must be supplied by the
    forthcoming coefficient-space implementation and its provenance/tests.
    """

    count = int(iterations)
    if count < 0 or count != iterations:
        raise ValueError("iterations must be a nonnegative integer")
    states = [reference]
    current = reference
    for _ in range(count):
        current = natural_picard_step(
            reference,
            current,
            certificate=certificate,
            remainder=remainder,
            algebra=algebra,
        )
        states.append(current)
    return tuple(states)
