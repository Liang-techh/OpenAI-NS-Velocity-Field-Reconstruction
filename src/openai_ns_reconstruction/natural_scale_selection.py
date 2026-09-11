"""Theorem-faithful scale/normalization selection for natural profiles.

This module materializes the selection algebra used by the pinned Lean chain
``AxisContraction.uniform_natural_fixedPoint`` ->
``AxisReference.exists_positive_scaled_profiles`` ->
``NaturalProfile.exists_natural_profiles``.

It intentionally does not manufacture the coefficient-space quantities that
enter those theorems. ``remainder_bound``, ``remainder_lipschitz``, and the
compact-set ``realPartSup`` must still be certified upstream from the actual
analytic coefficient family. Supplying numbers here therefore yields a
formal-structure witness, not a paper-exact profile by itself.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


# AxisEvaluation.jetBound(epsilon, 5, 0, 0) is independent of epsilon:
# sum_n (5/20)^n = 1/(1-1/4) = 4/3.
JET_VALUE_BOUND_AT_FIVE = 4.0 / 3.0

# AxisEvaluation.jetBound(epsilon, 5, 1, 0) is likewise independent of epsilon:
# sum_n (n+1)(5/20)^n = 1/(1-1/4)^2 = 16/9.
JET_RADIAL_DERIVATIVE_BOUND_AT_FIVE = 16.0 / 9.0


def _finite_nonnegative(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


def contraction_threshold(remainder_bound: float, remainder_lipschitz: float) -> float:
    """Pinned ``AxisContraction.contractionThreshold`` after B/L evaluation.

    Lean defines ``1 + remainderBound + remainderLip`` at the reference-pair
    ball. The caller is responsible for supplying those two actual evaluated
    quantities, not sampled or fitted substitutes.
    """

    bound = _finite_nonnegative(remainder_bound, "remainder_bound")
    lip = _finite_nonnegative(remainder_lipschitz, "remainder_lipschitz")
    value = 1.0 + bound + lip
    if not math.isfinite(value):
        raise ArithmeticError("contraction threshold overflowed")
    return value


def stability_scale(remainder_bound: float) -> float:
    """Pinned ``AxisReference.stabilityScale`` with its two jet sums closed.

    ``errorConstant`` in ``NaturalProfile`` is the same evaluated
    ``remainderBound`` used in the fixed-point error estimate. At ``R=5`` and
    parameter-derivative order ``m=0``, the pinned ``AxisEvaluation.majorant``
    sums exactly to ``4/3`` for ``k=0`` and ``16/9`` for ``k=1``. Therefore
    ``stabilityScale = 1 + (14000/9) * remainderBound``.
    """

    bound = _finite_nonnegative(remainder_bound, "remainder_bound")
    value = 1.0 + 500.0 * (
        JET_VALUE_BOUND_AT_FIVE + JET_RADIAL_DERIVATIVE_BOUND_AT_FIVE
    ) * bound
    if not math.isfinite(value):
        raise ArithmeticError("stability scale overflowed")
    return value


def normalization_threshold(Lambda: float, phase_real_part_sup: float) -> float:
    """Pinned ``AnalyticInputs.normalizationThreshold``.

    This is literally ``exp(Lambda * realPartSup(axisPhase, compactSet))``.
    ``phase_real_part_sup`` must therefore be the certified compact-set
    supremum from the actual ``AnalyticInputs``, not a sampled grid maximum.
    """

    Lambda = float(Lambda)
    if not math.isfinite(Lambda) or Lambda <= 0.0:
        raise ValueError("Lambda must be finite and positive")
    phase_sup = _finite_nonnegative(phase_real_part_sup, "phase_real_part_sup")
    exponent = Lambda * phase_sup
    if not math.isfinite(exponent):
        raise ArithmeticError("normalization exponent overflowed")
    try:
        value = math.exp(exponent)
    except OverflowError as exc:
        raise ArithmeticError(
            "normalization threshold is finite in the theorem but exceeds binary64 range"
        ) from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ArithmeticError("normalization threshold must remain finite and positive")
    return value


@dataclass(frozen=True)
class NaturalScaleSelection:
    """Executable witness for the theorem's deterministic admissible choices.

    ``Lambda`` is the smallest threshold produced by the two theorem layers
    once their actual coefficient-space remainder constants are supplied.
    ``C`` is then chosen at equality with
    ``AnalyticInputs.normalizationThreshold``, exactly as in
    ``NaturalProfile.exists_natural_profiles``.
    """

    remainder_bound: float
    remainder_lipschitz: float
    phase_real_part_sup: float
    contraction: float
    stability: float
    Lambda: float
    C: float

    def admits(self, Lambda: float, C: float) -> bool:
        """Check the two theorem inequalities for another candidate pair."""

        Lambda = float(Lambda)
        C = float(C)
        if not math.isfinite(Lambda) or not math.isfinite(C) or Lambda <= 0.0 or C <= 0.0:
            return False
        if Lambda < self.contraction or Lambda < self.stability:
            return False
        try:
            threshold = normalization_threshold(Lambda, self.phase_real_part_sup)
        except (ValueError, ArithmeticError):
            return False
        return C >= threshold


def select_natural_scale(
    *,
    remainder_bound: float,
    remainder_lipschitz: float,
    phase_real_part_sup: float,
) -> NaturalScaleSelection:
    """Select the exact theorem-side ``Lambda`` and ``C`` from certified inputs.

    The pinned proof first uses ``contractionThreshold = 1+B+L`` and then
    replaces it by ``max(contractionThreshold, stabilityScale)``. The final
    natural-profile theorem takes that threshold itself as ``Lambda`` and sets
    ``C = normalizationThreshold Lambda``.
    """

    bound = _finite_nonnegative(remainder_bound, "remainder_bound")
    lip = _finite_nonnegative(remainder_lipschitz, "remainder_lipschitz")
    phase_sup = _finite_nonnegative(phase_real_part_sup, "phase_real_part_sup")

    contraction = contraction_threshold(bound, lip)
    stability = stability_scale(bound)
    Lambda = max(contraction, stability)
    C = normalization_threshold(Lambda, phase_sup)
    return NaturalScaleSelection(
        remainder_bound=bound,
        remainder_lipschitz=lip,
        phase_real_part_sup=phase_sup,
        contraction=contraction,
        stability=stability,
        Lambda=Lambda,
        C=C,
    )
