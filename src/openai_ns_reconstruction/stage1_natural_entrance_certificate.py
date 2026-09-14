"""Fail-closed scalar premises for the Stage-1 natural entrance estimates.

The pinned formal chain proves two quantitative entrance facts for a genuine
coefficient-space solution near ``referencePair``:

* ``NaturalEntrance.coefficient_phi_lower`` gives ``phi > 1/8`` on the scaled
  entrance strip ``0 <= Y <= 41/10``; and
* ``AxisReference.log_slope_of_uniformMixedError`` gives
  ``23/10 < -8 * partialY(phi)(4, eta) / phi(4, eta)`` wherever the input
  ``chi`` coefficient is at least ``99/100``.

Both results consume the same scale hypothesis

``AxisReference.stabilityScale epsilon K <= Lambda``

and the fixed-point error bound ``||x-referencePair|| <= K/(2*Lambda)``.
For the current natural-profile chain, ``K`` is the same certified
``remainderBound`` used by ``AxisContraction.exists_unique_natural_fixedPoint``.
The pinned jet sums at radial radius five reduce the stability threshold to
``1 + (14000/9) K``; :mod:`natural_scale_selection_wide` already carries this
quantity in 96-digit upward-rounded ``Decimal`` arithmetic.

This module therefore certifies only those *scalar premises* and exact theorem
constants. It does not materialize ``AxisCoefficientSpace`` data, infer a
fixed point from tests, sample/fill coefficients, or claim the entrance
inequalities themselves before a genuine backend supplies the coefficient
state and its norm-error witness.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction

from .axis_fixed_point_picard import NaturalPicardContractionCertificate
from .natural_scale_selection_wide import (
    SymbolicExponentialThreshold,
    WideNaturalScaleSelection,
    contraction_threshold_wide,
    stability_scale_wide,
)
from .outgoing_tail import TailData
from .stage1_scale_chain_wide import diagnose_actual_schedule_scale_chain_wide


ENTRANCE_SCALED_RADIUS_MAX = Fraction(41, 10)
PHI_STRICT_LOWER = Fraction(1, 8)
LOG_SLOPE_SCALED_RADIUS = Fraction(4, 1)
CHI_LOG_SLOPE_MIN = Fraction(99, 100)
CHI_UPPER = Fraction(1, 1)
LOG_SLOPE_STRICT_LOWER = Fraction(23, 10)
_DECIMAL_PRECISION = 96


def _finite_nonnegative_decimal(value: Decimal, name: str) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite() or value < 0:
        raise ValueError(f"{name} must be a finite nonnegative Decimal")
    return value


def _finite_positive_decimal(value: Decimal, name: str) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite() or value <= 0:
        raise ValueError(f"{name} must be a finite positive Decimal")
    return value


def _exact_fraction(value: Fraction, name: str) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be an exact Fraction")
    return value


def _fixed_point_radius_down(bound: Decimal, Lambda: Decimal) -> Decimal:
    """Safe lower Decimal for the exact theorem radius ``bound/(2*Lambda)``.

    A backend supplies an *upper* bound for ``||x-referencePair||``. Admission
    therefore needs a lower approximation to the theorem's allowed radius,
    never the Picard certificate's upward-rounded product. The denominator is
    formed at double precision before the final division is rounded toward
    ``-infinity`` (the ratio is nonnegative).
    """

    bound = _finite_nonnegative_decimal(bound, "bound")
    Lambda = _finite_positive_decimal(Lambda, "Lambda")
    with localcontext() as ctx:
        ctx.prec = 2 * _DECIMAL_PRECISION
        denominator = Decimal(2) * Lambda
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.rounding = ROUND_FLOOR
        return +(bound / denominator)


@dataclass(frozen=True)
class NaturalEntranceScaleCertificate:
    """Exact scalar gate used before any profile-level entrance claim.

    ``fixed_point_error_radius_safe`` is a downward-rounded lower bound for the
    exact theorem radius ``K/(2*Lambda)``. Requiring a backend error upper bound
    to fit this value is conservative and cannot admit a state merely because
    an upward-rounded Decimal crossed the theorem threshold.
    """

    remainder_bound_upper: Decimal
    stability_threshold: Decimal
    Lambda: Decimal
    fixed_point_error_radius_safe: Decimal

    def __post_init__(self) -> None:
        bound = _finite_nonnegative_decimal(
            self.remainder_bound_upper, "remainder_bound_upper"
        )
        stability = _finite_positive_decimal(
            self.stability_threshold, "stability_threshold"
        )
        Lambda = _finite_positive_decimal(self.Lambda, "Lambda")
        radius = _finite_nonnegative_decimal(
            self.fixed_point_error_radius_safe, "fixed_point_error_radius_safe"
        )

        pinned = stability_scale_wide(bound)
        if stability != pinned:
            raise ValueError("stability threshold is not the pinned AxisReference value")
        if Lambda < stability:
            raise ValueError("Lambda does not dominate the pinned stability threshold")

        expected_safe = _fixed_point_radius_down(bound, Lambda)
        if radius != expected_safe:
            raise ValueError("fixed-point safe radius is not the pinned K/(2*Lambda) gate")

        # Cross-check directionality against the landed upward-rounded Picard
        # arithmetic. The safe admission radius must never exceed that upper
        # product; equality is allowed when the quotient is exactly representable.
        scalar_scale = WideNaturalScaleSelection(
            remainder_bound=bound,
            remainder_lipschitz=Decimal(0),
            phase_real_part_sup=Decimal(0),
            contraction=contraction_threshold_wide(bound, Decimal(0)),
            stability=stability,
            Lambda=Lambda,
            C=SymbolicExponentialThreshold(Decimal(0)),
        )
        picard_upper = NaturalPicardContractionCertificate.from_scale(
            scalar_scale
        ).one_step_radius_upper
        if radius > picard_upper:
            raise RuntimeError("safe entrance radius exceeds the Picard upper product")

    @classmethod
    def from_scale(
        cls,
        scale: WideNaturalScaleSelection,
    ) -> "NaturalEntranceScaleCertificate":
        """Bind the entrance gate to one already-certified wide scale object."""

        if not isinstance(scale, WideNaturalScaleSelection):
            raise TypeError("scale must be WideNaturalScaleSelection")

        pinned_stability = stability_scale_wide(scale.remainder_bound)
        if scale.stability != pinned_stability:
            raise ValueError("scale.stability does not equal the pinned stability threshold")
        if scale.Lambda < scale.stability:
            raise ValueError("scale.Lambda is below its pinned stability threshold")

        picard = NaturalPicardContractionCertificate.from_scale(scale)
        if picard.remainder_bound_upper != scale.remainder_bound:
            raise RuntimeError("Picard/entrance remainder-bound identity mismatch")
        if picard.Lambda != scale.Lambda:
            raise RuntimeError("Picard/entrance Lambda identity mismatch")

        safe_radius = _fixed_point_radius_down(scale.remainder_bound, scale.Lambda)
        if safe_radius > picard.one_step_radius_upper:
            raise RuntimeError("safe entrance radius exceeds actual Picard upper product")

        return cls(
            remainder_bound_upper=scale.remainder_bound,
            stability_threshold=scale.stability,
            Lambda=scale.Lambda,
            fixed_point_error_radius_safe=safe_radius,
        )

    def require_backend_norm_error(self, error_upper: Decimal) -> None:
        """Require a genuine backend error witness to fit the theorem radius.

        Passing this scalar gate is necessary for the formal entrance
        implications, but is deliberately not sufficient to assert that the
        caller actually owns the theorem's coefficient-space state.
        """

        error = _finite_nonnegative_decimal(error_upper, "error_upper")
        if error > self.fixed_point_error_radius_safe:
            raise ValueError("backend norm error exceeds safe K/(2*Lambda) radius")

    def require_log_slope_chi_bounds(
        self,
        *,
        chi_lower: Fraction,
        chi_upper: Fraction,
    ) -> None:
        """Check only the exact chi interval needed by the log-slope theorem.

        The caller must provide analytically/formally certified exact bounds at
        the parameter point. Float/Decimal approximations are rejected so a
        sampled ``0.99`` cannot masquerade as the theorem hypothesis.
        """

        lower = _exact_fraction(chi_lower, "chi_lower")
        upper = _exact_fraction(chi_upper, "chi_upper")
        if lower > upper:
            raise ValueError("chi lower bound exceeds chi upper bound")
        if lower < CHI_LOG_SLOPE_MIN:
            raise ValueError("chi lower bound is below 99/100")
        if upper > CHI_UPPER:
            raise ValueError("chi upper bound exceeds 1")

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def full_reconstruction(self) -> bool:
        return False


def actual_schedule_natural_entrance_scale_certificate(
    data: TailData,
    j: float,
) -> NaturalEntranceScaleCertificate:
    """Construct the entrance scalar gate from the actual SchedulePressure chain."""

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    diagnostic = diagnose_actual_schedule_scale_chain_wide(data, j)
    if diagnostic.upstream_obstruction is not None:
        raise RuntimeError(diagnostic.upstream_obstruction)
    if diagnostic.scale is None:
        raise RuntimeError("actual-schedule chain did not produce a wide natural scale")
    return NaturalEntranceScaleCertificate.from_scale(diagnostic.scale)
