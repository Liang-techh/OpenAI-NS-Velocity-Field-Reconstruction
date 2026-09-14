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
constants.  It does not materialize ``AxisCoefficientSpace`` data, infer a
fixed point from tests, sample/fill coefficients, or claim the entrance
inequalities themselves before a genuine backend supplies the coefficient
state and its norm-error witness.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction

from .axis_fixed_point_picard import NaturalPicardContractionCertificate
from .natural_scale_selection_wide import (
    WideNaturalScaleSelection,
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


@dataclass(frozen=True)
class NaturalEntranceScaleCertificate:
    """Exact scalar gate used before any profile-level entrance claim.

    ``fixed_point_error_radius_upper`` is the same upward-rounded quantity
    ``remainderBound/(2*Lambda)`` exposed by the landed Picard certificate.
    The certificate intentionally stores no coefficient values or profile
    samples.
    """

    remainder_bound_upper: Decimal
    stability_threshold: Decimal
    Lambda: Decimal
    fixed_point_error_radius_upper: Decimal

    def __post_init__(self) -> None:
        bound = _finite_nonnegative_decimal(
            self.remainder_bound_upper, "remainder_bound_upper"
        )
        stability = _finite_positive_decimal(
            self.stability_threshold, "stability_threshold"
        )
        Lambda = _finite_positive_decimal(self.Lambda, "Lambda")
        radius = _finite_nonnegative_decimal(
            self.fixed_point_error_radius_upper, "fixed_point_error_radius_upper"
        )

        pinned = stability_scale_wide(bound)
        if stability != pinned:
            raise ValueError("stability threshold is not the pinned AxisReference value")
        if Lambda < stability:
            raise ValueError("Lambda does not dominate the pinned stability threshold")

        # AxisContraction's fixed-point theorem gives exactly K/(2*Lambda).
        # Reuse the landed certificate rather than independently reimplementing
        # its upward-rounded Decimal division.
        scale = WideNaturalScaleSelection(
            remainder_bound=bound,
            remainder_lipschitz=Decimal(0),
            phase_real_part_sup=Decimal(0),
            contraction=Decimal(1) + bound,
            stability=stability,
            Lambda=Lambda,
            C=type(_zero_scale().C)(Decimal(0)),
        )
        expected = NaturalPicardContractionCertificate.from_scale(scale).one_step_radius_upper
        if radius != expected:
            raise ValueError("fixed-point error radius is not the pinned K/(2*Lambda) bound")

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

        return cls(
            remainder_bound_upper=scale.remainder_bound,
            stability_threshold=scale.stability,
            Lambda=scale.Lambda,
            fixed_point_error_radius_upper=picard.one_step_radius_upper,
        )

    def require_backend_norm_error(self, error_upper: Decimal) -> None:
        """Require a genuine backend error witness to fit the theorem radius.

        Passing this scalar gate is necessary for the formal entrance
        implications, but is deliberately not sufficient to assert that the
        caller actually owns the theorem's coefficient-space state.
        """

        error = _finite_nonnegative_decimal(error_upper, "error_upper")
        if error > self.fixed_point_error_radius_upper:
            raise ValueError("backend norm error exceeds K/(2*Lambda)")

    def require_log_slope_chi_bounds(
        self,
        *,
        chi_lower: Fraction,
        chi_upper: Fraction,
    ) -> None:
        """Check only the exact chi interval needed by the log-slope theorem.

        The caller must provide analytically/formally certified exact bounds at
        the parameter point.  Float/Decimal approximations are rejected so a
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


def _zero_scale() -> WideNaturalScaleSelection:
    """Internal constructor used only to obtain the symbolic-C dataclass type."""

    from .natural_scale_selection_wide import select_natural_scale_wide

    return select_natural_scale_wide(
        remainder_bound=Decimal(0),
        remainder_lipschitz=Decimal(0),
        phase_real_part_sup=Decimal(0),
    )


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
