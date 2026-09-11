"""Constructive norm inputs for the natural-axis fixed-point problem.

This module materializes the norm bookkeeping that is implicit in the pinned
OpenAI Lean construction ``NaturalAxisCoefficients`` + ``AxisResolvent``.
It deliberately stops *before* claiming a paper-exact fixed point.

For the canonical existence proof in ``exists_coefficientFamily_on_neighborhood``
we have

* coefficient radius ``epsilon = rho / 2``;
* ``radiusLoss(epsilon / rho) = radiusLoss(1/2) = 12``;
* every fixed coefficient-space field has norm at most ``12 * B`` when ``B``
  is a certified common complex-value bound on the compact neighborhood;
* the normalized amplitude has norm at most ``12``;
* ``||chi|| <= 12 * B`` therefore gives the resolvent factorial majorant
  ``K = 2560 * ||chi|| <= 30720 * B``.

The factorial-series upper bound is enclosed by a positive-term partial sum and
an explicit geometric tail after the successive-term ratio drops below 1/2.
Decimal arithmetic is rounded upward, then the result is converted to a binary64
upper value with one final upward correction.  This is executable conservative
bookkeeping, not an interval/Lean proof object.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, ROUND_CEILING, localcontext
import math

from .axis_remainder_bounds import AxisDataNormBounds, NaturalOperatorNormBounds


_RADIUS_LOSS_HALF = 12.0
_RESOLVENT_POWER_CONSTANT = 2560.0
_DECIMAL_PRECISION = 80
_MAX_SERIES_STEPS = 100_000


def _finite_positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _finite_nonnegative(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


def _mul_upper(left: float, right: float, name: str) -> float:
    """Return a binary64 upper rounding of a positive product."""

    left = _finite_nonnegative(left, f"{name} left factor")
    right = _finite_nonnegative(right, f"{name} right factor")
    product = left * right
    if not math.isfinite(product):
        raise ArithmeticError(f"{name} overflowed binary64")
    if product == 0.0:
        return 0.0
    return math.nextafter(product, math.inf)


def _decimal_to_float_upper(value: Decimal, name: str) -> float:
    candidate = float(value)
    if not math.isfinite(candidate):
        raise ArithmeticError(f"{name} exceeds binary64 range")
    if Decimal.from_float(candidate) < value:
        candidate = math.nextafter(candidate, math.inf)
    if not math.isfinite(candidate):
        raise ArithmeticError(f"{name} upward rounding overflowed binary64")
    return candidate


@dataclass(frozen=True)
class FactorialResolventMajorant:
    """Executable upper certificate for ``sum K^k/(k!(k+1)!)``."""

    majorant_parameter_upper: float
    partial_term_count: int
    tail_ratio_upper: float
    series_upper_decimal: Decimal

    def __post_init__(self) -> None:
        _finite_nonnegative(self.majorant_parameter_upper, "majorant_parameter_upper")
        if self.partial_term_count < 1:
            raise ValueError("partial_term_count must be positive")
        ratio = _finite_nonnegative(self.tail_ratio_upper, "tail_ratio_upper")
        if ratio > 0.5:
            raise ValueError("tail_ratio_upper must be at most 1/2")
        if not self.series_upper_decimal.is_finite() or self.series_upper_decimal < 0:
            raise ValueError("series_upper_decimal must be finite and nonnegative")

    @property
    def series_upper(self) -> float:
        """Binary64 upper rounding, failing closed if the real bound is too large."""

        return _decimal_to_float_upper(
            self.series_upper_decimal,
            "natural-resolvent factorial-series bound",
        )


def factorial_resolvent_majorant_upper(
    majorant_parameter_upper: float,
) -> FactorialResolventMajorant:
    """Bound the pinned ``AxisResolvent.factorialMajorant`` infinite sum.

    For ``a_k = K^k/(k!(k+1)!)`` the exact ratio is
    ``a_{k+1}/a_k = K / ((k+1)(k+2))`` and decreases with ``k``.  Once that
    ratio is at most 1/2, the untouched tail is bounded by the corresponding
    geometric series.  All Decimal operations in the positive recurrence are
    rounded toward +infinity.
    """

    K = _finite_nonnegative(majorant_parameter_upper, "majorant_parameter_upper")
    with localcontext() as ctx:
        ctx.prec = _DECIMAL_PRECISION
        ctx.rounding = ROUND_CEILING
        kd = Decimal.from_float(K)
        term = Decimal(1)  # a_0
        partial = Decimal(1)
        k = 0
        half = Decimal(1) / Decimal(2)
        one = Decimal(1)

        while k < _MAX_SERIES_STEPS:
            ratio = kd / Decimal((k + 1) * (k + 2))
            if ratio <= half:
                # partial currently contains a_0 + ... + a_k.  Since all later
                # ratios are <= ratio, tail <= a_k * ratio / (1-ratio).
                tail_upper = term * ratio / (one - ratio)
                series_upper = partial + tail_upper
                return FactorialResolventMajorant(
                    majorant_parameter_upper=K,
                    partial_term_count=k + 1,
                    tail_ratio_upper=_decimal_to_float_upper(ratio, "tail ratio"),
                    series_upper_decimal=+series_upper,
                )

            term = term * ratio
            partial = partial + term
            k += 1

    raise ArithmeticError(
        "factorial-series enclosure exceeded the deterministic step budget"
    )


@dataclass(frozen=True)
class ConstructiveAnalyticInputNormCertificate:
    """Norm ledger from the canonical ``rho/2`` coefficient-family choice.

    The two upstream inputs intentionally remain explicit:

    ``neighborhood_radius``
        A proved positive common analytic radius ``rho``.
    ``complex_field_sup_upper``
        A proved common value bound ``B`` for every fixed complex field on the
        compact neighborhood used by ``finite_family_bound``.

    No sampling routine in this module is allowed to manufacture either input.
    """

    neighborhood_radius: float
    complex_field_sup_upper: float
    epsilon: float
    radius_loss_half: float
    coefficient_family_norm_upper: float
    amplitude_norm_upper: float
    chi_norm_upper: float
    resolvent_majorant: FactorialResolventMajorant

    def __post_init__(self) -> None:
        _finite_positive(self.neighborhood_radius, "neighborhood_radius")
        _finite_positive(self.complex_field_sup_upper, "complex_field_sup_upper")
        _finite_positive(self.epsilon, "epsilon")
        if self.epsilon != self.neighborhood_radius / 2.0:
            raise ValueError("epsilon must be the canonical neighborhood_radius/2 choice")
        if self.radius_loss_half != _RADIUS_LOSS_HALF:
            raise ValueError("radius_loss_half must equal the exact value 12")
        _finite_positive(self.coefficient_family_norm_upper, "coefficient_family_norm_upper")
        _finite_positive(self.amplitude_norm_upper, "amplitude_norm_upper")
        _finite_positive(self.chi_norm_upper, "chi_norm_upper")
        if self.amplitude_norm_upper != _RADIUS_LOSS_HALF:
            raise ValueError("canonical amplitude norm upper must equal 12")
        if self.chi_norm_upper != self.coefficient_family_norm_upper:
            raise ValueError("chi norm upper must reuse the common coefficient-family bound")

    def operator_norm_bounds(self) -> NaturalOperatorNormBounds:
        """Instantiate the landed AxisOperators ledger with the resolvent bound."""

        return NaturalOperatorNormBounds.pinned_axis_operators(
            epsilon=self.epsilon,
            resolvent_norm_upper=self.resolvent_majorant.series_upper,
        )

    def axis_data_norm_bounds(self, *, h: float) -> AxisDataNormBounds:
        """Collapse ``CoefficientFamily.norm_le`` to the landed AxisData ledger.

        Every coefficient-space field in ``CoefficientFamily.axisData`` is one
        of ``elements k`` and is therefore bounded by the same family norm
        upper.  ``A``, ``D`` and ``h`` are scalar parameters rather than
        coefficient-space norms.
        """

        h = float(h)
        if not math.isfinite(h):
            raise ValueError("h must be finite")
        bound = self.coefficient_family_norm_upper
        return AxisDataNormBounds(
            A=0.5 + h,
            D=0.5 - h,
            h=h,
            one=bound,
            eta=bound,
            d=bound,
            inverse_l=bound,
            u_star=bound,
            u_star_eta=bound,
            w_star=bound,
            h_star=bound,
            normalized_gradient=bound,
            z_star=bound,
        )


def constructive_analytic_input_norm_certificate(
    *,
    neighborhood_radius: float,
    complex_field_sup_upper: float,
) -> ConstructiveAnalyticInputNormCertificate:
    """Materialize the canonical Lean norm choices from ``rho`` and ``B``.

    The identity ``radiusLoss(1/2)=12`` follows from

    ``sum_(m>=0) (m+1)^2 q^m = (1+q)/(1-q)^3``.

    Thus the canonical ``epsilon=rho/2`` construction gives both the common
    coefficient-family loss ``12*B`` and the normalized amplitude bound ``12``.
    The natural-resolvent majorant then uses
    ``K = 2560*||chi|| <= 2560*(12*B)``.
    """

    rho = _finite_positive(neighborhood_radius, "neighborhood_radius")
    B = _finite_positive(complex_field_sup_upper, "complex_field_sup_upper")
    epsilon = rho / 2.0
    if not math.isfinite(epsilon) or epsilon <= 0.0:
        raise ArithmeticError("canonical epsilon=rho/2 underflowed or overflowed")

    family_bound = _mul_upper(B, _RADIUS_LOSS_HALF, "coefficient-family norm bound")
    chi_bound = family_bound
    K = _mul_upper(
        _RESOLVENT_POWER_CONSTANT,
        chi_bound,
        "natural-resolvent majorant parameter",
    )
    resolvent = factorial_resolvent_majorant_upper(K)
    return ConstructiveAnalyticInputNormCertificate(
        neighborhood_radius=rho,
        complex_field_sup_upper=B,
        epsilon=epsilon,
        radius_loss_half=_RADIUS_LOSS_HALF,
        coefficient_family_norm_upper=family_bound,
        amplitude_norm_upper=_RADIUS_LOSS_HALF,
        chi_norm_upper=chi_bound,
        resolvent_majorant=resolvent,
    )
