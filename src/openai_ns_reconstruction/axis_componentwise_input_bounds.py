"""Componentwise analytic norm ledger for the natural-axis coefficient family.

The pinned Lean proof of ``exists_coefficientFamily_on_neighborhood`` first
chooses one common complex-value bound ``B`` for all eleven fixed fields.  That
is convenient for existence, but it is much looser than required by the
underlying ``boundedAxisElement_norm`` theorem, which is stated for one field
and one bound at a time.

For the already-certified actual SchedulePressure neighborhood we have an
individual complex-value upper bound for every fixed field.  This module keeps
those bounds separate while using the same common radius ``epsilon=rho/2`` and
the same exact Cauchy loss ``radiusLoss(1/2)=12``.  The resulting elements fit
the same ``CoefficientFamily`` structure by taking the maximum of the eleven
individual norm bounds as its common ``bound`` field.

The main practical consequence is theorem-faithful rather than heuristic:
``AxisResolvent`` depends on ``||chi||`` specifically.  A very large gradient
bound therefore no longer contaminates the resolvent parameter through the
unrelated common family bound.

This is still formal-structure bookkeeping.  It does not materialize the Lean
Banach-space elements or the fixed point and does not change any paper-exact
gate.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

from .axis_analytic_input_bounds import (
    FactorialResolventMajorant,
    factorial_resolvent_majorant_upper,
)
from .axis_remainder_bounds import AxisDataNormBounds, NaturalOperatorNormBounds

_RADIUS_LOSS_HALF = 12.0
_RESOLVENT_POWER_CONSTANT = 2560.0

_FIELD_ORDER = (
    "one",
    "eta",
    "d",
    "inverseL",
    "uStar",
    "uStarEta",
    "wStar",
    "hStar",
    "zStar",
    "chi",
    "gradient",
)


def _finite_positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _mul_upper(left: float, right: float, name: str) -> float:
    left = _finite_positive(left, f"{name} left factor")
    right = _finite_positive(right, f"{name} right factor")
    product = left * right
    if not math.isfinite(product):
        raise ArithmeticError(f"{name} overflowed binary64")
    return math.nextafter(product, math.inf)


@dataclass(frozen=True)
class ComponentwiseAnalyticInputNormCertificate:
    """Per-field Cauchy norm bounds at the canonical half-radius choice."""

    neighborhood_radius: float
    epsilon: float
    radius_loss_half: float
    field_value_sup_upper: Mapping[str, float]
    field_norm_upper: Mapping[str, float]
    coefficient_family_norm_upper: float
    amplitude_norm_upper: float
    chi_norm_upper: float
    resolvent_majorant: FactorialResolventMajorant

    def __post_init__(self) -> None:
        rho = _finite_positive(self.neighborhood_radius, "neighborhood_radius")
        eps = _finite_positive(self.epsilon, "epsilon")
        if eps != rho / 2.0:
            raise ValueError("epsilon must be the canonical neighborhood_radius/2 choice")
        if self.radius_loss_half != _RADIUS_LOSS_HALF:
            raise ValueError("radius_loss_half must equal the exact value 12")
        if tuple(self.field_value_sup_upper) != _FIELD_ORDER:
            raise ValueError("field_value_sup_upper must list the pinned fixed fields")
        if tuple(self.field_norm_upper) != _FIELD_ORDER:
            raise ValueError("field_norm_upper must list the pinned fixed fields")
        for name in _FIELD_ORDER:
            _finite_positive(self.field_value_sup_upper[name], f"field value bound {name}")
            _finite_positive(self.field_norm_upper[name], f"field norm bound {name}")
        common = _finite_positive(
            self.coefficient_family_norm_upper,
            "coefficient_family_norm_upper",
        )
        if common < max(self.field_norm_upper.values()):
            raise ValueError("coefficient family bound must dominate every field norm")
        if self.amplitude_norm_upper != _RADIUS_LOSS_HALF:
            raise ValueError("canonical amplitude norm upper must equal 12")
        chi = _finite_positive(self.chi_norm_upper, "chi_norm_upper")
        if chi != self.field_norm_upper["chi"]:
            raise ValueError("chi norm upper must use the componentwise chi bound")

    def operator_norm_bounds(self) -> NaturalOperatorNormBounds:
        """Instantiate the pinned AxisOperators ledger with this resolvent bound."""

        return NaturalOperatorNormBounds.pinned_axis_operators(
            epsilon=self.epsilon,
            resolvent_norm_upper=self.resolvent_majorant.series_upper,
        )

    def axis_data_norm_bounds(self, *, h: float) -> AxisDataNormBounds:
        """Map the individual fixed-field norms into ``CoefficientFamily.axisData``."""

        h = float(h)
        if not math.isfinite(h):
            raise ValueError("h must be finite")
        n = self.field_norm_upper
        return AxisDataNormBounds(
            A=0.5 + h,
            D=0.5 - h,
            h=h,
            one=n["one"],
            eta=n["eta"],
            d=n["d"],
            inverse_l=n["inverseL"],
            u_star=n["uStar"],
            u_star_eta=n["uStarEta"],
            w_star=n["wStar"],
            h_star=n["hStar"],
            normalized_gradient=n["gradient"],
            z_star=n["zStar"],
        )


def componentwise_analytic_input_norm_certificate(
    *,
    neighborhood_radius: float,
    field_value_sup_upper: Mapping[str, float],
) -> ComponentwiseAnalyticInputNormCertificate:
    """Apply ``boundedAxisElement_norm`` separately to the eleven fixed fields.

    Every field uses the same analytic radius ``rho`` and coefficient radius
    ``epsilon=rho/2``.  If ``B_k`` bounds one complex field on the common tube,
    the pinned Cauchy constructor gives ``||element_k|| <= 12 B_k``.  The
    ``CoefficientFamily.bound`` slot can then be the maximum of these eleven
    upper bounds, while ``AxisResolvent`` uses the much sharper chi-specific
    norm instead of that maximum.
    """

    rho = _finite_positive(neighborhood_radius, "neighborhood_radius")
    if tuple(field_value_sup_upper) != _FIELD_ORDER:
        raise ValueError("field_value_sup_upper must list the pinned fixed fields")

    values = MappingProxyType(
        {
            name: _finite_positive(field_value_sup_upper[name], f"field value bound {name}")
            for name in _FIELD_ORDER
        }
    )
    epsilon = rho / 2.0
    if not math.isfinite(epsilon) or epsilon <= 0.0:
        raise ArithmeticError("canonical epsilon=rho/2 underflowed or overflowed")

    norms = MappingProxyType(
        {
            name: _mul_upper(values[name], _RADIUS_LOSS_HALF, f"{name} coefficient norm")
            for name in _FIELD_ORDER
        }
    )
    family_upper = max(norms.values())
    chi_upper = norms["chi"]
    K = _mul_upper(
        _RESOLVENT_POWER_CONSTANT,
        chi_upper,
        "natural-resolvent majorant parameter",
    )
    resolvent = factorial_resolvent_majorant_upper(K)

    return ComponentwiseAnalyticInputNormCertificate(
        neighborhood_radius=rho,
        epsilon=epsilon,
        radius_loss_half=_RADIUS_LOSS_HALF,
        field_value_sup_upper=values,
        field_norm_upper=norms,
        coefficient_family_norm_upper=family_upper,
        amplitude_norm_upper=_RADIUS_LOSS_HALF,
        chi_norm_upper=chi_upper,
        resolvent_majorant=resolvent,
    )
