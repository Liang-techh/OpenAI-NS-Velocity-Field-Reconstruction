"""Clean-room exact replay of the compact radial moment/stress inverse seam.

This module intentionally implements only the finite algebra in the public R6--R12
formulas of the corrected Kokuno reader.  It is not a constructor for the paper's
admissible stress and it does not prove support/regularity/existence hypotheses.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

Exact = Fraction
Polynomial = tuple[Exact, ...]  # ascending powers of r


def _exact(value: object, *, name: str) -> Exact:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an int or Fraction")
    return Fraction(value)


def _poly(values: Sequence[object], *, name: str) -> Polynomial:
    if not values:
        return (Fraction(0),)
    coeffs = tuple(_exact(v, name=f"{name}[{i}]") for i, v in enumerate(values))
    end = len(coeffs)
    while end > 1 and coeffs[end - 1] == 0:
        end -= 1
    return coeffs[:end]


def poly_eval(coeffs: Sequence[object], r: object) -> Exact:
    p = _poly(coeffs, name="coeffs")
    x = _exact(r, name="r")
    acc = Fraction(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc


def poly_scale(coeffs: Sequence[object], factor: object) -> Polynomial:
    p = _poly(coeffs, name="coeffs")
    a = _exact(factor, name="factor")
    return _poly(tuple(a * c for c in p), name="scaled")


def poly_sub(left: Sequence[object], right: Sequence[object]) -> Polynomial:
    a = _poly(left, name="left")
    b = _poly(right, name="right")
    n = max(len(a), len(b))
    out = []
    for i in range(n):
        out.append((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0))
    return _poly(out, name="difference")


def weighted_moment(
    coeffs: Sequence[object], *, e: int, lower: object, upper: object
) -> Exact:
    """Exact integral of r**e * polynomial(r) on [lower, upper]."""
    if e not in (1, 2):
        raise ValueError("e must be 1 or 2")
    lo = _exact(lower, name="lower")
    hi = _exact(upper, name="upper")
    if not lo > 0 or not hi > lo:
        raise ValueError("require 0 < lower < upper")
    p = _poly(coeffs, name="coeffs")
    total = Fraction(0)
    for degree, coeff in enumerate(p):
        power = degree + e + 1
        total += coeff * (hi**power - lo**power) / power
    return total


def weighted_primitive(
    coeffs: Sequence[object], *, e: int, lower: object, radius: object
) -> Exact:
    """Exact integral of s**e * polynomial(s) from lower to radius."""
    if e not in (1, 2):
        raise ValueError("e must be 1 or 2")
    lo = _exact(lower, name="lower")
    r = _exact(radius, name="radius")
    if not lo > 0 or not r >= lo:
        raise ValueError("require radius >= lower > 0")
    p = _poly(coeffs, name="coeffs")
    total = Fraction(0)
    for degree, coeff in enumerate(p):
        power = degree + e + 1
        total += coeff * (r**power - lo**power) / power
    return total


@dataclass(frozen=True)
class RadialMomentInverseResult:
    e: int
    source_moment: Exact
    profile_moment: Exact
    projected_moment: Exact
    radius: Exact
    projected_source_value: Exact
    weighted_primitive_value: Exact
    stress_value: Exact
    stress_radial_derivative: Exact
    differential_value: Exact
    expected_differential_value: Exact
    inverse_residual: Exact
    upper_weighted_primitive: Exact

    @property
    def verified(self) -> bool:
        return (
            self.profile_moment == 1
            and self.projected_moment == 0
            and self.inverse_residual == 0
            and self.upper_weighted_primitive == 0
        )


def replay_radial_moment_inverse(
    *,
    e: int,
    source: Sequence[object],
    normalized_profile: Sequence[object],
    lower: object,
    upper: object,
    radius: object,
) -> RadialMomentInverseResult:
    """Replay the R6--R12 moment complement and radial differential identity.

    The inputs are exact polynomial fixtures on one positive radial interval.
    `normalized_profile` must have weighted moment exactly one.  This is an
    algebraic audit interface only: it does not assert smooth zero-extension of
    the polynomial fixture or materialize a paper stage/stress witness.
    """
    if e not in (1, 2):
        raise ValueError("e must be 1 or 2")
    lo = _exact(lower, name="lower")
    hi = _exact(upper, name="upper")
    r = _exact(radius, name="radius")
    if not lo > 0 or not hi > lo:
        raise ValueError("require 0 < lower < upper")
    if not (lo <= r <= hi):
        raise ValueError("radius must lie in [lower, upper]")

    f = _poly(source, name="source")
    b = _poly(normalized_profile, name="normalized_profile")
    source_moment = weighted_moment(f, e=e, lower=lo, upper=hi)
    profile_moment = weighted_moment(b, e=e, lower=lo, upper=hi)
    if profile_moment != 1:
        raise ValueError("normalized_profile must have exact weighted moment 1")

    projected = poly_sub(f, poly_scale(b, source_moment))
    projected_moment = weighted_moment(projected, e=e, lower=lo, upper=hi)
    h_value = poly_eval(projected, r)
    primitive = weighted_primitive(projected, e=e, lower=lo, radius=r)
    stress = -(r ** (-e)) * primitive
    stress_r = e * (r ** (-e - 1)) * primitive - h_value
    differential = stress_r + Fraction(e, 1) * stress / r
    expected = -poly_eval(f, r) + poly_eval(b, r) * source_moment
    upper_primitive = weighted_primitive(projected, e=e, lower=lo, radius=hi)

    return RadialMomentInverseResult(
        e=e,
        source_moment=source_moment,
        profile_moment=profile_moment,
        projected_moment=projected_moment,
        radius=r,
        projected_source_value=h_value,
        weighted_primitive_value=primitive,
        stress_value=stress,
        stress_radial_derivative=stress_r,
        differential_value=differential,
        expected_differential_value=expected,
        inverse_residual=differential - expected,
        upper_weighted_primitive=upper_primitive,
    )


def normalize_profile(
    profile: Sequence[object], *, e: int, lower: object, upper: object
) -> Polynomial:
    """Normalize an exact polynomial fixture to weighted moment one."""
    p = _poly(profile, name="profile")
    moment = weighted_moment(p, e=e, lower=lower, upper=upper)
    if moment == 0:
        raise ValueError("profile weighted moment must be nonzero")
    return poly_scale(p, Fraction(1, 1) / moment)
