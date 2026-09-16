"""Clean-room exact replay of the public profile-assembly PA.14 density increments.

This module intentionally contains no Kokuno source implementation.  It encodes only
short public mathematical formulas using exact rational arithmetic so that downstream
typed assembly code can fail closed on missing quadratic/cross terms.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from typing import Any

Exact = Fraction


def _exact(name: str, value: Any) -> Exact:
    """Accept theorem-facing rational inputs only; reject approximate/nonnumeric data."""
    if isinstance(value, bool) or isinstance(value, (float, Decimal)):
        raise TypeError(f"{name} must be an exact rational, not {type(value).__name__}")
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, Fraction):
        return value
    raise TypeError(f"{name} must be int or Fraction")


@dataclass(frozen=True)
class FiveMomentDensity:
    """Pointwise density tuple for the five profile assembly constraints."""

    M: Exact
    I: Exact
    J: Exact
    S: Exact
    C_p: Exact

    def __sub__(self, other: "FiveMomentDensity") -> "FiveMomentDensity":
        return FiveMomentDensity(
            self.M - other.M,
            self.I - other.I,
            self.J - other.J,
            self.S - other.S,
            self.C_p - other.C_p,
        )


@dataclass(frozen=True)
class ProfileAssemblyInput:
    """Exact local data sufficient for the finite PA.14 algebraic seam.

    ``radial_factor`` represents the displayed ``sqrt(2 X)`` factor.  Requiring
    ``radial_factor**2 == 2*X`` keeps the replay exact and avoids silently replacing
    an irrational square root by a floating approximation.
    """

    X: Exact
    radial_factor: Exact
    U: Exact
    E: Exact
    u: Exact
    e: Exact

    @classmethod
    def exact(
        cls,
        *,
        X: Any,
        radial_factor: Any,
        U: Any,
        E: Any,
        u: Any,
        e: Any,
    ) -> "ProfileAssemblyInput":
        values = cls(
            _exact("X", X),
            _exact("radial_factor", radial_factor),
            _exact("U", U),
            _exact("E", E),
            _exact("u", u),
            _exact("e", e),
        )
        if values.X <= 0:
            raise ValueError("X must be strictly positive for the C_p density")
        if values.radial_factor <= 0:
            raise ValueError("radial_factor must be the positive square root of 2*X")
        if values.radial_factor * values.radial_factor != 2 * values.X:
            raise ValueError("radial_factor must satisfy radial_factor**2 == 2*X exactly")
        return values


def five_moment_density(*, X: Any, radial_factor: Any, U: Any, E: Any) -> FiveMomentDensity:
    """Return the five unperturbed pointwise densities used by the PA.14 difference."""
    data = ProfileAssemblyInput.exact(
        X=X, radial_factor=radial_factor, U=U, E=E, u=0, e=0
    )
    R = data.radial_factor
    return FiveMomentDensity(
        M=data.U,
        I=R * data.E,
        J=data.U * R * data.E,
        S=data.U * data.U - data.E * data.E / 2,
        C_p=data.E * data.E / (2 * data.X),
    )


def pa14_increment_formula(data: ProfileAssemblyInput) -> FiveMomentDensity:
    """Evaluate the exact five perturbation-density increments from public PA.14."""
    X, R, U, E, u, e = (
        data.X,
        data.radial_factor,
        data.U,
        data.E,
        data.u,
        data.e,
    )
    H = R * E
    return FiveMomentDensity(
        M=u,
        I=R * e,
        J=H * u + U * R * e + R * u * e,
        S=2 * U * u + u * u - E * e - e * e / 2,
        C_p=E * e / X + e * e / (2 * X),
    )


def direct_density_difference(data: ProfileAssemblyInput) -> FiveMomentDensity:
    """Compute the same increment independently as new density minus old density."""
    before = five_moment_density(
        X=data.X, radial_factor=data.radial_factor, U=data.U, E=data.E
    )
    after = five_moment_density(
        X=data.X,
        radial_factor=data.radial_factor,
        U=data.U + data.u,
        E=data.E + data.e,
    )
    return after - before


def pa14_residual(data: ProfileAssemblyInput) -> FiveMomentDensity:
    """Exact formula-minus-direct residual; every field must be zero for a replay PASS."""
    return pa14_increment_formula(data) - direct_density_difference(data)


def all_zero(value: FiveMomentDensity) -> bool:
    return all(
        item == 0
        for item in (value.M, value.I, value.J, value.S, value.C_p)
    )
