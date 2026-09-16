"""Clean-room exact replay of the public Kokuno MC29 pressure-moment sign seam.

This module intentionally implements only short public mathematical identities.
It does not copy the source checker/proof implementation and does not establish
full mean-correction or stage convergence.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Sequence

Exact = Fraction | int
Poly = tuple[Fraction, ...]


def _fraction(value: Exact) -> Fraction:
    if isinstance(value, bool):
        raise TypeError("bool is not an exact polynomial coefficient")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    raise TypeError("coefficients must be int or Fraction")


def polynomial(coefficients: Iterable[Exact]) -> Poly:
    values = tuple(_fraction(value) for value in coefficients)
    if not values:
        return (Fraction(0),)
    end = len(values)
    while end > 1 and values[end - 1] == 0:
        end -= 1
    return values[:end]


def add(left: Sequence[Exact], right: Sequence[Exact]) -> Poly:
    a = polynomial(left)
    b = polynomial(right)
    n = max(len(a), len(b))
    return polynomial(
        (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
        for i in range(n)
    )


def scale(poly: Sequence[Exact], scalar: Exact) -> Poly:
    factor = _fraction(scalar)
    return polynomial(factor * coefficient for coefficient in polynomial(poly))


def multiply(left: Sequence[Exact], right: Sequence[Exact]) -> Poly:
    a = polynomial(left)
    b = polynomial(right)
    out = [Fraction(0) for _ in range(len(a) + len(b) - 1)]
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return polynomial(out)


def radial_shift(poly: Sequence[Exact], power: int = 1) -> Poly:
    if isinstance(power, bool) or not isinstance(power, int) or power < 0:
        raise TypeError("power must be a nonnegative integer")
    return polynomial((Fraction(0),) * power + polynomial(poly))


def integrate_unit_interval(poly: Sequence[Exact]) -> Fraction:
    """Return the exact integral over R in [0,1]."""
    p = polynomial(poly)
    return sum((coefficient / Fraction(index + 1) for index, coefficient in enumerate(p)), Fraction(0))


def pressure_source_linear_moment(v_over_r: Sequence[Exact], delta_v: Sequence[Exact]) -> Fraction:
    """Replay -1/2 ∫ R^2 (2 V/R) Δv dR exactly.

    ``v_over_r`` is V/R, so this avoids introducing a spurious division at R=0.
    """
    source = scale(multiply(v_over_r, delta_v), 2)
    return -Fraction(1, 2) * integrate_unit_interval(radial_shift(source, 2))


def mc24_fifth_row_linear(
    axial_base: Sequence[Exact],
    axial_correction: Sequence[Exact],
    v_over_r: Sequence[Exact],
    delta_v: Sequence[Exact],
) -> Fraction:
    """Return ∫(2 R G d - R V u)dR with V=R(V/R)."""
    gd = scale(multiply(axial_base, axial_correction), 2)
    first = integrate_unit_interval(radial_shift(gd, 1))
    v_delta = multiply(v_over_r, delta_v)
    second = integrate_unit_interval(radial_shift(v_delta, 2))
    return first - second


def mc29_linear_jz(
    axial_base: Sequence[Exact],
    axial_correction: Sequence[Exact],
    v_over_r: Sequence[Exact],
    delta_v: Sequence[Exact],
) -> Fraction:
    """Linear part of MC29's new J_z moment."""
    gd = scale(multiply(axial_base, axial_correction), 2)
    axial_linear = integrate_unit_interval(radial_shift(gd, 1))
    return axial_linear + pressure_source_linear_moment(v_over_r, delta_v)


def mc29_quadratic_defect(
    axial_correction: Sequence[Exact],
    rg_nonlinear: Sequence[Exact],
) -> Fraction:
    """The terms retained beyond the five-row linear cancellation.

    This is ∫R d^2 dR - 1/2 ∫R^2 R_g,nonlinear dR.
    """
    d2 = multiply(axial_correction, axial_correction)
    first = integrate_unit_interval(radial_shift(d2, 1))
    second = Fraction(1, 2) * integrate_unit_interval(radial_shift(rg_nonlinear, 2))
    return first - second


def mc29_full_jz(
    axial_base: Sequence[Exact],
    axial_correction: Sequence[Exact],
    v_over_r: Sequence[Exact],
    delta_v: Sequence[Exact],
    rg_nonlinear: Sequence[Exact],
) -> Fraction:
    return mc29_linear_jz(axial_base, axial_correction, v_over_r, delta_v) + mc29_quadratic_defect(
        axial_correction, rg_nonlinear
    )


def verify_pressure_moment_sign(
    axial_base: Sequence[Exact],
    axial_correction: Sequence[Exact],
    v_over_r: Sequence[Exact],
    delta_v: Sequence[Exact],
) -> bool:
    return mc29_linear_jz(axial_base, axial_correction, v_over_r, delta_v) == mc24_fifth_row_linear(
        axial_base, axial_correction, v_over_r, delta_v
    )


paper_exact_velocity_available = False
full_reconstruction = False
