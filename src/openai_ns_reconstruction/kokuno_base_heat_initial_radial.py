"""Exact clean-room checks for the NS base-heat radial normalization seam.

This module implements only short public formulas needed to regression-test the
chain rule from ``s=r^2`` and the ``H(0)=1`` normalization.  It does not encode
or prove existence of the heat profile H.
"""

from __future__ import annotations

from fractions import Fraction
from typing import TypeAlias

Exact: TypeAlias = int | Fraction


def _q(value: Exact, *, name: str) -> Fraction:
    """Coerce an exact rational input and reject inexact/boolean inputs."""
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an int or Fraction")
    return Fraction(value)


def heat_s_radial_core(
    A: Exact, Z: Exact, H: Exact, H_prime: Exact
) -> Fraction:
    """Return the normalized ``s * dK/ds`` core: ``-(A H + Z H')``."""
    A_q = _q(A, name="A")
    Z_q = _q(Z, name="Z")
    H_q = _q(H, name="H")
    Hp_q = _q(H_prime, name="H_prime")
    return -(A_q * H_q + Z_q * Hp_q)


def heat_r_radial_core(
    A: Exact, Z: Exact, H: Exact, H_prime: Exact
) -> Fraction:
    """Return ``r*d_r K/(c_inf*s^-A)`` for ``s=r^2`` and ``Z=2*tau/s``."""
    return 2 * heat_s_radial_core(A, Z, H, H_prime)


def base_swirl_radial_core(h: Exact) -> Fraction:
    """Return ``r*d_r b/b`` for ``b proportional to r^-(1+2h)``."""
    h_q = _q(h, name="h")
    return -(1 + 2 * h_q)


def initial_radial_seam_residual(h: Exact, H0: Exact = 1) -> Fraction:
    """Compare the heat radial core at ``Z=0`` with base-swirl homogeneity.

    A is fixed by the public relation ``A=1/2+h``.  The seam closes exactly
    when the imported heat-profile normalization ``H(0)=1`` is supplied.
    """
    h_q = _q(h, name="h")
    H0_q = _q(H0, name="H0")
    A_q = Fraction(1, 2) + h_q
    heat_core_at_zero = heat_r_radial_core(A_q, 0, H0_q, 0)
    return heat_core_at_zero - base_swirl_radial_core(h_q)
