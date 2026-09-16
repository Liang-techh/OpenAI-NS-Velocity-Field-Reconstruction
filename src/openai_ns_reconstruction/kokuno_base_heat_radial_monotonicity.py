"""Exact clean-room checks for the public NS base-heat radial sign identity.

The public heat representation has a positive kernel and, after differentiating
with respect to its similarity variable Z, the combination ``A H + Z H'``
has a pointwise positive factor.  This module verifies only that algebraic
factor and its exact positivity margin.  It does not prove existence of the
imported heat profile or the analytic steps needed to justify its integral.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import TypeAlias

Exact: TypeAlias = int | Fraction


def _q(value: Exact, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an int or Fraction")
    return Fraction(value)


def _source_inputs(h: Exact, z: Exact, v: Exact) -> tuple[Fraction, Fraction, Fraction]:
    h_q = _q(h, name="h")
    z_q = _q(z, name="z")
    v_q = _q(v, name="v")
    if h_q <= 0:
        raise ValueError("h must be positive in the source heat-profile regime")
    if z_q < 0:
        raise ValueError("z must be nonnegative")
    if v_q < 0:
        raise ValueError("v must be nonnegative")
    return h_q, z_q, v_q


def heat_similarity_exponent(h: Exact) -> Fraction:
    """Return the public exponent ``A = 1/2 + h`` exactly."""
    h_q = _q(h, name="h")
    if h_q <= 0:
        raise ValueError("h must be positive in the source heat-profile regime")
    return Fraction(1, 2) + h_q


def local_log_slope(h: Exact, z: Exact, v: Exact) -> Fraction:
    """Return ``h Z v / (1 + Z v)`` for one positive heat-kernel node."""
    h_q, z_q, v_q = _source_inputs(h, z, v)
    x = z_q * v_q
    return h_q * x / (1 + x)


def combined_factor_direct(h: Exact, z: Exact, v: Exact) -> Fraction:
    """Return the factor ``A - h Z v/(1+Zv)`` in ``A H + Z H'``."""
    h_q, z_q, v_q = _source_inputs(h, z, v)
    return heat_similarity_exponent(h_q) - local_log_slope(h_q, z_q, v_q)


def combined_factor_positive_form(h: Exact, z: Exact, v: Exact) -> Fraction:
    """Return the equivalent manifestly positive factor.

    With ``A=1/2+h``, exact algebra gives
      A - hZv/(1+Zv) = (A + Zv/2)/(1+Zv)
                        = 1/2 + h/(1+Zv).
    """
    h_q, z_q, v_q = _source_inputs(h, z, v)
    x = z_q * v_q
    return Fraction(1, 2) + h_q / (1 + x)


def combined_factor_identity_residual(h: Exact, z: Exact, v: Exact) -> Fraction:
    """Return direct-minus-positive-form residual; exact inputs should give 0."""
    return combined_factor_direct(h, z, v) - combined_factor_positive_form(h, z, v)


@dataclass(frozen=True)
class RadialSignCertificate:
    """Exact pointwise sign data used by the public radial-monotonicity argument."""

    factor: Fraction
    half_margin: Fraction
    outward_bracket: Fraction


def radial_sign_certificate(h: Exact, z: Exact, v: Exact) -> RadialSignCertificate:
    """Certify the pointwise sign before integration against the positive kernel.

    ``factor = A - hZv/(1+Zv)`` is the relative integrand in ``A H + ZH'``.
    It is strictly larger than 1/2 by ``h/(1+Zv)``.  Therefore the relative
    integrand in ``-A H - ZH'`` is strictly negative.
    """
    h_q, z_q, v_q = _source_inputs(h, z, v)
    factor = combined_factor_positive_form(h_q, z_q, v_q)
    margin = factor - Fraction(1, 2)
    if margin <= 0:
        raise AssertionError("positive-kernel margin must be strictly positive")
    return RadialSignCertificate(
        factor=factor,
        half_margin=margin,
        outward_bracket=-factor,
    )
