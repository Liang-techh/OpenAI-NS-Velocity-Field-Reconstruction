"""Clean-room exact audit of the public terminal endpoint scaling adapters.

This module intentionally implements only the short scaling formulas exposed in
Kokuno's public corrected reader.  It does not import or copy source-bundle
checker/proof code and it does not certify endpoint existence or full closure.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class EndpointScalingSignature:
    """Multiplicative bookkeeping for one endpoint adapter.

    ``x_argument`` and ``t_argument`` are the factors multiplying the input
    arguments.  ``period`` and ``singular_time`` are factors on the represented
    solution's period and singular time.  ``energy_sq`` is the factor on
    ``||u||_2^2``.
    """

    velocity: Fraction
    pressure: Fraction
    force: Fraction
    x_argument: Fraction
    t_argument: Fraction
    momentum_residual: Fraction
    divergence: Fraction
    energy_sq: Fraction
    period: Fraction
    singular_time: Fraction


def _exact_positive(value: Fraction, *, name: str) -> Fraction:
    if type(value) is not Fraction:
        raise TypeError(f"{name} must be fractions.Fraction")
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def viscosity_spatial_signature(*, nu: Fraction, sqrt_nu: Fraction) -> EndpointScalingSignature:
    """Return the exact public ``A_nu`` scaling signature.

    The caller supplies an exact rational square root; irrational square roots
    are outside this finite exact audit rather than approximated.
    """

    nu = _exact_positive(nu, name="nu")
    sqrt_nu = _exact_positive(sqrt_nu, name="sqrt_nu")
    if sqrt_nu * sqrt_nu != nu:
        raise ValueError("sqrt_nu must square exactly to nu")
    return EndpointScalingSignature(
        velocity=sqrt_nu,
        pressure=nu,
        force=sqrt_nu,
        x_argument=Fraction(1, 1) / sqrt_nu,
        t_argument=Fraction(1, 1),
        momentum_residual=sqrt_nu,
        divergence=Fraction(1, 1),
        energy_sq=nu * nu * sqrt_nu,
        period=sqrt_nu,
        singular_time=Fraction(1, 1),
    )


def parabolic_signature(*, lam: Fraction) -> EndpointScalingSignature:
    """Return the exact fixed-viscosity parabolic ``S_lambda`` signature."""

    lam = _exact_positive(lam, name="lam")
    return EndpointScalingSignature(
        velocity=lam,
        pressure=lam * lam,
        force=lam * lam * lam,
        x_argument=lam,
        t_argument=lam * lam,
        momentum_residual=lam * lam * lam,
        divergence=lam * lam,
        energy_sq=Fraction(1, 1) / lam,
        period=Fraction(1, 1) / lam,
        singular_time=Fraction(1, 1) / (lam * lam),
    )


def fixed_period_viscosity_signature(*, nu: Fraction) -> EndpointScalingSignature:
    """Return the exact public ``T_nu`` signature."""

    nu = _exact_positive(nu, name="nu")
    return EndpointScalingSignature(
        velocity=nu,
        pressure=nu * nu,
        force=nu * nu,
        x_argument=Fraction(1, 1),
        t_argument=nu,
        momentum_residual=nu * nu,
        divergence=nu,
        energy_sq=nu * nu,
        period=Fraction(1, 1),
        singular_time=Fraction(1, 1) / nu,
    )


def compose_signatures(first: EndpointScalingSignature, second: EndpointScalingSignature) -> EndpointScalingSignature:
    """Compose ``second`` after ``first`` for multiplicative bookkeeping."""

    if type(first) is not EndpointScalingSignature or type(second) is not EndpointScalingSignature:
        raise TypeError("first and second must be EndpointScalingSignature")
    return EndpointScalingSignature(
        velocity=first.velocity * second.velocity,
        pressure=first.pressure * second.pressure,
        force=first.force * second.force,
        x_argument=first.x_argument * second.x_argument,
        t_argument=first.t_argument * second.t_argument,
        momentum_residual=first.momentum_residual * second.momentum_residual,
        divergence=first.divergence * second.divergence,
        energy_sq=first.energy_sq * second.energy_sq,
        period=first.period * second.period,
        singular_time=first.singular_time * second.singular_time,
    )


def audit_fixed_period_identity(*, nu: Fraction, sqrt_nu: Fraction) -> tuple[EndpointScalingSignature, EndpointScalingSignature]:
    """Return ``S_sqrt(nu) o A_nu`` and ``T_nu`` after exact validation."""

    spatial = viscosity_spatial_signature(nu=nu, sqrt_nu=sqrt_nu)
    parabolic = parabolic_signature(lam=sqrt_nu)
    composed = compose_signatures(spatial, parabolic)
    direct = fixed_period_viscosity_signature(nu=nu)
    if composed != direct:
        raise AssertionError("public endpoint adapter composition does not close exactly")
    return composed, direct
