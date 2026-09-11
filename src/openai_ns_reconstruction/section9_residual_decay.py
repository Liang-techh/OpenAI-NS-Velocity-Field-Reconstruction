"""Exact exponent bookkeeping for the conditional Section 9 iteration.

Section 9 uses the schedule

    sigma_j = 1/5 + j/10

and gains ``g_j = h*j/10`` for the accumulated velocity corrections and
``h*sigma_j`` for the residual.  This module keeps those exponents in exact
rational arithmetic and computes the first iteration index at which a supplied
derivative loss is overcome.

This is deliberately only a ``formal-structure`` ledger.  In particular, it
does not prove the hypotheses behind Eqs. (9.17)-(9.18), does not manufacture
the constants ``C_{j,m}``, ``K_m`` or the flat remainders ``E_{j,m}``, and does
not construct the locally finite sums in Eq. (9.21).  A caller-supplied loss is
therefore a theorem datum, never a fitted residual estimate, and no result from
this module is a paper-exact velocity certificate.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math


SIGMA_0 = Fraction(1, 5)
SIGMA_STEP = Fraction(1, 10)

Scalar = int | float | Fraction


def _as_fraction(value: Scalar, name: str) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real number")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return Fraction(str(value))


def _nonnegative_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _positive_h(h: Scalar) -> Fraction:
    hq = _as_fraction(h, "h")
    if not Fraction(0) < hq < Fraction(1, 2):
        raise ValueError("h must satisfy 0 < h < 1/2 (similarity-chart domain)")
    return hq


def _ceil_fraction(value: Fraction) -> int:
    """Exact ceiling with no float conversion."""
    return -((-value.numerator) // value.denominator)


def section9_sigma(stage: int) -> Fraction:
    """The exact Section 9 schedule sigma_j = 1/5 + j/10."""
    stage = _nonnegative_integer(stage, "stage")
    return SIGMA_0 + stage * SIGMA_STEP


def section9_increment_gain(stage: int, h: Scalar) -> Fraction:
    """Return g_j = h*j/10 from Eq. (9.17)."""
    stage = _nonnegative_integer(stage, "stage")
    hq = _positive_h(h)
    return hq * stage * SIGMA_STEP


def section9_residual_gain(stage: int, h: Scalar) -> Fraction:
    """Return h*sigma_j, the leading residual gain in Eq. (9.18)."""
    return _positive_h(h) * section9_sigma(stage)


def section9_increment_derivative_loss(derivative_order: int, h: Scalar) -> Fraction:
    """Return one admissible ell_m used with Eq. (9.17).

    The Section 9 estimate permits

        ell_m = 2*A + (m+1)*(1 + 3*h/2),  A = 1/2 + h.
    """
    derivative_order = _nonnegative_integer(derivative_order, "derivative_order")
    hq = _positive_h(h)
    A = Fraction(1, 2) + hq
    return 2 * A + (derivative_order + 1) * (1 + Fraction(3, 2) * hq)


def first_stage_for_increment_power(
    h: Scalar,
    derivative_order: int,
    target_power: Scalar,
) -> int:
    """First j for which ``g_j - ell_m >= target_power``.

    This is only exponent arithmetic.  It is conditional on the actual
    Section 9 correction satisfying Eq. (9.17).
    """
    hq = _positive_h(h)
    target = _as_fraction(target_power, "target_power")
    if target < 0:
        raise ValueError("target_power must be nonnegative")
    loss = section9_increment_derivative_loss(derivative_order, hq)
    return max(0, _ceil_fraction(10 * (target + loss) / hq))


def first_stage_for_residual_power(
    h: Scalar,
    derivative_loss: Scalar,
    target_power: Scalar,
) -> int:
    """First j for which ``h*sigma_j - K_m >= target_power``.

    ``derivative_loss`` is the supplied ``K_m`` from Eq. (9.18).  The function
    never estimates or fits it from residual samples.
    """
    hq = _positive_h(h)
    loss = _as_fraction(derivative_loss, "derivative_loss")
    target = _as_fraction(target_power, "target_power")
    if loss < 0:
        raise ValueError("derivative_loss must be nonnegative")
    if target < 0:
        raise ValueError("target_power must be nonnegative")
    threshold = 10 * ((target + loss) / hq - SIGMA_0)
    return max(0, _ceil_fraction(threshold))


@dataclass(frozen=True)
class Section9ResidualPowerRecord:
    """Fail-closed record for one conditional Eq. (9.18) power threshold."""

    derivative_order: int
    h: Fraction
    derivative_loss: Fraction
    target_power: Fraction
    first_stage: int
    sigma: Fraction
    residual_gain: Fraction
    leading_exponent: Fraction
    status: str = "formal-structure"
    actual_residual_bound_verified: bool = False
    flat_remainder_verified: bool = False
    paper_exact_velocity_available: bool = False


def section9_residual_power_record(
    *,
    derivative_order: int,
    h: Scalar,
    derivative_loss: Scalar,
    target_power: Scalar,
) -> Section9ResidualPowerRecord:
    """Build an auditable conditional threshold record for Eq. (9.18)."""
    derivative_order = _nonnegative_integer(derivative_order, "derivative_order")
    hq = _positive_h(h)
    loss = _as_fraction(derivative_loss, "derivative_loss")
    target = _as_fraction(target_power, "target_power")
    stage = first_stage_for_residual_power(hq, loss, target)
    sigma = section9_sigma(stage)
    gain = hq * sigma
    return Section9ResidualPowerRecord(
        derivative_order=derivative_order,
        h=hq,
        derivative_loss=loss,
        target_power=target,
        first_stage=stage,
        sigma=sigma,
        residual_gain=gain,
        leading_exponent=gain - loss,
    )
