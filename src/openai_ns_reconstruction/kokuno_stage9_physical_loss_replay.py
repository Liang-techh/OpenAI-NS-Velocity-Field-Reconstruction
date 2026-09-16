"""Exact clean-room replay of one public Stage-9 exponent conversion.

This module does not copy the Kokuno checker.  It encodes only the short
public formulas needed to audit the physical residual loss and stage gain.
All theorem inputs are exact ``fractions.Fraction`` values.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


class Stage9ReplayMismatch(ValueError):
    """Raised when a supplied claimed source constant is not exact."""


def _require_fraction(name: str, value: Fraction) -> Fraction:
    if type(value) is not Fraction:
        raise TypeError(f"{name} must be fractions.Fraction, got {type(value).__name__}")
    return value


def _require_index(name: str, value: int) -> int:
    if type(value) is not int or value < 0:
        raise TypeError(f"{name} must be a nonnegative int")
    return value


@dataclass(frozen=True)
class Stage9PhysicalLossReceipt:
    """Exact exponent bookkeeping for one derivative order and finite stage."""

    h: Fraction
    m: int
    j: int
    A: Fraction
    physical_residual_loss: Fraction
    spatial_derivative_loss: Fraction
    temporal_derivative_loss: Fraction
    fixed_loss_Km: Fraction
    sigma_j: Fraction
    stage_gain: Fraction
    next_stage_gain: Fraction
    gain_increment: Fraction
    net_power: Fraction


def replay_stage9_physical_loss(
    h: Fraction,
    m: int,
    j: int,
    *,
    claimed_A: Fraction | None = None,
    claimed_spatial_loss: Fraction | None = None,
    claimed_temporal_loss: Fraction | None = None,
) -> Stage9PhysicalLossReceipt:
    """Replay the displayed physical-loss/gain arithmetic with zero tolerance.

    The public Stage-9 formulas use ``A=1/2+h``, the physical residual factor
    ``Q^(-2A-1/2)``, spatial derivative cost ``1/2+h/2``, time derivative
    cost ``1+3h/2``, and ``sigma_j=1/5+j/10``.  Since the time cost dominates
    the spatial cost for positive ``h``, a sufficient total-order-``m`` fixed
    loss is ``K_m=2A+1/2+m(1+3h/2)``, which is independent of ``j``.

    Optional claimed constants are accepted only when they equal the public
    values exactly; this makes transcription drift fail closed.
    """

    h = _require_fraction("h", h)
    m = _require_index("m", m)
    j = _require_index("j", j)
    if not (Fraction(0) < h < Fraction(1, 100)):
        raise ValueError("source construction requires 0 < h < 1/100")

    A = Fraction(1, 2) + h
    spatial_loss = Fraction(1, 2) + h / 2
    temporal_loss = Fraction(1) + 3 * h / 2

    claims = (
        ("A", claimed_A, A),
        ("spatial derivative loss", claimed_spatial_loss, spatial_loss),
        ("temporal derivative loss", claimed_temporal_loss, temporal_loss),
    )
    for name, claimed, expected in claims:
        if claimed is None:
            continue
        claimed = _require_fraction(f"claimed {name}", claimed)
        if claimed != expected:
            raise Stage9ReplayMismatch(
                f"{name} mismatch: claimed {claimed}, expected {expected}, "
                f"residual {claimed - expected}"
            )

    # Exact dominance needed for the common total-order loss.
    if temporal_loss < spatial_loss:
        raise AssertionError("public time-derivative loss must dominate spatial loss")

    physical_residual_loss = 2 * A + Fraction(1, 2)
    fixed_loss = physical_residual_loss + m * temporal_loss
    sigma_j = Fraction(1, 5) + Fraction(j, 10)
    stage_gain = h * sigma_j
    next_stage_gain = h * (Fraction(1, 5) + Fraction(j + 1, 10))

    return Stage9PhysicalLossReceipt(
        h=h,
        m=m,
        j=j,
        A=A,
        physical_residual_loss=physical_residual_loss,
        spatial_derivative_loss=spatial_loss,
        temporal_derivative_loss=temporal_loss,
        fixed_loss_Km=fixed_loss,
        sigma_j=sigma_j,
        stage_gain=stage_gain,
        next_stage_gain=next_stage_gain,
        gain_increment=next_stage_gain - stage_gain,
        net_power=stage_gain - fixed_loss,
    )


def underestimated_time_loss_defect(h: Fraction, m: int) -> Fraction:
    """Return the exact deficit from replacing ``1+3h/2`` by ``1+h``."""

    h = _require_fraction("h", h)
    m = _require_index("m", m)
    if h <= 0:
        raise ValueError("h must be positive")
    return m * ((Fraction(1) + 3 * h / 2) - (Fraction(1) + h))
