"""Exact chi-threshold bridge from the actual SchedulePressure low-|Z| margin.

The pinned ``NaturalAxisRange.exists_sigma`` proof chooses

    sigma = sqrt(m) / 20

after establishing ``m <= H^2`` on ``|Z| <= j/10``.  Therefore, on that
same region,

    chi = H^2 / (H^2 + sigma^2) >= 400/401 > 99/100.

The repository already materializes the actual-schedule analytic low-|Z|
margin and the theorem-side cutoff choice.  This module closes only the
remaining exact-arithmetic handoff: it binds that actual witness to the exact
rational chi lower bound needed downstream by the NaturalEntrance log-slope
gate.  It does not sample eta, infer a margin numerically, or materialize any
coefficient-space fixed point/profile field.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math

from .outgoing_tail import TailData
from .schedule_axis_margin import (
    ScheduleLowZMarginWitness,
    certify_schedule_low_Z_margin,
)


SIGMA_SQUARED_MARGIN_RATIO = Fraction(1, 400)
CHI_EXACT_LOWER = Fraction(400, 401)
CHI_REQUIRED_LOWER = Fraction(99, 100)
CHI_STRICT_SLACK = CHI_EXACT_LOWER - CHI_REQUIRED_LOWER


@dataclass(frozen=True)
class ActualScheduleChiCutoffCertificate:
    """Fail-closed exact chi implication for the actual schedule low-|Z| set."""

    margin_witness: ScheduleLowZMarginWitness
    sigma_squared_margin_ratio: Fraction
    chi_lower_exact: Fraction
    chi_required_lower: Fraction

    def __post_init__(self) -> None:
        witness = self.margin_witness
        if not isinstance(witness, ScheduleLowZMarginWitness):
            raise TypeError("margin_witness must be ScheduleLowZMarginWitness")
        if not math.isfinite(witness.margin) or witness.margin <= 0.0:
            raise ValueError("actual-schedule H^2 margin must be finite and positive")
        if witness.cutoff.margin != witness.margin:
            raise ValueError("cutoff/margin witness identity mismatch")
        if witness.cutoff.delta != witness.parameters.cutoff_delta:
            raise ValueError("cutoff delta is not the pinned j/10 choice")
        expected_sigma = math.sqrt(witness.margin) / 20.0
        if witness.cutoff.sigma != expected_sigma:
            raise ValueError("cutoff sigma is not the pinned sqrt(margin)/20 choice")
        if witness.cutoff.guaranteed_chi_lower_bound != 400.0 / 401.0:
            raise ValueError("runtime cutoff chi lower metadata drifted from 400/401")

        if not isinstance(self.sigma_squared_margin_ratio, Fraction):
            raise TypeError("sigma_squared_margin_ratio must be an exact Fraction")
        if not isinstance(self.chi_lower_exact, Fraction):
            raise TypeError("chi_lower_exact must be an exact Fraction")
        if not isinstance(self.chi_required_lower, Fraction):
            raise TypeError("chi_required_lower must be an exact Fraction")
        if self.sigma_squared_margin_ratio != SIGMA_SQUARED_MARGIN_RATIO:
            raise ValueError("sigma^2/margin ratio must equal the pinned 1/400")
        if self.chi_lower_exact != CHI_EXACT_LOWER:
            raise ValueError("chi exact lower bound must equal 400/401")
        if self.chi_required_lower != CHI_REQUIRED_LOWER:
            raise ValueError("required chi threshold must equal 99/100")
        if not self.chi_required_lower < self.chi_lower_exact:
            raise ValueError("actual-schedule chi bound must strictly exceed 99/100")

    @property
    def low_Z_delta(self) -> float:
        """The actual theorem region width, exactly the witness' ``j/10`` choice."""

        return self.margin_witness.cutoff.delta

    @property
    def sigma(self) -> float:
        """The actual schedule cutoff sigma selected from the analytic margin."""

        return self.margin_witness.cutoff.sigma

    @property
    def strict_slack(self) -> Fraction:
        """Exact rational slack ``400/401 - 99/100 = 301/40100``."""

        return self.chi_lower_exact - self.chi_required_lower

    @property
    def implication(self) -> str:
        """Human-readable theorem implication carried by this certificate."""

        return "|Z(eta)| <= j/10 -> 99/100 < chi(h,j,sigma,eta)"

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def full_reconstruction(self) -> bool:
        return False

    def require_exact_downstream_bound(self, bound: Fraction) -> None:
        """Admit only an exact rational lower bound justified by this certificate."""

        if not isinstance(bound, Fraction):
            raise TypeError("downstream chi lower bound must be an exact Fraction")
        if bound < self.chi_required_lower:
            raise ValueError("downstream chi lower bound is below the pinned 99/100 threshold")
        if bound > self.chi_lower_exact:
            raise ValueError("downstream chi lower bound exceeds the certified 400/401 bound")


def actual_schedule_chi_cutoff_certificate(
    data: TailData,
    j: float,
) -> ActualScheduleChiCutoffCertificate:
    """Construct the exact chi threshold bridge from actual SchedulePressure data.

    No margin, sigma, delta, or chi bound is caller supplied.  The actual
    schedule margin constructor owns those choices; this factory only attaches
    exact rational algebra to that already-certified analytic witness.
    """

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")
    witness = certify_schedule_low_Z_margin(data, j)
    return ActualScheduleChiCutoffCertificate(
        margin_witness=witness,
        sigma_squared_margin_ratio=SIGMA_SQUARED_MARGIN_RATIO,
        chi_lower_exact=CHI_EXACT_LOWER,
        chi_required_lower=CHI_REQUIRED_LOWER,
    )
