"""Clean-room exact replay of one public finite-stage gain calculation.

The public Kokuno workbench does not expose a repository-level reuse license.
This module therefore does not copy ``proof_sources/stage9/exact_checks.py``.
It independently evaluates the short Stage-9 exponent formulas recorded in the
public reader with exact rational arithmetic.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

SOURCE_RECORD_COMMIT = "e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f"
SOURCE_WORKBENCH_COMMIT = "fab69fdc4ac197159b8e6ae8d73a82bde2b20d55"
SOURCE_WORKBENCH_BLOB_SHA1 = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_BUNDLE_SHA256 = "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"
STAGE9_BODY_SHA256 = "bce636fb3ed7348ed9c73912183f81a75557cc975418ce9e162c4be994f82cb0"
STAGE9_CHECKER_PATH = "proof_sources/stage9/exact_checks.py"
STAGE9_CHECKER_SHA256 = "228935772385483e92c5dcac118290243135ba467541abba5e4f345b9fc513c1"
STAGE9_RESULT_SHA256 = "a54a985c47a5da769f68f7dc04d1ff59ee46025371cd6f6a852369d8898d349f"
SOURCE_DOI = "10.5281/zenodo.22678406"

SOURCE_KAPPA_S = Fraction(1, 100_000)
COMMON_GAIN = Fraction(1, 10)
RADIAL_GAIN_FLOOR = Fraction(17, 100)


@dataclass(frozen=True)
class Stage9GainMarginReplay:
    stage: int
    kappa_s: Fraction
    sigma_j: Fraction
    b_j: Fraction
    c_j_star: Fraction
    h_j: Fraction
    stage_step: Fraction
    transport_radial_gain: Fraction
    curl_radial_gain: Fraction
    self_product_radial_gain: Fraction
    mean_to_wave_gain: Fraction
    auxiliary_gain: Fraction
    pressure_gain: Fraction
    high_gain: Fraction
    tolerance: Fraction

    @property
    def published_constants_exact(self) -> bool:
        if self.kappa_s != SOURCE_KAPPA_S:
            return False
        return (
            self.transport_radial_gain == Fraction(8_999, 50_000)
            and self.curl_radial_gain == Fraction(49_997, 100_000)
            and self.mean_to_wave_gain == Fraction(24_999, 50_000)
            and self.auxiliary_gain == Fraction(39_999, 100_000)
            and self.pressure_gain == Fraction(24_999, 25_000)
            and self.high_gain == Fraction(22_499, 25_000)
        )

    @property
    def stage_relations_exact(self) -> bool:
        return (
            self.sigma_j == Fraction(1, 5) + Fraction(self.stage, 10)
            and self.b_j == Fraction(1, 2) + self.sigma_j
            and self.c_j_star == 1 + self.sigma_j
            and self.c_j_star == self.b_j + Fraction(1, 2)
            and self.h_j == self.c_j_star - 2 * self.kappa_s
            and self.stage_step == Fraction(1, 10)
        )

    @property
    def radial_margin_floor_holds(self) -> bool:
        return (
            self.transport_radial_gain > RADIAL_GAIN_FLOOR
            and self.curl_radial_gain > RADIAL_GAIN_FLOOR
            and self.self_product_radial_gain > RADIAL_GAIN_FLOOR
        )

    @property
    def common_gain_holds(self) -> bool:
        return all(
            margin > COMMON_GAIN
            for margin in (
                self.transport_radial_gain,
                self.curl_radial_gain,
                self.self_product_radial_gain,
                self.mean_to_wave_gain,
                self.auxiliary_gain,
                self.pressure_gain,
                self.high_gain,
            )
        )

    @property
    def passed(self) -> bool:
        return (
            self.tolerance == 0
            and self.published_constants_exact
            and self.stage_relations_exact
            and self.radial_margin_floor_holds
            and self.common_gain_holds
        )


def replay_stage9_gain_margins(
    stage: int = 0,
    *,
    kappa_s: Fraction = SOURCE_KAPPA_S,
) -> Stage9GainMarginReplay:
    """Evaluate the displayed Stage-9 gain arithmetic exactly.

    ``stage`` is the finite correction index ``j``.  The public reader fixes
    ``kappa_s = 10^-5`` and ``sigma_j = 1/5 + j/10``.  A custom exact
    ``kappa_s`` is accepted only for negative/mutation tests; ``passed`` stays
    false unless the published value is used.
    """

    if isinstance(stage, bool) or not isinstance(stage, int) or stage < 0:
        raise TypeError("stage must be a nonnegative integer")
    if not isinstance(kappa_s, Fraction):
        raise TypeError("kappa_s must be fractions.Fraction")
    if kappa_s <= 0:
        raise ValueError("kappa_s must be positive")

    sigma_j = Fraction(1, 5) + Fraction(stage, 10)
    b_j = Fraction(1, 2) + sigma_j
    c_j_star = 1 + sigma_j
    h_j = c_j_star - 2 * kappa_s
    sigma_next = Fraction(1, 5) + Fraction(stage + 1, 10)
    b_next = Fraction(1, 2) + sigma_next

    return Stage9GainMarginReplay(
        stage=stage,
        kappa_s=kappa_s,
        sigma_j=sigma_j,
        b_j=b_j,
        c_j_star=c_j_star,
        h_j=h_j,
        stage_step=b_next - b_j,
        transport_radial_gain=Fraction(9, 50) - 2 * kappa_s,
        curl_radial_gain=Fraction(1, 2) - 3 * kappa_s,
        self_product_radial_gain=sigma_j - 3 * kappa_s,
        mean_to_wave_gain=h_j - b_j,
        auxiliary_gain=Fraction(2, 5) - kappa_s,
        pressure_gain=1 - 4 * kappa_s,
        high_gain=Fraction(9, 10) - 4 * kappa_s,
        tolerance=Fraction(0),
    )
