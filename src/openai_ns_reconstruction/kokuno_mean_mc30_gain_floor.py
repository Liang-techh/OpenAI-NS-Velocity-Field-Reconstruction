"""Clean-room exact replay of the exposed MC30 gain-floor arithmetic.

This module intentionally implements only short public exponent formulas.  It
is not a copy of the Kokuno checker and it is not a complete mean-correction
verification.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Mapping

SOURCE_KAPPA_S = Fraction(1, 100000)
SOURCE_ALPHA_MIN = Fraction(9, 10)

_EXPOSED_KEYS = (
    "rg_time_epsilon_dT_B",
    "rg_radial_2bB",
    "rg_radial_2betaB",
    "rg_radial_B2",
    "rg_axial_bd",
    "mc29_product_alpha_plus_09",
    "mc29_product_2alpha",
)


def _require_fraction(name: str, value: Fraction) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Fraction):
        raise TypeError(f"{name} must be fractions.Fraction")
    return value


@dataclass(frozen=True)
class MC30ExposedGainFloor:
    alpha: Fraction
    kappa_s: Fraction
    target_exponent: Fraction
    term_exponents: tuple[tuple[str, Fraction], ...]
    margins: tuple[tuple[str, Fraction], ...]

    @property
    def minimum_margin(self) -> Fraction:
        return min(margin for _, margin in self.margins)


def source_exposed_mc30_exponents(alpha: Fraction, kappa_s: Fraction = SOURCE_KAPPA_S) -> dict[str, Fraction]:
    """Return the public exposed MC30/MC29 exponent rows exactly.

    The source construction fixes kappa_s=10^-5 and assumes alpha>=0.9.
    """

    alpha = _require_fraction("alpha", alpha)
    kappa_s = _require_fraction("kappa_s", kappa_s)
    if kappa_s != SOURCE_KAPPA_S:
        raise ValueError("kappa_s must equal the published source value 10^-5")
    if alpha < SOURCE_ALPHA_MIN:
        raise ValueError("MC30 source hypothesis requires alpha >= 9/10")

    return {
        "rg_time_epsilon_dT_B": alpha + 2,
        "rg_radial_2bB": alpha + 2 - kappa_s,
        "rg_radial_2betaB": alpha + Fraction(29, 10) - kappa_s,
        "rg_radial_B2": 2 * alpha + 2 - kappa_s,
        "rg_axial_bd": alpha + 2,
        "mc29_product_alpha_plus_09": alpha + Fraction(9, 10),
        "mc29_product_2alpha": 2 * alpha,
    }


def verify_exposed_mc30_gain_floor(
    alpha: Fraction,
    term_exponents: Mapping[str, Fraction],
    kappa_s: Fraction = SOURCE_KAPPA_S,
) -> MC30ExposedGainFloor:
    """Verify the exposed finite term ledger against the MC30 target floor.

    Every named exponent is supplied explicitly so a lost power/loss factor is
    detectable.  The required target is alpha + 0.9 - 2*kappa_s.
    """

    alpha = _require_fraction("alpha", alpha)
    kappa_s = _require_fraction("kappa_s", kappa_s)
    if kappa_s != SOURCE_KAPPA_S:
        raise ValueError("kappa_s must equal the published source value 10^-5")
    if alpha < SOURCE_ALPHA_MIN:
        raise ValueError("MC30 source hypothesis requires alpha >= 9/10")
    if set(term_exponents) != set(_EXPOSED_KEYS):
        missing = sorted(set(_EXPOSED_KEYS) - set(term_exponents))
        extra = sorted(set(term_exponents) - set(_EXPOSED_KEYS))
        raise ValueError(f"term ledger keys mismatch; missing={missing}, extra={extra}")

    checked: list[tuple[str, Fraction]] = []
    for name in _EXPOSED_KEYS:
        checked.append((name, _require_fraction(name, term_exponents[name])))

    target = alpha + Fraction(9, 10) - 2 * kappa_s
    margins = tuple((name, exponent - target) for name, exponent in checked)
    failed = [(name, margin) for name, margin in margins if margin < 0]
    if failed:
        raise ValueError(f"MC30 target exponent is not met: {failed}")

    return MC30ExposedGainFloor(
        alpha=alpha,
        kappa_s=kappa_s,
        target_exponent=target,
        term_exponents=tuple(checked),
        margins=margins,
    )


def replay_exposed_mc30_gain_floor(alpha: Fraction, kappa_s: Fraction = SOURCE_KAPPA_S) -> MC30ExposedGainFloor:
    """Compute and verify the short public source formulas in one exact replay."""

    rows = source_exposed_mc30_exponents(alpha, kappa_s)
    return verify_exposed_mc30_gain_floor(alpha, rows, kappa_s)
