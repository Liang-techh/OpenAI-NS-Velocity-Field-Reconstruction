"""Independent exact replay of one public Kokuno finite-stage identity.

The Kokuno workbench does not expose a clear repository reuse license.  This
module therefore does not copy ``proof_sources/stage9/exact_checks.py``.  It
independently reconstructs the two endpoint scaling identities displayed for
the ``l = -1`` hold in the public workbench.
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
SOURCE_DOI = "10.5281/zenodo.22678406"


@dataclass(frozen=True)
class Stage9HoldReplay:
    """Exact exponent bookkeeping for the stage-9 ``l=-1`` hold."""

    l_value: Fraction
    hold_length_log_inverse_h_coefficient: Fraction
    xe2_log_rate: Fraction
    e_log_rate: Fraction
    xe2_endpoint_h_power: Fraction
    e_endpoint_h_power: Fraction
    tolerance: Fraction

    @property
    def rates_exact(self) -> bool:
        return (
            self.xe2_log_rate == 2 * self.l_value
            and self.e_log_rate == self.l_value - Fraction(1, 2)
        )

    @property
    def endpoint_powers_exact(self) -> bool:
        return (
            self.xe2_endpoint_h_power == Fraction(8)
            and self.e_endpoint_h_power == Fraction(6)
        )

    @property
    def passed(self) -> bool:
        return self.tolerance == 0 and self.rates_exact and self.endpoint_powers_exact


def replay_stage9_l_minus_one_hold(
    *, hold_length_log_inverse_h_coefficient: Fraction = Fraction(4)
) -> Stage9HoldReplay:
    """Return the exact endpoint powers for the public ``l=-1`` hold.

    The workbench specifies a hold of length ``4 log(1/h)`` and records

    ``d/dy log(X E^2) = 2 l = -2``
    ``d/dy log(E) = l - 1/2 = -3/2``.

    If a quantity has log-rate ``r`` over ``c log(1/h)``, its endpoint ratio
    is ``h**(-r*c)``.  The exponent multiplication below is exact rational
    arithmetic; no logarithm or floating-point approximation is evaluated.
    """

    coefficient = Fraction(hold_length_log_inverse_h_coefficient)
    l_value = Fraction(-1)
    xe2_rate = 2 * l_value
    e_rate = l_value - Fraction(1, 2)

    return Stage9HoldReplay(
        l_value=l_value,
        hold_length_log_inverse_h_coefficient=coefficient,
        xe2_log_rate=xe2_rate,
        e_log_rate=e_rate,
        xe2_endpoint_h_power=-xe2_rate * coefficient,
        e_endpoint_h_power=-e_rate * coefficient,
        tolerance=Fraction(0),
    )


def concrete_endpoint_ratios(h: Fraction) -> tuple[Fraction, Fraction]:
    """Evaluate the replayed endpoint ratios for an exact rational ``0<h<1``."""

    h = Fraction(h)
    if not 0 < h < 1:
        raise ValueError("h must satisfy 0 < h < 1")
    replay = replay_stage9_l_minus_one_hold()
    if replay.xe2_endpoint_h_power.denominator != 1 or replay.e_endpoint_h_power.denominator != 1:
        raise AssertionError("expected integral endpoint powers")
    return (
        h ** replay.xe2_endpoint_h_power.numerator,
        h ** replay.e_endpoint_h_power.numerator,
    )
