from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from typing import TypeAlias

Exact: TypeAlias = int | Fraction


def _q(value: Exact) -> Fraction:
    if isinstance(value, bool) or isinstance(value, (float, Decimal)):
        raise TypeError("exact checker inputs must be int or Fraction")
    if not isinstance(value, (int, Fraction)):
        raise TypeError("exact checker inputs must be int or Fraction")
    return Fraction(value)


@dataclass(frozen=True)
class HeatJet:
    h: Fraction
    Z: Fraction
    H: Fraction
    H1: Fraction
    H2: Fraction

    @classmethod
    def exact(cls, h: Exact, Z: Exact, H: Exact, H1: Exact, H2: Exact) -> "HeatJet":
        return cls(*map(_q, (h, Z, H, H1, H2)))

    @property
    def A(self) -> Fraction:
        return Fraction(1, 2) + self.h


def heat_ode_residual(jet: HeatJet) -> Fraction:
    h, Z, H, H1, H2 = jet.h, jet.Z, jet.H, jet.H1, jet.H2
    return Z * Z * H2 + (((2 + 2 * h) * Z) + 1) * H1 + h * (1 + h) * H


def time_core(jet: HeatJet) -> Fraction:
    return -2 * jet.H1


def cylindrical_swirl_core(jet: HeatJet) -> Fraction:
    A, Z, H, H1, H2 = jet.A, jet.Z, jet.H, jet.H1, jet.H2
    return (
        (2 * A * A - Fraction(1, 2)) * H
        + (4 * A + 2) * Z * H1
        + 2 * Z * Z * H2
    )


def heat_equation_residual(jet: HeatJet) -> Fraction:
    return cylindrical_swirl_core(jet) - time_core(jet)


def replay_identity(jet: HeatJet) -> dict[str, Fraction | bool]:
    ode = heat_ode_residual(jet)
    residual = heat_equation_residual(jet)
    return {
        "ode_residual": ode,
        "time_core": time_core(jet),
        "cylindrical_swirl_core": cylindrical_swirl_core(jet),
        "heat_equation_residual": residual,
        "bridge_residual": residual - 2 * ode,
        "passes": ode == 0 and residual == 0,
    }


PAPER_EXACT_VELOCITY_AVAILABLE = False
FULL_RECONSTRUCTION = False
