"""Exact clean-room replay of the public MC15/MC27 chart-scaling seam."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


def _fraction(value: Fraction, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Fraction):
        raise TypeError(f"{name} must be fractions.Fraction")
    return value


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be int")
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


@dataclass(frozen=True)
class MeanChartScaling:
    h: Fraction
    A: Fraction
    D: Fraction
    weighted_moment_exponent: Fraction
    radial_primitive_exponent: Fraction
    velocity_exponent: Fraction
    residual_exponent: Fraction
    pressure_exponent: Fraction
    potential_exponent: Fraction
    normalized_negative_dz_potential_exponent: Fraction
    normalized_negative_dz_potential_sign: int
    mc27_amplitude_ratio_exponent: Fraction
    mc27_argument_ratio_exponent: Fraction


def mean_chart_scaling(*, h: Fraction, moment_a: Fraction, moment_e: int) -> MeanChartScaling:
    """Return the exact exponent ledger behind public formulas MC15 and MC27.

    No source implementation is copied.  The formulas replayed are:
      integral R^e f_* dR = Q^(a-(e+1)/2) integral r^e f_phys dr,
      T_e f_* = Q^(a-1/2) ...,
      u_* = Q^A u_phys, g_* = Q^(2A+1/2) g_phys,
      p_* = Q^(2A) p_phys, Psi_* = Q^(A-1/2) Psi_phys,
    with A=1/2+h, D=1/2-h, and the normalized -dz Psi term carrying
    exponent 1/2-D=h with the minus sign preserved.

    MC27 uses q^(-A) physical corrections.  In a fixed-Q chart this gives
    an amplitude factor (Q/q)^A and profile argument
    r/sqrt(q)=sqrt(Q/q) R.
    """
    h = _fraction(h, "h")
    moment_a = _fraction(moment_a, "moment_a")
    moment_e = _index(moment_e, "moment_e")
    if not (Fraction(0) < h < Fraction(1, 100)):
        raise ValueError("source construction requires 0 < h < 1/100")

    A = Fraction(1, 2) + h
    D = Fraction(1, 2) - h
    weighted = moment_a - Fraction(moment_e + 1, 2)
    primitive = moment_a - Fraction(1, 2)

    return MeanChartScaling(
        h=h,
        A=A,
        D=D,
        weighted_moment_exponent=weighted,
        radial_primitive_exponent=primitive,
        velocity_exponent=A,
        residual_exponent=2 * A + Fraction(1, 2),
        pressure_exponent=2 * A,
        potential_exponent=A - Fraction(1, 2),
        normalized_negative_dz_potential_exponent=Fraction(1, 2) - D,
        normalized_negative_dz_potential_sign=-1,
        mc27_amplitude_ratio_exponent=A,
        mc27_argument_ratio_exponent=Fraction(1, 2),
    )


def source_mean_moment_exponents(h: Fraction) -> tuple[Fraction, Fraction]:
    """Return the public chart/physical exponents for M_theta and M_z."""
    scaling_theta = mean_chart_scaling(h=h, moment_a=Fraction(1, 2) + h, moment_e=2)
    scaling_z = mean_chart_scaling(h=h, moment_a=Fraction(1, 2) + h, moment_e=1)
    return scaling_theta.weighted_moment_exponent, scaling_z.weighted_moment_exponent


def verify_negative_dz_potential(*, h: Fraction, reported_sign: int, reported_exponent: Fraction) -> bool:
    """Fail closed unless the normalized -dz Psi sign and exponent are exact."""
    if isinstance(reported_sign, bool) or not isinstance(reported_sign, int):
        raise TypeError("reported_sign must be int")
    reported_exponent = _fraction(reported_exponent, "reported_exponent")
    scaling = mean_chart_scaling(h=h, moment_a=Fraction(1), moment_e=0)
    return (
        reported_sign == scaling.normalized_negative_dz_potential_sign
        and reported_exponent == scaling.normalized_negative_dz_potential_exponent
    )


def verify_mc27_ratio_exponents(
    *, h: Fraction, amplitude_exponent: Fraction, argument_exponent: Fraction
) -> bool:
    """Check the fixed-Q chart exponents induced by the public MC27 physical field."""
    amplitude_exponent = _fraction(amplitude_exponent, "amplitude_exponent")
    argument_exponent = _fraction(argument_exponent, "argument_exponent")
    scaling = mean_chart_scaling(h=h, moment_a=Fraction(1), moment_e=0)
    return (
        amplitude_exponent == scaling.mc27_amplitude_ratio_exponent
        and argument_exponent == scaling.mc27_argument_ratio_exponent
    )
