from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


def _fraction(name: str, value: Fraction) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Fraction):
        raise TypeError(f"{name} must be fractions.Fraction")
    return value


def _positive(name: str, value: Fraction) -> Fraction:
    value = _fraction(name, value)
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _exact_nth_root(value: int, n: int) -> int:
    if value < 0 or n < 1:
        raise ValueError("exact root requires value >= 0 and n >= 1")
    if value in (0, 1) or n == 1:
        return value
    lo, hi = 0, 1
    while hi**n < value:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        power = mid**n
        if power < value:
            lo = mid
        else:
            hi = mid
    if hi**n != value:
        raise ValueError(f"{value} is not an exact {n}-th power")
    return hi


def exact_fraction_power(base: Fraction, exponent: Fraction) -> Fraction:
    base = _positive("base", base)
    exponent = _fraction("exponent", exponent)
    root_degree = exponent.denominator
    numerator_root = _exact_nth_root(base.numerator, root_degree)
    denominator_root = _exact_nth_root(base.denominator, root_degree)
    rooted = Fraction(numerator_root, denominator_root)
    return rooted ** exponent.numerator


@dataclass(frozen=True)
class BaseHeatCoordinateBridge:
    h: Fraction
    A: Fraction
    d: Fraction
    tau: Fraction
    s: Fraction
    z_profile: Fraction
    z_physical: Fraction
    profile_value: Fraction
    physical_value: Fraction
    z_residual: Fraction
    value_residual: Fraction


def verify_base_heat_coordinate_bridge(
    *,
    q: Fraction,
    X: Fraction,
    eta: Fraction,
    d: Fraction,
    tau: Fraction,
    s: Fraction,
    h: Fraction,
    c_inf: Fraction,
    H_value: Fraction,
) -> BaseHeatCoordinateBridge:
    q = _positive("q", q)
    X = _positive("X", X)
    eta = _fraction("eta", eta)
    d = _fraction("d", d)
    tau = _positive("tau", tau)
    s = _positive("s", s)
    h = _fraction("h", h)
    c_inf = _fraction("c_inf", c_inf)
    H_value = _fraction("H_value", H_value)

    if not Fraction(0) < h < Fraction(1, 100):
        raise ValueError("source construction requires 0 < h < 1/100")
    if abs(eta) > 1:
        raise ValueError("source profile endpoint domain requires |eta| <= 1")

    expected_d = 1 - eta * eta
    if d != expected_d:
        raise ValueError("d must equal 1-eta^2 exactly")
    if tau != q * d:
        raise ValueError("tau must equal q*d exactly")
    if s != q * X:
        raise ValueError("source chart requires s=q*X exactly, with s=r^2/2")

    A = Fraction(1, 2) + h
    z_profile = 2 * d / X
    z_physical = 2 * tau / s
    z_residual = z_profile - z_physical
    if z_residual != 0:
        raise ArithmeticError("heat arguments disagree")

    profile_power = exact_fraction_power(q, -A) * exact_fraction_power(X, -A)
    physical_power = exact_fraction_power(s, -A)
    profile_value = c_inf * profile_power * H_value
    physical_value = c_inf * physical_power * H_value
    value_residual = profile_value - physical_value
    if value_residual != 0:
        raise ArithmeticError("heat amplitudes disagree")

    return BaseHeatCoordinateBridge(
        h=h,
        A=A,
        d=d,
        tau=tau,
        s=s,
        z_profile=z_profile,
        z_physical=z_physical,
        profile_value=profile_value,
        physical_value=physical_value,
        z_residual=z_residual,
        value_residual=value_residual,
    )
