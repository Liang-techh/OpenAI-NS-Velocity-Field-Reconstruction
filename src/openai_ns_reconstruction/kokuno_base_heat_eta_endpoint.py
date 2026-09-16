"""Exact clean-room endpoint-parameter identities for the Kokuno base-heat factor.

This module intentionally implements only short public formulas from the
Kokuno Navier--Stokes reader.  It does not construct the imported heat profile
or prove differentiation under the defining integral.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


def _fraction(value: Fraction, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Fraction):
        raise TypeError(f"{name} must be fractions.Fraction")
    return value


def _positive(value: Fraction, name: str) -> Fraction:
    value = _fraction(value, name)
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _nonnegative_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be int")
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def rising_factorial(x: Fraction, order: int) -> Fraction:
    """Return (x)_order exactly."""
    x = _fraction(x, "x")
    order = _nonnegative_int(order, "order")
    out = Fraction(1)
    for k in range(order):
        out *= x + k
    return out


def heat_origin_derivative(h: Fraction, order: int) -> Fraction:
    """Public origin jet H^(m)(0)=(-1)^m (h)_m (1+h)_m."""
    h = _positive(h, "h")
    order = _nonnegative_int(order, "order")
    if h >= Fraction(1, 100):
        raise ValueError("source construction requires 0 < h < 1/100")
    sign = -1 if order % 2 else 1
    return sign * rising_factorial(h, order) * rising_factorial(1 + h, order)


@dataclass(frozen=True)
class HeatEtaEndpointJet:
    """One-sided eta jet of H(2(1-eta^2)/X) at eta=+/-1."""

    h: Fraction
    X: Fraction
    eta_endpoint: int
    Z: Fraction
    Z_eta: Fraction
    Z_etaeta: Fraction
    H: Fraction
    H_eta: Fraction
    H_etaeta: Fraction


def heat_eta_endpoint_jet(h: Fraction, X: Fraction, eta_endpoint: int) -> HeatEtaEndpointJet:
    """Compute the exact order-two endpoint eta jet from the public formulas."""
    h = _positive(h, "h")
    X = _positive(X, "X")
    if h >= Fraction(1, 100):
        raise ValueError("source construction requires 0 < h < 1/100")
    if isinstance(eta_endpoint, bool) or eta_endpoint not in (-1, 1):
        raise ValueError("eta_endpoint must be -1 or +1")

    eta = Fraction(eta_endpoint)
    Z = Fraction(0)
    Z_eta = -4 * eta / X
    Z_etaeta = -4 / X
    H0 = Fraction(1)
    H1 = heat_origin_derivative(h, 1)
    H2 = heat_origin_derivative(h, 2)
    H_eta = H1 * Z_eta
    H_etaeta = H2 * Z_eta * Z_eta + H1 * Z_etaeta
    return HeatEtaEndpointJet(
        h=h,
        X=X,
        eta_endpoint=eta_endpoint,
        Z=Z,
        Z_eta=Z_eta,
        Z_etaeta=Z_etaeta,
        H=H0,
        H_eta=H_eta,
        H_etaeta=H_etaeta,
    )


def verify_endpoint_pair(h: Fraction, X: Fraction) -> tuple[HeatEtaEndpointJet, HeatEtaEndpointJet]:
    """Fail closed unless the +/- endpoint jets obey exact parity identities."""
    minus = heat_eta_endpoint_jet(h, X, -1)
    plus = heat_eta_endpoint_jet(h, X, 1)
    if minus.H != plus.H:
        raise AssertionError("endpoint values must agree")
    if minus.H_eta != -plus.H_eta:
        raise AssertionError("first eta derivatives must be odd across endpoints")
    if minus.H_etaeta != plus.H_etaeta:
        raise AssertionError("second eta derivatives must agree")
    if plus.Z_eta != -minus.Z_eta:
        raise AssertionError("Z_eta parity mismatch")
    if plus.Z_etaeta != minus.Z_etaeta:
        raise AssertionError("Z_etaeta parity mismatch")
    return minus, plus
