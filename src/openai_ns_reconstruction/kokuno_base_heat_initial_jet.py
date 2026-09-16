"""Exact clean-room checks for the public NS base-heat initial derivative jet.

Only the short public Gamma/Laplace formula and its algebraic Gamma-moment
consequence are encoded here.  This module does not prove that the imported
heat profile exists in the source construction or that it has the other
regularity/support properties required by the full reconstruction.
"""

from __future__ import annotations

from fractions import Fraction
from typing import TypeAlias

Exact: TypeAlias = int | Fraction


def _q(value: Exact, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an int or Fraction")
    return Fraction(value)


def _nonnegative_int(value: int, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an int")
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def rising_factorial(start: Exact, order: int) -> Fraction:
    """Return ``(start)_order`` using exact rational arithmetic."""
    start_q = _q(start, name="start")
    order_i = _nonnegative_int(order, name="order")
    product = Fraction(1)
    for j in range(order_i):
        product *= start_q + j
    return product


def heat_initial_derivative(h: Exact, order: int) -> Fraction:
    """Return the public integral's exact derivative ``H^(order)(0)``.

    From
      H(Z) = Gamma(1+h)^(-1) int exp(-v) v^h (1+Zv)^(-h) dv,
    differentiation followed by the Gamma moment at Z=0 gives
      H^(k)(0) = (-1)^k (h)_k (1+h)_k.

    The source regime uses positive h; this routine therefore fails closed
    outside that regime rather than pretending to certify a broader one.
    """
    h_q = _q(h, name="h")
    order_i = _nonnegative_int(order, name="order")
    if h_q <= 0:
        raise ValueError("h must be positive in the source heat-profile regime")
    sign = -1 if order_i % 2 else 1
    return sign * rising_factorial(h_q, order_i) * rising_factorial(1 + h_q, order_i)


def heat_initial_jet(h: Exact, order: int) -> tuple[Fraction, ...]:
    """Return ``(H(0), H'(0), ..., H^(order)(0))`` exactly."""
    order_i = _nonnegative_int(order, name="order")
    return tuple(heat_initial_derivative(h, k) for k in range(order_i + 1))


def ode_taylor_recurrence_residual(h: Exact, order: int) -> Fraction:
    """Check the coefficient recurrence induced by the public base-heat ODE.

    For
      Z^2 H'' + (1 + 2(1+h)Z) H' + h(1+h) H = 0,
    the Z^k Taylor coefficient gives
      H^(k+1)(0) + (k+h)(k+h+1) H^k(0) = 0.
    """
    h_q = _q(h, name="h")
    order_i = _nonnegative_int(order, name="order")
    current = heat_initial_derivative(h_q, order_i)
    following = heat_initial_derivative(h_q, order_i + 1)
    return following + (order_i + h_q) * (order_i + h_q + 1) * current


def supplied_recurrence_residual(
    h: Exact, order: int, current: Exact, following: Exact
) -> Fraction:
    """Evaluate the same recurrence on caller-supplied jet entries."""
    h_q = _q(h, name="h")
    order_i = _nonnegative_int(order, name="order")
    if h_q <= 0:
        raise ValueError("h must be positive in the source heat-profile regime")
    current_q = _q(current, name="current")
    following_q = _q(following, name="following")
    return following_q + (order_i + h_q) * (order_i + h_q + 1) * current_q
