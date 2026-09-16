"""Clean-room exact replay of the public NS-profiles stress primitives.

The public Kokuno workbench does not expose a repository-level reuse license.
This module therefore copies no source checker implementation.  It independently
re-evaluates the short integrating-factor identities displayed in the public
reader, using exact rational polynomial arithmetic.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable

SOURCE_RECORD_COMMIT = "e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f"
SOURCE_WORKBENCH_COMMIT = "fab69fdc4ac197159b8e6ae8d73a82bde2b20d55"
SOURCE_WORKBENCH_BLOB_SHA1 = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_BUNDLE_SHA256 = "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"
PROFILES_BODY_SHA256 = "63327f4a6d339d39de096230af1810cb83e7c427ba0b0cc78c576b2750437a95"
PROFILES_CHECKER_PATH = "proof_sources/profiles/exact_checks.py"
PROFILES_CHECKER_SHA256 = "a7ee77024a92d06c1b69055a00b10a7d5fd600ca424fe084e2f4399310cf56e7"
PROFILES_RESULT_SHA256 = "5ef2b5f7cee750088143acda8d4a5ee4f0b12f36f0344ab991f5c9333a842e24"
SOURCE_DOI = "10.5281/zenodo.22678406"

Polynomial = tuple[Fraction, ...]


def _as_poly(coefficients: Iterable[Fraction], *, name: str) -> Polynomial:
    out = tuple(coefficients)
    if not out:
        raise ValueError(f"{name} must contain at least one coefficient")
    if any(not isinstance(c, Fraction) for c in out):
        raise TypeError(f"{name} coefficients must be fractions.Fraction")
    return out


def _eval(poly: Polynomial, x: Fraction) -> Fraction:
    acc = Fraction(0)
    for coefficient in reversed(poly):
        acc = acc * x + coefficient
    return acc


def _derivative(poly: Polynomial) -> Polynomial:
    if len(poly) == 1:
        return (Fraction(0),)
    return tuple(Fraction(i) * poly[i] for i in range(1, len(poly)))


def _multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return tuple(out)


def _integral_from_zero(poly: Polynomial) -> Polynomial:
    return (Fraction(0),) + tuple(poly[i] / Fraction(i + 1) for i in range(len(poly)))


@dataclass(frozen=True)
class ProfileStressPrimitiveReplay:
    x: Fraction
    h_value: Fraction
    h_prime: Fraction
    s_q: Fraction
    s_n: Fraction
    l_value: Fraction
    q_s: Fraction
    d_x_q_s: Fraction
    q_ode_residual: Fraction
    n_s: Fraction
    d_x_n_s: Fraction
    n_ode_residual: Fraction
    q_axis_limit: Fraction
    q_axis_expected: Fraction
    n_axis_limit: Fraction
    n_axis_expected: Fraction
    q_integration_constant: Fraction
    n_integration_constant: Fraction
    tolerance: Fraction

    @property
    def regular_axis_constants(self) -> bool:
        return self.q_integration_constant == 0 and self.n_integration_constant == 0

    @property
    def axis_limits_exact(self) -> bool:
        return self.q_axis_limit == self.q_axis_expected and self.n_axis_limit == self.n_axis_expected

    @property
    def ode_identities_exact(self) -> bool:
        return self.q_ode_residual == 0 and self.n_ode_residual == 0

    @property
    def passed(self) -> bool:
        return (
            self.tolerance == 0
            and self.ode_identities_exact
            and self.axis_limits_exact
            and self.regular_axis_constants
        )


def replay_profile_stress_primitives(
    x: Fraction,
    h_coefficients: Iterable[Fraction],
    s_q_coefficients: Iterable[Fraction],
    s_n_coefficients: Iterable[Fraction],
    *,
    q_integration_constant: Fraction = Fraction(0),
    n_integration_constant: Fraction = Fraction(0),
) -> ProfileStressPrimitiveReplay:
    """Replay the displayed regular stress primitives at one exact ``X>0``.

    The reader defines ``D_X = X d/dX`` and gives

      ``D_X Q_s + (1+l) Q_s = S_q``,  ``l = X H'/H``,
      ``Q_s = (int_0^X H S_q dx)/(X H)``,

    together with

      ``D_X N_s + N_s = S_n``,
      ``N_s = (int_0^X S_n dx)/X``.

    A nonzero integration constant still solves the positive-X ODE but produces
    the forbidden axis singularity.  ``passed`` therefore requires both exact
    ODE residuals and the regular zero-constant choice.
    """

    if not isinstance(x, Fraction):
        raise TypeError("x must be fractions.Fraction")
    if x <= 0:
        raise ValueError("x must be positive")
    if not isinstance(q_integration_constant, Fraction) or not isinstance(
        n_integration_constant, Fraction
    ):
        raise TypeError("integration constants must be fractions.Fraction")

    h_poly = _as_poly(h_coefficients, name="h_coefficients")
    sq_poly = _as_poly(s_q_coefficients, name="s_q_coefficients")
    sn_poly = _as_poly(s_n_coefficients, name="s_n_coefficients")

    if len(h_poly) < 2 or h_poly[0] != 0 or h_poly[1] == 0:
        raise ValueError("H must have an axis-compatible simple zero: H(0)=0 and H_X(0)!=0")

    hp_poly = _derivative(h_poly)
    h_value = _eval(h_poly, x)
    if h_value == 0:
        raise ValueError("H(x) must be nonzero")
    h_prime = _eval(hp_poly, x)
    s_q = _eval(sq_poly, x)
    s_n = _eval(sn_poly, x)
    l_value = x * h_prime / h_value

    q_integral_poly = _integral_from_zero(_multiply(h_poly, sq_poly))
    q_integral = _eval(q_integral_poly, x) + q_integration_constant
    q_integral_prime = h_value * s_q
    q_denominator = x * h_value
    q_denominator_prime = h_value + x * h_prime
    q_s = q_integral / q_denominator
    d_x_q_s = (
        q_integral_prime * q_denominator - q_integral * q_denominator_prime
    ) / (q_denominator * q_denominator)
    q_ode_residual = x * d_x_q_s + (1 + l_value) * q_s - s_q

    n_integral_poly = _integral_from_zero(sn_poly)
    n_integral = _eval(n_integral_poly, x) + n_integration_constant
    n_s = n_integral / x
    d_x_n_s = (s_n * x - n_integral) / (x * x)
    n_ode_residual = x * d_x_n_s + n_s - s_n

    # These limits come from the independently integrated coefficient series.
    # H = h1 X + O(X^2), so int(H S_q)=h1*S_q(0) X^2/2+O(X^3)
    # while XH=h1 X^2+O(X^3).
    h1 = h_poly[1]
    q_integral_x2 = _multiply(h_poly, sq_poly)[1] / 2
    q_axis_limit = q_integral_x2 / h1
    q_axis_expected = sq_poly[0] / 2
    n_axis_limit = sn_poly[0]
    n_axis_expected = sn_poly[0]

    return ProfileStressPrimitiveReplay(
        x=x,
        h_value=h_value,
        h_prime=h_prime,
        s_q=s_q,
        s_n=s_n,
        l_value=l_value,
        q_s=q_s,
        d_x_q_s=d_x_q_s,
        q_ode_residual=q_ode_residual,
        n_s=n_s,
        d_x_n_s=d_x_n_s,
        n_ode_residual=n_ode_residual,
        q_axis_limit=q_axis_limit,
        q_axis_expected=q_axis_expected,
        n_axis_limit=n_axis_limit,
        n_axis_expected=n_axis_expected,
        q_integration_constant=q_integration_constant,
        n_integration_constant=n_integration_constant,
        tolerance=Fraction(0),
    )
