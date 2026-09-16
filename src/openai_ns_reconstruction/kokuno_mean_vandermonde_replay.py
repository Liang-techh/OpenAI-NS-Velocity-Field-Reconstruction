"""Exact clean-room regression for one mean-correction interpolation seam.

This module covers only the MC23--MC26 algebraic interface: translated-bump
moment scaling, the three-node angular interpolation, its nondegeneracy
condition, and the two-node axial solve.  It does not construct the complete
mean correction or a Navier--Stokes witness.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Tuple

Exact = Fraction
Quadratic = Tuple[Exact, Exact, Exact]

SOURCE_COMPONENT = "NS-mean_corrections"
SOURCE_ARCHIVE_PATH = "proof_sources/mean_corrections/mean_corrections_body.tex"
SOURCE_COMPONENT_SHA256 = (
    "c494cebcaf476a01c5aae17d38f84dfc6b0092040c03b9aecbad1f5c3831ea29"
)
SOURCE_ASSOCIATED_CHECKER = "proof_sources/mean_corrections/verify_exact.py"
SOURCE_CHECKER_SHA256 = (
    "5d54f54000e0cf03e70cfde1d298398bd29d091addf07887f70660af7c98a366"
)
SOURCE_RESULT_SHA256 = (
    "5b81365d43a7108b16d11eba8cdf1dc267266226a0d90ff221ff67827f37495f"
)


def _exact(value: Exact | int) -> Exact:
    if isinstance(value, bool):
        raise TypeError("boolean inputs are not accepted as exact scalars")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    raise TypeError("expected fractions.Fraction or int; floating inputs fail closed")


def translated_moment(mu_p: Exact | int, t: Exact | int, j: int) -> Exact:
    """MC23 in exact algebraic form, with ``t`` standing for exp(d*p).

    The public source formula is integral x^p eta_j(x) dx = mu_p * exp(j*d*p).
    Supplying ``t = exp(d*p)`` isolates the exact finite-dimensional identity
    without pretending that transcendental exponentials are rational data.
    """

    if isinstance(j, bool) or not isinstance(j, int) or j < 0:
        raise ValueError("translation index j must be a nonnegative integer")
    return _exact(mu_p) * _exact(t) ** j


def _mul_linear(left: Tuple[Exact, ...], root: Exact) -> Tuple[Exact, ...]:
    out = [Fraction(0, 1)] * (len(left) + 1)
    for degree, coefficient in enumerate(left):
        out[degree] -= root * coefficient
        out[degree + 1] += coefficient
    return tuple(out)


def _scaled_quadratic_from_roots(root_a: Exact, root_b: Exact, scale: Exact) -> Quadratic:
    coeffs = _mul_linear((Fraction(1, 1),), root_a)
    coeffs = _mul_linear(coeffs, root_b)
    return tuple(scale * coefficient for coefficient in coeffs)  # type: ignore[return-value]


def angular_interpolant(
    t0: Exact | int,
    t1: Exact | int,
    t2: Exact | int,
    b1: Exact | int,
    b2: Exact | int,
) -> Quadratic:
    """Return coefficients of the unique degree <=2 MC25 interpolant.

    It enforces C(t0)=0, C(t1)=b1, C(t2)=b2 and fails closed when the
    three nodes are not pairwise distinct.
    """

    x0, x1, x2 = (_exact(t0), _exact(t1), _exact(t2))
    y1, y2 = _exact(b1), _exact(b2)
    if len({x0, x1, x2}) != 3:
        raise ValueError("MC25 requires three pairwise-distinct nodes")

    den1 = (x1 - x0) * (x1 - x2)
    den2 = (x2 - x0) * (x2 - x1)
    term1 = _scaled_quadratic_from_roots(x0, x2, y1 / den1)
    term2 = _scaled_quadratic_from_roots(x0, x1, y2 / den2)
    return tuple(a + b for a, b in zip(term1, term2))  # type: ignore[return-value]


def evaluate_quadratic(coefficients: Iterable[Exact | int], s: Exact | int) -> Exact:
    coeffs = tuple(_exact(value) for value in coefficients)
    if len(coeffs) != 3:
        raise ValueError("expected exactly three quadratic coefficients")
    x = _exact(s)
    return coeffs[0] + coeffs[1] * x + coeffs[2] * x * x


def angular_targets(
    *,
    pressure_target: Exact | int,
    axial_target: Exact | int,
    a: Exact | int,
    mu_p1: Exact | int,
    mu_p2: Exact | int,
) -> Tuple[Exact, Exact]:
    """Return the two nonzero MC24/MC25 interpolation values b1,b2."""

    aa = _exact(a)
    m1 = _exact(mu_p1)
    m2 = _exact(mu_p2)
    if aa == 0 or m1 == 0 or m2 == 0:
        raise ValueError("a and the two moment factors must be nonzero")
    return (
        -_exact(pressure_target) / (2 * aa * m1),
        _exact(axial_target) / (aa * m2),
    )


def angular_source_determinant(
    *,
    a: Exact | int,
    mu_p0: Exact | int,
    mu_p1: Exact | int,
    mu_p2: Exact | int,
    t0: Exact | int,
    t1: Exact | int,
    t2: Exact | int,
) -> Exact:
    """The published determinant factor for the three angular source rows."""

    aa = _exact(a)
    m0, m1, m2 = _exact(mu_p0), _exact(mu_p1), _exact(mu_p2)
    x0, x1, x2 = _exact(t0), _exact(t1), _exact(t2)
    return -2 * aa * aa * m0 * m1 * m2 * (x1 - x0) * (x2 - x0) * (x2 - x1)


def axial_coefficients(
    *,
    t_b: Exact | int,
    t_a: Exact | int,
    j_theta: Exact | int,
    a: Exact | int,
    mu_weighted: Exact | int,
) -> Tuple[Exact, Exact]:
    """MC26 exact two-node axial solve, returning (s0,s1)."""

    tb, ta = _exact(t_b), _exact(t_a)
    aa, mu = _exact(a), _exact(mu_weighted)
    if ta == tb:
        raise ValueError("MC26 degenerates when t_a == t_b (including lambda=0)")
    if aa == 0 or mu == 0:
        raise ValueError("a and mu_weighted must be nonzero")
    s1 = -_exact(j_theta) / (aa * mu * (ta - tb))
    s0 = -tb * s1
    return s0, s1


def axial_residuals(
    *,
    s0: Exact | int,
    s1: Exact | int,
    t_b: Exact | int,
    t_a: Exact | int,
    j_theta: Exact | int,
    a: Exact | int,
    mu_zero: Exact | int,
    mu_weighted: Exact | int,
) -> Tuple[Exact, Exact]:
    """Return residuals of the zero-moment and target equations in MC26."""

    c0, c1 = _exact(s0), _exact(s1)
    tb, ta = _exact(t_b), _exact(t_a)
    aa = _exact(a)
    zero_residual = _exact(mu_zero) * (c0 + tb * c1)
    target_residual = aa * _exact(mu_weighted) * (c0 + ta * c1) + _exact(j_theta)
    return zero_residual, target_residual
