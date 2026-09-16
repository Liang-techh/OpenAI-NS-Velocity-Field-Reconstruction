"""Clean-room exact checker for the corrected oscillatory radial transport product.

This module reimplements the finite Leibniz bookkeeping recorded in the public
corrected reader for the NS-oscillations transport defect.  It intentionally
does not import or copy the source workbench checker.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb
from typing import Mapping, Sequence, Tuple

MultiIndex = Tuple[int, ...]

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_REPO_COMMIT = "fab69fdc4ac197159b8e6ae8d73a82bde2b20d55"
RESEARCH_STATE_SOURCE_COMMIT = "e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f"
SOURCE_COMPONENT = "proof_sources/oscillations/oscillations_body.tex"
SOURCE_COMPONENT_SHA256 = (
    "3fd5c61de39b6f65a581ead5b29c741c2e0f46fb30f62e582dc24d93bbe83615"
)
SOURCE_PAGES = "62–87,157–165"
ASSOCIATED_CHECKER = "proof_sources/oscillations/exact_checks.py"
BUNDLE_SHA256 = "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"
DOI = "10.5281/zenodo.22678406"

PAPER_EXACT_VELOCITY_AVAILABLE = False
FULL_RECONSTRUCTION = False


def _fraction(value: Fraction | int, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, int)):
        raise TypeError(f"{name} must be an exact int or Fraction")
    return Fraction(value)


def _multiindex(index: Sequence[int], name: str = "index") -> MultiIndex:
    result = tuple(index)
    if not result:
        raise ValueError(f"{name} must be nonempty")
    if any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in result):
        raise TypeError(f"{name} entries must be nonnegative integers")
    return result


def subindices(index: Sequence[int]) -> tuple[MultiIndex, ...]:
    """Enumerate J <= I componentwise for one finite multiindex I."""
    idx = _multiindex(index)
    return tuple(tuple(j) for j in product(*(range(i + 1) for i in idx)))


def subtract_index(left: Sequence[int], right: Sequence[int]) -> MultiIndex:
    left_i = _multiindex(left, "left")
    right_i = _multiindex(right, "right")
    if len(left_i) != len(right_i):
        raise ValueError("multiindices must have the same dimension")
    if any(r > l for l, r in zip(left_i, right_i)):
        raise ValueError("right must be <= left componentwise")
    return tuple(l - r for l, r in zip(left_i, right_i))


def multi_binomial(index: Sequence[int], subindex: Sequence[int]) -> int:
    """Return prod_k binom(I_k, J_k)."""
    idx = _multiindex(index)
    sub = _multiindex(subindex, "subindex")
    if len(idx) != len(sub):
        raise ValueError("multiindices must have the same dimension")
    if any(j > i for i, j in zip(idx, sub)):
        raise ValueError("subindex must be <= index componentwise")
    result = 1
    for i, j in zip(idx, sub):
        result *= comb(i, j)
    return result


def exact_product_derivative(
    index: Sequence[int],
    b_derivatives: Mapping[MultiIndex, Fraction | int],
    n_derivatives: Mapping[MultiIndex, Fraction | int],
) -> Fraction:
    """Exact finite Leibniz sum for D^I(b n_{Phi,r})."""
    idx = _multiindex(index)
    total = Fraction(0)
    for sub in subindices(idx):
        complement = subtract_index(idx, sub)
        try:
            b_val = _fraction(b_derivatives[sub], f"D^{sub}b")
            n_val = _fraction(n_derivatives[complement], f"D^{complement}n")
        except KeyError as exc:
            raise KeyError(f"missing derivative datum for Leibniz term {sub}") from exc
        total += multi_binomial(idx, sub) * b_val * n_val
    return total


def corrected_polynomial_bound(
    index: Sequence[int],
    epsilon: Fraction | int,
    s_star: Fraction | int,
    b_constants: Mapping[MultiIndex, Fraction | int],
    n_constants: Mapping[MultiIndex, Fraction | int],
    n_degrees: Mapping[MultiIndex, int],
) -> tuple[Fraction, Fraction, int]:
    """Return the exact AS40 sum and its max-degree S_* envelope.

    Inputs encode |D^J b| <= C_b,J epsilon and
    |D^K n_{Phi,r}| <= C_n,K S_*^{d_K}.  For S_* >= 1 the second
    returned quantity is the finite max-degree envelope.
    """
    idx = _multiindex(index)
    eps = _fraction(epsilon, "epsilon")
    s = _fraction(s_star, "s_star")
    if eps < 0:
        raise ValueError("epsilon must be nonnegative")
    if s < 1:
        raise ValueError("s_star must be at least 1")

    exact = Fraction(0)
    coefficient_sum = Fraction(0)
    max_degree = 0
    for sub in subindices(idx):
        complement = subtract_index(idx, sub)
        try:
            cb = _fraction(b_constants[sub], f"C_b,{sub}")
            cn = _fraction(n_constants[complement], f"C_n,{complement}")
            degree = n_degrees[complement]
        except KeyError as exc:
            raise KeyError(f"missing bound datum for Leibniz term {sub}") from exc
        if cb < 0 or cn < 0:
            raise ValueError("bound constants must be nonnegative")
        if isinstance(degree, bool) or not isinstance(degree, int) or degree < 0:
            raise TypeError("n_degrees must be nonnegative integers")
        weight = multi_binomial(idx, sub) * cb * cn
        exact += eps * weight * (s ** degree)
        coefficient_sum += weight
        max_degree = max(max_degree, degree)

    envelope = eps * (s ** max_degree) * coefficient_sum
    return exact, envelope, max_degree


def transport_defect(
    epsilon: Fraction | int,
    v: Fraction | int,
    h_t: Fraction | int,
    g: Fraction | int,
    h_z: Fraction | int,
    b: Fraction | int,
    n_r: Fraction | int,
) -> Fraction:
    """Exact scalar form epsilon*v*(H_T-G H_Z) + b*n_{Phi,r}."""
    eps = _fraction(epsilon, "epsilon")
    return (
        eps * _fraction(v, "v") * (_fraction(h_t, "h_t") - _fraction(g, "g") * _fraction(h_z, "h_z"))
        + _fraction(b, "b") * _fraction(n_r, "n_r")
    )
