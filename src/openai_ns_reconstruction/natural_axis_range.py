"""Constructive adapters from the pinned Lean ``NaturalAxisRange`` layer.

This module implements only choices that are explicit in the official
formalization at the repository-pinned commit.  In particular:

* ``0 < h <= 1/100`` and ``0 < j <= 1/20`` are the full printed axis range;
* the unique root of ``H(h,j,eta)`` lies in ``(-j/4,-j/5)``;
* the low-|Z| cutoff width is chosen explicitly as ``delta = j/10``;
* once a certified positive lower margin ``m <= H^2`` on that low-|Z| set is
  available, the Lean proof chooses ``sigma = sqrt(m)/20``.

The missing step is important: this file does NOT manufacture the margin ``m``
or the pressure datum needed to prove it.  Therefore the returned objects are
constructive witnesses conditional on caller-supplied certified margin data,
not a completed Theorem 4.6 profile.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .natural_axis import H, chi


@dataclass(frozen=True)
class NaturalAxisRangeParameters:
    """The complete range proved in ``NaturalAxisRange.Parameters``."""

    h: float
    j: float

    def __post_init__(self) -> None:
        h, j = float(self.h), float(self.j)
        if not math.isfinite(h) or not 0.0 < h <= 1.0 / 100.0:
            raise ValueError("h must satisfy 0 < h <= 1/100")
        if not math.isfinite(j) or not 0.0 < j <= 1.0 / 20.0:
            raise ValueError("j must satisfy 0 < j <= 1/20")
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "j", j)

    @property
    def root_bracket(self) -> tuple[float, float]:
        """Open bracket from ``exists_unique_root``: ``(-j/4,-j/5)``."""

        return (-self.j / 4.0, -self.j / 5.0)

    @property
    def cutoff_delta(self) -> float:
        """The explicit ``delta = j/10`` used by ``exists_cutoff_parameters``."""

        return self.j / 10.0


@dataclass(frozen=True)
class NaturalAxisRootWitness:
    """Numerical realization of the formally unique root bracket.

    The bracket/uniqueness theorem is supplied by Lean; the floating-point
    ``eta0`` here is only a numerical location of that root.
    """

    eta0: float
    residual: float
    bracket: tuple[float, float]


def locate_unique_H_root(
    parameters: NaturalAxisRangeParameters,
    *,
    atol: float = 1e-14,
    max_iter: int = 200,
) -> NaturalAxisRootWitness:
    """Locate the unique root of ``H`` by bisection inside the Lean bracket."""

    atol = float(atol)
    if not math.isfinite(atol) or atol <= 0.0:
        raise ValueError("atol must be finite and positive")
    if isinstance(max_iter, bool) or not isinstance(max_iter, int) or max_iter < 1:
        raise ValueError("max_iter must be a positive integer")

    lo, hi = parameters.root_bracket
    flo = H(parameters.h, parameters.j, lo)
    fhi = H(parameters.h, parameters.j, hi)
    if not (flo < 0.0 < fhi):
        raise ArithmeticError("formal H-root bracket lost its sign change numerically")

    for _ in range(max_iter):
        mid = lo + (hi - lo) / 2.0
        fm = H(parameters.h, parameters.j, mid)
        if abs(fm) <= atol or hi - lo <= atol:
            return NaturalAxisRootWitness(mid, fm, parameters.root_bracket)
        if fm <= 0.0:
            lo = mid
        else:
            hi = mid
    raise RuntimeError("H-root bisection did not converge")


@dataclass(frozen=True)
class CutoffParametersFromMargin:
    """Explicit delta/sigma choices conditional on a certified H^2 margin."""

    delta: float
    sigma: float
    margin: float
    guaranteed_chi_lower_bound: float


def cutoff_parameters_from_margin(
    parameters: NaturalAxisRangeParameters,
    margin: float,
) -> CutoffParametersFromMargin:
    """Instantiate the explicit Lean choice ``sigma = sqrt(m)/20``.

    If the upstream pressure argument has established ``H(eta)^2 >= m > 0``
    whenever ``|Z(eta)| <= delta``, then with ``delta=j/10`` the cutoff obeys

        chi = H^2 / (H^2 + sigma^2) >= 400/401 > 99/100.

    This function deliberately accepts ``m`` rather than trying to infer it
    from samples; a sampled minimum is not a proof of the required uniform
    margin.
    """

    margin = float(margin)
    if not math.isfinite(margin) or margin <= 0.0:
        raise ValueError("margin must be finite and positive")
    sigma = math.sqrt(margin) / 20.0
    lower = 400.0 / 401.0
    return CutoffParametersFromMargin(
        delta=parameters.cutoff_delta,
        sigma=sigma,
        margin=margin,
        guaranteed_chi_lower_bound=lower,
    )


def verify_chi_margin_point(
    parameters: NaturalAxisRangeParameters,
    witness: CutoffParametersFromMargin,
    eta: float,
) -> float:
    """Check the pointwise algebra once ``H^2 >= m`` is independently known.

    This is a diagnostic consistency check for a supplied margin witness, not a
    replacement for proving the uniform low-|Z| hypothesis.
    """

    eta = float(eta)
    if not math.isfinite(eta) or abs(eta) > 1.0:
        raise ValueError("eta must be finite and lie in [-1,1]")
    h2 = H(parameters.h, parameters.j, eta) ** 2
    if h2 < witness.margin:
        raise ValueError("point does not satisfy the supplied H^2 margin")
    value = chi(parameters.h, parameters.j, witness.sigma, eta)
    if value + 1e-15 < witness.guaranteed_chi_lower_bound:
        raise ArithmeticError("chi lower-bound algebra failed numerically")
    return value
