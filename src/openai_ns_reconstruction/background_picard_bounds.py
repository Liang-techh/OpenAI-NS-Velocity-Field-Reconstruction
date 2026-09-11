"""Lemma 5.1 / Eq. (5.8) Picard-majorant utilities.

Section 5 writes the positive-order inner problem as Eq. (5.7) and bounds the
``k``-th Picard term on a smaller complex parameter strip by

    ||K^k G f_n||_{S_{rho'}}
      <= (C_n^(k+1) a^(k+1) / (k+1)!)
         * max(1, p_k / Delta)^p_k,

where ``p_k = ceil(k/2)`` and ``Delta = rho-rho'``.  The sparse block form of
``A1`` is what limits a nonzero composition of ``k`` factors to at most
``p_k`` parameter derivatives.

This module records that theorem-side majorant and a conservative tail bound.
It does *not* infer ``C_n``, ``rho`` or the matrices/sources from samples.  In
particular, caller-supplied constants remain formal-structure inputs until the
materialized Section-4/5 hierarchy supplies certified analytic bounds.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral
import sys


_LOG_MAX_FLOAT = math.log(sys.float_info.max)


def _order(k: int, *, name: str = "k") -> int:
    if isinstance(k, bool) or not isinstance(k, Integral) or k < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(k)


def _positive_finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _strip_loss(rho: float, rho_prime: float) -> tuple[float, float, float]:
    rho = _positive_finite(rho, "rho")
    rho_prime = float(rho_prime)
    if not math.isfinite(rho_prime) or not 0.0 < rho_prime < rho:
        raise ValueError("rho_prime must be finite with 0 < rho_prime < rho")
    return rho, rho_prime, rho - rho_prime


def picard_parameter_derivative_count(k: int) -> int:
    """Return ``p_k = ceil(k/2)`` from the sparse-``A1`` proof of Eq. (5.8)."""

    k = _order(k)
    return (k + 1) // 2


def eq_5_8_log_majorant(
    k: int,
    C_n: float,
    radial_extent: float,
    rho: float,
    rho_prime: float,
) -> float:
    """Return the logarithm of the displayed Eq. (5.8) term majorant.

    Working in logarithms prevents the factorial and the Cauchy-loss factor
    from overflowing separately at moderately large Picard order.
    """

    k = _order(k)
    C_n = _positive_finite(C_n, "C_n")
    radial_extent = _positive_finite(radial_extent, "radial_extent")
    _, _, delta = _strip_loss(rho, rho_prime)

    p_k = picard_parameter_derivative_count(k)
    q = C_n * radial_extent
    if not math.isfinite(q):
        raise OverflowError("C_n * radial_extent is nonfinite")

    cauchy = max(1.0, p_k / delta)
    return (
        (k + 1) * math.log(q)
        - math.lgamma(k + 2)
        + p_k * math.log(cauchy)
    )


def eq_5_8_majorant(
    k: int,
    C_n: float,
    radial_extent: float,
    rho: float,
    rho_prime: float,
) -> float:
    """Evaluate the Eq. (5.8) upper-bound expression in binary64.

    ``inf`` is returned rather than wrapping if the finite analytic expression
    lies outside binary64 range.  Use :func:`eq_5_8_log_majorant` for stable
    comparisons in that regime.
    """

    log_bound = eq_5_8_log_majorant(k, C_n, radial_extent, rho, rho_prime)
    if log_bound > _LOG_MAX_FLOAT:
        return math.inf
    value = math.exp(log_bound)
    # Move one representable float upward.  This protects the public numeric
    # value from a final downward rounding of exp; it is not a substitute for
    # interval arithmetic for the preceding transcendental evaluations.
    return math.nextafter(value, math.inf)


def _two_step_ratio_upper_log(
    k: int,
    C_n: float,
    radial_extent: float,
    delta: float,
) -> float:
    """Logarithm of a uniform upper bound for ``B_{k+2}/B_k``.

    Once ``p_k >= max(1, Delta)``, the Cauchy factor is on its power branch.
    With ``p=p_k`` we use

        (1 + 1/p)^p < 3,     p+1 <= (k+3)/2,

    which gives

        B_{k+2}/B_k <= 3 (C_n a)^2 / (2 Delta (k+2)).

    The right-hand side decreases with ``k``, so the same ratio controls every
    later two-step term of the same parity.
    """

    p_k = picard_parameter_derivative_count(k)
    if p_k < max(1.0, delta):
        raise ValueError(
            "tail start is before the certified Cauchy-factor power regime"
        )
    q = C_n * radial_extent
    if not math.isfinite(q):
        raise OverflowError("C_n * radial_extent is nonfinite")
    return math.log(1.5) + 2.0 * math.log(q) - math.log(delta) - math.log(k + 2.0)


@dataclass(frozen=True)
class PicardTailCertificate:
    """A finite-order analytic tail enclosure derived from Eq. (5.8).

    ``start_order=N`` means the retained Picard prefix is ``k=0,...,N-1`` and
    this object bounds ``sum_{k>=N} ||K^k G f_n||`` under the supplied Eq. (5.8)
    hypotheses.
    """

    start_order: int
    delta: float
    derivative_count_at_start: int
    two_step_ratio_upper: float
    first_log_majorant: float
    second_log_majorant: float
    log_tail_upper_bound: float
    tail_upper_bound: float

    def certifies_tolerance(self, tolerance: float) -> bool:
        tolerance = _positive_finite(tolerance, "tolerance")
        return self.log_tail_upper_bound <= math.log(tolerance)


def picard_tail_certificate(
    start_order: int,
    C_n: float,
    radial_extent: float,
    rho: float,
    rho_prime: float,
) -> PicardTailCertificate:
    """Bound the complete Eq. (5.8) tail starting at ``start_order``.

    The proof splits the tail into even/odd subsequences.  If ``r`` is the
    two-step ratio upper bound at the first omitted order, monotonicity gives

        tail <= (B_N + B_{N+1}) / (1-r).

    The function fails closed until the Cauchy factor is on its power branch
    and ``r < 1``; no finite prefix is extrapolated outside that regime.
    """

    start_order = _order(start_order, name="start_order")
    C_n = _positive_finite(C_n, "C_n")
    radial_extent = _positive_finite(radial_extent, "radial_extent")
    rho, rho_prime, delta = _strip_loss(rho, rho_prime)

    log_ratio = _two_step_ratio_upper_log(
        start_order, C_n, radial_extent, delta
    )
    if log_ratio >= 0.0:
        raise ValueError(
            "tail start is too early: the certified two-step ratio is not < 1"
        )
    ratio = math.exp(log_ratio)

    first = eq_5_8_log_majorant(
        start_order, C_n, radial_extent, rho, rho_prime
    )
    second = eq_5_8_log_majorant(
        start_order + 1, C_n, radial_extent, rho, rho_prime
    )
    high, low = (first, second) if first >= second else (second, first)
    log_numerator = high + math.log1p(math.exp(low - high))
    log_tail = log_numerator - math.log1p(-ratio)

    if log_tail > _LOG_MAX_FLOAT:
        tail = math.inf
    else:
        tail = math.nextafter(math.exp(log_tail), math.inf)

    return PicardTailCertificate(
        start_order=start_order,
        delta=delta,
        derivative_count_at_start=picard_parameter_derivative_count(start_order),
        two_step_ratio_upper=ratio,
        first_log_majorant=first,
        second_log_majorant=second,
        log_tail_upper_bound=log_tail,
        tail_upper_bound=tail,
    )


def find_picard_truncation(
    tolerance: float,
    C_n: float,
    radial_extent: float,
    rho: float,
    rho_prime: float,
    *,
    max_start_order: int = 100_000,
) -> PicardTailCertificate:
    """Find the first certified Picard-tail start meeting ``tolerance``.

    This chooses an order from the manuscript's analytic bound only.  It does
    not fit a truncation order to sampled residuals or to a numerical Picard
    trajectory.
    """

    tolerance = _positive_finite(tolerance, "tolerance")
    C_n = _positive_finite(C_n, "C_n")
    radial_extent = _positive_finite(radial_extent, "radial_extent")
    _, _, delta = _strip_loss(rho, rho_prime)
    max_start_order = _order(max_start_order, name="max_start_order")

    # First enter the p_k >= max(1,Delta) regime used by the two-step proof.
    p_required = max(1, math.ceil(delta))
    start = max(0, 2 * p_required - 1)

    # Also jump directly to the point where the conservative ratio can be < 1:
    # 3(Ca)^2 / (2 Delta (k+2)) < 1.
    q = C_n * radial_extent
    if not math.isfinite(q):
        raise OverflowError("C_n * radial_extent is nonfinite")
    ratio_threshold = 1.5 * q * q / delta
    if not math.isfinite(ratio_threshold):
        raise OverflowError("Eq. (5.8) two-step ratio threshold is nonfinite")
    start = max(start, math.floor(ratio_threshold - 2.0) + 1)
    start = max(start, 0)

    log_tolerance = math.log(tolerance)
    for candidate in range(start, max_start_order + 1):
        try:
            certificate = picard_tail_certificate(
                candidate, C_n, radial_extent, rho, rho_prime
            )
        except ValueError:
            continue
        if certificate.log_tail_upper_bound <= log_tolerance:
            return certificate

    raise ValueError(
        "no Eq. (5.8) Picard truncation certificate met the tolerance "
        f"through start_order={max_start_order}"
    )
