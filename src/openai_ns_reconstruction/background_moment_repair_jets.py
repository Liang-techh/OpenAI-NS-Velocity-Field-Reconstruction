"""Eta-jet lift of the Lemma 5.2 compact moment repair.

The landed :mod:`background_moment_repair` module solves Eqs. (5.14)-(5.16)
pointwise in ``eta``.  For recursive use at the next Section-5 order one also
needs parameter derivatives of the repaired coefficient functions.  Because
the two moment matrices are fixed in ``eta``, those derivatives can be
propagated algebraically from jets of the five unrepaired moments and the
patch factor ``p(eta)=e_* f(eta)``.

This module performs that finite-jet propagation.  Its inputs remain caller-
supplied until the actual positive-order hierarchy is materialized, so the
result is solver infrastructure / formal structure rather than a paper-exact
background coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import comb
import math

import numpy as np

from .background_moment_repair import Lemma52MomentRepair


@dataclass(frozen=True)
class Lemma52RepairJetCoefficients:
    """Ordinary eta derivatives of the five-bump repair at one parameter value."""

    alpha_derivatives: np.ndarray
    beta_derivatives: np.ndarray
    corrected_moment_derivatives: np.ndarray

    def __post_init__(self) -> None:
        alpha = np.asarray(self.alpha_derivatives, dtype=float)
        beta = np.asarray(self.beta_derivatives, dtype=float)
        corrected = np.asarray(self.corrected_moment_derivatives, dtype=float)
        if alpha.ndim != 2 or alpha.shape[0] < 1 or alpha.shape[1] != 2:
            raise ValueError("alpha_derivatives must have shape (J+1, 2)")
        if beta.shape != (alpha.shape[0], 3):
            raise ValueError("beta_derivatives must have shape (J+1, 3)")
        if corrected.shape != (alpha.shape[0], 5):
            raise ValueError("corrected_moment_derivatives must have shape (J+1, 5)")
        if not all(np.all(np.isfinite(a)) for a in (alpha, beta, corrected)):
            raise ValueError("repair jet arrays must be finite")
        alpha = alpha.copy()
        beta = beta.copy()
        corrected = corrected.copy()
        alpha.setflags(write=False)
        beta.setflags(write=False)
        corrected.setflags(write=False)
        object.__setattr__(self, "alpha_derivatives", alpha)
        object.__setattr__(self, "beta_derivatives", beta)
        object.__setattr__(self, "corrected_moment_derivatives", corrected)

    @property
    def max_order(self) -> int:
        return int(self.alpha_derivatives.shape[0] - 1)


def _moment_derivative_array(value: object) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.ndim != 2 or out.shape[0] < 1 or out.shape[1] != 5:
        raise ValueError("base_moment_derivatives must have shape (J+1, 5)")
    if not np.all(np.isfinite(out)):
        raise ValueError("base_moment_derivatives must be finite")
    return out


def _patch_derivative_array(value: object, count: int) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (count,) or not np.all(np.isfinite(out)):
        raise ValueError("patch_factor_derivatives must have shape (J+1,) and be finite")
    if out[0] == 0.0:
        raise ValueError("the center patch factor e_* f(eta) must be nonzero")
    return out


def _quotient_derivatives(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    """Raw derivatives of ``numerator / denominator`` from Leibniz recursion."""

    count = numerator.shape[0]
    result = np.empty(count, dtype=float)
    p0 = float(denominator[0])
    for order in range(count):
        lower = math.fsum(
            float(comb(order, j)) * float(denominator[j]) * float(result[order - j])
            for j in range(1, order + 1)
        )
        value = (float(numerator[order]) - lower) / p0
        if not math.isfinite(value):
            raise OverflowError("Lemma 5.2 quotient jet is outside floating-point range")
        result[order] = value
    return result


def _product_derivatives(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Raw derivatives of a product using the exact Leibniz coefficients."""

    if left.shape != right.shape:
        raise ValueError("product derivative arrays must have matching shapes")
    count = left.shape[0]
    result = np.empty(count, dtype=float)
    for order in range(count):
        value = math.fsum(
            float(comb(order, j)) * float(left[j]) * float(right[order - j])
            for j in range(order + 1)
        )
        if not math.isfinite(value):
            raise OverflowError("Lemma 5.2 product jet is outside floating-point range")
        result[order] = value
    return result


def solve_lemma52_eta_jets(
    repair: Lemma52MomentRepair,
    base_moment_derivatives: object,
    patch_factor_derivatives: object,
) -> Lemma52RepairJetCoefficients:
    """Propagate a finite ``eta`` jet through Eqs. (5.14)-(5.16).

    Arrays contain ordinary derivatives, not factorial-normalized Taylor
    coefficients.  If ``m_i^(k)`` and ``p^(k)`` are supplied for
    ``k=0,...,J``, the routine differentiates the quotient data in Eq. (5.16)
    via

    ``(m/p)^(n) = [m^(n)-sum_{j=1}^n C(n,j)p^(j)(m/p)^(n-j)]/p``

    and then solves the two *constant* moment matrices at every derivative
    order.  The returned five corrected-moment derivatives are recomputed by
    the Leibniz rule and checked to vanish to a condition-number-scaled
    floating-point tolerance.

    This is only an implication from supplied jets.  It does not certify that
    ``p`` stays nonzero on an eta interval, that the input jets come from the
    manuscript hierarchy, or that the resulting functions satisfy uniform
    analytic bounds.
    """

    if not isinstance(repair, Lemma52MomentRepair):
        raise TypeError("repair must be a Lemma52MomentRepair")
    moments = _moment_derivative_array(base_moment_derivatives)
    patch = _patch_derivative_array(patch_factor_derivatives, moments.shape[0])

    ratio_m3 = _quotient_derivatives(moments[:, 3], patch)
    ratio_m2 = _quotient_derivatives(moments[:, 2], patch)
    ratio_m4 = _quotient_derivatives(moments[:, 4], patch)

    count = moments.shape[0]
    alpha = np.empty((count, 2), dtype=float)
    beta = np.empty((count, 3), dtype=float)
    for order in range(count):
        d_u = np.array([moments[order, 0], ratio_m3[order]], dtype=float)
        d_e = np.array(
            [moments[order, 1], 0.5 * ratio_m2[order], -ratio_m4[order]],
            dtype=float,
        )
        alpha[order] = -np.linalg.solve(repair.u_matrix, d_u)
        beta[order] = -np.linalg.solve(repair.e_matrix, d_e)
    if not np.all(np.isfinite(alpha)) or not np.all(np.isfinite(beta)):
        raise OverflowError("Lemma 5.2 repair coefficient jet is nonfinite")

    # Rows of the constant moment matrices applied to alpha/beta jets.
    u_first = alpha @ repair.u_matrix[0]
    u_weighted = alpha @ repair.u_matrix[1]
    e_first = beta @ repair.e_matrix[0]
    e_weighted = beta @ repair.e_matrix[1]
    e_last = beta @ repair.e_matrix[2]

    patch_u = _product_derivatives(patch, u_weighted)
    patch_e_weighted = _product_derivatives(patch, e_weighted)
    patch_e_last = _product_derivatives(patch, e_last)

    corrected = moments.copy()
    corrected[:, 0] += u_first
    corrected[:, 1] += e_first
    corrected[:, 2] += 2.0 * patch_e_weighted
    corrected[:, 3] += patch_u
    corrected[:, 4] -= patch_e_last
    if not np.all(np.isfinite(corrected)):
        raise OverflowError("corrected moment derivative jet is nonfinite")

    # This tolerance checks only the floating-point algebraic solve.  It is not
    # used to turn a nonzero mathematical recurrence defect into an exact zero.
    condition = max(
        1.0,
        float(np.linalg.cond(repair.u_matrix)),
        float(np.linalg.cond(repair.e_matrix)),
    )
    eps = np.finfo(float).eps
    for order in range(count):
        correction_scale = max(
            abs(float(u_first[order])),
            abs(float(e_first[order])),
            abs(float(2.0 * patch_e_weighted[order])),
            abs(float(patch_u[order])),
            abs(float(patch_e_last[order])),
        )
        scale = max(1.0, float(np.max(np.abs(moments[order]))), correction_scale)
        tolerance = 4096.0 * eps * condition * scale
        if float(np.max(np.abs(corrected[order]))) > tolerance:
            raise RuntimeError(
                f"differentiated Eq. (5.16) failed to cancel moments at order {order}"
            )

    return Lemma52RepairJetCoefficients(alpha, beta, corrected)
