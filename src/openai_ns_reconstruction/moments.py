"""Pointwise finite-dimensional primitives for Appendix A, Lemmas A.1-A.2.

The caller must still construct the paper's moment data B(eta), Q_eta and d(eta).
These floating-point checks do NOT establish uniform-in-eta invertibility,
parameter smoothness, derivative bounds, or the complete leading profile.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral
import numpy as np
from .quadrature import unit_rule


@dataclass(frozen=True)
class MomentSolution:
    coefficients: np.ndarray
    iterations: int
    residual_norm: float
    inverse_norm_estimate: float
    bilinear_norm_upper_bound: float
    smallness_estimate: float
    ball_radius_estimate: float
    contraction_estimate: float


def _norm(value: np.ndarray) -> float:
    """Scaled Euclidean norm without spurious square underflow/overflow."""
    scale = float(np.max(np.abs(value)))
    return 0.0 if scale == 0 else scale*float(np.linalg.norm(value/scale))


def _array(value: object, shape: tuple[int, ...], name: str) -> np.ndarray:
    value = np.asarray(value, dtype=float)
    if value.shape != shape or not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must be a finite array of shape {shape}")
    return value


def solve_quadratic_moments(B: np.ndarray, Q: np.ndarray, d: np.ndarray, *,
                            rtol: float = 1e-12, atol: float = 0.0,
                            max_iter: int = 200) -> MomentSolution:
    """Solve B c + Q(c,c) = d by the zero-initialized iteration of Lemma A.2.

    Q[i,j,k] represents output i, inputs j,k (symmetry is not required).
    Use the Euclidean matrix inverse norm and a conservative Frobenius bound
    on the bilinear norm. Reject data outside 8*beta^2*kappa*||d|| <= 1;
    this is a sufficient floating-point screen, NOT an interval certificate.
    For Q=0 there is no mathematical smallness requirement.
    """
    d = np.asarray(d, dtype=float)
    if d.ndim != 1 or not len(d):
        raise ValueError("d must be a nonempty vector")
    m = len(d)
    d, B, Q = _array(d, (m,), "d"), _array(B, (m,m), "B"), _array(Q, (m,m,m), "Q")
    if not math.isfinite(rtol) or not 4*np.finfo(float).eps <= rtol < 1:
        raise ValueError("rtol must be in [4*machine_epsilon, 1)")
    if not math.isfinite(atol) or atol < 0:
        raise ValueError("atol must be finite and nonnegative")
    if isinstance(max_iter,bool) or not isinstance(max_iter,Integral) or max_iter < 1:
        raise ValueError("max_iter must be a positive integer")
    singular_values = np.linalg.svd(B, compute_uv=False)
    if singular_values[-1] <= np.finfo(float).eps*m*singular_values[0]:
        raise ValueError("B is singular or too ill-conditioned for this float64 solver")
    beta = float(1/singular_values[-1])
    kappa = _norm(Q.ravel())
    dnorm = _norm(d)
    radius = 2*beta*dnorm
    smallness = 0.0 if kappa == 0 else 8*beta**2*kappa*dnorm
    contraction = 2*beta*kappa*radius
    if not all(math.isfinite(v) for v in (beta,kappa,dnorm,radius,smallness,contraction)):
        raise OverflowError("rescale the moment system: norm bounds overflowed")
    if smallness > 1:
        raise ValueError("Lemma A.2 sufficient smallness test failed (estimate > 1)")
    c = np.zeros(m)
    for iteration in range(1, int(max_iter)+1):
        quadratic = np.einsum('ijk,j,k->i', Q, c, c)
        new = np.linalg.solve(B, d-quadratic)
        if not np.all(np.isfinite(new)):
            raise OverflowError("moment iteration produced a nonfinite coefficient")
        quadratic_new = np.einsum('ijk,j,k->i', Q, new, new)
        residual = _norm(B@new+quadratic_new-d)
        difference = _norm(new-c)
        cnorm = _norm(new)
        scale = max(dnorm, float(singular_values[0])*cnorm, kappa*cnorm**2)
        if (residual <= atol+rtol*scale and
                contraction*difference/(1-contraction) <= atol+rtol*cnorm):
            new.setflags(write=False)
            return MomentSolution(new,iteration,residual,beta,kappa,smallness,radius,contraction)
        c = new
    raise RuntimeError("moment iteration did not converge; increase max_iter or rescale")


def power_moment_matrix(exponents: np.ndarray, intervals: np.ndarray, *, n: int = 64) -> np.ndarray:
    """Numerically integrate x**alpha_i against normalized disjoint C-infinity bumps.

    Positive, ordered disjoint intervals and distinct exponents implement the
    finite moment-matrix layout of Lemma A.1. Quadrature is not a proof of
    nonzero determinant; inspect singular values before solving a given matrix.
    """
    exponents = np.asarray(exponents, dtype=float)
    if (exponents.ndim != 1 or not len(exponents) or not np.all(np.isfinite(exponents))
            or len(np.unique(exponents)) != len(exponents)):
        raise ValueError("exponents must be a finite nonempty vector of distinct values")
    m = len(exponents)
    intervals = _array(intervals, (m,2), "intervals")
    if (np.any(intervals[:,0] <= 0) or np.any(intervals[:,1] <= intervals[:,0])
            or np.any(intervals[1:,0] <= intervals[:-1,1])):
        raise ValueError("intervals must be positive, ordered, and disjoint")
    x,w = unit_rule(n)
    bump = np.exp(4-1/(x*(1-x)))
    weights = w*bump
    weights = weights/np.sum(weights)
    physical = intervals[:,0,None]+(intervals[:,1]-intervals[:,0])[:,None]*x
    with np.errstate(over='raise', invalid='raise'):
        try:
            matrix = np.sum(physical[None,:,:]**exponents[:,None,None]*weights, axis=2)
        except FloatingPointError as exc:
            raise OverflowError("power moments overflowed; rescale intervals/exponents") from exc
    return matrix
