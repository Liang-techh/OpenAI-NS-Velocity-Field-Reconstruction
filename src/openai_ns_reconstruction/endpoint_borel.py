"""Executable Taylor--Borel right-extension infrastructure for the t=1 endpoint.

This mirrors the algebraic/local-finiteness part of the pinned Lean construction:
``BorelExtension.term``, ``DiagonalScale.doublingEnvelope`` and
``SpatialBorelExtension.rightExtension``.  It does *not* certify the analytic
bounds used by Lean to prove all-order smoothness.  Until those bounds are
materialized from the actual residual jets, this module is formal-structure
infrastructure rather than a paper-exact force constructor.
"""
from __future__ import annotations

from collections.abc import Callable
import math
import operator

import numpy as np

from .coordinates import _finite
from .cutoffs import standard_cutoff
from .verify import finite_vector

SpatialJetFamily = Callable[[int, float, float, float], np.ndarray]
LocalScale = Callable[[int], int]
VectorField = Callable[[float, float, float, float], np.ndarray]


def _natural(value: object, name: str) -> int:
    """Return an exact nonnegative integer, rejecting lossy float coercions."""
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a nonnegative integer")
    try:
        out = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be a nonnegative integer") from exc
    if out < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(out)


class DoublingEnvelope:
    """Executable form of ``DiagonalScale.doublingEnvelope``.

    For local natural-number scales ``b_j`` this returns
    ``a_0=max(1,b_0)`` and ``a_{j+1}=max(b_{j+1},2 a_j)``.  Consequently every
    returned scale is positive and every step at least doubles.  The latter is
    the local-finiteness mechanism for the Borel sum away from the endpoint.
    """

    def __init__(self, local_scale: LocalScale):
        if not callable(local_scale):
            raise ValueError("local_scale must be callable")
        self._local_scale = local_scale
        self._cache: list[int] = []

    def __call__(self, degree: int) -> int:
        degree = _natural(degree, "degree")
        while len(self._cache) <= degree:
            j = len(self._cache)
            b = _natural(self._local_scale(j), f"local_scale({j})")
            scale = max(1, b) if j == 0 else max(b, 2 * self._cache[-1])
            self._cache.append(scale)
        return self._cache[degree]


def borel_monomial_coefficient(s: float, degree: int) -> float:
    """Return ``s**degree / degree!`` without constructing a huge factorial."""
    s = _finite(s, "s")
    degree = _natural(degree, "degree")
    value = 1.0
    for j in range(1, degree + 1):
        value *= s / j
    return value


def _cutoff_argument(scale: int, s: float) -> float:
    """Compute scale*s, treating integer-to-float overflow as outside support."""
    try:
        return float(scale) * s
    except OverflowError:
        return math.inf


class BorelRightExtension:
    """Pointwise evaluator for the pinned Taylor--Borel future-branch series.

    ``jet(j,x,y,z)`` supplies the j-th normal endpoint jet and ``local_scale``
    supplies the pre-envelope integer scales.  For ``t>T`` the doubling
    envelope makes the nominal infinite sum pointwise finite: once
    ``a_j*(t-T) >= 1``, that and every later cutoff is zero.

    This class intentionally does not expose a ``paper_exact`` flag.  A
    caller-supplied jet family and local scales are not evidence of the
    template derivative bounds or of locally uniform limits of the actual
    Navier--Stokes residual derivatives required by ``CandidateFromLimits``.
    """

    def __init__(self, jet: SpatialJetFamily, local_scale: LocalScale, *, endpoint: float = 1.0):
        if not callable(jet):
            raise ValueError("jet must be callable")
        self.jet = jet
        self.endpoint = _finite(endpoint, "endpoint")
        self.scales = DoublingEnvelope(local_scale)

    def scale(self, degree: int) -> int:
        return self.scales(degree)

    def term(self, degree: int, x: float, y: float, z: float, t: float) -> np.ndarray:
        """Evaluate one cutoff monomial term on the future side ``t>=T``."""
        degree = _natural(degree, "degree")
        t = _finite(t, "t")
        s = t - self.endpoint
        if s < 0.0:
            raise ValueError("Borel right-extension terms are defined for t>=endpoint")
        argument = _cutoff_argument(self.scale(degree), s)
        cutoff = standard_cutoff(argument)
        if cutoff == 0.0:
            return np.zeros(3)
        coefficient = borel_monomial_coefficient(s, degree)
        if coefficient == 0.0:
            return np.zeros(3)
        return cutoff * coefficient * finite_vector(self.jet(degree, x, y, z))

    def __call__(self, x: float, y: float, z: float, t: float) -> np.ndarray:
        """Evaluate the locally finite future series, including its value at T."""
        x, y, z, t = (_finite(v, name) for v, name in zip(
            (x, y, z, t), ("x", "y", "z", "t")
        ))
        s = t - self.endpoint
        if s < 0.0:
            raise ValueError("Borel right extension may only be evaluated for t>=endpoint")
        if s == 0.0:
            return finite_vector(self.jet(0, x, y, z))

        out = np.zeros(3)
        coefficient = 1.0
        degree = 0
        # A finite IEEE-754 s>0 needs fewer than ~1100 doublings from scale>=1.
        # 4096 is a fail-closed guard against a broken scale implementation.
        while degree < 4096:
            scale = self.scale(degree)
            argument = _cutoff_argument(scale, s)
            if argument >= 1.0:
                return finite_vector(out)
            cutoff = standard_cutoff(argument)
            if cutoff != 0.0 and coefficient != 0.0:
                out += cutoff * coefficient * finite_vector(self.jet(degree, x, y, z))
            degree += 1
            coefficient *= s / degree
        raise RuntimeError("doubling-envelope local-finiteness guard was exceeded")

    def glue_to_past(self, past: VectorField) -> VectorField:
        """Return the algebraic t=T glue, without asserting C-infinity matching."""
        if not callable(past):
            raise ValueError("past must be callable")

        def glued(x: float, y: float, z: float, t: float) -> np.ndarray:
            t = _finite(t, "t")
            if t <= self.endpoint:
                return finite_vector(past(x, y, z, t))
            return self(x, y, z, t)

        return glued
