"""Compact Lemma 5.2 moment repair for Section 5 background coefficients.

For each positive order n, Lemma 5.2 extends the inner solution and adds five
compactly supported radial bumps in the reserved positive-order patch:

    U_n = U_tilde_n + sum_{j=1}^2 alpha_{n,j} b^U_j(R),
    E_n = E_tilde_n + sum_{j=1}^3 beta_{n,j} b^E_j(R),       (5.14)

with R = sqrt(2X).  On that patch the leading profile has
E_0 = e_* f(eta) R^(-1-2 lambda), so the five moment equations split into
the two constant systems in (5.16).

This module materializes that finite repair with actual C-infinity compact
bumps and fixed support geometry.  The base moments m^0_n(eta) and the factor
e_* f(eta) remain explicit inputs until the upstream paper-exact profiles are
available.  Numerical quadrature/inversion is verification infrastructure,
not a formal proof of Lemma A.1 or a paper-exact completed coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from numbers import Integral
from typing import Sequence
import math

import numpy as np

from .quadrature import integrate


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


@dataclass(frozen=True)
class CompactMomentBump:
    """Unit-mass C-infinity bump supported on one positive R interval."""

    left: float
    right: float
    quadrature_points: int = 96

    def __post_init__(self) -> None:
        left = _finite(self.left, "left")
        right = _finite(self.right, "right")
        if not 0.0 < left < right:
            raise ValueError("bump support must satisfy 0 < left < right")
        if (
            isinstance(self.quadrature_points, bool)
            or not isinstance(self.quadrature_points, Integral)
            or not 8 <= int(self.quadrature_points) <= 2048
        ):
            raise ValueError("quadrature_points must be an integer in [8, 2048]")

    def _raw(self, radius: float) -> float:
        if radius <= self.left or radius >= self.right:
            return 0.0
        s = (radius - self.left) / (self.right - self.left)
        return math.exp(4.0 - 1.0 / (s * (1.0 - s)))

    @cached_property
    def normalization(self) -> float:
        value = integrate(
            self._raw,
            self.left,
            self.right,
            n=int(self.quadrature_points),
        )
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError("bump normalization must be finite and positive")
        return value

    def __call__(self, radius: float) -> float:
        radius = _finite(radius, "radius")
        if radius <= self.left or radius >= self.right:
            return 0.0
        return self._raw(radius) / self.normalization

    def moment(self, power: float) -> float:
        """Return integral R**power b(R) dR over this compact support."""
        power = _finite(power, "power")
        return integrate(
            lambda radius: radius**power * self(radius),
            self.left,
            self.right,
            n=int(self.quadrature_points),
        )


@dataclass(frozen=True)
class Lemma52RepairCoefficients:
    """The five coefficient functions evaluated at one eta."""

    alpha: np.ndarray
    beta: np.ndarray
    corrected_moments: np.ndarray

    def __post_init__(self) -> None:
        alpha = np.asarray(self.alpha, dtype=float)
        beta = np.asarray(self.beta, dtype=float)
        moments = np.asarray(self.corrected_moments, dtype=float)
        if alpha.shape != (2,) or not np.all(np.isfinite(alpha)):
            raise ValueError("alpha must be a finite vector of shape (2,)")
        if beta.shape != (3,) or not np.all(np.isfinite(beta)):
            raise ValueError("beta must be a finite vector of shape (3,)")
        if moments.shape != (5,) or not np.all(np.isfinite(moments)):
            raise ValueError("corrected_moments must be a finite vector of shape (5,)")
        alpha = alpha.copy()
        beta = beta.copy()
        moments = moments.copy()
        alpha.setflags(write=False)
        beta.setflags(write=False)
        moments.setflags(write=False)
        object.__setattr__(self, "alpha", alpha)
        object.__setattr__(self, "beta", beta)
        object.__setattr__(self, "corrected_moments", moments)


@dataclass(frozen=True)
class Lemma52MomentRepair:
    """Fixed compact-support repair geometry for Eqs. (5.14)-(5.16)."""

    lambda_exponent: float
    u_bumps: tuple[CompactMomentBump, CompactMomentBump]
    e_bumps: tuple[CompactMomentBump, CompactMomentBump, CompactMomentBump]

    def __post_init__(self) -> None:
        lam = _finite(self.lambda_exponent, "lambda_exponent")
        if lam <= 0.0:
            raise ValueError("lambda_exponent must be positive")
        if len(self.u_bumps) != 2 or len(self.e_bumps) != 3:
            raise ValueError("Lemma 5.2 requires two U bumps and three E bumps")
        bumps = tuple(self.u_bumps) + tuple(self.e_bumps)
        if not all(isinstance(bump, CompactMomentBump) for bump in bumps):
            raise TypeError("all bump entries must be CompactMomentBump instances")
        intervals = sorted((b.left, b.right) for b in bumps)
        if any(intervals[i + 1][0] <= intervals[i][1] for i in range(4)):
            raise ValueError("the five bump supports must be pairwise disjoint")

    @classmethod
    def from_intervals(
        cls,
        lambda_exponent: float,
        u_intervals: Sequence[Sequence[float]],
        e_intervals: Sequence[Sequence[float]],
        *,
        quadrature_points: int = 96,
    ) -> "Lemma52MomentRepair":
        """Construct the five fixed unit-mass bumps from positive R intervals."""
        if len(u_intervals) != 2 or len(e_intervals) != 3:
            raise ValueError("expected two U intervals and three E intervals")
        u_bumps = tuple(
            CompactMomentBump(float(pair[0]), float(pair[1]), quadrature_points)
            for pair in u_intervals
        )
        e_bumps = tuple(
            CompactMomentBump(float(pair[0]), float(pair[1]), quadrature_points)
            for pair in e_intervals
        )
        return cls(lambda_exponent, u_bumps, e_bumps)  # type: ignore[arg-type]

    @cached_property
    def u_matrix(self) -> np.ndarray:
        """The 2x2 B_U matrix in Eq. (5.16)."""
        powers = (1.0, 1.0 - 2.0 * self.lambda_exponent)
        matrix = np.array(
            [[bump.moment(power) for bump in self.u_bumps] for power in powers],
            dtype=float,
        )
        self._check_matrix(matrix, "B_U")
        matrix.setflags(write=False)
        return matrix

    @cached_property
    def e_matrix(self) -> np.ndarray:
        """The 3x3 B_E matrix in Eq. (5.16)."""
        lam = self.lambda_exponent
        powers = (2.0, -2.0 - 2.0 * lam, -2.0 * lam)
        matrix = np.array(
            [[bump.moment(power) for bump in self.e_bumps] for power in powers],
            dtype=float,
        )
        self._check_matrix(matrix, "B_E")
        matrix.setflags(write=False)
        return matrix

    @staticmethod
    def _check_matrix(matrix: np.ndarray, name: str) -> None:
        if not np.all(np.isfinite(matrix)):
            raise ValueError(f"{name} contains nonfinite entries")
        singular_values = np.linalg.svd(matrix, compute_uv=False)
        threshold = np.finfo(float).eps * matrix.shape[0] * singular_values[0]
        if singular_values[-1] <= threshold:
            raise ValueError(
                f"{name} is singular or too ill-conditioned for float64 inversion"
            )

    def solve(
        self,
        base_moments: Sequence[float],
        patch_factor: float,
    ) -> Lemma52RepairCoefficients:
        """Solve Eq. (5.16) from the five unrepaired moments m^0_n(eta).

        ``patch_factor`` is the manuscript quantity ``e_* f(eta)``. It must be
        nonzero; a paper-exact caller must additionally certify the stated
        uniform lower bound for f on |eta| <= 1.
        """
        moments = np.asarray(base_moments, dtype=float)
        if moments.shape != (5,) or not np.all(np.isfinite(moments)):
            raise ValueError("base_moments must be a finite vector of shape (5,)")
        patch_factor = _finite(patch_factor, "patch_factor")
        if patch_factor == 0.0:
            raise ValueError("patch_factor=e_* f(eta) must be nonzero")

        d_u = np.array([moments[0], moments[3] / patch_factor], dtype=float)
        d_e = np.array(
            [
                moments[1],
                moments[2] / (2.0 * patch_factor),
                -moments[4] / patch_factor,
            ],
            dtype=float,
        )
        alpha = -np.linalg.solve(self.u_matrix, d_u)
        beta = -np.linalg.solve(self.e_matrix, d_e)
        if not np.all(np.isfinite(alpha)) or not np.all(np.isfinite(beta)):
            raise OverflowError("Lemma 5.2 moment coefficients are nonfinite")

        corrected = moments.copy()
        corrected[0] += self.u_matrix[0] @ alpha
        corrected[3] += patch_factor * (self.u_matrix[1] @ alpha)
        corrected[1] += self.e_matrix[0] @ beta
        corrected[2] += 2.0 * patch_factor * (self.e_matrix[1] @ beta)
        corrected[4] -= patch_factor * (self.e_matrix[2] @ beta)
        if not np.all(np.isfinite(corrected)):
            raise OverflowError("corrected moment vector is nonfinite")

        scale = max(1.0, float(np.max(np.abs(moments))))
        tolerance = 512.0 * np.finfo(float).eps * max(
            scale,
            float(np.linalg.cond(self.u_matrix)) * scale,
            float(np.linalg.cond(self.e_matrix)) * scale,
        )
        if float(np.max(np.abs(corrected))) > tolerance:
            raise RuntimeError("Eq. (5.16) solve failed to cancel the five moments")

        return Lemma52RepairCoefficients(alpha, beta, corrected)

    def corrected_u_value(
        self,
        X: float,
        base_value: float,
        coefficients: Lemma52RepairCoefficients,
    ) -> float:
        """Evaluate the compactly repaired U_n value from Eq. (5.14)."""
        X = _finite(X, "X")
        if X < 0.0:
            raise ValueError("X must be nonnegative")
        base_value = _finite(base_value, "base_value")
        radius = math.sqrt(2.0 * X)
        value = base_value + sum(
            float(a) * bump(radius)
            for a, bump in zip(coefficients.alpha, self.u_bumps)
        )
        return _finite(value, "corrected U_n")

    def corrected_e_value(
        self,
        X: float,
        base_value: float,
        coefficients: Lemma52RepairCoefficients,
    ) -> float:
        """Evaluate the compactly repaired E_n value from Eq. (5.14)."""
        X = _finite(X, "X")
        if X < 0.0:
            raise ValueError("X must be nonnegative")
        base_value = _finite(base_value, "base_value")
        radius = math.sqrt(2.0 * X)
        value = base_value + sum(
            float(b) * bump(radius)
            for b, bump in zip(coefficients.beta, self.e_bumps)
        )
        return _finite(value, "corrected E_n")
