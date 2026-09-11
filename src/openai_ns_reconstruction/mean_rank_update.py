"""Section 8 compact five-row mean correction, at a fail-closed execution boundary.

This module ports the finite-dimensional/scaling algebra pinned in
``NavierStokes/FiveRowRank.lean`` and ``NavierStokes/MeanRankUpdate.lean``:

* the three debts are ordered ``(P, J_theta, J_z)``;
* the angular powers are ``(2, -2-2*lam, -2*lam)`` and the axial powers are
  ``(1, 1-2*lam)``;
* three angular and two axial correction cells are the exact separated cells
  used by ``FiveRowRank.cellLower/cellUpper``;
* debt normalization/scaling uses the physical ``ell`` and ``U`` powers from
  ``MeanRankUpdate``.

The official Lean construction uses a noncomputable Mathlib ``ContDiffBump``.
For executable evaluation only, this module reuses the repository's explicit
C-infinity :class:`background_moment_repair.CompactMomentBump` on the *same
proved inner support intervals*.  Its interior values are not definitionally
the Mathlib bump and its moment matrices are inverted numerically.  Therefore
this adapter is ``formal-structure`` only.  It must not be used to promote a
wave, mean correction, or velocity field to paper-exact status.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from numbers import Integral
from typing import Sequence
import math

import numpy as np

from .background_moment_repair import CompactMomentBump


STATUS = "formal-structure"
PAPER_EXACT = False


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _count(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("count must be a positive integer")
    return int(value)


def cell_intervals(a: float, b: float, count: int) -> tuple[tuple[float, float], ...]:
    """Return the pinned ``FiveRowRank.cellLower/cellUpper`` intervals."""
    a = _finite(a, "a")
    b = _finite(b, "b")
    count = _count(count)
    if not 0.0 < a < b:
        raise ValueError("cell interval must satisfy 0 < a < b")
    step = (b - a) / (2.0 * count + 1.0)
    return tuple(
        (
            a + (2.0 * j + 1.0) * step,
            a + (2.0 * j + 2.0) * step,
        )
        for j in range(count)
    )


def inner_support(interval: Sequence[float]) -> tuple[float, float]:
    """Support interval of ``LocalizedMomentRepair.bump`` inside one cell."""
    if len(interval) != 2:
        raise ValueError("interval must contain exactly two endpoints")
    left = _finite(interval[0], "left")
    right = _finite(interval[1], "right")
    if not 0.0 < left < right:
        raise ValueError("interval must satisfy 0 < left < right")
    return (3.0 * left + right) / 4.0, (left + 3.0 * right) / 4.0


@dataclass(frozen=True)
class MeanRankUpdatePlan:
    """Parameterized Section 8 five-row correction plan.

    ``debt`` is the physical three-vector ``(P, J_theta, J_z)``.  Upstream
    Sections 6--8 must eventually supply this debt and the actual patch
    parameters.  This class never infers them from samples.
    """

    lam: float
    C: float
    a: float
    b: float
    ell: float
    U: float
    debt: tuple[float, float, float]
    quadrature_points: int = 128

    def __post_init__(self) -> None:
        lam = _finite(self.lam, "lam")
        C = _finite(self.C, "C")
        a = _finite(self.a, "a")
        b = _finite(self.b, "b")
        ell = _finite(self.ell, "ell")
        U = _finite(self.U, "U")
        if lam <= 0.0:
            raise ValueError("lam must be positive")
        if C == 0.0:
            raise ValueError("C must be nonzero")
        if not 0.0 < a < b:
            raise ValueError("repair patch must satisfy 0 < a < b")
        if ell <= 0.0:
            raise ValueError("ell must be positive")
        if U == 0.0:
            raise ValueError("U must be nonzero")
        if len(self.debt) != 3:
            raise ValueError("debt must contain (P, J_theta, J_z)")
        debt = tuple(_finite(v, f"debt[{i}]") for i, v in enumerate(self.debt))
        if (
            isinstance(self.quadrature_points, bool)
            or not isinstance(self.quadrature_points, Integral)
            or not 8 <= int(self.quadrature_points) <= 2048
        ):
            raise ValueError("quadrature_points must be an integer in [8, 2048]")
        object.__setattr__(self, "debt", debt)

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def angular_powers(self) -> tuple[float, float, float]:
        return 2.0, -2.0 - 2.0 * self.lam, -2.0 * self.lam

    @property
    def axial_powers(self) -> tuple[float, float]:
        return 1.0, 1.0 - 2.0 * self.lam

    @property
    def normalized_debt(self) -> np.ndarray:
        d0, d1, d2 = self.debt
        result = np.array(
            [
                d0 / self.U**2,
                d1 / (self.ell**3 * self.U**2),
                d2 / (self.ell**2 * self.U**2),
            ],
            dtype=float,
        )
        result.setflags(write=False)
        return result

    @property
    def angular_targets(self) -> np.ndarray:
        d = self.normalized_debt
        result = np.array([0.0, -d[0] / (2.0 * self.C), d[2] / self.C])
        result.setflags(write=False)
        return result

    @property
    def axial_targets(self) -> np.ndarray:
        d = self.normalized_debt
        result = np.array([0.0, -d[1] / self.C])
        result.setflags(write=False)
        return result

    @cached_property
    def angular_cells(self) -> tuple[tuple[float, float], ...]:
        return cell_intervals(self.a, self.b, 3)

    @cached_property
    def axial_cells(self) -> tuple[tuple[float, float], ...]:
        return cell_intervals(self.a, self.b, 2)

    @cached_property
    def angular_supports(self) -> tuple[tuple[float, float], ...]:
        return tuple(inner_support(cell) for cell in self.angular_cells)

    @cached_property
    def axial_supports(self) -> tuple[tuple[float, float], ...]:
        return tuple(inner_support(cell) for cell in self.axial_cells)

    @cached_property
    def angular_bumps(self) -> tuple[CompactMomentBump, ...]:
        return tuple(
            CompactMomentBump(left, right, int(self.quadrature_points))
            for left, right in self.angular_supports
        )

    @cached_property
    def axial_bumps(self) -> tuple[CompactMomentBump, ...]:
        return tuple(
            CompactMomentBump(left, right, int(self.quadrature_points))
            for left, right in self.axial_supports
        )

    @staticmethod
    def _moment_matrix(
        powers: Sequence[float], bumps: Sequence[CompactMomentBump], name: str
    ) -> np.ndarray:
        matrix = np.array(
            [[bump.moment(power) for bump in bumps] for power in powers], dtype=float
        )
        if not np.all(np.isfinite(matrix)):
            raise ValueError(f"{name} contains nonfinite entries")
        singular = np.linalg.svd(matrix, compute_uv=False)
        threshold = np.finfo(float).eps * matrix.shape[0] * singular[0]
        if singular[-1] <= threshold:
            raise ValueError(f"{name} is singular or too ill-conditioned in float64")
        matrix.setflags(write=False)
        return matrix

    @cached_property
    def angular_matrix(self) -> np.ndarray:
        return self._moment_matrix(self.angular_powers, self.angular_bumps, "angular matrix")

    @cached_property
    def axial_matrix(self) -> np.ndarray:
        return self._moment_matrix(self.axial_powers, self.axial_bumps, "axial matrix")

    def solve(self) -> "MeanRankCorrection":
        angular = np.linalg.solve(self.angular_matrix, self.angular_targets)
        axial = np.linalg.solve(self.axial_matrix, self.axial_targets)
        return MeanRankCorrection(self, angular, axial)

    @property
    def physical_row_targets(self) -> np.ndarray:
        d0, d1, d2 = self.debt
        result = np.array([0.0, 0.0, -d0, -d1, -d2], dtype=float)
        result.setflags(write=False)
        return result

    def rows_from_normalized_moments(
        self, angular: Sequence[float], axial: Sequence[float]
    ) -> np.ndarray:
        """Apply the exact MeanRankUpdate scaling identities to five moments.

        The inputs are ordered by :attr:`angular_powers` and
        :attr:`axial_powers`.  This is finite-dimensional theorem algebra, not
        a claim that a supplied function actually has those moments.
        """
        am = np.asarray(angular, dtype=float)
        ax = np.asarray(axial, dtype=float)
        if am.shape != (3,) or ax.shape != (2,) or not (
            np.all(np.isfinite(am)) and np.all(np.isfinite(ax))
        ):
            raise ValueError("angular/axial moments must have shapes (3,) and (2,)")
        result = np.array(
            [
                self.ell**3 * self.U * am[0],
                self.ell**2 * self.U * ax[0],
                self.U**2 * (2.0 * self.C * am[1]),
                self.ell**3 * self.U**2 * (self.C * ax[1]),
                self.ell**2 * self.U**2 * (-self.C * am[2]),
            ],
            dtype=float,
        )
        result.setflags(write=False)
        return result


@dataclass(frozen=True)
class MeanRankCorrection:
    """Executable compact correction for one supplied five-row debt."""

    plan: MeanRankUpdatePlan
    angular_coefficients: np.ndarray
    axial_coefficients: np.ndarray

    def __post_init__(self) -> None:
        angular = np.asarray(self.angular_coefficients, dtype=float)
        axial = np.asarray(self.axial_coefficients, dtype=float)
        if angular.shape != (3,) or axial.shape != (2,) or not (
            np.all(np.isfinite(angular)) and np.all(np.isfinite(axial))
        ):
            raise ValueError("correction coefficients must have shapes (3,) and (2,)")
        angular = angular.copy()
        axial = axial.copy()
        angular.setflags(write=False)
        axial.setflags(write=False)
        object.__setattr__(self, "angular_coefficients", angular)
        object.__setattr__(self, "axial_coefficients", axial)

    @property
    def paper_exact(self) -> bool:
        return False

    def _normalized_angular(self, R: float) -> float:
        R = _finite(R, "R")
        return float(
            sum(c * bump(R) for c, bump in zip(self.angular_coefficients, self.plan.angular_bumps))
        )

    def _normalized_axial(self, R: float) -> float:
        R = _finite(R, "R")
        return float(
            sum(c * bump(R) for c, bump in zip(self.axial_coefficients, self.plan.axial_bumps))
        )

    def angular_increment(self, radius: float) -> float:
        radius = _finite(radius, "radius")
        return self.plan.U * self._normalized_angular(radius / self.plan.ell)

    def desired_axial_increment(self, radius: float) -> float:
        radius = _finite(radius, "radius")
        return self.plan.U * self._normalized_axial(radius / self.plan.ell)

    def background(self, radius: float) -> float:
        """Physical power-law patch ``scaleField ell U (C R^(-1-2lam))``."""
        radius = _finite(radius, "radius")
        if radius <= 0.0:
            raise ValueError("background is evaluated only at positive radius")
        R = radius / self.plan.ell
        return self.plan.U * self.plan.C * R ** (-1.0 - 2.0 * self.plan.lam)

    @property
    def physical_angular_supports(self) -> tuple[tuple[float, float], ...]:
        return tuple(
            (self.plan.ell * left, self.plan.ell * right)
            for left, right in self.plan.angular_supports
        )

    @property
    def physical_axial_supports(self) -> tuple[tuple[float, float], ...]:
        return tuple(
            (self.plan.ell * left, self.plan.ell * right)
            for left, right in self.plan.axial_supports
        )
