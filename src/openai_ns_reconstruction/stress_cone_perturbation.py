"""Section 7 actual-vs-reference covariance perturbation certificate.

This module is a theorem-shaped, caller-parameterized bridge for the finite-
dimensional perturbation step behind Eqs. (7.28)-(7.29). It starts from the
landed exact signed reference cone and from *certified normalized column error
vectors* ``e_-`` and ``e_+`` satisfying ``|e_sigma| <= delta``. The resulting
columns are

``scale_- * ((-a,-b) + e_-)`` and
``scale_+ * ((-a,+b) + e_+)``.

The implementation then derives conservative determinant, inverse-norm and
coefficient-positivity margins without identifying caller data with the paper's
actual pulse covariance. It remains ``formal-structure`` until upstream code
constructs the true Eq. (7.27) integrals and proves the Eq. (7.28) error bound.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np

from .stress_cone import PositiveSignedStressDecomposition


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _nonnegative(value: float, name: str) -> float:
    value = _finite(value, name)
    if value < 0.0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _vector2(values: Iterable[float], name: str) -> tuple[float, float]:
    arr = np.asarray(tuple(values), dtype=float)
    if arr.shape != (2,):
        raise ValueError(f"{name} must have shape (2,)")
    if not np.isfinite(arr).all():
        raise ValueError(f"{name} must contain finite values")
    return float(arr[0]), float(arr[1])


def _up(value: float, name: str) -> float:
    value = _finite(value, name)
    out = math.nextafter(value, math.inf)
    if not math.isfinite(out):
        raise OverflowError(f"{name} overflowed")
    return out


def _down(value: float, name: str) -> float:
    value = _finite(value, name)
    return math.nextafter(value, -math.inf)


def _mul_up(x: float, y: float, name: str) -> float:
    if x < 0.0 or y < 0.0:
        raise ValueError(f"{name} expects nonnegative factors")
    return _up(x * y, name)


def _mul_down(x: float, y: float, name: str) -> float:
    if x < 0.0 or y < 0.0:
        raise ValueError(f"{name} expects nonnegative factors")
    return _down(x * y, name)


def _add_up(x: float, y: float, name: str) -> float:
    if x < 0.0 or y < 0.0:
        raise ValueError(f"{name} expects nonnegative summands")
    return _up(x + y, name)


def _sub_down(x: float, y: float, name: str) -> float:
    if x < 0.0 or y < 0.0:
        raise ValueError(f"{name} expects nonnegative terms")
    return _down(x - y, name)


def _div_up(x: float, y: float, name: str) -> float:
    if x < 0.0 or y <= 0.0:
        raise ValueError(f"{name} expects nonnegative numerator and positive denominator")
    return _up(x / y, name)


def _explicit_solve_2x2(matrix: np.ndarray, target: np.ndarray) -> np.ndarray:
    h00, h01 = float(matrix[0, 0]), float(matrix[0, 1])
    h10, h11 = float(matrix[1, 0]), float(matrix[1, 1])
    det = h00 * h11 - h01 * h10
    if not math.isfinite(det) or det == 0.0:
        raise FloatingPointError("actual covariance matrix is singular in binary64")
    y0 = (h11 * float(target[0]) - h01 * float(target[1])) / det
    y1 = (-h10 * float(target[0]) + h00 * float(target[1])) / det
    out = np.array([y0, y1], dtype=float)
    if not np.isfinite(out).all():
        raise FloatingPointError("actual squared-amplitude solve is not finite")
    return out


@dataclass(frozen=True)
class CovariancePerturbationCertificate:
    """Fail-closed finite-dimensional implication for Eqs. (7.28)-(7.29).

    ``reference`` supplies the strict signed reference cone. ``error_minus``
    and ``error_plus`` are the normalized vector errors ``e_sigma`` from the
    paper-shaped relation

    ``H_sigma = scale_sigma * (reference_direction_sigma + e_sigma)``.

    ``normalized_error_bound`` must be a previously certified uniform bound on
    both Euclidean error norms. Supplying a number here is not itself a proof
    of Eq. (7.28); the constructor checks the concrete vectors against it and
    only derives the downstream finite-dimensional consequences.
    """

    reference: PositiveSignedStressDecomposition
    error_minus: tuple[float, float]
    error_plus: tuple[float, float]
    normalized_error_bound: float

    def __post_init__(self) -> None:
        if not isinstance(self.reference, PositiveSignedStressDecomposition):
            raise TypeError("reference must be a PositiveSignedStressDecomposition")
        em = _vector2(self.error_minus, "error_minus")
        ep = _vector2(self.error_plus, "error_plus")
        delta = _nonnegative(self.normalized_error_bound, "normalized_error_bound")
        object.__setattr__(self, "error_minus", em)
        object.__setattr__(self, "error_plus", ep)
        object.__setattr__(self, "normalized_error_bound", delta)

        em_norm_up = _up(math.hypot(*em), "error_minus norm")
        ep_norm_up = _up(math.hypot(*ep), "error_plus norm")
        if em_norm_up > delta or ep_norm_up > delta:
            raise ValueError("normalized column error exceeds certified bound")

        if self.normalized_determinant_margin <= 0.0:
            raise ValueError("column-error bound does not preserve determinant nondegeneracy")
        if self.squared_amplitude_lower_bound <= 0.0:
            raise ValueError("column-error bound does not preserve positive squared amplitudes")

        actual = self.actual_squared_amplitudes
        if not np.all(actual > 0.0):
            raise FloatingPointError(
                "analytic positivity margin passed but binary64 actual solve lost positivity"
            )

    @property
    def reference_direction_norm_upper(self) -> float:
        return _up(math.hypot(self.reference.a, self.reference.b), "reference direction norm")

    @property
    def normalized_determinant_error_upper(self) -> float:
        """Upper bound for ``|det(V)-det(V0)|`` before column scales.

        For two error vectors of norm at most ``delta`` and reference column
        norm ``r=sqrt(a^2+b^2)``, multilinearity gives
        ``2*r*delta + delta^2``.
        """
        r = self.reference_direction_norm_upper
        delta = self.normalized_error_bound
        linear = _mul_up(_mul_up(2.0, r, "det linear factor"), delta, "det linear error")
        quadratic = _mul_up(delta, delta, "det quadratic error")
        return _add_up(linear, quadratic, "normalized determinant error")

    @property
    def normalized_reference_determinant_lower(self) -> float:
        return _mul_down(
            _mul_down(2.0, self.reference.a, "reference determinant 2a"),
            self.reference.b,
            "reference determinant 2ab",
        )

    @property
    def normalized_determinant_margin(self) -> float:
        return _sub_down(
            self.normalized_reference_determinant_lower,
            self.normalized_determinant_error_upper,
            "normalized determinant margin",
        )

    @property
    def determinant_abs_lower_bound(self) -> float:
        scale_product = _mul_down(
            self.reference.scale_minus,
            self.reference.scale_plus,
            "column scale product",
        )
        return _mul_down(
            scale_product,
            self.normalized_determinant_margin,
            "actual determinant lower bound",
        )

    @property
    def actual_matrix(self) -> np.ndarray:
        a, b = self.reference.a, self.reference.b
        em = np.asarray(self.error_minus, dtype=float)
        ep = np.asarray(self.error_plus, dtype=float)
        col_minus = self.reference.scale_minus * (np.array([-a, -b]) + em)
        col_plus = self.reference.scale_plus * (np.array([-a, b]) + ep)
        out = np.column_stack((col_minus, col_plus))
        if not np.isfinite(out).all():
            raise OverflowError("actual covariance columns overflowed")
        return out

    @property
    def actual_determinant(self) -> float:
        matrix = self.actual_matrix
        out = float(matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0])
        if not math.isfinite(out):
            raise OverflowError("actual determinant overflowed")
        return out

    @property
    def matrix_perturbation_norm_upper(self) -> float:
        scales = _up(
            math.hypot(self.reference.scale_minus, self.reference.scale_plus),
            "column-scale Euclidean norm",
        )
        return _mul_up(
            self.normalized_error_bound,
            scales,
            "matrix perturbation norm",
        )

    @property
    def actual_matrix_frobenius_upper(self) -> float:
        direction = _add_up(
            self.reference_direction_norm_upper,
            self.normalized_error_bound,
            "perturbed direction norm",
        )
        scales = _up(
            math.hypot(self.reference.scale_minus, self.reference.scale_plus),
            "column-scale Euclidean norm",
        )
        return _mul_up(direction, scales, "actual Frobenius norm bound")

    @property
    def inverse_two_norm_upper(self) -> float:
        """Conservative spectral inverse bound from ``sigma_max/|det|``."""
        return _div_up(
            self.actual_matrix_frobenius_upper,
            self.determinant_abs_lower_bound,
            "inverse two-norm upper bound",
        )

    @property
    def reference_squared_amplitudes(self) -> np.ndarray:
        return self.reference.squared_amplitudes.copy()

    @property
    def coefficient_error_two_norm_upper(self) -> float:
        y0_norm = _up(
            math.hypot(*map(float, self.reference_squared_amplitudes)),
            "reference coefficient norm",
        )
        first = _mul_up(
            self.inverse_two_norm_upper,
            self.matrix_perturbation_norm_upper,
            "coefficient error operator factor",
        )
        return _mul_up(first, y0_norm, "coefficient error bound")

    @property
    def squared_amplitude_lower_bound(self) -> float:
        reference_min = min(map(float, self.reference_squared_amplitudes))
        reference_min_lower = math.nextafter(reference_min, 0.0)
        return _sub_down(
            reference_min_lower,
            self.coefficient_error_two_norm_upper,
            "positive coefficient lower bound",
        )

    @property
    def actual_squared_amplitudes(self) -> np.ndarray:
        return _explicit_solve_2x2(self.actual_matrix, self.reference.target)

    @property
    def actual_amplitudes(self) -> np.ndarray:
        out = np.sqrt(self.actual_squared_amplitudes)
        if not np.isfinite(out).all() or not np.all(out > 0.0):
            raise FloatingPointError("actual positive amplitude square root failed")
        return out

    @property
    def reconstructed_target(self) -> np.ndarray:
        return self.actual_matrix @ self.actual_squared_amplitudes
