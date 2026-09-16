from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

Vector = tuple[Fraction, ...]
Matrix = tuple[tuple[Fraction, ...], ...]


def _require_fraction(value: Fraction, name: str) -> None:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be fractions.Fraction")


def _validate_vector(vector: Vector, dimension: int, name: str) -> None:
    if len(vector) != dimension:
        raise ValueError(f"{name} has the wrong dimension")
    for index, value in enumerate(vector):
        _require_fraction(value, f"{name}[{index}]")


def _zero_matrix(dimension: int) -> list[list[Fraction]]:
    return [[Fraction(0) for _ in range(dimension)] for _ in range(dimension)]


def _freeze(matrix: list[list[Fraction]]) -> Matrix:
    return tuple(tuple(row) for row in matrix)


def matrix_add(*matrices: Matrix) -> Matrix:
    if not matrices:
        raise ValueError("at least one matrix is required")
    dimension = len(matrices[0])
    if dimension == 0 or any(len(row) != dimension for row in matrices[0]):
        raise ValueError("matrices must be nonempty and square")
    out = _zero_matrix(dimension)
    for matrix in matrices:
        if len(matrix) != dimension or any(len(row) != dimension for row in matrix):
            raise ValueError("all matrices must have the same square shape")
        for i in range(dimension):
            for j in range(dimension):
                _require_fraction(matrix[i][j], f"matrix[{i}][{j}]")
                out[i][j] += matrix[i][j]
    return _freeze(out)


def matrix_sub(left: Matrix, right: Matrix) -> Matrix:
    dimension = len(left)
    if dimension == 0 or len(right) != dimension:
        raise ValueError("matrices must have the same nonzero dimension")
    out = _zero_matrix(dimension)
    for i in range(dimension):
        if len(left[i]) != dimension or len(right[i]) != dimension:
            raise ValueError("matrices must be square")
        for j in range(dimension):
            _require_fraction(left[i][j], f"left[{i}][{j}]")
            _require_fraction(right[i][j], f"right[{i}][{j}]")
            out[i][j] = left[i][j] - right[i][j]
    return _freeze(out)


@dataclass(frozen=True)
class ExactEnsemble:
    weights: tuple[Fraction, ...]
    samples: tuple[Vector, ...]

    def __post_init__(self) -> None:
        if not self.weights or len(self.weights) != len(self.samples):
            raise ValueError("weights and samples must be nonempty and aligned")
        dimension = len(self.samples[0])
        if dimension == 0:
            raise ValueError("sample dimension must be positive")
        total = Fraction(0)
        for index, weight in enumerate(self.weights):
            _require_fraction(weight, f"weights[{index}]")
            if weight <= 0:
                raise ValueError("weights must be positive")
            total += weight
            _validate_vector(self.samples[index], dimension, f"samples[{index}]")
        if total != 1:
            raise ValueError("weights must sum exactly to one")

    @property
    def dimension(self) -> int:
        return len(self.samples[0])

    def add(self, other: ExactEnsemble) -> ExactEnsemble:
        _require_compatible(self, other)
        samples = tuple(
            tuple(left + right for left, right in zip(a, b))
            for a, b in zip(self.samples, other.samples)
        )
        return ExactEnsemble(self.weights, samples)


def _require_compatible(left: ExactEnsemble, right: ExactEnsemble) -> None:
    if left.weights != right.weights:
        raise ValueError("ensembles must use identical exact averaging weights")
    if left.dimension != right.dimension or len(left.samples) != len(right.samples):
        raise ValueError("ensembles must have identical sample shape")


def covariance(field: ExactEnsemble) -> Matrix:
    out = _zero_matrix(field.dimension)
    for weight, vector in zip(field.weights, field.samples):
        for i in range(field.dimension):
            for j in range(field.dimension):
                out[i][j] += weight * vector[i] * vector[j]
    return _freeze(out)


def bilinear_covariance(left: ExactEnsemble, right: ExactEnsemble) -> Matrix:
    _require_compatible(left, right)
    out = _zero_matrix(left.dimension)
    for weight, x, y in zip(left.weights, left.samples, right.samples):
        for i in range(left.dimension):
            for j in range(left.dimension):
                out[i][j] += weight * (x[i] * y[j] + y[i] * x[j])
    return _freeze(out)


def covariance_increment_audit(
    w0: ExactEnsemble,
    error: ExactEnsemble,
    linear_inverse: ExactEnsemble,
    curl_remainder: ExactEnsemble,
    sigma: Matrix,
) -> dict[str, Matrix]:
    for field in (error, linear_inverse, curl_remainder):
        _require_compatible(w0, field)

    signed_inverse_response = bilinear_covariance(w0, linear_inverse)
    if signed_inverse_response != sigma:
        raise ValueError("the exact signed-inverse premise B(W0,L)=Sigma is not satisfied")

    u = w0.add(error)
    v = linear_inverse.add(curl_remainder)
    delta = matrix_sub(covariance(u.add(v)), covariance(u))

    error_linear = bilinear_covariance(error, linear_inverse)
    actual_remainder = bilinear_covariance(u, curl_remainder)
    quadratic = covariance(v)
    linear_quadratic = covariance(linear_inverse)
    cross_quadratic = bilinear_covariance(linear_inverse, curl_remainder)
    remainder_quadratic = covariance(curl_remainder)

    expanded = matrix_add(sigma, error_linear, actual_remainder, quadratic)
    quadratic_expanded = matrix_add(
        linear_quadratic, cross_quadratic, remainder_quadratic
    )
    if delta != expanded:
        raise ArithmeticError("covariance increment expansion did not close exactly")
    if quadratic != quadratic_expanded:
        raise ArithmeticError("quadratic covariance split did not close exactly")

    return {
        "sigma": sigma,
        "delta_covariance": delta,
        "error_linear_remainder": error_linear,
        "actual_curl_remainder": actual_remainder,
        "quadratic_remainder": quadratic,
        "linear_quadratic": linear_quadratic,
        "linear_curl_cross": cross_quadratic,
        "curl_quadratic": remainder_quadratic,
    }


def covariance_remainder_exponents(
    alpha: Fraction, rho: Fraction, kappa_s: Fraction
) -> dict[str, Fraction]:
    for value, name in ((alpha, "alpha"), (rho, "rho"), (kappa_s, "kappa_s")):
        _require_fraction(value, name)
    if kappa_s <= 0:
        raise ValueError("kappa_s must be positive")
    return {
        "B(E,L)": rho + alpha - Fraction(1, 2),
        "B(U,R)": alpha + Fraction(1, 2) - kappa_s,
        "C(L)": 2 * alpha - 1,
        "B(L,R)": 2 * alpha - Fraction(1, 2) - kappa_s,
        "C(R)": 2 * alpha - 2 * kappa_s,
    }
