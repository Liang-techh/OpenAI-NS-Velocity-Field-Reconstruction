"""Clean-room exact algebra for one Kokuno NS-stage_inputs covariance seam.

This module intentionally verifies only a finite two-component differential
right-inverse identity and squared-partition assembly.  It does not construct
an admissible stress, positive amplitudes, or an actual Section-9 stage.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Sequence

Exact = Fraction
Vector2 = tuple[Exact, Exact]
Matrix2 = tuple[Vector2, Vector2]


def _exact(value: int | Fraction, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an exact int/Fraction")
    return Fraction(value)


def _vector2(values: Sequence[int | Fraction], name: str) -> Vector2:
    if len(values) != 2:
        raise ValueError(f"{name} must have exactly two components")
    return (_exact(values[0], f"{name}[0]"), _exact(values[1], f"{name}[1]"))


def _matrix2(rows: Sequence[Sequence[int | Fraction]]) -> Matrix2:
    if len(rows) != 2:
        raise ValueError("H must have exactly two rows")
    return (_vector2(rows[0], "H[0]"), _vector2(rows[1], "H[1]"))


def _matvec(matrix: Matrix2, vector: Vector2) -> Vector2:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )


def invert_matrix2(rows: Sequence[Sequence[int | Fraction]]) -> Matrix2:
    """Return the exact inverse of a nonsingular 2x2 rational matrix."""

    matrix = _matrix2(rows)
    a, b = matrix[0]
    c, d = matrix[1]
    det = a * d - b * c
    if det == 0:
        raise ValueError("H must be invertible")
    return ((d / det, -b / det), (-c / det, a / det))


@dataclass(frozen=True)
class SignedDifferentialInverse:
    sigma: Vector2
    epsilon: Exact
    amplitudes: Vector2
    d_sigma: Vector2
    delta_amplitudes: Vector2
    component_variation: Vector2
    covariance_response: Vector2
    residual: Vector2

    @property
    def verified(self) -> bool:
        return self.residual == (Fraction(0), Fraction(0))


def signed_differential_inverse(
    H: Sequence[Sequence[int | Fraction]],
    sigma: Sequence[int | Fraction],
    epsilon: int | Fraction,
    amplitudes: Sequence[int | Fraction],
) -> SignedDifferentialInverse:
    """Reimplement the exact two-sign differential covariance inverse.

    Public source structure:
      d_Sigma = H^{-1}(Sigma / epsilon)
      delta a_sigma = (d_Sigma)_sigma / (2 a_sigma)

    The exact first covariance variation in H-coordinates is therefore
    2*epsilon*a_sigma*delta a_sigma = epsilon*(d_Sigma)_sigma, so applying H
    must recover Sigma exactly.
    """

    matrix = _matrix2(H)
    target = _vector2(sigma, "sigma")
    eps = _exact(epsilon, "epsilon")
    if eps <= 0:
        raise ValueError("epsilon must be strictly positive")
    amps = _vector2(amplitudes, "amplitudes")
    if amps[0] <= 0 or amps[1] <= 0:
        raise ValueError("both amplitudes must be strictly positive")

    inverse = invert_matrix2(matrix)
    scaled_target = (target[0] / eps, target[1] / eps)
    d_sigma = _matvec(inverse, scaled_target)
    delta_amplitudes = (
        d_sigma[0] / (2 * amps[0]),
        d_sigma[1] / (2 * amps[1]),
    )
    component_variation = (
        2 * eps * amps[0] * delta_amplitudes[0],
        2 * eps * amps[1] * delta_amplitudes[1],
    )
    response = _matvec(matrix, component_variation)
    residual = (response[0] - target[0], response[1] - target[1])
    return SignedDifferentialInverse(
        sigma=target,
        epsilon=eps,
        amplitudes=amps,
        d_sigma=d_sigma,
        delta_amplitudes=delta_amplitudes,
        component_variation=component_variation,
        covariance_response=response,
        residual=residual,
    )


def covariance_response_from_deltas(
    H: Sequence[Sequence[int | Fraction]],
    epsilon: int | Fraction,
    amplitudes: Sequence[int | Fraction],
    delta_amplitudes: Sequence[int | Fraction],
) -> Vector2:
    """Evaluate the exact linear covariance response for supplied deltas."""

    matrix = _matrix2(H)
    eps = _exact(epsilon, "epsilon")
    if eps <= 0:
        raise ValueError("epsilon must be strictly positive")
    amps = _vector2(amplitudes, "amplitudes")
    deltas = _vector2(delta_amplitudes, "delta_amplitudes")
    if amps[0] <= 0 or amps[1] <= 0:
        raise ValueError("both amplitudes must be strictly positive")
    variation = (
        2 * eps * amps[0] * deltas[0],
        2 * eps * amps[1] * deltas[1],
    )
    return _matvec(matrix, variation)


@dataclass(frozen=True)
class PartitionedCovarianceAssembly:
    sigma: Vector2
    square_sum: Exact
    assembled: Vector2
    residual: Vector2

    @property
    def verified(self) -> bool:
        return self.square_sum == 1 and self.residual == (Fraction(0), Fraction(0))


def assemble_squared_partition(
    sigma: Sequence[int | Fraction],
    weights: Iterable[int | Fraction],
) -> PartitionedCovarianceAssembly:
    """Assemble identical local covariance responses with squared weights.

    This finite algebra is the target-side seam for the public identity
    sum_beta eta_beta^2 = 1.  It fails closed unless that equality is exact.
    """

    target = _vector2(sigma, "sigma")
    exact_weights = tuple(_exact(w, f"weights[{i}]") for i, w in enumerate(weights))
    if not exact_weights:
        raise ValueError("at least one partition weight is required")
    square_sum = sum((w * w for w in exact_weights), Fraction(0))
    assembled = (square_sum * target[0], square_sum * target[1])
    residual = (assembled[0] - target[0], assembled[1] - target[1])
    return PartitionedCovarianceAssembly(
        sigma=target,
        square_sum=square_sum,
        assembled=assembled,
        residual=residual,
    )
