from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

Q = Fraction
Vec3 = tuple[Q, Q, Q]
Mat2 = tuple[tuple[Q, Q], tuple[Q, Q]]
Mat23 = tuple[Vec3, Vec3]
Mat32 = tuple[tuple[Q, Q], tuple[Q, Q], tuple[Q, Q]]
Mat33 = tuple[Vec3, Vec3, Vec3]


def _q(value: Q) -> Q:
    if not isinstance(value, Fraction):
        raise TypeError("exact replay inputs must be fractions.Fraction")
    return value


def _matmul(a: Sequence[Sequence[Q]], b: Sequence[Sequence[Q]]) -> tuple[tuple[Q, ...], ...]:
    if not a or not b:
        raise ValueError("matrices must be nonempty")
    width = len(a[0])
    if any(len(row) != width for row in a):
        raise ValueError("left matrix is ragged")
    b_width = len(b[0])
    if any(len(row) != b_width for row in b):
        raise ValueError("right matrix is ragged")
    if width != len(b):
        raise ValueError("matrix dimensions do not match")
    return tuple(
        tuple(sum((a[i][k] * b[k][j] for k in range(width)), Q(0)) for j in range(b_width))
        for i in range(len(a))
    )


def _sub(a: Sequence[Sequence[Q]], b: Sequence[Sequence[Q]]) -> tuple[tuple[Q, ...], ...]:
    if len(a) != len(b) or any(len(ar) != len(br) for ar, br in zip(a, b)):
        raise ValueError("matrix dimensions do not match")
    return tuple(tuple(x - y for x, y in zip(ar, br)) for ar, br in zip(a, b))


def _identity(n: int) -> tuple[tuple[Q, ...], ...]:
    return tuple(tuple(Q(1) if i == j else Q(0) for j in range(n)) for i in range(n))


def _transpose(a: Sequence[Sequence[Q]]) -> tuple[tuple[Q, ...], ...]:
    if not a:
        raise ValueError("matrix must be nonempty")
    width = len(a[0])
    if any(len(row) != width for row in a):
        raise ValueError("matrix is ragged")
    return tuple(tuple(a[i][j] for i in range(len(a))) for j in range(width))


# Defined separately instead of relying on numpy so the replay is exact and dependency-free.
def outer3(a: Vec3, b: Vec3) -> Mat33:
    return (
        (a[0] * b[0], a[0] * b[1], a[0] * b[2]),
        (a[1] * b[0], a[1] * b[1], a[1] * b[2]),
        (a[2] * b[0], a[2] * b[1], a[2] * b[2]),
    )


def max_abs_entry(a: Sequence[Sequence[Q]]) -> Q:
    values = [abs(x) for row in a for x in row]
    return max(values, default=Q(0))


@dataclass(frozen=True)
class MovingPlaneFixture:
    s_a: Q
    k_theta: Q
    k_z: Q
    tangential_norm: Q
    gamma: Q

    def __post_init__(self) -> None:
        for value in (self.s_a, self.k_theta, self.k_z, self.tangential_norm, self.gamma):
            _q(value)
        if self.k_theta * self.k_theta + self.k_z * self.k_z != 1:
            raise ValueError("K_a must be an exact unit tangential vector")
        if self.tangential_norm <= 0:
            raise ValueError("|n_{Phi,tan}| must be positive")
        if self.gamma == 0:
            raise ValueError("the two-column coordinate matrix must be invertible")

    @property
    def e_r(self) -> Vec3:
        return (Q(1), Q(0), Q(0))

    @property
    def k_a(self) -> Vec3:
        return (Q(0), self.k_theta, self.k_z)

    @property
    def n_a(self) -> Vec3:
        return (Q(0), self.k_z, -self.k_theta)

    @property
    def n_phi(self) -> Vec3:
        n = self.tangential_norm
        return (n * self.s_a, n * self.k_theta, n * self.k_z)

    @property
    def u(self) -> Mat32:
        k = self.k_a
        n = self.n_a
        return (
            (Q(1), Q(0)),
            (-self.s_a * k[1], n[1]),
            (-self.s_a * k[2], n[2]),
        )

    @property
    def m(self) -> Mat2:
        g = self.gamma
        return ((Q(1), Q(1)), (g, -g))

    @property
    def m_inverse(self) -> Mat2:
        g = self.gamma
        return ((Q(1, 2), Q(1, 2) / g), (Q(1, 2), -Q(1, 2) / g))

    @property
    def w(self) -> Mat23:
        return (self.e_r, self.n_a)

    @property
    def b(self) -> Mat32:
        product = _matmul(self.u, self.m)
        return tuple(tuple(x for x in row) for row in product)  # type: ignore[return-value]

    @property
    def b_left(self) -> Mat23:
        product = _matmul(self.m_inverse, self.w)
        return tuple(tuple(x for x in row) for row in product)  # type: ignore[return-value]

    def ambient_projector(self, *, n_phi: Vec3 | None = None) -> Mat33:
        n = self.n_phi if n_phi is None else n_phi
        if any(not isinstance(x, Fraction) for x in n):
            raise TypeError("normal components must be fractions.Fraction")
        scaled_outer = outer3(self.k_a, n)
        identity = _identity(3)
        return tuple(
            tuple(identity[i][j] - scaled_outer[i][j] / self.tangential_norm for j in range(3))
            for i in range(3)
        )  # type: ignore[return-value]


@dataclass(frozen=True)
class PlaneInverseReplay:
    left_inverse_residual: tuple[tuple[Q, ...], ...]
    ambient_formula_residual: tuple[tuple[Q, ...], ...]
    projector_idempotence_residual: tuple[tuple[Q, ...], ...]
    normal_range_residual: tuple[tuple[Q, ...], ...]
    kernel_residual: tuple[tuple[Q, ...], ...]

    @property
    def passed(self) -> bool:
        return all(
            max_abs_entry(matrix) == 0
            for matrix in (
                self.left_inverse_residual,
                self.ambient_formula_residual,
                self.projector_idempotence_residual,
                self.normal_range_residual,
                self.kernel_residual,
            )
        )


def replay_plane_inverse(fixture: MovingPlaneFixture, *, projector_normal: Vec3 | None = None) -> PlaneInverseReplay:
    b = fixture.b
    b_left = fixture.b_left
    b_left_b = _matmul(b_left, b)
    bb_left = _matmul(b, b_left)
    projector = fixture.ambient_projector(n_phi=projector_normal)
    projector_squared = _matmul(projector, projector)
    normal_row = (fixture.n_phi,)
    normal_range = _matmul(normal_row, projector)
    kernel = _matmul(projector, _transpose((fixture.k_a,)))
    return PlaneInverseReplay(
        left_inverse_residual=_sub(b_left_b, _identity(2)),
        ambient_formula_residual=_sub(bb_left, projector),
        projector_idempotence_residual=_sub(projector_squared, projector),
        normal_range_residual=normal_range,
        kernel_residual=kernel,
    )


PAPER_EXACT_VELOCITY_AVAILABLE = False
FULL_RECONSTRUCTION = False
