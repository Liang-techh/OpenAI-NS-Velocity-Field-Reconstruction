from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence


def _exact_fraction(value: Fraction) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Fraction):
        raise TypeError("PA.16 theorem data must be fractions.Fraction")
    return value


def _ordered_positive(points: Sequence[Fraction], expected: int) -> tuple[Fraction, ...]:
    if len(points) != expected:
        raise ValueError(f"expected exactly {expected} support coordinates")
    exact = tuple(_exact_fraction(x) for x in points)
    if any(x <= 0 for x in exact):
        raise ValueError("support coordinates must be positive")
    if any(a >= b for a, b in zip(exact, exact[1:])):
        raise ValueError("support coordinates must be strictly increasing")
    return exact


def _det3(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


@dataclass(frozen=True)
class PA16PowerBlockAudit:
    """Exact normalized point-evaluation audit for the public PA.16 power blocks."""

    u_points_y: tuple[Fraction, Fraction]
    e_points_y: tuple[Fraction, Fraction, Fraction]
    u_determinant: Fraction
    e_determinant: Fraction

    @property
    def invertible(self) -> bool:
        return self.u_determinant != 0 and self.e_determinant != 0

    @property
    def u_sign(self) -> int:
        return (self.u_determinant > 0) - (self.u_determinant < 0)

    @property
    def e_sign(self) -> int:
        return (self.e_determinant > 0) - (self.e_determinant < 0)


def audit_pa16_power_blocks(
    *,
    u_points_y: Sequence[Fraction],
    e_points_y: Sequence[Fraction],
) -> PA16PowerBlockAudit:
    """Audit PA.16 after the exact substitution y=x^(1/10).

    Nonzero source constants multiply rows and therefore do not affect
    invertibility. The normalized u block has powers (0, 6); the normalized
    e block has signed powers (5, 1, -9), matching the displayed PA.16 rows.
    """

    u = _ordered_positive(u_points_y, 2)
    e = _ordered_positive(e_points_y, 3)

    u_det = u[1] ** 6 - u[0] ** 6
    e_matrix = (
        tuple(y**5 for y in e),
        tuple(-y for y in e),
        tuple(y**-9 for y in e),
    )
    e_det = _det3(e_matrix)
    return PA16PowerBlockAudit(
        u_points_y=(u[0], u[1]),
        e_points_y=(e[0], e[1], e[2]),
        u_determinant=u_det,
        e_determinant=e_det,
    )
