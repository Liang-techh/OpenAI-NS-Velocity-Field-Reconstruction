"""Exact clean-room checks for one oscillatory cylindrical-curl seam.

This module intentionally covers only a local algebraic interface: transverse
phase-polarization recovery and the non-phase part of the normalized
cylindrical curl.  It does not construct the full oscillatory witness.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Tuple

Exact = Fraction
Vector3 = Tuple[Exact, Exact, Exact]

SOURCE_COMPONENT = "NS-oscillations"
SOURCE_ARCHIVE_PATH = "proof_sources/oscillations/oscillations_body.tex"
SOURCE_COMPONENT_SHA256 = (
    "3fd5c61de39b6f65a581ead5b29c741c2e0f46fb30f62e582dc24d93bbe83615"
)
SOURCE_ASSOCIATED_CHECKER = "proof_sources/oscillations/exact_checks.py"
SOURCE_CHECKER_SHA256 = (
    "a44ba49990c3b604c883d5b4e1726cb33068867cd85b23f3a2c22847e4ce8184"
)
SOURCE_RESULT_SHA256 = (
    "30d7a4d934e984a77cfcb04188d72c35391c46e38bcbf25a2024e38bcb783a12"
)


def _exact(value: Exact | int) -> Exact:
    if isinstance(value, bool):
        raise TypeError("boolean inputs are not accepted as exact scalars")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    raise TypeError("expected fractions.Fraction or int; floating inputs fail closed")


def _vector3(values: Iterable[Exact | int]) -> Vector3:
    data = tuple(_exact(value) for value in values)
    if len(data) != 3:
        raise ValueError("expected exactly three cylindrical/vector components")
    return data  # type: ignore[return-value]


def dot(left: Iterable[Exact | int], right: Iterable[Exact | int]) -> Exact:
    a = _vector3(left)
    b = _vector3(right)
    return sum((x * y for x, y in zip(a, b)), Fraction(0, 1))


def cross(left: Iterable[Exact | int], right: Iterable[Exact | int]) -> Vector3:
    a = _vector3(left)
    b = _vector3(right)
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def transverse_phase_recovery(
    normal: Iterable[Exact | int], target: Iterable[Exact | int]
) -> Vector3:
    """Return -n x (n x t) / |n|^2 using exact rational arithmetic.

    For transverse target data (n dot t = 0), this is exactly t.  For data
    with a normal component it returns the tangential projection instead,
    which lets tests detect a missing transversality hypothesis.
    """

    n = _vector3(normal)
    t = _vector3(target)
    norm_sq = dot(n, n)
    if norm_sq == 0:
        raise ValueError("phase normal must be nonzero")
    nested = cross(n, cross(n, t))
    return tuple(-component / norm_sq for component in nested)  # type: ignore[return-value]


def vector_residual(actual: Iterable[Exact | int], expected: Iterable[Exact | int]) -> Vector3:
    a = _vector3(actual)
    b = _vector3(expected)
    return tuple(x - y for x, y in zip(a, b))  # type: ignore[return-value]


def normalized_cylindrical_curl_remainder(
    *,
    epsilon: Exact | int,
    radius: Exact | int,
    d_z_c_theta: Exact | int,
    d_theta_c_z: Exact | int,
    d_z_c_r: Exact | int,
    d_r_c_z: Exact | int,
    c_theta: Exact | int,
    d_r_c_theta: Exact | int,
    d_theta_c_r: Exact | int,
) -> Vector3:
    """Non-phase normalized cylindrical-curl remainder in (r, theta, z).

    The z component expands d_r(R*C_theta) exactly as
    C_theta + R*d_r(C_theta); retaining C_theta is the cylindrical product
    term that a Cartesianized transcription would miss.
    """

    eps = _exact(epsilon)
    r = _exact(radius)
    if r == 0:
        raise ValueError("cylindrical radius must be nonzero")

    dz_ct = _exact(d_z_c_theta)
    dt_cz = _exact(d_theta_c_z)
    dz_cr = _exact(d_z_c_r)
    dr_cz = _exact(d_r_c_z)
    ct = _exact(c_theta)
    dr_ct = _exact(d_r_c_theta)
    dt_cr = _exact(d_theta_c_r)

    return (
        -dz_ct + (eps / r) * dt_cz,
        dz_cr - eps * dr_cz,
        (eps / r) * (ct + r * dr_ct) - (eps / r) * dt_cr,
    )


def cartesianized_z_mutation(
    *,
    epsilon: Exact | int,
    radius: Exact | int,
    d_r_c_theta: Exact | int,
    d_theta_c_r: Exact | int,
) -> Exact:
    """Intentionally wrong z remainder with the cylindrical C_theta term dropped."""

    eps = _exact(epsilon)
    r = _exact(radius)
    if r == 0:
        raise ValueError("cylindrical radius must be nonzero")
    dr_ct = _exact(d_r_c_theta)
    dt_cr = _exact(d_theta_c_r)
    return eps * dr_ct - (eps / r) * dt_cr
