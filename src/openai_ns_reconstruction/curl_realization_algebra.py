"""Coefficient-level curl-realization algebra for Section 7 oscillatory waves.

This module transcribes the finite-dimensional identities used by the pinned
official Lean modules ``NavierStokes/CurlClassBounds.lean`` and
``NavierStokes/LocalizedCurlRealization.lean`` before any spatial derivative
or support/gluing theorem is invoked.

For a nonzero phase normal ``n`` and complex raw coefficient ``a`` define

``B = |n|^-2 (n x a)``.

Then the vector triple-product identity gives

``n x B = n (n.a)/|n|^2 - a``.

Consequently the principal harmonic curl factor
``-(n x B)`` is exactly the tangential projection of ``a`` and equals ``a``
when the upstream construction has proved ``n.a = 0``.  With carrier
frequency ``K != 0``, the paper/Lean potential uses ``i/K`` and the remaining
coefficient derivative enters as ``(i/K) curl(B)``.

This is deliberately only ``formal-structure``.  A caller-supplied value for
``curl(B)`` is not evidence that it is the cylindrical curl of the paper's
coefficient.  This module does not construct the actual Proposition 5.5 base
field, amplitude ODE, phase patch, cutoff, cylindrical geometry, derivatives,
or supported physical potential, and it does not itself certify a
paper-exact divergence-free velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np


def _real3(value: Iterable[float], name: str) -> np.ndarray:
    out = np.asarray(tuple(value), dtype=float)
    if out.shape != (3,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite real three-vector")
    return out


def _complex3(value: Iterable[complex], name: str) -> np.ndarray:
    out = np.asarray(tuple(value), dtype=complex)
    if out.shape != (3,) or not np.all(np.isfinite(out.real)) or not np.all(np.isfinite(out.imag)):
        raise ValueError(f"{name} must be a finite complex three-vector")
    return out


def _nonzero_frequency(value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or value == 0.0:
        raise ValueError("frequency must be finite and nonzero")
    return value


def normal_dot(normal: Iterable[float], amplitude: Iterable[complex]) -> complex:
    """Lean ``normalDot`` for a real phase normal and complex coefficient."""
    n = _real3(normal, "normal")
    a = _complex3(amplitude, "amplitude")
    return complex(np.dot(n.astype(complex), a))


def normal_cross(normal: Iterable[float], amplitude: Iterable[complex]) -> np.ndarray:
    """Lean ``normalCross``: the ordinary cross product ``n x a``."""
    n = _real3(normal, "normal")
    a = _complex3(amplitude, "amplitude")
    return np.cross(n.astype(complex), a)


def normal_coefficient(normal: Iterable[float], amplitude: Iterable[complex]) -> np.ndarray:
    """Return ``|n|^-2 (n x a)`` from Formula (30)'s vector potential."""
    n = _real3(normal, "normal")
    a = _complex3(amplitude, "amplitude")
    norm2 = float(np.dot(n, n))
    if norm2 == 0.0:
        raise ValueError("phase normal must be nonzero")
    return np.cross(n.astype(complex), a) / norm2


def tangential_projection(normal: Iterable[float], amplitude: Iterable[complex]) -> np.ndarray:
    """Orthogonally remove the complex normal component of ``amplitude``."""
    n = _real3(normal, "normal")
    a = _complex3(amplitude, "amplitude")
    norm2 = float(np.dot(n, n))
    if norm2 == 0.0:
        raise ValueError("phase normal must be nonzero")
    return a - n.astype(complex) * (np.dot(n.astype(complex), a) / norm2)


def principal_curl_coefficient(normal: Iterable[float], amplitude: Iterable[complex]) -> np.ndarray:
    """Principal harmonic curl coefficient ``-n x normalCoefficient(n,a)``.

    The triple-product identity makes this equal to the tangential projection
    for arbitrary ``a``.  It equals ``a`` only after exact tangency is supplied
    by the upstream theorem.
    """
    n = _real3(normal, "normal")
    a = _complex3(amplitude, "amplitude")
    return -np.cross(n.astype(complex), normal_coefficient(n, a))


def inverse_carrier(frequency: float) -> complex:
    """Pinned Lean ``inverseCarrier K = i/K``."""
    K = _nonzero_frequency(frequency)
    return 1j / K


def phase_factor(frequency: float) -> complex:
    """Harmonic derivative factor ``i K`` used by the pinned Lean identity."""
    K = _nonzero_frequency(frequency)
    return 1j * K


def curl_remainder_from_coefficient_curl(
    frequency: float,
    coefficient_curl: Iterable[complex],
) -> np.ndarray:
    """Return the displayed derivative remainder ``(i/K) curl(B)``.

    ``coefficient_curl`` must come from an independently justified cylindrical
    derivative calculation.  This function only applies the exact scalar
    factor; it does not compute or certify that derivative.
    """
    curl_b = _complex3(coefficient_curl, "coefficient_curl")
    return inverse_carrier(frequency) * curl_b


@dataclass(frozen=True)
class TangentCurlRealization:
    """Fail-closed adapter for the tangent specialization of Formula (30).

    Exact binary arithmetic equality ``n.a == 0`` is intentionally required
    here.  Numerical near-tangency is not promoted to theorem evidence.  For
    approximate data use :func:`principal_curl_coefficient` directly and keep
    the normal-projection defect explicit.
    """

    frequency: float
    normal: tuple[float, float, float]
    amplitude: tuple[complex, complex, complex]

    def __post_init__(self) -> None:
        K = _nonzero_frequency(self.frequency)
        n = _real3(self.normal, "normal")
        a = _complex3(self.amplitude, "amplitude")
        if float(np.dot(n, n)) == 0.0:
            raise ValueError("phase normal must be nonzero")
        tangency = complex(np.dot(n.astype(complex), a))
        if tangency != 0j:
            raise ValueError(
                "exact tangency n.a=0 is required; numerical near-zero is not a certificate"
            )
        object.__setattr__(self, "frequency", K)

    @property
    def coefficient(self) -> np.ndarray:
        return normal_coefficient(self.normal, self.amplitude)

    @property
    def principal(self) -> np.ndarray:
        """Principal harmonic curl coefficient, exactly the supplied tangent ``a``."""
        return principal_curl_coefficient(self.normal, self.amplitude)

    @property
    def carrier_identity(self) -> complex:
        """Return ``inverseCarrier(K) * phaseFactor(K) = -1``."""
        return inverse_carrier(self.frequency) * phase_factor(self.frequency)

    def realized_coefficient(self, coefficient_curl: Iterable[complex]) -> np.ndarray:
        """Return ``a + (i/K) curl(B)`` under the exact tangency gate.

        The derivative value remains a caller-supplied formal input.  The
        actual divergence-free conclusion additionally needs the genuine
        cylindrical curl identity and smooth geometry from the upstream Lean
        hypotheses.
        """
        a = _complex3(self.amplitude, "amplitude")
        return a + curl_remainder_from_coefficient_curl(self.frequency, coefficient_curl)
