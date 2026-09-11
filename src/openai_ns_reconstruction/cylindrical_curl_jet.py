"""Structured cylindrical coefficient derivatives for Section 7 curl realization.

The landed :mod:`curl_realization_algebra` module implements Formula (30)'s
finite-dimensional split but deliberately accepts ``curl(B)`` as a supplied
complex three-vector.  This module moves one step closer to the genuine
cylindrical construction: callers provide the coordinate derivative jet of the
coefficient ``B`` and the cylindrical curl is assembled here, including the
moving-frame ``B_theta/r`` term.

This remains ``formal-structure``.  Supplying derivative vectors is not proof
that they are derivatives of the paper's localized coefficient, and this
module does not establish support, smoothness, phase-patch compatibility or
``div(curl)=0`` for an actual paper wave.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np

from .curl_realization_algebra import TangentCurlRealization, inverse_carrier


def _complex3(value: Iterable[complex], name: str) -> np.ndarray:
    out = np.asarray(tuple(value), dtype=complex)
    if (
        out.shape != (3,)
        or not np.all(np.isfinite(out.real))
        or not np.all(np.isfinite(out.imag))
    ):
        raise ValueError(f"{name} must be a finite complex three-vector")
    return out


def cylindrical_coefficient_curl(
    *,
    radius: float,
    coefficient: Iterable[complex],
    radial_derivative: Iterable[complex],
    angular_derivative: Iterable[complex],
    axial_derivative: Iterable[complex],
) -> np.ndarray:
    """Return the cylindrical curl of a complex coefficient at one point.

    ``coefficient`` contains the physical components
    ``(B_r, B_theta, B_z)`` in the orthonormal cylindrical frame.  The three
    derivative arguments are coordinate derivatives of those *components*
    with respect to ``r``, ``theta`` and ``z``.  For ``r>0`` the standard
    cylindrical formula is

    ``curl(B)_r     = (1/r) d_theta B_z - d_z B_theta``
    ``curl(B)_theta = d_z B_r - d_r B_z``
    ``curl(B)_z     = d_r B_theta + B_theta/r - (1/r) d_theta B_r``.

    The ``B_theta/r`` term is intentionally explicit.  Omitting it would treat
    cylindrical components as if their basis were Cartesian and would not be
    the derivative used by the paper/Lean cylindrical curl realization.
    """
    r = float(radius)
    if not math.isfinite(r) or r <= 0.0:
        raise ValueError("cylindrical coefficient curl requires finite radius>0")

    b = _complex3(coefficient, "coefficient")
    d_r = _complex3(radial_derivative, "radial_derivative")
    d_theta = _complex3(angular_derivative, "angular_derivative")
    d_z = _complex3(axial_derivative, "axial_derivative")

    return np.array(
        [
            d_theta[2] / r - d_z[1],
            d_z[0] - d_r[2],
            d_r[1] + b[1] / r - d_theta[0] / r,
        ],
        dtype=complex,
    )


@dataclass(frozen=True)
class TangentCylindricalCurlJet:
    """Formula (30) with a structured cylindrical derivative jet.

    Exact tangency and the pinned coefficient value
    ``B = |n|^-2 (n cross a)`` are inherited from
    :class:`~openai_ns_reconstruction.curl_realization_algebra.TangentCurlRealization`.
    The caller supplies only the three coordinate derivative vectors of that
    coefficient.  Thus an arbitrary caller-supplied ``curl(B)`` can no longer
    be inserted at this layer.

    The derivative jet itself is still an upstream hypothesis.  Until it is
    produced from the actual localized amplitude/cutoff/phase data, this class
    is not a paper-exact wave constructor.
    """

    frequency: float
    normal: tuple[float, float, float]
    amplitude: tuple[complex, complex, complex]
    radius: float
    radial_derivative: tuple[complex, complex, complex]
    angular_derivative: tuple[complex, complex, complex]
    axial_derivative: tuple[complex, complex, complex]

    def __post_init__(self) -> None:
        tangent = TangentCurlRealization(self.frequency, self.normal, self.amplitude)
        r = float(self.radius)
        if not math.isfinite(r) or r <= 0.0:
            raise ValueError("cylindrical coefficient curl requires finite radius>0")
        _complex3(self.radial_derivative, "radial_derivative")
        _complex3(self.angular_derivative, "angular_derivative")
        _complex3(self.axial_derivative, "axial_derivative")
        object.__setattr__(self, "frequency", tangent.frequency)
        object.__setattr__(self, "radius", r)

    @property
    def coefficient(self) -> np.ndarray:
        """Pinned normal coefficient ``|n|^-2 (n cross a)``."""
        return TangentCurlRealization(
            self.frequency, self.normal, self.amplitude
        ).coefficient

    @property
    def coefficient_curl(self) -> np.ndarray:
        """Cylindrical ``curl(B)`` assembled from the supplied derivative jet."""
        return cylindrical_coefficient_curl(
            radius=self.radius,
            coefficient=self.coefficient,
            radial_derivative=self.radial_derivative,
            angular_derivative=self.angular_derivative,
            axial_derivative=self.axial_derivative,
        )

    @property
    def realized_coefficient(self) -> np.ndarray:
        """Return the tangent Formula (30) coefficient ``a + (i/K) curl(B)``."""
        amplitude = _complex3(self.amplitude, "amplitude")
        return amplitude + inverse_carrier(self.frequency) * self.coefficient_curl
