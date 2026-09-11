"""Section 7 signed two-slot stress-cone algebra.

This module implements the exact *reference* covariance solve used to expose the
strict positive cone before the paper perturbs to the actual integrated pulse
covariances.  The formulas are cross-checked against
``NavierStokes/Covariance.lean`` at the pinned upstream commit
``f9e8bc5b38b6e212696e8a30e3e91517af887bbd``.

The reference columns in the orthonormal normal/transverse coordinates are

``(-a, -b) * scale_minus`` and ``(-a, +b) * scale_plus``,

and the target is ``(-m, t)``.  The strict cone ``|a*t| < b*m`` makes the two
explicit squared amplitudes positive.  Their positive square roots therefore
realize the target covariance exactly.

This is deliberately only ``formal-structure`` for the reconstruction.  It does
**not** identify the paper's actual pulse-integrated matrix ``H`` with this
reference matrix, prove the Section 7 perturbation/inverse estimates, produce
the paper-exact stress target, or construct any oscillatory wave.  Caller data
must never be promoted to paper-exact merely because this finite-dimensional
solve succeeds.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


PINNED_LEAN_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _positive(value: float, name: str) -> float:
    value = _finite(value, name)
    if value <= 0.0:
        raise ValueError(f"{name} must be positive")
    return value


def normal_magnitude(c_star: float, u_star: float) -> float:
    """Reference normal magnitude ``-c_* sqrt(1+u_*^2)``.

    This is the quantity called ``normalMagnitude`` in the pinned
    ``Covariance.lean`` model.  The theorem-side signed convention requires
    ``c_* < 0`` and ``u_* > 0``.
    """
    c_star = _finite(c_star, "c_star")
    u_star = _positive(u_star, "u_star")
    if c_star >= 0.0:
        raise ValueError("reference signed-slot convention requires c_star<0")
    out = -c_star * math.sqrt(1.0 + u_star * u_star)
    if not math.isfinite(out) or out <= 0.0:
        raise OverflowError("normal magnitude is not a finite positive number")
    return out


def paper_ratio_margin(*, c_star: float, u_star: float, m: float, t: float) -> float:
    """Return the strict Section-7 reference ratio margin.

    The pinned Lean implication uses

    ``|c_* t / m| < u_* / sqrt(1+u_*^2)``.

    A positive return value is exactly the strict margin of that inequality.
    """
    c_star = _finite(c_star, "c_star")
    u_star = _positive(u_star, "u_star")
    m = _positive(m, "m")
    t = _finite(t, "t")
    if c_star >= 0.0:
        raise ValueError("reference signed-slot convention requires c_star<0")
    rhs = u_star / math.sqrt(1.0 + u_star * u_star)
    lhs = abs(c_star * t / m)
    margin = rhs - lhs
    if not math.isfinite(margin):
        raise OverflowError("paper ratio margin is not finite")
    return margin


def signed_covariance_matrix(
    *, a: float, b: float, scale_minus: float, scale_plus: float
) -> np.ndarray:
    """Exact two-column reference covariance matrix from ``Covariance.lean``."""
    a = _positive(a, "a")
    b = _positive(b, "b")
    scale_minus = _positive(scale_minus, "scale_minus")
    scale_plus = _positive(scale_plus, "scale_plus")
    matrix = np.array(
        [
            [-a * scale_minus, -a * scale_plus],
            [-b * scale_minus, b * scale_plus],
        ],
        dtype=float,
    )
    if not np.isfinite(matrix).all():
        raise OverflowError("reference covariance matrix overflowed")
    return matrix


def stress_target(*, m: float, t: float) -> np.ndarray:
    """Target vector ``(-m,t)`` in the reference normal/transverse frame."""
    m = _finite(m, "m")
    t = _finite(t, "t")
    target = np.array([-m, t], dtype=float)
    if not np.isfinite(target).all():
        raise OverflowError("stress target overflowed")
    return target


def squared_amplitude_coefficients(
    *,
    a: float,
    b: float,
    scale_minus: float,
    scale_plus: float,
    m: float,
    t: float,
) -> np.ndarray:
    """Explicit inverse solve for the two squared amplitudes.

    This function evaluates the pinned closed formulas.  It does not require the
    strict cone and therefore may return non-positive entries.  Use
    :class:`PositiveSignedStressDecomposition` when positivity is part of the
    certificate.
    """
    a = _positive(a, "a")
    b = _positive(b, "b")
    scale_minus = _positive(scale_minus, "scale_minus")
    scale_plus = _positive(scale_plus, "scale_plus")
    m = _finite(m, "m")
    t = _finite(t, "t")

    denom_minus = 2.0 * a * b * scale_minus
    denom_plus = 2.0 * a * b * scale_plus
    coeff = np.array(
        [
            (b * m - a * t) / denom_minus,
            (b * m + a * t) / denom_plus,
        ],
        dtype=float,
    )
    if not np.isfinite(coeff).all():
        raise OverflowError("squared-amplitude solve overflowed")
    return coeff


@dataclass(frozen=True)
class PositiveSignedStressDecomposition:
    """Fail-closed certificate for the exact signed reference cone solve.

    Parameters describe only the reference two-slot covariance.  In particular,
    ``scale_minus`` and ``scale_plus`` must already be certified positive by an
    upstream pulse/covariance argument; arbitrary caller values do not certify
    the paper's actual integrated covariance matrix.
    """

    a: float
    b: float
    scale_minus: float
    scale_plus: float
    m: float
    t: float

    def __post_init__(self) -> None:
        a = _positive(self.a, "a")
        b = _positive(self.b, "b")
        scale_minus = _positive(self.scale_minus, "scale_minus")
        scale_plus = _positive(self.scale_plus, "scale_plus")
        m = _positive(self.m, "m")
        t = _finite(self.t, "t")

        margin = b * m - abs(a * t)
        if not math.isfinite(margin):
            raise OverflowError("strict-cone margin overflowed")
        if margin <= 0.0:
            raise ValueError("strict cone requires |a*t| < b*m")

        # Evaluate the theorem formula now and fail closed if binary64 cannot
        # retain the theorem's strict positivity even though the scalar margin
        # passed.  This guards downstream square roots against silent underflow.
        coeff = squared_amplitude_coefficients(
            a=a,
            b=b,
            scale_minus=scale_minus,
            scale_plus=scale_plus,
            m=m,
            t=t,
        )
        if not np.all(coeff > 0.0):
            raise FloatingPointError(
                "strict cone was positive but binary64 lost positive squared amplitudes"
            )

    @property
    def determinant(self) -> float:
        """Pinned determinant ``-2*a*b*scale_minus*scale_plus``."""
        out = -2.0 * self.a * self.b * self.scale_minus * self.scale_plus
        if not math.isfinite(out):
            raise OverflowError("reference determinant overflowed")
        return out

    @property
    def cone_margin(self) -> float:
        return self.b * self.m - abs(self.a * self.t)

    @property
    def matrix(self) -> np.ndarray:
        return signed_covariance_matrix(
            a=self.a,
            b=self.b,
            scale_minus=self.scale_minus,
            scale_plus=self.scale_plus,
        )

    @property
    def target(self) -> np.ndarray:
        return stress_target(m=self.m, t=self.t)

    @property
    def squared_amplitudes(self) -> np.ndarray:
        return squared_amplitude_coefficients(
            a=self.a,
            b=self.b,
            scale_minus=self.scale_minus,
            scale_plus=self.scale_plus,
            m=self.m,
            t=self.t,
        )

    @property
    def amplitudes(self) -> np.ndarray:
        coeff = self.squared_amplitudes
        out = np.sqrt(coeff)
        if not np.isfinite(out).all() or not np.all(out > 0.0):
            raise FloatingPointError("positive amplitude square root failed")
        return out

    @property
    def reconstructed_target(self) -> np.ndarray:
        """Reference covariance reconstructed from the positive squared amplitudes."""
        return self.matrix @ np.square(self.amplitudes)

    def scaled_wave_amplitudes(self, *, epsilon: float, mask: float) -> np.ndarray:
        """Return ``sqrt(epsilon) * amplitude * mask`` for each signed slot.

        This mirrors ``Covariance.scaled_primary_covariance``.  ``mask`` may be
        any finite real because the covariance depends on ``mask^2``; this
        routine does not assert that it is an actual paper partition mask.
        """
        epsilon = _finite(epsilon, "epsilon")
        mask = _finite(mask, "mask")
        if epsilon < 0.0:
            raise ValueError("epsilon must be nonnegative")
        out = math.sqrt(epsilon) * self.amplitudes * mask
        if not np.isfinite(out).all():
            raise OverflowError("scaled wave amplitudes overflowed")
        return out

    def scaled_reconstructed_target(self, *, epsilon: float, mask: float) -> np.ndarray:
        """Covariance of the scaled signed amplitudes.

        Under exact arithmetic this equals ``epsilon * mask^2 * target``.
        """
        waves = self.scaled_wave_amplitudes(epsilon=epsilon, mask=mask)
        return self.matrix @ np.square(waves)

    @classmethod
    def from_paper_ratio(
        cls,
        *,
        c_star: float,
        u_star: float,
        scale_minus: float,
        scale_plus: float,
        m: float,
        t: float,
    ) -> "PositiveSignedStressDecomposition":
        """Instantiate the reference cone from the pinned manuscript ratio test.

        The strict ratio is checked independently before it is converted to the
        equivalent signed-matrix cone.  This is still only the reference model:
        the actual Section 7 pulse-integrated columns require separate analytic
        perturbation and derivative estimates.
        """
        margin = paper_ratio_margin(c_star=c_star, u_star=u_star, m=m, t=t)
        if margin <= 0.0:
            raise ValueError(
                "paper reference ratio requires |c_star*t/m| < "
                "u_star/sqrt(1+u_star^2)"
            )
        return cls(
            a=normal_magnitude(c_star, u_star),
            b=_positive(u_star, "u_star"),
            scale_minus=scale_minus,
            scale_plus=scale_plus,
            m=m,
            t=t,
        )
