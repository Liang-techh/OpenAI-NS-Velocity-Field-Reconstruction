"""Section 7 frame/damping implication bounds from the pinned Lean source.

This module covers the scalar implication layer after a uniform phase-normal
comparison has already been certified.  It transcribes the constants in
``BasePhaseGeometry.frame_errors_of_normal_close`` and
``BasePhaseGeometry.damping_error_of_normal_close`` from
``openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd``.

Nothing here manufactures the missing paper-exact Proposition 5.5 base fields
or proves the vector hypotheses of those theorems.  In particular, callers
must separately certify the actual normal/slot derivative bounds and actual
base-field comparison bounds.  These routines only turn such certified inputs
into the quantitative frame and viscosity-error envelopes used downstream.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .phase_estimates import phase_constant as phase_estimates_constant


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _finite_nonnegative(value: float, name: str) -> float:
    value = _finite(value, name)
    if value < 0.0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def damping_denominator(u: float) -> float:
    """Pinned ``BasePhaseGeometry.dampingDenominator``."""
    u = _finite(u, "u")
    radius_sq = 1.0 + u * u
    return radius_sq * math.sqrt(radius_sq)


def normal_lower(M: float, u: float) -> float:
    """Pinned ``BasePhaseGeometry.normalLower M u``.

    The formal theorem only uses this quantity after assuming ``M >= 1``.
    """
    M = _finite(M, "M")
    if M < 1.0:
        raise ValueError("normalLower specialization requires M>=1")
    return math.sqrt((1.0 / M) / (4.0 * damping_denominator(u)))


def base_phase_constant(M: float) -> float:
    """Pinned ``BasePhaseGeometry.phaseConstant M``.

    This is *not* the same constant as ``PhaseEstimates.phaseConstant``:
    ``normalConstant M = 8 * PhaseEstimates.phaseConstant (2*M)`` and the
    family-level name ``phaseConstant`` abbreviates that normal constant.
    """
    M = _finite(M, "M")
    if M < 1.0:
        raise ValueError("BasePhaseGeometry phase constant requires M>=1")
    return 8.0 * phase_estimates_constant(2.0 * M)


def frame_coordinate_error_bound(
    *, M: float, A: float, b: float, delta: float, eta: float,
) -> float:
    """Conclusion bound from ``frame_errors_of_normal_close``.

    Under the theorem's geometric/vector hypotheses, each of the three changing
    tangent-frame coefficient errors is bounded by

    ``16 * G^2 * (1+G) * E``

    with ``G=M+2+2A`` and ``E=eta+8(1+A)delta/b``.
    """
    M = _finite(M, "M")
    A = _finite_nonnegative(A, "A")
    b = _finite(b, "b")
    delta = _finite_nonnegative(delta, "delta")
    eta = _finite_nonnegative(eta, "eta")
    if M < 1.0:
        raise ValueError("frame implication requires M>=1")
    if b <= 0.0:
        raise ValueError("frame implication requires b>0")
    G = M + 2.0 + 2.0 * A
    E = eta + 8.0 * (1.0 + A) * delta / b
    out = 16.0 * G * G * (1.0 + G) * E
    if not math.isfinite(out):
        raise OverflowError("frame-coordinate error bound overflowed")
    return out


def damping_error_bound(*, M: float, A: float, delta: float) -> float:
    """Conclusion bound from ``damping_error_of_normal_close``.

    Provided the actual viscosity coefficient is in ``[0,4]`` and the theorem's
    normal-closeness hypotheses hold, the two-sided damping discrepancy is at
    most ``4*M*(2*A+5)*delta``.
    """
    M = _finite(M, "M")
    A = _finite_nonnegative(A, "A")
    delta = _finite_nonnegative(delta, "delta")
    if M < 0.0:
        raise ValueError("damping implication requires M>=0")
    out = 4.0 * M * (2.0 * A + 5.0) * delta
    if not math.isfinite(out):
        raise OverflowError("damping error bound overflowed")
    return out


def coordinate_constant(M: float, u: float) -> float:
    """Pinned family-level ``BasePhaseGeometry.coordinateConstant M u``."""
    M = _finite(M, "M")
    if M < 1.0:
        raise ValueError("coordinate constant requires M>=1")
    A = 3.0 * M
    G = M + 2.0 + 2.0 * A
    D = base_phase_constant(M)
    b = normal_lower(M, u)
    out = 16.0 * G * G * (1.0 + G) * (
        16.0 * M * M + 8.0 * (1.0 + A) * D / b
    )
    if not math.isfinite(out):
        raise OverflowError("coordinate constant overflowed")
    return out


def damping_constant(M: float) -> float:
    """Pinned family-level ``BasePhaseGeometry.dampingConstant M``."""
    M = _finite(M, "M")
    if M < 1.0:
        raise ValueError("damping constant requires M>=1")
    out = 4.0 * M * (6.0 * M + 5.0) * base_phase_constant(M)
    if not math.isfinite(out):
        raise OverflowError("damping constant overflowed")
    return out


@dataclass(frozen=True)
class FrameDampingEnvelope:
    """Fail-closed scalar gate for the Section 7 frame/damping implications.

    ``delta`` is a certified common upper bound for both the actual normal error
    and its slot derivative; ``eta`` is a certified bound for the actual/base
    coefficient comparison entering ``frame_errors_of_normal_close``.

    This class checks only scalar hypotheses that can be checked without the
    unresolved paper-exact background.  It does **not** prove the required
    vector identities, unit/orthogonality facts, differentiability, or the
    actual inequalities represented by ``delta`` and ``eta``.
    """

    M: float
    A: float
    b: float
    B: float
    delta: float
    eta: float
    viscosity: float

    def __post_init__(self) -> None:
        M = _finite(self.M, "M")
        A = _finite_nonnegative(self.A, "A")
        b = _finite(self.b, "b")
        B = _finite(self.B, "B")
        delta = _finite_nonnegative(self.delta, "delta")
        _finite_nonnegative(self.eta, "eta")
        viscosity = _finite(self.viscosity, "viscosity")
        if M < 1.0:
            raise ValueError("frame/damping envelope requires M>=1")
        if b <= 0.0 or B <= 0.0:
            raise ValueError("frame/damping envelope requires b>0 and B>0")
        if b > B:
            raise ValueError("theorem requires b<=B")
        if B > M:
            raise ValueError("damping theorem requires B<=M")
        if delta > B / 2.0:
            raise ValueError("normal error is not certified <=B/2")
        if not 0.0 <= viscosity <= 4.0:
            raise ValueError("damping theorem requires viscosity in [0,4]")
        # Force overflow detection at construction time rather than later.
        _ = frame_coordinate_error_bound(M=M, A=A, b=b, delta=delta, eta=self.eta)
        _ = damping_error_bound(M=M, A=A, delta=delta)

    @property
    def frame_error(self) -> float:
        return frame_coordinate_error_bound(
            M=self.M, A=self.A, b=self.b, delta=self.delta, eta=self.eta
        )

    @property
    def damping_error(self) -> float:
        return damping_error_bound(M=self.M, A=self.A, delta=self.delta)

    @classmethod
    def from_large_band(
        cls, *, M: float, S: float, u: float, B: float, viscosity: float,
    ) -> "FrameDampingEnvelope":
        """Instantiate the exact family-level specialization used in Lean.

        The pinned construction uses ``A=3M``, normal/normal-velocity error
        ``delta=phaseConstant(M)/S``, base comparison ``eta=16 M^2/S``, and
        ``b=normalLower(M,u)``.  The constructor still requires the actual
        reference normal scale ``B`` and actual viscosity coefficient so it can
        fail closed on the scalar hypotheses ``b<=B``, ``delta<=B/2``,
        ``B<=M`` and ``0<=viscosity<=4``.
        """
        M = _finite(M, "M")
        S = _finite(S, "S")
        if M < 1.0:
            raise ValueError("large-band specialization requires M>=1")
        if S < 1.0:
            raise ValueError("large-band specialization requires S>=1")
        D = base_phase_constant(M)
        return cls(
            M=M,
            A=3.0 * M,
            b=normal_lower(M, u),
            B=B,
            delta=D / S,
            eta=16.0 * M * M / S,
            viscosity=viscosity,
        )

    def specialized_constant_checks(self, *, u: float, S: float, rel_tol: float = 1e-12) -> dict[str, bool]:
        """Cross-check the family specialization against the named Lean constants."""
        S = _finite(S, "S")
        rel_tol = _finite_nonnegative(rel_tol, "rel_tol")
        if S <= 0.0:
            raise ValueError("S must be positive")
        return {
            "frame_equals_coordinateConstant_over_S": math.isclose(
                self.frame_error, coordinate_constant(self.M, u) / S,
                rel_tol=rel_tol, abs_tol=0.0,
            ),
            "damping_equals_dampingConstant_over_S": math.isclose(
                self.damping_error, damping_constant(self.M) / S,
                rel_tol=rel_tol, abs_tol=0.0,
            ),
        }
