"""Compose the exact large-band phase gate with Section 7 frame/damping bounds.

``phase_large_band_scale`` certifies the asymptotic scalar hypotheses behind
``PhaseEstimates.phaseError_le_four_div`` without ever forming the underflow-
prone chart scale ``Q=2^-ell``.  ``phase_frame_bounds`` separately transcribes
the scalar consequences of the pinned ``BasePhaseGeometry`` frame and damping
theorems.  This module connects those two landed pieces at one integer band.

The bridge is intentionally fail-closed and strictly ``formal-structure``.  It
checks that the same ``ell``/``S_*=ell^2`` has crossed both the exact phase-
scale gate and the family-level scalar smallness requirements
``b<=B``, ``delta<=B/2``, ``B<=M`` and ``0<=nu<=4``.  It then exposes the
corresponding normal, frame and damping envelopes used in the route to Eqs.
(7.9)-(7.11).

Nothing here proves the still-missing vector/local-base hypotheses on an active
slow box, constructs the Proposition 5.5 background, or manufactures a wave.
In particular, the smaller ``4*PhaseEstimates.phaseConstant(M)/S_*`` bound and
the larger family-level ``BasePhaseGeometry.phaseConstant(M)/S_*`` delta are
kept as distinct quantities rather than silently identifying two different
Lean constants.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math

from .phase_frame_bounds import FrameDampingEnvelope
from .phase_large_band_scale import LargeBandPhaseScaleCertificate


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be finite")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


@dataclass(frozen=True)
class LargeBandPhaseFrameCertificate:
    """One-band scalar admission gate toward uniform Eqs. (7.9)-(7.11).

    ``scale`` must already be a successful
    :class:`LargeBandPhaseScaleCertificate`, so the exact integer phase-scale
    hypotheses are inherited.  The remaining inputs are theorem-side scalar
    data used by ``BasePhaseGeometry``; they are not fitted from sampled fields.

    Passing this object means only that the scalar implication chain is ready
    at this band.  Actual LocalBase/vector hypotheses and paper-exact base
    fields remain external obligations.
    """

    scale: LargeBandPhaseScaleCertificate
    u: float
    B: float
    viscosity: float
    envelope: FrameDampingEnvelope = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.scale, LargeBandPhaseScaleCertificate):
            raise TypeError("scale must be LargeBandPhaseScaleCertificate")
        u = _finite(self.u, "u")
        B = _finite(self.B, "B")
        viscosity = _finite(self.viscosity, "viscosity")
        object.__setattr__(self, "u", u)
        object.__setattr__(self, "B", B)
        object.__setattr__(self, "viscosity", viscosity)

        try:
            S = float(self.scale.S_star)
        except OverflowError as exc:
            raise OverflowError("S_*=ell^2 is not representable for frame arithmetic") from exc
        if not math.isfinite(S):
            raise OverflowError("S_*=ell^2 is not representable for frame arithmetic")

        envelope = FrameDampingEnvelope.from_large_band(
            M=self.scale.M,
            S=S,
            u=u,
            B=B,
            viscosity=viscosity,
        )
        object.__setattr__(self, "envelope", envelope)

    @property
    def ell(self) -> int:
        return self.scale.ell

    @property
    def S_star(self) -> int:
        return self.scale.S_star

    @property
    def local_phase_normal_bound(self) -> float:
        """Pinned rounded-normal consequence ``4*C_phase(M)/S_*``.

        This is the consequence of the exact large-band scale gate alone; it
        does not assert the missing LocalBase hypotheses needed to apply the
        rounded-normal theorem on an actual box.
        """

        return self.scale.simplified_rounded_normal_bound

    @property
    def family_normal_delta(self) -> float:
        """Family-level ``BasePhaseGeometry.phaseConstant(M)/S_*`` input."""

        return self.envelope.delta

    @property
    def frame_error_bound(self) -> float:
        return self.envelope.frame_error

    @property
    def damping_error_bound(self) -> float:
        return self.envelope.damping_error

    def scalar_checks(self) -> dict[str, bool]:
        """Expose the scalar gates without claiming any vector-field proof."""

        return {
            "exact_large_band_phase_scale": all(self.scale.scale_hypotheses().values()),
            "normal_lower_le_B": self.envelope.b <= self.envelope.B,
            "family_normal_delta_le_B_over_2": self.envelope.delta <= self.envelope.B / 2.0,
            "B_le_M": self.envelope.B <= self.envelope.M,
            "viscosity_in_0_4": 0.0 <= self.envelope.viscosity <= 4.0,
        }

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_base_fields_verified(self) -> bool:
        return False

    @property
    def local_base_vector_hypotheses_verified(self) -> bool:
        return False

    @property
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
