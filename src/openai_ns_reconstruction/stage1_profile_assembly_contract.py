"""Fail-closed Stage-1 identity contract for the fixed-point/profile handoff.

This module owns no coefficient-space operations. Agent/backend code that
materializes the genuine ``AxisCoefficientSpace`` fixed point can hand its
identity metadata to this contract before the resulting ``phi/u`` fields are
allowed to flow into :class:`natural_axis.NaturalProfileAssembly`.

The contract is intentionally strict: schedule parameters, theorem-selected
``sigma``/``epsilon``, wide ``Lambda``/``log C``, the coefficient window, and
the amplitude convention must agree *exactly*. No tolerance, sampled fit,
default coefficient, or binary64 projection is used to repair a mismatch.

Official structure pinned by ``references/provenance_manifest.json``:

* ``AxisContraction.lean`` / ``AxisCoefficientSpace.lean`` -- fixed-point state;
* ``NaturalProfile.lean`` -- rescaling from the solved state to E/U/Pi;
* ``SchedulePressure.lean`` -- actual outgoing pressure datum.

This is a downstream assembly guard only. It does not materialize the fixed
point, average, pressure correction, support/moment/cone certificates, or a
paper-exact velocity.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
import math
from typing import Final

from .axis_coefficient_reference_state import (
    WINDOW_LEFT,
    WINDOW_RIGHT,
    ActualScheduleReferenceAxisState,
    actual_schedule_reference_axis_state,
)
from .axis_fixed_point_picard import NaturalPicardContractionCertificate
from .natural_scale_selection_wide import WideNaturalScaleSelection
from .outgoing_tail import TailData
from .stage1_scale_chain_wide import diagnose_actual_schedule_scale_chain_wide


NATURAL_AMPLITUDE_CONVENTION: Final[str] = "exp(Lambda*realPhase)/C"


def _finite_positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _finite_decimal(value: Decimal, name: str, *, positive: bool) -> Decimal:
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError(f"{name} must be a finite Decimal")
    if positive and value <= 0:
        raise ValueError(f"{name} must be positive")
    if not positive and value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


@dataclass(frozen=True)
class Stage1ProfileBackendIdentity:
    """Identity metadata a genuine coefficient backend must preserve.

    ``schedule_key`` is ``(P, m, lam, wait, h)`` from the exact ``TailData``
    instance. Keeping those values in the handoff prevents a fixed point from
    a different admissible schedule from being silently combined with this
    profile's pressure datum or theorem scale.
    """

    schedule_key: tuple[float, float, float, float, float]
    j: float
    sigma: float
    epsilon: float
    Lambda: Decimal
    log_C: Decimal
    window_left: float
    window_right: float
    amplitude_convention: str = NATURAL_AMPLITUDE_CONVENTION

    def __post_init__(self) -> None:
        if len(self.schedule_key) != 5:
            raise ValueError("schedule_key must be (P,m,lam,wait,h)")
        key = tuple(float(value) for value in self.schedule_key)
        if not all(math.isfinite(value) for value in key):
            raise ValueError("schedule_key values must be finite")
        if key[0] <= 0.0 or key[1] <= 0.0 or not 0.0 < key[2] < 0.1:
            raise ValueError("schedule_key violates OutgoingCoreParameters constraints")
        if key[3] <= 25.0 or key[4] <= 0.0 or 2.0 * key[4] >= key[2]:
            raise ValueError("schedule_key violates TailData constraints")
        object.__setattr__(self, "schedule_key", key)
        object.__setattr__(self, "j", _finite_positive(self.j, "j"))
        object.__setattr__(self, "sigma", _finite_positive(self.sigma, "sigma"))
        object.__setattr__(self, "epsilon", _finite_positive(self.epsilon, "epsilon"))
        _finite_decimal(self.Lambda, "Lambda", positive=True)
        _finite_decimal(self.log_C, "log_C", positive=False)

        left = float(self.window_left)
        right = float(self.window_right)
        if not math.isfinite(left) or not math.isfinite(right) or not left < right:
            raise ValueError("coefficient window must be finite and ordered")
        object.__setattr__(self, "window_left", left)
        object.__setattr__(self, "window_right", right)

        if self.amplitude_convention != NATURAL_AMPLITUDE_CONVENTION:
            raise ValueError("unrecognized natural-profile amplitude convention")


@dataclass(frozen=True)
class ActualScheduleProfileAssemblyContract:
    """Exact identity handshake from current Stage-1 data to a future backend."""

    reference: ActualScheduleReferenceAxisState
    scale: WideNaturalScaleSelection
    picard: NaturalPicardContractionCertificate
    identity: Stage1ProfileBackendIdentity

    def __post_init__(self) -> None:
        if self.reference.epsilon != self.identity.epsilon:
            raise ValueError("reference/assembly epsilon mismatch")
        if self.reference.reference.sigma != self.identity.sigma:
            raise ValueError("reference/assembly sigma mismatch")
        if self.scale.Lambda != self.identity.Lambda:
            raise ValueError("scale/assembly Lambda mismatch")
        if self.scale.C.exponent_upper != self.identity.log_C:
            raise ValueError("scale/assembly log C mismatch")
        if self.picard.Lambda != self.identity.Lambda:
            raise ValueError("Picard/assembly Lambda mismatch")
        if (
            self.identity.window_left != WINDOW_LEFT
            or self.identity.window_right != WINDOW_RIGHT
        ):
            raise ValueError("assembly identity must use the pinned coefficient window")

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def full_reconstruction(self) -> bool:
        return False

    def require_backend_identity(self, candidate: Stage1ProfileBackendIdentity) -> None:
        """Reject any backend state produced from a different theorem identity."""

        if not isinstance(candidate, Stage1ProfileBackendIdentity):
            raise TypeError("candidate must be Stage1ProfileBackendIdentity")

        names = (
            "schedule_key",
            "j",
            "sigma",
            "epsilon",
            "Lambda",
            "log_C",
            "window_left",
            "window_right",
            "amplitude_convention",
        )
        mismatches = [
            name
            for name in names
            if getattr(candidate, name) != getattr(self.identity, name)
        ]
        if mismatches:
            raise ValueError(
                "fixed-point/profile assembly identity mismatch: " + ", ".join(mismatches)
            )


def actual_schedule_profile_assembly_contract(
    data: TailData,
    j: float,
) -> ActualScheduleProfileAssemblyContract:
    """Construct the downstream handoff contract from actual SchedulePressure data.

    The factory accepts no independent ``sigma``, ``epsilon``, ``Lambda``,
    normalization constant, pressure datum, or coefficient window. It
    reconstructs those identities from the same current-main providers that
    own the reference state and wide contraction scale.
    """

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")

    reference = actual_schedule_reference_axis_state(data, j)
    diagnostic = diagnose_actual_schedule_scale_chain_wide(data, j)
    if diagnostic.upstream_obstruction is not None:
        raise RuntimeError(diagnostic.upstream_obstruction)
    if diagnostic.scale is None:
        raise RuntimeError("actual-schedule chain did not produce the pinned wide scale")

    scale = diagnostic.scale
    picard = NaturalPicardContractionCertificate.from_scale(scale)
    if picard.Lambda != scale.Lambda:
        raise RuntimeError("wide scale and Picard certificate disagree on Lambda")

    identity = Stage1ProfileBackendIdentity(
        schedule_key=(
            data.core.P,
            data.core.m,
            data.core.lam,
            data.core.wait,
            data.h,
        ),
        j=float(j),
        sigma=reference.reference.sigma,
        epsilon=reference.epsilon,
        Lambda=scale.Lambda,
        log_C=scale.C.exponent_upper,
        window_left=WINDOW_LEFT,
        window_right=WINDOW_RIGHT,
    )
    return ActualScheduleProfileAssemblyContract(
        reference=reference,
        scale=scale,
        picard=picard,
        identity=identity,
    )