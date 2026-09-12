"""Typed mixed-scale ``naturalRemainder(x0)`` pair for the actual schedule.

The angular and axial halves of the pinned natural-axis remainder are already
materialized separately at the contraction centre ``x0 = referencePair``.  The
theorem-selected scale is too wide for a faithful binary64 collapse:

* the angular half keeps ``slow1 / Lambda`` as an explicit Decimal ratio;
* the axial half keeps ``slow2 / Lambda`` likewise and the genuine ``a^2``
  pressure contribution in signed-log form.

This module closes only the *pair assembly* seam.  It verifies that both landed
halves come from the same actual SchedulePressure datum, ``j``, certified
``sigma``/``epsilon``, and theorem-selected ``Lambda``, then exposes them as one
typed value representing the complete pinned ``naturalRemainder(x0)``.  It does
not narrow any wide scale and deliberately does not yet apply the outer
fixed-point factor ``1/(2*Lambda)``.

No caller-supplied sigma, epsilon, Lambda, C, amplitude, coefficient table, or
cutoff is accepted by the production factory.  This remains
``formal-structure`` rather than a paper-exact fixed point or velocity field.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal

from .axis_coefficient_reference_state import ActualScheduleReferenceAxisState
from .axis_coefficient_wide_angular_remainder import (
    ActualScheduleReferenceWideAngularRemainderState,
    MixedScaleAngularCoefficientJet,
    actual_schedule_reference_wide_angular_remainder_state,
)
from .axis_coefficient_wide_axial_remainder import (
    ActualScheduleReferenceWideAxialRemainderState,
    MixedScaleAxialCoefficientJet,
    actual_schedule_reference_wide_axial_remainder_state,
)
from .outgoing_tail import TailData


@dataclass(frozen=True)
class ActualScheduleReferenceWideNaturalRemainderState:
    """Complete mixed-scale pair ``naturalRemainder(referencePair)``.

    The object is intentionally a structural pairing of the two independently
    landed theorem-grounded halves.  Its validation is strict so a branch from
    one SchedulePressure construction cannot be silently paired with another.
    """

    angular: ActualScheduleReferenceWideAngularRemainderState
    axial: ActualScheduleReferenceWideAxialRemainderState

    def __post_init__(self) -> None:
        if not isinstance(self.angular, ActualScheduleReferenceWideAngularRemainderState):
            raise TypeError("angular must be ActualScheduleReferenceWideAngularRemainderState")
        if not isinstance(self.axial, ActualScheduleReferenceWideAxialRemainderState):
            raise TypeError("axial must be ActualScheduleReferenceWideAxialRemainderState")

        a_ref = self.angular.reference
        u_ref = self.axial.reference
        if self.angular.epsilon != self.axial.epsilon:
            raise ValueError("angular/axial epsilon mismatch")
        if self.angular.Lambda != self.axial.Lambda:
            raise ValueError("angular/axial theorem-selected Lambda mismatch")
        if a_ref.reference.data != u_ref.reference.data:
            raise ValueError("angular/axial SchedulePressure TailData mismatch")
        if a_ref.reference.j != u_ref.reference.j:
            raise ValueError("angular/axial schedule j mismatch")
        if a_ref.reference.sigma != u_ref.reference.sigma:
            raise ValueError("angular/axial certified sigma mismatch")

        if not self.angular.mixed_scale_angular_remainder_complete:
            raise ValueError("angular mixed-scale remainder is incomplete")
        if not self.axial.mixed_scale_axial_remainder_complete:
            raise ValueError("axial mixed-scale remainder is incomplete")

    @property
    def reference(self) -> ActualScheduleReferenceAxisState:
        """Return the shared actual-schedule reference-state view."""

        return self.angular.reference

    @property
    def epsilon(self) -> float:
        return self.angular.epsilon

    @property
    def Lambda(self) -> Decimal:
        return self.angular.Lambda

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def mixed_scale_natural_remainder_x0_complete(self) -> bool:
        """The two pinned remainder branches are assembled, before Picard scaling."""

        return True

    @property
    def picard_x1_materialized(self) -> bool:
        """Fail closed until the outer ``1/(2*Lambda)`` scale is implemented."""

        return False

    def jet_pair(
        self,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[MixedScaleAngularCoefficientJet, MixedScaleAxialCoefficientJet]:
        """Return the angular/axial mixed-scale jets at one coefficient coordinate."""

        return self.angular.jet(n, m, eta), self.axial.jet(n, m, eta)


def actual_schedule_reference_wide_natural_remainder_state(
    data: TailData,
    j: float,
    *,
    phase_samples: int = 4001,
) -> ActualScheduleReferenceWideNaturalRemainderState:
    """Materialize the complete mixed-scale ``naturalRemainder(x0)`` pair.

    Both branches are independently rebuilt from the same caller-supplied
    schedule datum and ``j`` through their existing actual SchedulePressure
    theorem chains.  The returned type then checks exact compatibility of all
    theorem-selected scales before exposing the pair.
    """

    if not isinstance(data, TailData):
        raise TypeError("data must be TailData")

    angular = actual_schedule_reference_wide_angular_remainder_state(
        data,
        j,
        phase_samples=phase_samples,
    )
    axial = actual_schedule_reference_wide_axial_remainder_state(
        data,
        j,
        phase_samples=phase_samples,
    )
    return ActualScheduleReferenceWideNaturalRemainderState(
        angular=angular,
        axial=axial,
    )
