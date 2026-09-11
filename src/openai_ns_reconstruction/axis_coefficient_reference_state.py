"""Lazy coefficient-state view of the actual natural-axis reference pair.

The pinned ``AxisCoefficientSpace`` does not store unrelated samples.  Its raw
coordinate ``(n, m, eta)`` represents the *actual* ``m``-th parameter
derivative of radial coefficient ``n``, divided by the positive axis weight,
and compatibility is the adjacent-jet fundamental-theorem-of-calculus (FTC)
identity on ``NaturalAxisCoefficients.window = [-11/10, 11/10]``.

PRs #87 and #92 already materialized arbitrary finite actual eta-jets for the
angular and axial halves of the schedule-derived ``referencePair``.  This
module closes the representation seam between those two families: it exposes a
single pair of lazy coefficient states with the same window, weight
normalization, and unnormalized actual-jet semantics used by
``AxisCoefficientSpace.ofJetFamily``.

This is intentionally *not* a claim that the complete Banach-space backend has
landed.  In particular, this module does not certify the global supremum over
all radial/parameter indices, does not implement ``AxisOperators`` or
``coefficientOperators``, and does not evaluate ``naturalRemainder``.  It
therefore cannot produce a genuine Picard iterate and always keeps
``paper_exact`` false.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable
import math

from .axis_reference_angular_jets import (
    ActualScheduleAngularReferenceJets,
    actual_schedule_angular_reference_jets,
    axis_weight,
)
from .axis_reference_axial_jets import (
    ActualScheduleAxialReferenceJets,
    actual_schedule_axial_reference_jets,
)
from .axis_reference_pair import ActualScheduleReferencePair
from .outgoing_tail import TailData

WINDOW_LEFT = -11.0 / 10.0
WINDOW_RIGHT = 11.0 / 10.0

JetProvider = Callable[[int, int, float], float]


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _eta_in_window(value: float) -> float:
    eta = float(value)
    if not math.isfinite(eta) or not WINDOW_LEFT <= eta <= WINDOW_RIGHT:
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]")
    return eta


@dataclass(frozen=True)
class AxisCoefficientJetState:
    """Lazy actual-derivative coordinates in the pinned coefficient-space shape.

    ``jet(n, m, eta)`` is the unnormalized actual derivative, matching Lean's
    ``AxisCoefficientSpace.jet``.  ``normalized_coordinate`` divides by the
    pinned positive weight, matching the stored ``RawJets`` coordinate.

    Instances in the public construction below are backed only by the landed
    analytic derivative families.  The class deliberately exposes no norm or
    ``AxisSpace`` membership certificate; a future backend must add the global
    all-index weighted bound required by ``ofJetFamily``.
    """

    epsilon: float
    origin: str
    _jet_provider: JetProvider = field(repr=False, compare=False)

    def __post_init__(self) -> None:
        epsilon = float(self.epsilon)
        if not math.isfinite(epsilon) or epsilon <= 0.0:
            raise ValueError("epsilon must be finite and positive")
        if not isinstance(self.origin, str) or not self.origin.strip():
            raise ValueError("origin must be a nonempty provenance label")
        if not callable(self._jet_provider):
            raise TypeError("_jet_provider must be callable")
        object.__setattr__(self, "epsilon", epsilon)

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        """Fail-closed marker for the still-missing all-index supremum proof."""

        return False

    def jet(self, n: int, m: int, eta: float) -> float:
        """Return the unnormalized actual parameter derivative ``J[n,m](eta)``."""

        n = _index(n, "n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)
        value = float(self._jet_provider(n, m, eta))
        if not math.isfinite(value):
            raise ArithmeticError("coefficient jet must remain finite")
        return value

    def coefficient(self, n: int, eta: float) -> float:
        """Return radial coefficient ``n`` (the zeroth actual parameter jet)."""

        return self.jet(n, 0, eta)

    def normalized_coordinate(self, n: int, m: int, eta: float) -> float:
        """Return the pinned ``RawJets`` coordinate ``J[n,m]/weight(epsilon,n,m)``."""

        n = _index(n, "n")
        m = _index(m, "m")
        return self.jet(n, m, eta) / axis_weight(self.epsilon, n, m)


@dataclass(frozen=True)
class ActualScheduleReferenceAxisState:
    """Both actual-schedule halves of ``AxisContraction.referencePair``.

    The two component states share exactly the theorem-selected ``epsilon`` and
    the same schedule-derived reference data.  This object is a representation
    adapter only; it is not yet the complete ``AxisSpace × AxisSpace`` fixed
    point state used by the contraction theorem.
    """

    reference: ActualScheduleReferencePair
    angular_jets: ActualScheduleAngularReferenceJets = field(repr=False)
    axial_jets: ActualScheduleAxialReferenceJets = field(repr=False)
    phi: AxisCoefficientJetState
    u: AxisCoefficientJetState

    def __post_init__(self) -> None:
        if self.angular_jets.epsilon != self.axial_jets.epsilon:
            raise ValueError("angular and axial reference jets must use the same epsilon")
        if self.phi.epsilon != self.u.epsilon or self.phi.epsilon != self.angular_jets.epsilon:
            raise ValueError("reference coefficient states must use the same epsilon")
        if self.angular_jets.reference.sigma != self.axial_jets.reference.sigma:
            raise ValueError("angular and axial reference jets must use the same certified sigma")
        if self.reference.sigma != self.angular_jets.reference.sigma:
            raise ValueError("reference state must use the same certified sigma")
        if self.reference.data != self.angular_jets.reference.data or self.reference.data != self.axial_jets.reference.data:
            raise ValueError("reference state components must use the same TailData")
        if self.reference.j != self.angular_jets.reference.j or self.reference.j != self.axial_jets.reference.j:
            raise ValueError("reference state components must use the same j")

    @property
    def epsilon(self) -> float:
        return self.phi.epsilon

    @property
    def paper_exact(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    def jet_pair(self, n: int, m: int, eta: float) -> tuple[float, float]:
        """Return ``(phi0[n]^(m), u0[n]^(m))`` at ``eta``."""

        return self.phi.jet(n, m, eta), self.u.jet(n, m, eta)

    def coefficient_pair(self, n: int, eta: float) -> tuple[float, float]:
        """Return the zeroth-jet reference coefficient pair."""

        return self.phi.coefficient(n, eta), self.u.coefficient(n, eta)


def actual_schedule_reference_axis_state(
    data: TailData,
    j: float,
) -> ActualScheduleReferenceAxisState:
    """Assemble the landed actual angular/axial jet families into one state view.

    No caller-supplied pressure, sigma, rho, epsilon, or coefficient table is
    accepted.  Both families independently traverse the actual SchedulePressure
    certificate chain; their theorem-selected scales are checked for exact
    agreement before the state is returned.
    """

    angular = actual_schedule_angular_reference_jets(data, j)
    axial = actual_schedule_axial_reference_jets(data, j)
    reference = angular.reference
    phi = AxisCoefficientJetState(
        epsilon=angular.epsilon,
        origin="actual-schedule angular referencePair eta-jets",
        _jet_provider=angular.parameter_jet,
    )
    u = AxisCoefficientJetState(
        epsilon=axial.epsilon,
        origin="actual-schedule axial referencePair eta-jets",
        _jet_provider=axial.parameter_jet,
    )
    return ActualScheduleReferenceAxisState(
        reference=reference,
        angular_jets=angular,
        axial_jets=axial,
        phi=phi,
        u=u,
    )
