"""Theorem-shaped coefficient evaluator for the natural-axis reference pair.

This module materializes the *zeroth parameter jets* of the actual reference
pair ``x0 = referencePair O d S`` used by the pinned natural-axis contraction.
It is deliberately narrower than a complete ``AxisCoefficientSpace.AxisSpace``
backend: parameter-derivative jets, generic coefficient operators, and the
nonlinear ``naturalRemainder`` still have to be implemented before a real
Picard step can be executed.

For the angular component, ``AxisReference.reference_coefficient`` proves that
when ``S`` is the constructed natural resolvent,

    phi0[n](eta) = (-chi(eta)/2)^n / (n! (n+1)!).

For the axial component, ``AxisContraction.referencePair`` gives

    u0 = -(1/2) J_1(inverseL * zStar),

and ``AxisWeightEstimates.regularInverseJet`` gives a zero constant coefficient
and shifts a radial coefficient by one degree with divisor
``radialDivisor(r,n)=(n+1)(n+r)``.  Because ``inverseL`` and ``zStar`` are
radially constant fixed data, the axial reference therefore has exactly one
nonzero radial coefficient:

    u0[1](eta) = -(1/2) inverseL(eta) zStar(eta).

The factory :func:`actual_schedule_reference_pair` obtains ``sigma`` from the
already-landed analytic low-|Z| margin for the actual ``SchedulePressure``
datum.  It does not accept an independent cutoff scale or pressure function.
The pressure values themselves retain the documented numerical-quadrature
caveat of :mod:`schedule_axis_pressure`.

Nothing in this module is promoted to ``paper-exact``.  In particular, an
actual coefficient-space fixed point ``(phi,u)`` is not yet materialized.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .natural_axis import L, chi
from .natural_axis_range import CutoffParametersFromMargin
from .outgoing_tail import TailData
from .schedule_axis_margin import (
    ScheduleLowZMarginWitness,
    certify_schedule_low_Z_margin,
    schedule_natural_axis_Z,
)

_WINDOW_LEFT = -11.0 / 10.0
_WINDOW_RIGHT = 11.0 / 10.0


def _radial_index(value: int, name: str = "n") -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _eta_in_window(eta: float) -> float:
    eta = float(eta)
    if not math.isfinite(eta) or not _WINDOW_LEFT <= eta <= _WINDOW_RIGHT:
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]")
    return eta


def radial_divisor(r: int, n: int) -> float:
    """Pinned ``AxisWeightEstimates.radialDivisor r n``.

    Only positive ``r`` occurs in the reference-pair construction.  Keeping
    this tiny definition executable avoids hiding the exact regular-inverse
    normalization behind a numerical fit.
    """

    r = _radial_index(r, "r")
    n = _radial_index(n)
    if r < 1:
        raise ValueError("r must be a positive integer")
    return float((n + 1) * (n + r))


@dataclass(frozen=True)
class ActualScheduleReferencePair:
    """Coefficient-level realization of the pinned ``referencePair``.

    ``margin_witness`` fixes the theorem-side ``sigma=sqrt(m)/20`` from the
    actual schedule.  The object exposes coefficient functions, not arbitrary
    independent coefficient arrays; consequently the radial recurrence is the
    one proved for the official resolvent/reference pair.

    This remains **formal-structure** because only the zeroth parameter jets
    are executable here and no nonlinear Picard iterate is claimed.
    """

    data: TailData
    j: float
    margin_witness: ScheduleLowZMarginWitness

    def __post_init__(self) -> None:
        if not isinstance(self.data, TailData):
            raise TypeError("data must be TailData")
        j = float(self.j)
        if not math.isfinite(j):
            raise ValueError("j must be finite")
        p = self.margin_witness.parameters
        if p.h != self.data.h or p.j != j:
            raise ValueError("margin witness must belong to the same schedule h and j")
        sigma = float(self.margin_witness.cutoff.sigma)
        if not math.isfinite(sigma) or sigma <= 0.0:
            raise ValueError("certified sigma must be finite and positive")
        object.__setattr__(self, "j", j)

    @property
    def cutoff(self) -> CutoffParametersFromMargin:
        return self.margin_witness.cutoff

    @property
    def sigma(self) -> float:
        return float(self.cutoff.sigma)

    @property
    def paper_exact(self) -> bool:
        """Fail-closed publication gate for this partial construction."""

        return False

    def chi0(self, eta: float) -> float:
        eta = _eta_in_window(eta)
        value = chi(self.data.h, self.j, self.sigma, eta)
        if not math.isfinite(value):
            raise ArithmeticError("reference chi coefficient must remain finite")
        return value

    def inverse_l(self, eta: float) -> float:
        eta = _eta_in_window(eta)
        denominator = L(self.data.h, eta)
        if not math.isfinite(denominator) or denominator == 0.0:
            raise ArithmeticError("NaturalAxisData.L vanished on the pinned coefficient window")
        value = 1.0 / denominator
        if not math.isfinite(value):
            raise ArithmeticError("inverseL coefficient must remain finite")
        return value

    def z_star(self, eta: float) -> float:
        eta = _eta_in_window(eta)
        return schedule_natural_axis_Z(self.data, self.j, eta)

    def phi_coefficient(self, n: int, eta: float) -> float:
        """Return the actual reference angular radial coefficient ``phi0[n]``.

        The recurrence is exactly ``AxisReference.reference_coefficient_succ``.
        Iterating it is numerically preferable to converting large factorials
        to binary64 and is algebraically identical to the proved closed form.
        """

        n = _radial_index(n)
        c = -0.5 * self.chi0(eta)
        value = 1.0
        for k in range(n):
            value *= c / radial_divisor(2, k)
        if not math.isfinite(value):
            raise ArithmeticError("reference angular coefficient must remain finite")
        return value

    def u_coefficient(self, n: int, eta: float) -> float:
        """Return the actual reference axial radial coefficient ``u0[n]``.

        ``inverseL*zStar`` is radially constant.  The pinned regular inverse
        ``J_1`` therefore puts it only in degree one, with divisor
        ``radialDivisor(1,0)=1``; ``referencePair`` contributes ``-1/2``.
        """

        n = _radial_index(n)
        eta = _eta_in_window(eta)
        if n != 1:
            return 0.0
        value = -0.5 * self.inverse_l(eta) * self.z_star(eta) / radial_divisor(1, 0)
        if not math.isfinite(value):
            raise ArithmeticError("reference axial coefficient must remain finite")
        return value

    def coefficient_pair(self, n: int, eta: float) -> tuple[float, float]:
        """Return ``(phi0[n](eta), u0[n](eta))`` for the pinned reference pair."""

        return self.phi_coefficient(n, eta), self.u_coefficient(n, eta)


def actual_schedule_reference_pair(data: TailData, j: float) -> ActualScheduleReferencePair:
    """Construct ``referencePair`` coefficient data from the actual schedule.

    No caller-supplied ``sigma``, pressure, ``chi`` or ``zStar`` is accepted.
    The already-landed no-sampling low-|Z| witness determines ``sigma`` and the
    actual SchedulePressure evaluator determines ``zStar``.
    """

    witness = certify_schedule_low_Z_margin(data, j)
    return ActualScheduleReferencePair(data=data, j=float(j), margin_witness=witness)
