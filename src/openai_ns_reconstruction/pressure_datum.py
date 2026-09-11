"""Executable ideal-prefix pressure datum from the pinned Lean pressure layer.

This module closes one *specific* Stage-1 existential input without inventing a
profile.  ``NaturalAxisRange.lean`` accepts any ``PressureDatum.Admissible``
weight/exponent pair satisfying the ideal prefix

    g(y) = B^2 exp(y/5),  a(y) = 1  for y <= 0,  B >= 2.

A minimal admissible extension is therefore available explicitly: keep that
prefix, set the weight to zero for ``y > 0``, and keep ``a = 1`` everywhere.
This is not the manuscript's later ``SchedulePressure.axisPressure`` (whose
positive-y tail is built from the outgoing schedule); it is a genuine witness
for the abstract NaturalAxisRange pressure hypotheses and is kept labelled as
such.

For this witness the pressure integral is elementary,

    P(eta) = -(5/2) B^2 (1 + eta^2)^(-2),

so the pressure and its derivative can be evaluated without numerical
quadrature.  Choosing ``B = 2`` is the smallest amplitude explicitly allowed by
the pinned ``pressureData_of_ideal_prefix`` / ``ideal_prefix_*`` theorems.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .natural_axis import A, H, axis_U, d


@dataclass(frozen=True)
class IdealPrefixPressureDatum:
    """Concrete ``PressureDatum.Admissible`` ideal-prefix witness.

    ``B`` is restricted to the exact theorem-side range ``B >= 2``.  The
    positive-clock continuation is the zero weight, which is sufficient for
    the abstract pressure-datum theorem but deliberately does *not* claim to
    reproduce the later outgoing-schedule pressure.
    """

    B: float = 2.0

    def __post_init__(self) -> None:
        B = float(self.B)
        if not math.isfinite(B) or B < 2.0:
            raise ValueError("B must be finite and satisfy the Lean hypothesis B >= 2")
        object.__setattr__(self, "B", B)

    @classmethod
    def minimal(cls) -> "IdealPrefixPressureDatum":
        """Return the smallest amplitude allowed by the pinned theorem."""

        return cls(B=2.0)

    @property
    def cap(self) -> float:
        """Exponent cap ``A=1`` for the admissibility record."""

        return 1.0

    def weight(self, y: float) -> float:
        """Nonnegative integrable clock weight ``g``.

        The value at ``y=0`` belongs to the required ideal prefix.  A jump just
        to the right of zero is harmless: ``PressureDatum.Admissible`` requires
        integrability/nonnegativity of ``g``, not smoothness of ``g``.
        """

        y = float(y)
        if not math.isfinite(y):
            raise ValueError("y must be finite")
        if y <= 0.0:
            return self.B * self.B * math.exp(y / 5.0)
        return 0.0

    def exponent(self, y: float) -> float:
        """Measurable exponent ``a=1`` with ``0 <= a <= cap``."""

        y = float(y)
        if not math.isfinite(y):
            raise ValueError("y must be finite")
        return 1.0

    @property
    def total_mass(self) -> float:
        """Exact integral ``integral g = 5 B^2``."""

        return 5.0 * self.B * self.B

    @property
    def positive_exponent_mass(self) -> float:
        """Exact integral ``integral g*a = 5 B^2 > 0``."""

        return self.total_mass

    def kernel(self, eta: float) -> float:
        """``PressureDatum.kernel 1 eta = (1+eta^2)^(-2)``."""

        eta = float(eta)
        if not math.isfinite(eta):
            raise ValueError("eta must be finite")
        base = 1.0 + eta * eta
        return 1.0 / (base * base)

    def pressure(self, eta: float) -> float:
        """Closed-form pressure integral for this admissible witness."""

        return -0.5 * self.total_mass * self.kernel(eta)

    def pressure_derivative(self, eta: float) -> float:
        """Exact derivative ``10 B^2 eta / (1+eta^2)^3``."""

        eta = float(eta)
        if not math.isfinite(eta):
            raise ValueError("eta must be finite")
        base = 1.0 + eta * eta
        return 10.0 * self.B * self.B * eta / (base * base * base)

    @property
    def least_negative_pressure_on_unit_interval(self) -> float:
        """Exact maximum of ``P`` on ``[-1,1]``: ``-5 B^2 / 8``.

        Since ``B >= 2``, this is at most ``-5/2 < -1`` and therefore
        discharges the negativity hypothesis used by ``NaturalAxisData``.
        """

        return -5.0 * self.B * self.B / 8.0

    def natural_axis_Z(self, h: float, j: float, eta: float) -> float:
        """Evaluate the exact ``NaturalAxisData.Z`` formula for this datum.

        This is an executable cross-check of the pressure input.  A numerical
        value of ``Z`` is not promoted to a formal uniform certificate.
        """

        eta = float(eta)
        U = axis_U(j, eta)
        P = self.pressure(eta)
        dP = self.pressure_derivative(eta)
        return (
            -A(h) * (1.0 - 2.0 * eta * U) * U
            - 4.0 * H(h, j, eta)
            - d(eta) * dP
            + 4.0 * A(h) * eta * P
        )

    def pressure_data_point(self, eta: float) -> tuple[float, float]:
        """Return ``(P,P')`` after checking the exact unit-interval signs.

        The checks mirror the real hypotheses consumed by the natural-axis
        layer.  They are guards around closed-form identities, not a sampled
        proof of any later compact-set minimum.
        """

        eta = float(eta)
        if not math.isfinite(eta) or not -1.0 <= eta <= 1.0:
            raise ValueError("eta must be finite and lie in [-1,1]")
        P = self.pressure(eta)
        dP = self.pressure_derivative(eta)
        if P > -1.0:
            raise ArithmeticError("closed-form pressure lost the required P <= -1 bound")
        if eta < 0.0 and dP > 0.0:
            raise ArithmeticError("closed-form pressure derivative lost the eta<0 sign")
        if eta > 0.0 and dP < 0.0:
            raise ArithmeticError("closed-form pressure derivative lost the eta>0 sign")
        return P, dP
