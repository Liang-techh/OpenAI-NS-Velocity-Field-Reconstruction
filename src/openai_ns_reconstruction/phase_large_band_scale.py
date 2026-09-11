"""Exact large-band scale gate for the Section 7 phase estimates.

The existing :mod:`charts` implementation intentionally evaluates ``Q=2^-ell``
in binary64 and therefore stops before deep asymptotic bands.  The pinned Lean
phase estimate, however, is an eventual statement in the integer band index and
does not require evaluating ``Q`` numerically.  This module keeps those two
roles separate.

For the manuscript choices

    S = ell^2,
    epsilon = Q^h = 2^(-ell*h),
    k = ceil(epsilon^(-1/2)),

the pinned theorem ``PhaseEstimates.phaseError_le_four_div`` needs

    S^2 epsilon^2 <= 1,
    S^2 / k <= 1,
    epsilon S^2 <= 1.

When ``h=a/b`` is positive rational, the first and third conditions are exactly
implied by the integer inequalities

    ell^(4b) <= 2^(2 ell a),
    ell^(4b) <= 2^(ell a).

For the carrier condition, ``k >= epsilon^(-1/2)`` gives the rigorous
sufficient condition

    ell^4 <= 2^(ell h / 2),

equivalently ``ell^(8b) <= 2^(ell a)``.  The latter is sufficient rather than
necessary because the ceiling in ``k`` can only make ``S^2/k`` smaller.

All comparisons below use Python integers; neither ``Q`` nor ``epsilon`` is
formed as binary64.  Passing this gate certifies only the scalar large-band
hypotheses and the consequent ``phaseError <= 4/S`` envelope.  It does not
materialize Proposition 5.5 base fields, prove LocalBase/vector hypotheses on a
slow box, or construct an oscillatory wave.  The status is therefore strictly
``formal-structure`` and ``paper_exact_velocity_available`` remains false.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from numbers import Integral, Rational, Real

from .phase_estimates import phase_constant


def _positive_band_index(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("ell must be a positive integer")
    return int(value)


def _rational_h(value: object) -> Fraction:
    if isinstance(value, bool):
        raise ValueError("h must be a finite rationalizable real with 0<h<1/2")
    if isinstance(value, Fraction):
        out = value
    elif isinstance(value, Integral):
        out = Fraction(int(value), 1)
    elif isinstance(value, Rational):
        out = Fraction(value)
    elif isinstance(value, Real):
        x = float(value)
        if not math.isfinite(x):
            raise ValueError("h must be finite")
        # Decimal-string conversion intentionally recovers manuscript values
        # such as 0.005 as exactly 1/200 rather than their binary64 fraction.
        out = Fraction(str(x))
    else:
        raise TypeError("h must be a rationalizable real number")
    if not Fraction(0, 1) < out < Fraction(1, 2):
        raise ValueError("h must satisfy 0<h<1/2")
    return out


def _finite_M(value: object) -> float:
    if isinstance(value, bool):
        raise ValueError("M must be finite with M>=1")
    out = float(value)
    if not math.isfinite(out) or out < 1.0:
        raise ValueError("M must be finite with M>=1")
    return out


def exact_large_band_hypotheses(
    ell: int,
    h: Fraction | float | int,
) -> dict[str, bool]:
    """Return exact integer checks implying the pinned phase-scale hypotheses.

    The middle entry is deliberately named ``*_via_carrier_lower_bound``:
    unlike the other two integer comparisons it is a sufficient implication
    using ``ceil(x)>=x``, not an equivalence with the ceiling-rounded condition.
    """

    n = _positive_band_index(ell)
    hq = _rational_h(h)
    a, b = hq.numerator, hq.denominator

    return {
        "S2_epsilon2_le_1": pow(n, 4 * b) <= (1 << (2 * n * a)),
        "S2_over_carrier_le_1_via_carrier_lower_bound": (
            pow(n, 8 * b) <= (1 << (n * a))
        ),
        "epsilon_S2_le_1": pow(n, 4 * b) <= (1 << (n * a)),
    }


@dataclass(frozen=True)
class LargeBandPhaseScaleCertificate:
    """Fail-closed scalar certificate for the actual dyadic large-band scale.

    This intentionally does not construct a :class:`charts.DyadicChart`; deep
    bands may have ``Q=2^-ell`` below binary64.  The exact integer band index and
    rational ``h`` are enough for the scalar theorem gate.
    """

    ell: int
    h: Fraction | float | int
    M: float

    def __post_init__(self) -> None:
        ell = _positive_band_index(self.ell)
        h = _rational_h(self.h)
        M = _finite_M(self.M)
        object.__setattr__(self, "ell", ell)
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "M", M)

        checks = exact_large_band_hypotheses(ell, h)
        failed = [name for name, ok in checks.items() if not ok]
        if failed:
            raise ValueError(
                "uncertified large-band phase hypotheses: " + ", ".join(failed)
            )

    @property
    def h_fraction(self) -> Fraction:
        return self.h if isinstance(self.h, Fraction) else _rational_h(self.h)

    @property
    def S_star(self) -> int:
        return self.ell * self.ell

    @property
    def epsilon_log2(self) -> Fraction:
        """Exact ``log2(epsilon)=-ell*h`` without forming epsilon."""

        return -self.ell * self.h_fraction

    @property
    def inverse_sqrt_epsilon_log2(self) -> Fraction:
        """Exact ``log2(epsilon^(-1/2))=ell*h/2``."""

        return self.ell * self.h_fraction / 2

    def scale_hypotheses(self) -> dict[str, bool]:
        return exact_large_band_hypotheses(self.ell, self.h_fraction)

    @property
    def simplified_phase_error_bound(self) -> float:
        """Pinned consequence ``phaseError <= 4/S`` after this scalar gate."""

        return 4.0 / float(self.S_star)

    @property
    def simplified_rounded_normal_bound(self) -> float:
        """``phaseConstant(M) * 4/S``; base-field hypotheses remain external."""

        return phase_constant(self.M) * self.simplified_phase_error_bound

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_base_fields_verified(self) -> bool:
        return False

    @property
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
