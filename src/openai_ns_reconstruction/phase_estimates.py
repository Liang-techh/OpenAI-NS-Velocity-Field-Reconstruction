"""Scale-level Section 7 phase-error certificates from the pinned Lean source.

This module transcribes the quantitative scalar bounds proved in
``openai/NavierStokesAndEuler@f9e8bc5.../NavierStokes/PhaseEstimates.lean``.
It is deliberately independent of the still-missing paper-exact base fields:
it certifies only the *scale hypotheses* used by the rounded-normal estimates,
not the local C2/base-field hypotheses themselves.

In particular, a successful ``PhaseScaleCertificate`` must not be interpreted
as a paper-exact wave or as proof that Eqs. (7.9)-(7.11) hold on every active
box.  It is a reusable fail-closed gate for the scale part of that argument.
"""
from __future__ import annotations

from dataclasses import dataclass
import math


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def phase_error(S: float, epsilon: float, k: float) -> float:
    """Lean ``phaseError S ε k = 1/S + S ε² + S/k + ε S``."""
    S = _finite(S, "S")
    epsilon = _finite(epsilon, "epsilon")
    k = _finite(k, "k")
    if S <= 0.0 or epsilon < 0.0 or k <= 0.0:
        raise ValueError("phase-error parameters require S>0, epsilon>=0, k>0")
    return 1.0 / S + S * epsilon**2 + S / k + epsilon * S


def phase_constant(M: float) -> float:
    """Lean ``phaseConstant M = 8 M^3 + 2 M^4``."""
    M = _finite(M, "M")
    if M < 0.0:
        raise ValueError("M must be nonnegative")
    return 8.0 * M**3 + 2.0 * M**4


@dataclass(frozen=True)
class PhaseScaleCertificate:
    """Fail-closed executable form of ``phaseError_le_four_div`` hypotheses.

    The pinned Lean theorem proves ``phaseError S ε k <= 4/S`` when

    ``S > 0``, ``S^2 ε^2 <= 1``, ``S^2/k <= 1``, ``ε S^2 <= 1``.

    We additionally require ``M>=1, S>=1, k>=1`` because those are the scale
    assumptions used by the downstream ``rounded_normal_estimates`` theorem.
    This class does not manufacture the C2/base-field bounds needed there.
    """

    M: float
    S: float
    epsilon: float
    k: float

    def __post_init__(self) -> None:
        M = _finite(self.M, "M")
        S = _finite(self.S, "S")
        epsilon = _finite(self.epsilon, "epsilon")
        k = _finite(self.k, "k")
        if M < 1.0 or S < 1.0 or epsilon < 0.0 or k < 1.0:
            raise ValueError("rounded phase estimates require M>=1, S>=1, epsilon>=0, k>=1")
        checks = self.scale_hypotheses()
        failed = [name for name, ok in checks.items() if not ok]
        if failed:
            raise ValueError("uncertified phase scale hypotheses: " + ", ".join(failed))

    def scale_hypotheses(self) -> dict[str, bool]:
        S, eps, k = float(self.S), float(self.epsilon), float(self.k)
        return {
            "S2_epsilon2_le_1": S * S * eps * eps <= 1.0,
            "S2_over_k_le_1": S * S / k <= 1.0,
            "epsilon_S2_le_1": eps * S * S <= 1.0,
        }

    @property
    def error(self) -> float:
        return phase_error(self.S, self.epsilon, self.k)

    @property
    def four_over_S(self) -> float:
        return 4.0 / float(self.S)

    @property
    def rounded_normal_bound(self) -> float:
        """``phaseConstant(M) * phaseError`` from ``rounded_normal_estimates``."""
        return phase_constant(self.M) * self.error

    @property
    def simplified_normal_bound(self) -> float:
        """Consequence after ``phaseError_le_four_div``: ``4*C(M)/S``."""
        return phase_constant(self.M) * self.four_over_S

    def verify_scalar_inequality(self, *, slack: float = 1e-12) -> bool:
        """Numerically cross-check the proved scalar inequality without weakening it."""
        slack = _finite(slack, "slack")
        if slack < 0.0:
            raise ValueError("slack must be nonnegative")
        return self.error <= self.four_over_S * (1.0 + slack)
