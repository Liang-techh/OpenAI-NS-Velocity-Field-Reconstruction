"""Uniform LocalBase numeric consequences for Section 7 phase estimates.

This module implements the *numeric implication* used between the pinned Lean
``PhaseEstimates.LocalBaseBounds`` hypotheses and
``PhaseEstimates.rounded_normal_estimates``.  It deliberately does not try to
manufacture the still-missing paper-exact base functions ``F,G,F0,G0`` or a
proof of convexity/differentiability from samples.

An upstream analytic/formal argument must certify the sup bounds recorded in
:class:`LocalBaseNumericEnvelope`.  Once those bounds are available, the
pinned theorem ``localBase_derivative_errors`` gives, uniformly for points in a
convex slow box of diameter ``d``,

``|F_R(q)-F0_R(q0)|, |G_R(q)-G0_R(q0)| <= M * (d + epsilon**2)``.

If ``d <= S**-3``, this is exactly the pair of reference-derivative hypotheses
required by ``rounded_normal_estimates``.  The adapter therefore closes a
previous interface gap without promoting the construction beyond
``formal-structure``.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .phase_estimates import PhaseScaleCertificate


def _finite_nonnegative(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


@dataclass(frozen=True)
class LocalBaseNumericEnvelope:
    """Certified numeric part of pinned ``LocalBaseBounds``.

    These values are *supremum certificates supplied by an upstream analytic
    or formal argument*, not estimates inferred from sampled field values.
    The qualitative Lean hypotheses (convex domain and differentiability) are
    intentionally outside this record and remain separate proof obligations.

    The field names mirror the quantitative fields of the pinned structure:
    second-derivative bounds for ``F0,G0``; a first-derivative bound for ``F0``;
    C0/C1 actual-to-reference errors; and the actual directional bounds used by
    ``rounded_normal_estimates``.
    """

    M: float
    epsilon: float
    second_F0_bound: float
    second_G0_bound: float
    first_F0_bound: float
    value_F_error_bound: float
    first_F_error_bound: float
    first_G_error_bound: float
    radial_F_bound: float
    axial_F_bound: float
    axial_G_bound: float

    def __post_init__(self) -> None:
        M = _finite_nonnegative(self.M, "M")
        eps = _finite_nonnegative(self.epsilon, "epsilon")
        values = {
            "second_F0_bound": self.second_F0_bound,
            "second_G0_bound": self.second_G0_bound,
            "first_F0_bound": self.first_F0_bound,
            "value_F_error_bound": self.value_F_error_bound,
            "first_F_error_bound": self.first_F_error_bound,
            "first_G_error_bound": self.first_G_error_bound,
            "radial_F_bound": self.radial_F_bound,
            "axial_F_bound": self.axial_F_bound,
            "axial_G_bound": self.axial_G_bound,
        }
        checked = {name: _finite_nonnegative(value, name) for name, value in values.items()}
        eps2_threshold = M * eps * eps
        direct_M = (
            "second_F0_bound",
            "second_G0_bound",
            "first_F0_bound",
            "radial_F_bound",
            "axial_F_bound",
            "axial_G_bound",
        )
        failed = [name for name in direct_M if checked[name] > M]
        failed += [
            name for name in ("value_F_error_bound", "first_F_error_bound", "first_G_error_bound")
            if checked[name] > eps2_threshold
        ]
        if failed:
            raise ValueError(
                "numeric envelope exceeds pinned LocalBaseBounds thresholds: " + ", ".join(failed)
            )

    @property
    def epsilon_squared_error_threshold(self) -> float:
        return float(self.M) * float(self.epsilon) ** 2

    def representative_derivative_error(self, diameter: float) -> float:
        """Pinned ``localBase_derivative_errors`` envelope ``M(d+epsilon^2)``."""
        diameter = _finite_nonnegative(diameter, "diameter")
        return float(self.M) * (diameter + float(self.epsilon) ** 2)

    def representative_value_error(self, diameter: float) -> float:
        """Pinned ``localBase_value_error`` envelope for ``|F(q)-F0(q0)|``."""
        diameter = _finite_nonnegative(diameter, "diameter")
        return float(self.M) * (diameter + float(self.epsilon) ** 2)


@dataclass(frozen=True)
class UniformLocalBaseBridge:
    """Turn certified LocalBase sup bounds into uniform rounded-normal inputs.

    ``diameter`` must already be certified for the convex slow-box domain.  The
    bridge checks ``diameter <= S^-3`` and exact consistency with the landed
    :class:`~openai_ns_reconstruction.phase_estimates.PhaseScaleCertificate`.
    It then exposes the uniform bounds needed for the ``FR/GR`` reference-error
    hypotheses of ``rounded_normal_estimates`` together with the actual
    directional bounds already present in ``LocalBaseBounds``.

    This object is not a certificate that the paper-exact base fields exist;
    it is valid only conditional on the upstream analytic/formal sup bounds.
    """

    local: LocalBaseNumericEnvelope
    scale: PhaseScaleCertificate
    diameter: float

    def __post_init__(self) -> None:
        diameter = _finite_nonnegative(self.diameter, "diameter")
        if float(self.local.M) != float(self.scale.M):
            raise ValueError("LocalBase M must match the phase-scale M exactly")
        if float(self.local.epsilon) != float(self.scale.epsilon):
            raise ValueError("LocalBase epsilon must match the phase-scale epsilon exactly")
        if diameter > self.max_slow_box_diameter:
            raise ValueError("slow-box diameter is not certified <= S^-3")

    @property
    def max_slow_box_diameter(self) -> float:
        return 1.0 / float(self.scale.S) ** 3

    @property
    def derived_reference_error(self) -> float:
        """Sharper ``M(d+epsilon^2)`` error retained from the actual diameter."""
        return self.local.representative_derivative_error(self.diameter)

    @property
    def theorem_reference_error(self) -> float:
        """The exact looser input ``M(S^-3+epsilon^2)`` used downstream."""
        return float(self.scale.M) * (
            self.max_slow_box_diameter + float(self.scale.epsilon) ** 2
        )

    def implication_checks(self) -> dict[str, bool]:
        """Check only theorem implication arithmetic; no field sampling occurs."""
        derived = self.derived_reference_error
        theorem = self.theorem_reference_error
        M = float(self.scale.M)
        return {
            "diameter_le_S_minus_3": float(self.diameter) <= self.max_slow_box_diameter,
            "FR_reference_error_implied": derived <= theorem,
            "GR_reference_error_implied": derived <= theorem,
            "FR_uniform_bound_available": float(self.local.radial_F_bound) <= M,
            "FZ_uniform_bound_available": float(self.local.axial_F_bound) <= M,
            "GZ_uniform_bound_available": float(self.local.axial_G_bound) <= M,
        }

    def rounded_normal_local_bounds(self) -> dict[str, float]:
        """Uniform LocalBase-origin bounds consumed by ``rounded_normal_estimates``.

        The returned values are bounds, not sampled point values.  Representative
        frequency, inverse-radius, slot-length and other hypotheses remain the
        responsibility of the existing phase/box adapters.
        """
        if not all(self.implication_checks().values()):
            raise ArithmeticError("LocalBase numeric implication unexpectedly failed")
        return {
            "FR_abs": float(self.scale.M),
            "FZ_abs": float(self.scale.M),
            "GZ_abs": float(self.scale.M),
            "FR_minus_FR0_abs": self.theorem_reference_error,
            "GR_minus_GR0_abs": self.theorem_reference_error,
        }
