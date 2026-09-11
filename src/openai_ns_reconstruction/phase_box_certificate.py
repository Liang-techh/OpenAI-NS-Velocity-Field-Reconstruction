"""Fail-closed representative/point adapter for Section 7 phase estimates.

This module bridges the existing local phase algebra to the hypotheses of the
pinned official theorem ``PhaseEstimates.rounded_normal_estimates`` without
pretending that the still-missing paper-exact base fields have been built.

Two design choices are important:

* The representative shear vector, orthogonal unit direction, and frequency are
  constructed from one representative base jet, so the Lean equalities
  ``g = [R0*FR0, GR0]``, ``K ⟂ g`` and
  ``[target/R0,pz] = representativeFrequency ...`` hold by construction rather
  than by a floating tolerance check.
* The rounded frequency here follows the pinned Lean implementation exactly:
  ``floor(k*target)`` with the zero floor replaced by ``1``.  This is kept
  separate from the paper-level deterministic nearest-integer convention in
  ``phase.py`` because the quantitative Lean theorem is proved for this
  particular punctured-lattice rounding.

A successful point certificate checks the numerical inequalities appearing in
``rounded_normal_estimates`` and then independently evaluates the normal and
slot-derivative errors.  It is *not* a uniform-box proof: an upstream caller
must still certify that the local C1/C2 bounds used here hold at every point of
an active slow box.  Consequently this module is ``formal-structure`` only.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from .phase_estimates import PhaseScaleCertificate, phase_constant, phase_error


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def lean_nonzero_round(x: float) -> int:
    """Pinned Lean ``nonzeroRound``: floor, except a zero floor becomes one."""
    x = _finite(x, "x")
    n = math.floor(x)
    return 1 if n == 0 else int(n)


def lean_rounded_frequency(k: float, target: float) -> float:
    """Pinned Lean ``roundedFrequency k target`` with ``k>0``."""
    k = _finite(k, "k")
    target = _finite(target, "target")
    if k <= 0.0:
        raise ValueError("k must be positive")
    return lean_nonzero_round(k * target) / k


@dataclass(frozen=True)
class RepresentativePhaseData:
    """Representative data satisfying the exact algebraic Lean equalities.

    ``orientation`` selects one of the two unit quarter-turns of the shear
    vector.  The constructor derives ``g``, ``K``, ``target`` and ``pz`` from
    the representative slow derivatives instead of accepting redundant values.
    """

    R0: float
    FR0: float
    GR0: float
    B: float
    sigma: int
    u: float
    L: float
    orientation: int = 1

    def __post_init__(self) -> None:
        for name in ("R0", "FR0", "GR0", "B", "u", "L"):
            _finite(getattr(self, name), name)
        if self.R0 == 0.0:
            raise ValueError("R0 must be nonzero")
        if self.B <= 0.0:
            raise ValueError("B must be positive")
        if self.L == 0.0:
            raise ValueError("L must be nonzero")
        if self.sigma not in (-1, 1):
            raise ValueError("sigma must be +1 or -1")
        if self.orientation not in (-1, 1):
            raise ValueError("orientation must be +1 or -1")
        if self.g_norm == 0.0:
            raise ValueError("representative shear g must be nonzero")

    @property
    def g(self) -> np.ndarray:
        return np.array([self.R0 * self.FR0, self.GR0], dtype=float)

    @property
    def g_norm(self) -> float:
        return float(np.linalg.norm(np.array([self.R0 * self.FR0, self.GR0], dtype=float)))

    @property
    def K(self) -> np.ndarray:
        g = self.g
        out = self.orientation * np.array([-g[1], g[0]], dtype=float) / self.g_norm
        return out

    @property
    def representative_frequency(self) -> np.ndarray:
        g = self.g
        coeff = self.sigma * self.u / (self.L * self.g_norm**2)
        return self.B * (self.K - coeff * g)

    @property
    def target(self) -> float:
        return self.R0 * float(self.representative_frequency[0])

    @property
    def pz(self) -> float:
        return float(self.representative_frequency[1])

    def algebraic_residuals(self) -> dict[str, float]:
        """Floating cross-checks of identities that are true by construction."""
        g = self.g
        K = self.K
        freq = self.representative_frequency
        return {
            "K_norm_minus_1": float(np.linalg.norm(K) - 1.0),
            "K_dot_g": float(K @ g),
            "frequency_theta": float(freq[0] - self.target / self.R0),
            "frequency_z": float(freq[1] - self.pz),
        }


@dataclass(frozen=True)
class RoundedNormalPointCertificate:
    """Executable hypotheses of pinned ``rounded_normal_estimates`` at one point.

    The supplied derivatives are actual point values.  The reference values
    ``FR0,GR0`` and frequency data come from :class:`RepresentativePhaseData`.
    Uniformity over an entire slow box is deliberately *not* inferred from
    samples; the caller must provide that analytic/local-base step separately.
    """

    representative: RepresentativePhaseData
    scale: PhaseScaleCertificate
    R: float
    v: float
    FR: float
    GR: float
    FZ: float
    GZ: float

    def __post_init__(self) -> None:
        for name in ("R", "v", "FR", "GR", "FZ", "GZ"):
            _finite(getattr(self, name), name)
        failed = [name for name, ok in self.hypotheses().items() if not ok]
        if failed:
            raise ValueError("rounded-normal point hypotheses failed: " + ", ".join(failed))

    @property
    def p(self) -> float:
        return lean_rounded_frequency(self.scale.k, self.representative.target)

    @property
    def reference_normal(self) -> np.ndarray:
        r = self.representative
        signed_slot = r.sigma * (r.u / 2.0 + r.u * self.v / r.L)
        return np.array([r.B * signed_slot, *(r.B * r.K)], dtype=float)

    @property
    def explicit_normal(self) -> np.ndarray:
        r = self.representative
        eps = float(self.scale.epsilon)
        return np.array([
            r.sigma * r.B * r.u / 2.0 - self.v * (self.p * self.FR + r.pz * self.GR),
            self.p / self.R,
            r.pz - eps * self.v * (self.p * self.FZ + r.pz * self.GZ),
        ], dtype=float)

    @property
    def normal_velocity(self) -> np.ndarray:
        r = self.representative
        eps = float(self.scale.epsilon)
        return np.array([
            -(self.p * self.FR + r.pz * self.GR),
            0.0,
            -eps * (self.p * self.FZ + r.pz * self.GZ),
        ], dtype=float)

    @property
    def theorem_bound(self) -> float:
        s = self.scale
        return phase_constant(s.M) * phase_error(s.S, s.epsilon, s.k)

    @property
    def normal_error(self) -> float:
        return float(np.linalg.norm(self.explicit_normal - self.reference_normal))

    @property
    def normal_velocity_norm(self) -> float:
        return float(np.linalg.norm(self.normal_velocity))

    def hypotheses(self) -> dict[str, bool]:
        r = self.representative
        s = self.scale
        M, S, eps, k = map(float, (s.M, s.S, s.epsilon, s.k))
        base = M * (1.0 / S**3 + eps**2)
        return {
            "R_nonzero": self.R != 0.0,
            "R0_nonzero": r.R0 != 0.0,
            "g_nonzero": r.g_norm > 0.0,
            "M_ge_1": M >= 1.0,
            "S_ge_1": S >= 1.0,
            "k_ge_1": k >= 1.0,
            "epsilon_nonnegative": eps >= 0.0,
            "target_bound": abs(r.target) <= M,
            "pz_bound": abs(r.pz) <= M,
            "B_bound": abs(r.B) <= M,
            "sigma_unit": abs(r.sigma) == 1,
            "u_bound": abs(r.u) <= M,
            "inverse_g_bound": 1.0 / r.g_norm <= M,
            "inverse_L_bound": 1.0 / abs(r.L) <= M / S,
            "slot_bound": abs(self.v) <= M * S,
            "inverse_R_bound": self.R != 0.0 and abs(1.0 / self.R) <= M,
            "inverse_R0_bound": abs(1.0 / r.R0) <= M,
            "radius_difference": abs(self.R - r.R0) <= 1.0 / S**3,
            "FR_bound": abs(self.FR) <= M,
            "FZ_bound": abs(self.FZ) <= M,
            "GZ_bound": abs(self.GZ) <= M,
            "FR_reference_error": abs(self.FR - r.FR0) <= base,
            "GR_reference_error": abs(self.GR - r.GR0) <= base,
        }

    def numeric_crosscheck(self, *, slack: float = 1e-12) -> dict[str, bool]:
        """Independent floating evaluation of the two theorem conclusions."""
        slack = _finite(slack, "slack")
        if slack < 0.0:
            raise ValueError("slack must be nonnegative")
        bound = self.theorem_bound * (1.0 + slack)
        return {
            "normal_error_le_bound": self.normal_error <= bound,
            "normal_velocity_le_bound": self.normal_velocity_norm <= bound,
        }

    def geometry_smallness(self) -> bool:
        """The extra ``actual_phase_geometry`` gate ``delta <= B/2``."""
        return self.theorem_bound <= self.representative.B / 2.0

    def geometry_crosscheck(self) -> dict[str, float | bool]:
        """Evaluate lower/frame quantities used after the rounded-normal bound.

        This is a numerical cross-check of the consequences in
        ``actual_phase_geometry``.  It refuses to run unless the theorem-bound
        smallness condition holds.
        """
        r = self.representative
        if not self.geometry_smallness():
            raise ValueError("phase error is not small enough for actual_phase_geometry")
        n = self.explicit_normal
        tail = n[1:]
        beta = float(np.linalg.norm(tail))
        if beta == 0.0:
            raise ArithmeticError("punctured-lattice normal unexpectedly vanished")
        direction = tail / beta
        signed_slot = r.sigma * (r.u / 2.0 + r.u * self.v / r.L)
        radial_slope = float(n[0] / beta)
        delta = self.theorem_bound
        return {
            "normal_scale": beta,
            "normal_norm": float(np.linalg.norm(n)),
            "normal_scale_ge_B_over_2": beta >= r.B / 2.0,
            "normal_norm_ge_B_over_2": float(np.linalg.norm(n)) >= r.B / 2.0,
            "radial_slope_error": abs(radial_slope - signed_slot),
            "radial_slope_bound": 2.0 * (1.0 + abs(signed_slot)) * delta / r.B,
            "direction_error": float(np.linalg.norm(direction - r.K)),
            "direction_bound": 4.0 * delta / r.B,
        }
