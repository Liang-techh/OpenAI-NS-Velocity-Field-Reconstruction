"""Analytic second-jet bridge for the Lemma 5.2 compact moment repair.

The landed function-level repair materializes repaired ``U_n`` and ``E_n``
values from Eq. (5.14), while the strict-lower-history PositiveAxis source needs
second ``(X, eta)`` jets.  This module differentiates the *actual compact bump
corrections* analytically and combines them with caller-supplied second jets of
the unrepaired coefficient.

No finite-difference derivative is used in the production path.  The unrepaired
second jets, moment jets, and patch-factor jets remain explicit inputs, so this
is solver infrastructure / ``formal-structure`` only.  In particular it does
not identify ``E_n`` with the PositiveAxis ``phi_n`` variable, construct the
regular flux ``V_n/X``, or promote any coefficient to paper-exact status.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import math

from .background_lower_history_source import ProfileSecondJet
from .background_moment_repair import CompactMomentBump
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter

BaseSecondJetProvider = Callable[[float, float], ProfileSecondJet]


def _bump_x_second_jet(
    bump: CompactMomentBump,
    X: float,
) -> tuple[float, float, float]:
    """Return ``(b, d_X b, d_X^2 b)`` for ``b(sqrt(2X))`` analytically.

    The Lemma-5.2 bumps have support bounded away from the axis.  Hence the
    ``R=sqrt(2X)`` chain rule is nonsingular wherever the bump is nonzero, while
    outside (and at the flat support boundary) all three values are exactly 0.
    """

    X = float(X)
    if not math.isfinite(X) or X < 0.0:
        raise ValueError("X must be finite and nonnegative")
    if not isinstance(bump, CompactMomentBump):
        raise TypeError("bump must be a CompactMomentBump")

    radius = math.sqrt(2.0 * X)
    if radius <= bump.left or radius >= bump.right:
        return 0.0, 0.0, 0.0

    width = bump.right - bump.left
    s = (radius - bump.left) / width
    u = s * (1.0 - s)
    # g(s)=4-1/[s(1-s)] for the normalized bump b=exp(g)/normalization.
    g_s = (1.0 - 2.0 * s) / (u * u)
    g_ss = -2.0 / (u * u) - 2.0 * (1.0 - 2.0 * s) ** 2 / (u**3)

    value = float(bump(radius))
    d_radius = value * g_s / width
    d_radius2 = value * (g_s * g_s + g_ss) / (width * width)

    d_X = d_radius / radius
    d_X2 = d_radius2 / (radius * radius) - d_radius / (radius**3)
    if not all(math.isfinite(v) for v in (value, d_X, d_X2)):
        raise OverflowError("compact repair bump second jet is outside binary64 range")
    return value, d_X, d_X2


def _provider_jet(
    provider: BaseSecondJetProvider,
    X: float,
    eta: float,
    name: str,
) -> ProfileSecondJet:
    value = provider(X, eta)
    if not isinstance(value, ProfileSecondJet):
        raise TypeError(f"{name} must return ProfileSecondJet")
    return value


@dataclass(frozen=True)
class Lemma52RepairedSecondJetAdapter:
    """Lift the function-level Eq. (5.14) repair to second ``(X,eta)`` jets.

    ``base_U_second_jet`` and ``base_E_second_jet`` describe the unrepaired
    coefficient.  The compact correction is then differentiated exactly:

    ``d_X^a d_eta^b [c_j(eta) b_j(sqrt(2X))]``

    is the product of the corresponding coefficient derivative and analytic
    radial bump derivative for ``a,b <= 2`` with ``a+b <= 2``.
    """

    profile: Lemma52RepairedProfileAdapter
    base_U_second_jet: BaseSecondJetProvider
    base_E_second_jet: BaseSecondJetProvider
    provenance: str = (
        "Lemma 5.2 repaired U/E second-jet bridge from caller-supplied unrepaired "
        "second jets and hierarchy moment/patch jets; formal-structure only, not "
        "a paper-exact recursive coefficient."
    )

    def __post_init__(self) -> None:
        if not isinstance(self.profile, Lemma52RepairedProfileAdapter):
            raise TypeError("profile must be a Lemma52RepairedProfileAdapter")
        if not callable(self.base_U_second_jet):
            raise TypeError("base_U_second_jet must be callable")
        if not callable(self.base_E_second_jet):
            raise TypeError("base_E_second_jet must be callable")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty string")

    @staticmethod
    def _has_active_bump(X: float, bumps: tuple[CompactMomentBump, ...]) -> bool:
        radius = math.sqrt(2.0 * X)
        return any(bump.left < radius < bump.right for bump in bumps)

    @staticmethod
    def _combine(
        base: ProfileSecondJet,
        coefficient_derivatives: object,
        bumps: tuple[CompactMomentBump, ...],
        X: float,
    ) -> ProfileSecondJet:
        rows = coefficient_derivatives
        # ``repair_coefficients`` already returns a finite ndarray with rows
        # 0,1,2.  Keep the indexing explicit so eta-derivative order cannot be
        # silently confused with factorial-normalized Taylor coefficients.
        radial_jets = tuple(_bump_x_second_jet(bump, X) for bump in bumps)

        value_correction = math.fsum(
            float(rows[0, j]) * radial_jets[j][0] for j in range(len(bumps))
        )
        radial_correction = math.fsum(
            float(rows[0, j]) * radial_jets[j][1] for j in range(len(bumps))
        )
        radial2_correction = math.fsum(
            float(rows[0, j]) * radial_jets[j][2] for j in range(len(bumps))
        )
        parameter_correction = math.fsum(
            float(rows[1, j]) * radial_jets[j][0] for j in range(len(bumps))
        )
        radial_parameter_correction = math.fsum(
            float(rows[1, j]) * radial_jets[j][1] for j in range(len(bumps))
        )
        parameter2_correction = math.fsum(
            float(rows[2, j]) * radial_jets[j][0] for j in range(len(bumps))
        )

        return ProfileSecondJet(
            value=base.value + value_correction,
            radial=base.radial + radial_correction,
            radial2=base.radial2 + radial2_correction,
            parameter=base.parameter + parameter_correction,
            radial_parameter=base.radial_parameter + radial_parameter_correction,
            parameter2=base.parameter2 + parameter2_correction,
        )

    def U_second_jet(self, X: float, eta: float) -> ProfileSecondJet:
        """Return the repaired axial-coefficient second jet at ``(X,eta)``."""

        X, eta = self.profile.base_profile._point(X, eta)
        base = _provider_jet(self.base_U_second_jet, X, eta, "base_U_second_jet")
        bumps = tuple(self.profile.repair.u_bumps)
        if not self._has_active_bump(X, bumps):
            # The C-infinity bumps are flat off their compact supports, so no
            # repair coefficient derivatives are required at this point.
            return base
        coeff = self.profile.repair_coefficients(eta, max_order=2)
        return self._combine(base, coeff.alpha_derivatives, bumps, X)

    def E_second_jet(self, X: float, eta: float) -> ProfileSecondJet:
        """Return the repaired swirl-coefficient second jet at ``(X,eta)``."""

        X, eta = self.profile.base_profile._point(X, eta)
        base = _provider_jet(self.base_E_second_jet, X, eta, "base_E_second_jet")
        bumps = tuple(self.profile.repair.e_bumps)
        if not self._has_active_bump(X, bumps):
            return base
        coeff = self.profile.repair_coefficients(eta, max_order=2)
        return self._combine(base, coeff.beta_derivatives, bumps, X)
