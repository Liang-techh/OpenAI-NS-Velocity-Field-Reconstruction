"""Analytic fourth-mixed U-jet bridge for the Lemma 5.2 compact repair.

The landed Eq. (5.2) regular-flux third-mixed adapter needs the repaired axial
coefficient through exactly three additional total-order-four derivatives:
``U_XXetaeta``, ``U_Xetaetaeta`` and ``U_etaetaetaeta``.  The function-level
Lemma 5.2 repair already exposes arbitrary finite eta-jets of its repair
coefficients, while the landed radial bump bridge already provides the exact
``R=sqrt(2X)`` derivatives through ``d_X^2``.

This module composes those two pieces.  It does not finite-difference the repair
in production and it does not infer missing hierarchy data.  The unrepaired
fourth-mixed U jet plus the moment/patch eta-jets remain explicit upstream
inputs, so this is fail-closed ``formal-structure`` solver infrastructure rather
than a paper-exact recursive coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import math

from .background_moment_repair import CompactMomentBump
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter
from .background_moment_repair_second_jets import _bump_x_second_jet
from .background_omega_parameter_jet import RegularFluxThirdMixedJet
from .background_regular_flux_third_mixed_jets import (
    AxialFourthMixedJet,
    regular_flux_third_mixed_jet_eq_5_2,
)

BaseFourthMixedJetProvider = Callable[[float, float], AxialFourthMixedJet]


def _provider_jet(
    provider: BaseFourthMixedJetProvider,
    X: float,
    eta: float,
) -> AxialFourthMixedJet:
    value = provider(X, eta)
    if not isinstance(value, AxialFourthMixedJet):
        raise TypeError("base_U_fourth_mixed_jet must return AxialFourthMixedJet")
    return value


@dataclass(frozen=True)
class Lemma52RepairedFourthMixedJetAdapter:
    """Lift the compact Lemma 5.2 U repair to the fourth-mixed Eq. (5.2) jet.

    For a correction ``sum_j alpha_j(eta) b_j(sqrt(2X))``, the three new
    derivatives beyond the landed third-mixed adapter are exactly

    ``U_XXetaeta    += sum_j alpha_j''     d_X^2 b_j``,
    ``U_Xetaetaeta  += sum_j alpha_j'''    d_X b_j``, and
    ``U_etaetaetaeta += sum_j alpha_j''''   b_j``.

    All lower fields are corrected from the same coefficient jet, so the
    returned :class:`AxialFourthMixedJet` is coherent and can be passed directly
    to :func:`regular_flux_third_mixed_jet_eq_5_2`.
    """

    profile: Lemma52RepairedProfileAdapter
    base_U_fourth_mixed_jet: BaseFourthMixedJetProvider
    provenance: str = (
        "Lemma 5.2 repaired axial fourth-mixed jet from caller-supplied unrepaired "
        "fourth-mixed U data and hierarchy moment/patch eta-jets; formal-structure "
        "only, not a paper-exact recursive coefficient."
    )

    def __post_init__(self) -> None:
        if not isinstance(self.profile, Lemma52RepairedProfileAdapter):
            raise TypeError("profile must be a Lemma52RepairedProfileAdapter")
        if not callable(self.base_U_fourth_mixed_jet):
            raise TypeError("base_U_fourth_mixed_jet must be callable")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty string")

    @staticmethod
    def _has_active_bump(X: float, bumps: tuple[CompactMomentBump, ...]) -> bool:
        radius = math.sqrt(2.0 * X)
        return any(bump.left < radius < bump.right for bump in bumps)

    @staticmethod
    def _combine(
        base: AxialFourthMixedJet,
        coefficient_derivatives: object,
        bumps: tuple[CompactMomentBump, ...],
        X: float,
    ) -> AxialFourthMixedJet:
        rows = coefficient_derivatives
        radial_jets = tuple(_bump_x_second_jet(bump, X) for bump in bumps)

        def correction(eta_order: int, radial_order: int) -> float:
            return math.fsum(
                float(rows[eta_order, j]) * radial_jets[j][radial_order]
                for j in range(len(bumps))
            )

        return AxialFourthMixedJet(
            value=base.value + correction(0, 0),
            radial=base.radial + correction(0, 1),
            radial2=base.radial2 + correction(0, 2),
            parameter=base.parameter + correction(1, 0),
            radial_parameter=base.radial_parameter + correction(1, 1),
            parameter2=base.parameter2 + correction(2, 0),
            radial2_parameter=base.radial2_parameter + correction(1, 2),
            radial_parameter2=base.radial_parameter2 + correction(2, 1),
            parameter3=base.parameter3 + correction(3, 0),
            radial2_parameter2=base.radial2_parameter2 + correction(2, 2),
            radial_parameter3=base.radial_parameter3 + correction(3, 1),
            parameter4=base.parameter4 + correction(4, 0),
        )

    def U_fourth_mixed_jet(self, X: float, eta: float) -> AxialFourthMixedJet:
        """Return the repaired axial jet needed by the beta third-mixed adapter."""

        X, eta = self.profile.base_profile._point(X, eta)
        base = _provider_jet(self.base_U_fourth_mixed_jet, X, eta)
        bumps = tuple(self.profile.repair.u_bumps)
        if not self._has_active_bump(X, bumps):
            # The compact C-infinity repair is exactly flat off support, so no
            # fourth eta coefficient row is needed there.
            return base
        coeff = self.profile.repair_coefficients(eta, max_order=4)
        return self._combine(base, coeff.alpha_derivatives, bumps, X)

    def beta_third_mixed_jet(
        self,
        h: float,
        order: int,
        X: float,
        eta: float,
        *,
        quadrature_points: int = 32,
    ) -> RegularFluxThirdMixedJet:
        """Feed the repaired U fourth-mixed jet into the analytic Eq. (5.2) map."""

        return regular_flux_third_mixed_jet_eq_5_2(
            h,
            order,
            X,
            eta,
            self.U_fourth_mixed_jet,
            quadrature_points=quadrature_points,
        )
