"""Analytic third-mixed U-jet bridge for the Lemma 5.2 compact repair.

The landed Eq. (5.2) regular-flux adapter needs the repaired axial coefficient
through the three total-order-three derivatives ``U_XXeta``, ``U_Xetaeta`` and
``U_etaetaeta``.  The function-level Lemma 5.2 repair already exposes arbitrary
finite eta-jets of its repair coefficients, while the second-jet bridge already
contains the exact ``R=sqrt(2X)`` bump differentiation.  This module combines
those two landed pieces without introducing finite-difference production data.

The unrepaired third-mixed U jet and the hierarchy moment/patch eta-jets remain
explicit inputs.  Therefore this is solver infrastructure / ``formal-structure``
only; it does not claim a converged Section-5 coefficient or paper-exact status.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import math

from .background_moment_repair import CompactMomentBump
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter
from .background_moment_repair_second_jets import _bump_x_second_jet
from .background_regular_flux_second_jets import (
    AxialThirdMixedJet,
    regular_flux_second_jet_eq_5_2,
)
from .background_lower_history_source import ProfileSecondJet

BaseThirdMixedJetProvider = Callable[[float, float], AxialThirdMixedJet]


def _provider_jet(
    provider: BaseThirdMixedJetProvider,
    X: float,
    eta: float,
) -> AxialThirdMixedJet:
    value = provider(X, eta)
    if not isinstance(value, AxialThirdMixedJet):
        raise TypeError("base_U_third_mixed_jet must return AxialThirdMixedJet")
    return value


@dataclass(frozen=True)
class Lemma52RepairedThirdMixedJetAdapter:
    """Lift the compact Lemma 5.2 axial repair to the Eq. (5.2) jet interface.

    For a correction ``sum_j alpha_j(eta) b_j(sqrt(2X))``, the additional
    derivatives required by Eq. (5.2) are exactly

    ``U_XXeta    += sum_j alpha'_j   d_X^2 b_j``,
    ``U_Xetaeta  += sum_j alpha''_j  d_X b_j``, and
    ``U_etaetaeta += sum_j alpha'''_j b_j``.

    The first six fields are corrected consistently as well, so the returned
    :class:`AxialThirdMixedJet` is a single coherent repaired jet.
    """

    profile: Lemma52RepairedProfileAdapter
    base_U_third_mixed_jet: BaseThirdMixedJetProvider
    provenance: str = (
        "Lemma 5.2 repaired axial third-mixed jet from caller-supplied unrepaired "
        "third-mixed U data and hierarchy moment/patch eta-jets; formal-structure "
        "only, not a paper-exact recursive coefficient."
    )

    def __post_init__(self) -> None:
        if not isinstance(self.profile, Lemma52RepairedProfileAdapter):
            raise TypeError("profile must be a Lemma52RepairedProfileAdapter")
        if not callable(self.base_U_third_mixed_jet):
            raise TypeError("base_U_third_mixed_jet must be callable")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty string")

    @staticmethod
    def _has_active_bump(X: float, bumps: tuple[CompactMomentBump, ...]) -> bool:
        radius = math.sqrt(2.0 * X)
        return any(bump.left < radius < bump.right for bump in bumps)

    @staticmethod
    def _combine(
        base: AxialThirdMixedJet,
        coefficient_derivatives: object,
        bumps: tuple[CompactMomentBump, ...],
        X: float,
    ) -> AxialThirdMixedJet:
        rows = coefficient_derivatives
        radial_jets = tuple(_bump_x_second_jet(bump, X) for bump in bumps)

        def correction(eta_order: int, radial_order: int) -> float:
            return math.fsum(
                float(rows[eta_order, j]) * radial_jets[j][radial_order]
                for j in range(len(bumps))
            )

        return AxialThirdMixedJet(
            value=base.value + correction(0, 0),
            radial=base.radial + correction(0, 1),
            radial2=base.radial2 + correction(0, 2),
            parameter=base.parameter + correction(1, 0),
            radial_parameter=base.radial_parameter + correction(1, 1),
            parameter2=base.parameter2 + correction(2, 0),
            radial2_parameter=base.radial2_parameter + correction(1, 2),
            radial_parameter2=base.radial_parameter2 + correction(2, 1),
            parameter3=base.parameter3 + correction(3, 0),
        )

    def U_third_mixed_jet(self, X: float, eta: float) -> AxialThirdMixedJet:
        """Return the repaired axial jet needed by the Eq. (5.2) beta adapter."""

        X, eta = self.profile.base_profile._point(X, eta)
        base = _provider_jet(self.base_U_third_mixed_jet, X, eta)
        bumps = tuple(self.profile.repair.u_bumps)
        if not self._has_active_bump(X, bumps):
            # Every C-infinity repair bump is flat off its compact support, so
            # unavailable third eta coefficient data are genuinely unnecessary.
            return base
        coeff = self.profile.repair_coefficients(eta, max_order=3)
        return self._combine(base, coeff.alpha_derivatives, bumps, X)

    def beta_second_jet(
        self,
        h: float,
        order: int,
        X: float,
        eta: float,
        *,
        quadrature_points: int = 32,
    ) -> ProfileSecondJet:
        """Feed the repaired U jet directly into the landed exact Eq. (5.2) path."""

        return regular_flux_second_jet_eq_5_2(
            h,
            order,
            X,
            eta,
            self.U_third_mixed_jet,
            quadrature_points=quadrature_points,
        )
