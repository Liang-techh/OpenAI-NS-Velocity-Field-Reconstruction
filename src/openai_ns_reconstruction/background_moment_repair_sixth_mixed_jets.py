"""Analytic sixth-mixed U-jet bridge for the Lemma 5.2 compact repair.

The third-eta Eq. (5.6) ``Omega/X`` path ultimately needs a fifth-mixed
``beta=V/X`` jet. Through Eq. (5.2), that introduces exactly three new
total-order-six axial derivatives: ``U_XXetaetaetaeta``,
``U_Xetaetaetaetaeta`` and ``U_etaetaetaetaetaeta``.

The Lemma 5.2 compact repair already solves arbitrary finite eta jets of its
coefficient functions, and its compact radial bumps expose exact X derivatives
through second order. Lower fifth-mixed U data are delegated to the landed
adapter, so the new layer cannot silently fork lower semantics.

The unrepaired sixth-mixed U jet and moment/patch eta jets remain explicit
upstream inputs. This is fail-closed Stage-2 ``formal-structure``
infrastructure, not paper-exact coefficient data or an all-order result.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import math

from .background_moment_repair import CompactMomentBump
from .background_moment_repair_fifth_mixed_jets import Lemma52RepairedFifthMixedJetAdapter
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter
from .background_moment_repair_second_jets import _bump_x_second_jet
from .background_regular_flux_fifth_mixed_jets import (
    AxialSixthMixedJet,
    RegularFluxFifthMixedJet,
    regular_flux_fifth_mixed_jet_eq_5_2,
)

BaseSixthMixedJetProvider = Callable[[float, float], AxialSixthMixedJet]


def _provider_jet(provider: BaseSixthMixedJetProvider, X: float, eta: float) -> AxialSixthMixedJet:
    value = provider(X, eta)
    if not isinstance(value, AxialSixthMixedJet):
        raise TypeError("base_U_sixth_mixed_jet must return AxialSixthMixedJet")
    return value


@dataclass(frozen=True)
class Lemma52RepairedSixthMixedJetAdapter:
    """Lift one compact Lemma 5.2 U repair to the sixth-mixed Eq. (5.2) jet."""

    profile: Lemma52RepairedProfileAdapter
    base_U_sixth_mixed_jet: BaseSixthMixedJetProvider
    provenance: str = (
        "Lemma 5.2 repaired axial sixth-mixed jet from caller-supplied unrepaired "
        "sixth-mixed U data and moment/patch eta-jets; formal-structure only, "
        "not a paper-exact recursive coefficient."
    )

    def __post_init__(self) -> None:
        if not isinstance(self.profile, Lemma52RepairedProfileAdapter):
            raise TypeError("profile must be a Lemma52RepairedProfileAdapter")
        if not callable(self.base_U_sixth_mixed_jet):
            raise TypeError("base_U_sixth_mixed_jet must be callable")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty string")

    @staticmethod
    def _has_active_bump(X: float, bumps: tuple[CompactMomentBump, ...]) -> bool:
        radius = math.sqrt(2.0 * X)
        return any(bump.left < radius < bump.right for bump in bumps)

    def _fifth_bridge(self) -> Lemma52RepairedFifthMixedJetAdapter:
        return Lemma52RepairedFifthMixedJetAdapter(
            self.profile,
            lambda X, eta: _provider_jet(self.base_U_sixth_mixed_jet, X, eta).fifth(),
        )

    def U_sixth_mixed_jet(self, X: float, eta: float) -> AxialSixthMixedJet:
        """Return the compact-repaired axial sixth-mixed jet."""
        X, eta = self.profile.base_profile._point(X, eta)
        base = _provider_jet(self.base_U_sixth_mixed_jet, X, eta)
        bumps = tuple(self.profile.repair.u_bumps)
        if not self._has_active_bump(X, bumps):
            return base

        lower = self._fifth_bridge().U_fifth_mixed_jet(X, eta)
        coeff = self.profile.repair_coefficients(eta, max_order=6)
        rows = coeff.alpha_derivatives
        radial_jets = tuple(_bump_x_second_jet(bump, X) for bump in bumps)

        def correction(eta_order: int, radial_order: int) -> float:
            return math.fsum(
                float(rows[eta_order, j]) * radial_jets[j][radial_order]
                for j in range(len(bumps))
            )

        return AxialSixthMixedJet(
            value=lower.value, radial=lower.radial, radial2=lower.radial2,
            parameter=lower.parameter, radial_parameter=lower.radial_parameter,
            parameter2=lower.parameter2, radial2_parameter=lower.radial2_parameter,
            radial_parameter2=lower.radial_parameter2, parameter3=lower.parameter3,
            radial2_parameter2=lower.radial2_parameter2,
            radial_parameter3=lower.radial_parameter3, parameter4=lower.parameter4,
            radial2_parameter3=lower.radial2_parameter3,
            radial_parameter4=lower.radial_parameter4, parameter5=lower.parameter5,
            radial2_parameter4=base.radial2_parameter4 + correction(4, 2),
            radial_parameter5=base.radial_parameter5 + correction(5, 1),
            parameter6=base.parameter6 + correction(6, 0),
        )

    def beta_fifth_mixed_jet(
        self, h: float, order: int, X: float, eta: float, *, quadrature_points: int = 32,
    ) -> RegularFluxFifthMixedJet:
        """Feed repaired sixth-mixed U into the analytic Eq. (5.2) fifth jet."""
        return regular_flux_fifth_mixed_jet_eq_5_2(
            h, order, X, eta, self.U_sixth_mixed_jet,
            quadrature_points=quadrature_points,
        )
